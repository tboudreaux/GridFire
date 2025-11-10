#include "gridfire/reaction/weak/weak_interpolator.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/weak/weak.h"
#include "gridfire/utils/hashing.h"

#include <algorithm>
#include <set>
#include <unordered_map>
#include <vector>
#include <optional>
#include <expected>
#include <ranges>

#include "fourdst/atomic/species.h"

#include "quill/LogMacros.h"

namespace gridfire::rates::weak {

    WeakRateInterpolator::WeakRateInterpolator(const RowDataTable &raw_data) {
        // Group all raw data rows by their isotope ID.
        std::unordered_map<uint32_t, std::vector<const RateDataRow*>> grouped_rows;
        for (const auto& row : raw_data) {
            grouped_rows[utils::hash_atomic(row.A, row.Z)].push_back(&row);
        }

        // Process each isotope's data to build a simple 2D grid.
        for (auto const& [isotope_id, rows] : grouped_rows) {
            IsotopeGrid grid;

            // Establish the T9 and log(rho*Ye) axes
            std::set<float> unique_t9, unique_rhoYe;
            for (const auto* row : rows) {
                unique_t9.emplace(row->t9);
                unique_rhoYe.emplace(row->log_rhoye);
            }

            grid.t9_axis.assign(unique_t9.begin(), unique_t9.end());
            grid.rhoYe_axis.assign(unique_rhoYe.begin(), unique_rhoYe.end());

            const size_t nt9 = grid.t9_axis.size();
            const size_t nrhoYe = grid.rhoYe_axis.size();

            grid.data.resize(nt9 * nrhoYe);

            // Create reverse maps for efficient index lookups.
            std::unordered_map<double, size_t> t9_map, rhoYe_map;
            for (size_t i = 0; i < nt9; i++) { t9_map[grid.t9_axis[i]] = i; }
            for (size_t j = 0; j < nrhoYe; j++) { rhoYe_map[grid.rhoYe_axis[j]] = j; }

            // Populate the 2D grid.
            for (const auto* row: rows) {
                size_t i_t9 = t9_map.at(row->t9);
                size_t j_rhoYe = rhoYe_map.at(row->log_rhoye);

                size_t index = i_t9 * nrhoYe + j_rhoYe;
                grid.data[index] = WeakRatePayload{
                    row->log_beta_plus,
                    row->log_electron_capture,
                    row->log_neutrino_loss_ec,
                    row->log_beta_minus,
                    row->log_positron_capture,
                    row->log_antineutrino_loss_bd
                };
            }
            m_rate_table[isotope_id] = std::move(grid);
        }
    }



    std::vector<fourdst::atomic::Species> WeakRateInterpolator::available_isotopes() const {
        using namespace fourdst::atomic;
        std::vector<Species> isotopes;
        for (const auto &packed_id: m_rate_table | std::views::keys) {
            const auto A = static_cast<uint16_t>(packed_id >> 8);
            const auto Z = static_cast<uint8_t>(packed_id & 0xFF);
            std::expected<Species, SpeciesErrorType> result = az_to_species(A, Z);
            if (!result.has_value()) {
                std::string msg = "Could not convert A=" + std::to_string(A) + ", Z=" + std::to_string(Z) + " to Species: ";
                msg += (result.error() == SpeciesErrorType::ELEMENT_SYMBOL_NOT_FOUND) ? "Unknown element (Z out of range)." : "Invalid isotope (A < Z or A out of range).";
                LOG_TRACE_L3(m_logger, "{}", msg);
            } else {
                isotopes.emplace_back(result.value());
            }
        }
        return isotopes;
    }

        std::expected<WeakRatePayload, InterpolationError> WeakRateInterpolator::get_rates(
        const uint16_t A,
        const uint8_t Z,
        const double t9,
        const double log_rhoYe
    ) const {
        const auto it = m_rate_table.find(utils::hash_atomic(A, Z));
        if (it == m_rate_table.end()) {
            return std::unexpected(InterpolationError{InterpolationErrorType::UNKNOWN_SPECIES_ERROR});
        }
        const auto& grid = it->second;
        const auto& t9_axis = grid.t9_axis;
        const auto& rhoYe_axis = grid.rhoYe_axis;

        // Find bracketing indices for the 2D (t9, rhoYe) grid
        auto find_lower_index = [](const std::vector<double>& axis, const double value) -> std::optional<size_t> {
            const auto upperBoundIterator = std::ranges::upper_bound(axis, value);
            if (upperBoundIterator == axis.begin() || upperBoundIterator == axis.end()) {
                return std::nullopt; // Out of bounds
            }
            return std::distance(axis.begin(), upperBoundIterator) - 1;
        };

        const auto i_t9_opt = find_lower_index(t9_axis, t9);
        const auto j_rhoYe_opt = find_lower_index(rhoYe_axis, log_rhoYe);

        // Handle bounds errors for the 2D grid
        if (!i_t9_opt || !j_rhoYe_opt) {
            std::unordered_map<TableAxes, BoundsErrorInfo> boundsInfo;
            if (!i_t9_opt.has_value()) {
                boundsInfo[TableAxes::T9] = BoundsErrorInfo{TableAxes::T9, t9_axis.front(), t9_axis.back(), t9};
            }
            if (!j_rhoYe_opt.has_value()) {
                boundsInfo[TableAxes::LOG_RHOYE] = BoundsErrorInfo{TableAxes::LOG_RHOYE, rhoYe_axis.front(), rhoYe_axis.back(), log_rhoYe};
            }
            return std::unexpected(InterpolationError{InterpolationErrorType::BOUNDS_ERROR, boundsInfo});
        }

        const size_t i = i_t9_opt.value();
        const size_t j = j_rhoYe_opt.value();
        const size_t nrhoYe = rhoYe_axis.size();

        // Get the four corner payloads for the bilinear interpolation
        const auto& p00 = grid.data[(i * nrhoYe) + j];
        const auto& p01 = grid.data[(i * nrhoYe) + j + 1];
        const auto& p10 = grid.data[((i + 1) * nrhoYe) + j];
        const auto& p11 = grid.data[((i + 1) * nrhoYe) + j + 1];

        // Fractional distances for the 2D bilinear interpolation
        const double td = (t9 - t9_axis[i]) / (t9_axis[i + 1] - t9_axis[i]);
        const double rd = (log_rhoYe - rhoYe_axis[j]) / (rhoYe_axis[j + 1] - rhoYe_axis[j]);

        // Helper lambda to linearly interpolate between two full payloads
        auto lerp_payload = [](const WeakRatePayload& p0, const WeakRatePayload& p1, double t) {
            return WeakRatePayload{
                .log_beta_plus = std::lerp(p0.log_beta_plus, p1.log_beta_plus, t),
                .log_electron_capture = std::lerp(p0.log_electron_capture, p1.log_electron_capture, t),
                .log_neutrino_loss_ec = std::lerp(p0.log_neutrino_loss_ec, p1.log_neutrino_loss_ec, t),
                .log_beta_minus = std::lerp(p0.log_beta_minus, p1.log_beta_minus, t),
                .log_positron_capture = std::lerp(p0.log_positron_capture, p1.log_positron_capture, t),
                .log_antineutrino_loss_bd = std::lerp(p0.log_antineutrino_loss_bd, p1.log_antineutrino_loss_bd, t),
            };
        };

        // Perform the bilinear interpolation
        const WeakRatePayload p0 = lerp_payload(p00, p01, rd);
        const WeakRatePayload p1 = lerp_payload(p10, p11, rd);

        return lerp_payload(p0, p1, td);
    }

    std::expected<WeakRateDerivatives, InterpolationError> WeakRateInterpolator::get_rate_derivatives(
        uint16_t A,
        uint8_t Z,
        double t9,
        double log_rhoYe
    ) const {
        WeakRateDerivatives result{};
        //TODO: Make this perturbation scale aware
        constexpr double eps = 1e-6; // Small perturbation for finite difference

        // Perturbations for finite difference
        const double h_t9 = (t9 > 1e-9) ? t9 * eps : eps;
        const auto payload_plus_t9 = get_rates(A, Z, t9 + h_t9, log_rhoYe);
        const auto payload_minus_t9 = get_rates(A, Z, t9 - h_t9, log_rhoYe);

        const double h_rhoYe = (std::abs(log_rhoYe) > 1e-9) ? std::abs(log_rhoYe) * eps : eps;
        const auto payload_plus_rhoYe = get_rates(A, Z, t9, log_rhoYe + h_rhoYe);
        const auto payload_minus_rhoYe = get_rates(A, Z, t9, log_rhoYe - h_rhoYe);


        if (!payload_plus_t9 || !payload_minus_t9 || !payload_plus_rhoYe || !payload_minus_rhoYe) {
            // Determine which perturbation failed and return a consolidated error
            auto first_error = !payload_plus_t9 ? payload_plus_t9.error() :
                               !payload_minus_t9 ? payload_minus_t9.error() :
                               !payload_plus_rhoYe ? payload_plus_rhoYe.error() :
                               payload_minus_rhoYe.error();
            return std::unexpected(first_error);
        }

        // Derivatives wrt. T9
        const double t9_denominator = 2 * h_t9;
        result.d_log_beta_plus[0] = (payload_plus_t9->log_beta_plus - payload_minus_t9->log_beta_plus) / t9_denominator;
        result.d_log_beta_minus[0] = (payload_plus_t9->log_beta_minus - payload_minus_t9->log_beta_minus) / t9_denominator;
        result.d_log_electron_capture[0] = (payload_plus_t9->log_electron_capture - payload_minus_t9->log_electron_capture) / t9_denominator;
        result.d_log_neutrino_loss_ec[0] = (payload_plus_t9->log_neutrino_loss_ec - payload_minus_t9->log_neutrino_loss_ec) / t9_denominator;
        result.d_log_positron_capture[0] = (payload_plus_t9->log_positron_capture - payload_minus_t9->log_positron_capture) / t9_denominator;
        result.d_log_antineutrino_loss_bd[0] = (payload_plus_t9->log_antineutrino_loss_bd - payload_minus_t9->log_antineutrino_loss_bd) / t9_denominator;

        // Derivatives wrt. logRhoYe
        const double rhoYe_denominator = 2 * h_rhoYe;
        result.d_log_beta_plus[1] = (payload_plus_rhoYe->log_beta_plus - payload_minus_rhoYe->log_beta_plus) / rhoYe_denominator;
        result.d_log_beta_minus[1] = (payload_plus_rhoYe->log_beta_minus - payload_minus_rhoYe->log_beta_minus) / rhoYe_denominator;
        result.d_log_electron_capture[1] = (payload_plus_rhoYe->log_electron_capture - payload_minus_rhoYe->log_electron_capture) / rhoYe_denominator;
        result.d_log_neutrino_loss_ec[1] = (payload_plus_rhoYe->log_neutrino_loss_ec - payload_minus_rhoYe->log_neutrino_loss_ec) / rhoYe_denominator;
        result.d_log_positron_capture[1] = (payload_plus_rhoYe->log_positron_capture - payload_minus_rhoYe->log_positron_capture) / rhoYe_denominator;
        result.d_log_antineutrino_loss_bd[1] = (payload_plus_rhoYe->log_antineutrino_loss_bd - payload_minus_rhoYe->log_antineutrino_loss_bd) / rhoYe_denominator;

        return result;
    }

}
