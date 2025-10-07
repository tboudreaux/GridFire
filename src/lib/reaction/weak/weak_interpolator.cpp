#include "gridfire/reaction/weak/weak_interpolator.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/weak/weak.h"

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <optional>
#include <expected>
#include <ranges>

#include "fourdst/composition/species.h"

namespace gridfire::rates::weak {

    WeakRateInterpolator::WeakRateInterpolator(const RowDataTable &raw_data) {
        std::map<uint32_t, std::vector<const RateDataRow*>> grouped_rows;
        for (const auto& row : raw_data) {
            grouped_rows[pack_isotope_id(row.A, row.Z)].push_back(&row);
        }

        for (auto const& [isotope_id, rows] : grouped_rows) {
            IsotopeGrid grid;

            std::set<float> unique_t9, unique_rhoYe, unique_mue;
            for (const auto* row : rows) {
                unique_t9.emplace(row->t9);
                unique_rhoYe.emplace(row->log_rhoye);
                unique_mue.emplace(row->mu_e);
            }

            grid.t9_axis.reserve(unique_t9.size());
            grid.rhoYe_axis.reserve(unique_rhoYe.size());
            grid.mue_axis.reserve(unique_mue.size());

            grid.t9_axis.insert(grid.t9_axis.begin(), unique_t9.begin(), unique_t9.end());
            grid.rhoYe_axis.insert(grid.rhoYe_axis.begin(), unique_rhoYe.begin(), unique_rhoYe.end());
            grid.mue_axis.insert(grid.mue_axis.begin(), unique_mue.begin(), unique_mue.end());

            std::ranges::sort(grid.t9_axis);
            std::ranges::sort(grid.rhoYe_axis);
            std::ranges::sort(grid.mue_axis);

            const size_t nt9 = grid.t9_axis.size();
            const size_t nrhoYe = grid.rhoYe_axis.size();
            const size_t nmue = grid.mue_axis.size();

            grid.data.resize(nt9 * nrhoYe * nmue);

            // Reverse map for quick index lookup
            std::unordered_map<float, size_t> t9_map, rhoYe_map, mue_map;
            for (size_t i = 0; i < nt9; i++) { t9_map[grid.t9_axis[i]] = i; }
            for (size_t j = 0; j < nrhoYe; j++) { rhoYe_map[grid.rhoYe_axis[j]] = j; }
            for (size_t k = 0; k < nmue; k++) { mue_map[grid.mue_axis[k]] = k; }

            for (const auto* row: rows) {
                size_t i_t9 = t9_map.at(row->t9);
                size_t j_rhoYe = rhoYe_map.at(row->log_rhoye);
                size_t k_mue = mue_map.at(row->mu_e);

                size_t index = (i_t9 * nrhoYe + j_rhoYe) * nmue + k_mue;
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
        std::vector<fourdst::atomic::Species> isotopes;
        for (const auto &packed_id: m_rate_table | std::views::keys) {
            const uint16_t A = static_cast<uint16_t>(packed_id >> 8);
            const uint8_t Z = static_cast<uint8_t>(packed_id & 0xFF);
            try {
                fourdst::atomic::Species species = fourdst::atomic::az_to_species(A, Z);
                isotopes.push_back(species);
            } catch (const std::exception& e) {
                throw std::runtime_error("Error converting A=" + std::to_string(A) + ", Z=" + std::to_string(Z) + " to Species: " + e.what());
            }
        }
        return isotopes;
    }

    std::expected<WeakRatePayload, InterpolationError> WeakRateInterpolator::get_rates(
        const uint16_t A,
        const uint8_t Z,
        const double t9,
        const double log_rhoYe,
        const double mu_e
    ) const {
        const auto it = m_rate_table.find(pack_isotope_id(A, Z));
        if (it == m_rate_table.end()) {
            return std::unexpected(InterpolationError{InterpolationErrorType::UNKNOWN_SPECIES_ERROR});
        }
        const auto&[t9_axis, rhoYe_axis, mue_axis, data] = it->second;

        // Now find the bracketing indices for t9, log_rhoYe, and mu_e
        auto find_lower_index = [](const std::vector<double>& axis, const double value) -> std::optional<size_t> {
            const auto upperBoundIterator = std::ranges::upper_bound(axis, value);
            if (upperBoundIterator == axis.begin() || upperBoundIterator == axis.end()) {
                return std::nullopt; // Out of bounds
            }
            return std::distance(axis.begin(), upperBoundIterator) - 1;
        };

        const auto i_t9_opt = find_lower_index(t9_axis, t9);
        const auto j_rhoYe_opt = find_lower_index(rhoYe_axis, log_rhoYe);
        const auto k_mue_opt = find_lower_index(mue_axis, mu_e);

        if (!i_t9_opt || !j_rhoYe_opt || !k_mue_opt) {
            std::unordered_map<TableAxes, BoundsErrorInfo> boundsInfo;
            if (!i_t9_opt) {
                boundsInfo[TableAxes::T9] = BoundsErrorInfo{
                    TableAxes::T9,
                    t9_axis.front(),
                    t9_axis.back(),
                    t9
                };
            }
            if (!j_rhoYe_opt) {
                boundsInfo[TableAxes::LOG_RHOYE] = BoundsErrorInfo{
                    TableAxes::LOG_RHOYE,
                    rhoYe_axis.front(),
                    rhoYe_axis.back(),
                    log_rhoYe
                };
            }
            if (!k_mue_opt) {
                boundsInfo[TableAxes::MUE] = BoundsErrorInfo{
                    TableAxes::MUE,
                    mue_axis.front(),
                    mue_axis.back(),
                    mu_e
                };
            }
            return std::unexpected(
                InterpolationError{
                    InterpolationErrorType::BOUNDS_ERROR,
                    boundsInfo
                }
            );
        }

        const size_t i = i_t9_opt.value();
        const size_t j = j_rhoYe_opt.value();
        const size_t k = k_mue_opt.value();

        // Coordinates of the bounding cube
        const double t1 = t9_axis[i];
        const double t2 = t9_axis[i + 1];
        const double r1 = rhoYe_axis[j];
        const double r2 = rhoYe_axis[j + 1];
        const double m1 = mue_axis[k];
        const double m2 = mue_axis[k + 1];

        const double td = (t9 - t1) / (t2 - t1);
        const double rd = (log_rhoYe - r1) / (r2 - r1);
        const double md = (mu_e - m1) / (m2 - m1);

        auto lerp = [](const double v0, const double v1, const double t) {
          return v0 * (1 - t) + v1 * t;
        };

        auto interpolationField = [&](auto field_accessor) {
            const size_t nrhoYe = rhoYe_axis.size();
            const size_t nmue = mue_axis.size();

            auto get_val = [&](const size_t i_t, const size_t j_r, const size_t k_m) {
                return field_accessor(data[(i_t * nrhoYe + j_r) * nmue + k_m]);
            };

            const double c000 = get_val(i, j, k);
            const double c001 = get_val(i, j, k + 1);
            const double c010 = get_val(i, j + 1, k);
            const double c011 = get_val(i, j + 1, k + 1);
            const double c100 = get_val(i + 1, j, k);
            const double c101 = get_val(i + 1, j, k + 1);
            const double c110 = get_val(i + 1, j + 1, k);
            const double c111 = get_val(i + 1, j + 1, k + 1);

            const double c00 = lerp(c000, c001, md);
            const double c01 = lerp(c010, c011, md);
            const double c10 = lerp(c100, c101, md);
            const double c11 = lerp(c110, c111, md);

            const double c0 = lerp(c00, c01, rd);
            const double c1 = lerp(c10, c11, rd);

            return lerp(c0, c1, td);

        };

        WeakRatePayload result;

        result.log_beta_plus = interpolationField([](const WeakRatePayload& p) { return p.log_beta_plus; });
        result.log_electron_capture = interpolationField([](const WeakRatePayload& p) { return p.log_electron_capture; });
        result.log_neutrino_loss_ec = interpolationField([](const WeakRatePayload& p) { return p.log_neutrino_loss_ec; });
        result.log_beta_minus = interpolationField([](const WeakRatePayload& p) { return p.log_beta_minus; });
        result.log_positron_capture = interpolationField([](const WeakRatePayload& p) { return p.log_positron_capture; });
        result.log_antineutrino_loss_bd = interpolationField([](const WeakRatePayload& p) { return p.log_antineutrino_loss_bd; });
        return result;
    }

    std::expected<WeakRateDerivatives, InterpolationError> WeakRateInterpolator::get_rate_derivatives(
        uint16_t A,
        uint8_t Z,
        double t9,
        double log_rhoYe,
        double mu_e
    ) const {
        WeakRateDerivatives result;
        constexpr double eps = 1e-6; // Small perturbation for finite difference

        // Perturbations for finite difference
        const double h_t9 = (t9 > 1e-9) ? t9 * eps : eps;
        const auto payload_plus_t9 = get_rates(A, Z, t9 + h_t9, log_rhoYe, mu_e);
        const auto payload_minus_t9 = get_rates(A, Z, t9 - h_t9, log_rhoYe, mu_e);

        const double h_rhoYe = (std::abs(log_rhoYe) > 1e-9) ? std::abs(log_rhoYe) * eps : eps;
        const auto payload_plus_rhoYe = get_rates(A, Z, t9, log_rhoYe + h_rhoYe, mu_e);
        const auto payload_minus_rhoYe = get_rates(A, Z, t9, log_rhoYe - h_rhoYe, mu_e);

        const double h_mue = (std::abs(mu_e) > 1e-9) ? std::abs(mu_e) * eps : eps;
        const auto payload_plus_mue = get_rates(A, Z, t9, log_rhoYe, mu_e + h_mue);
        const auto payload_minus_mue = get_rates(A, Z, t9, log_rhoYe, mu_e - h_mue);

        if (!payload_plus_t9 || !payload_minus_t9 || !payload_plus_rhoYe || !payload_minus_rhoYe || !payload_plus_mue || !payload_minus_mue) {
            const auto it = m_rate_table.find(pack_isotope_id(A, Z));
            if (it == m_rate_table.end()) {
                return std::unexpected(InterpolationError{InterpolationErrorType::UNKNOWN_SPECIES_ERROR});
            }

            const IsotopeGrid& grid = it->second;
            InterpolationError error;
            std::unordered_map<TableAxes, BoundsErrorInfo> boundsInfo;
            if (!payload_minus_t9 || !payload_plus_t9) {
                boundsInfo[TableAxes::T9] = BoundsErrorInfo{
                    TableAxes::T9,
                    grid.t9_axis.front(),
                    grid.t9_axis.back(),
                    t9
                };
            }
            if (!payload_minus_rhoYe || !payload_plus_rhoYe) {
                boundsInfo[TableAxes::LOG_RHOYE] = BoundsErrorInfo{
                    TableAxes::LOG_RHOYE,
                    grid.rhoYe_axis.front(),
                    grid.rhoYe_axis.back(),
                    log_rhoYe
                };
            }
            if (!payload_minus_mue || !payload_plus_mue) {
                boundsInfo[TableAxes::MUE] = BoundsErrorInfo{
                    TableAxes::MUE,
                    grid.mue_axis.front(),
                    grid.mue_axis.back(),
                    mu_e
                };
            }
            error.type = InterpolationErrorType::BOUNDS_ERROR;
            error.boundsErrorInfo = boundsInfo;
            return std::unexpected(error);
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

        // Derivatives wrt. MuE
        const double mue_denominator = 2 * h_mue;
        result.d_log_beta_plus[2] = (payload_plus_mue->log_beta_plus - payload_minus_mue->log_beta_plus) / mue_denominator;
        result.d_log_beta_minus[2] = (payload_plus_mue->log_beta_minus - payload_minus_mue->log_beta_minus) / mue_denominator;
        result.d_log_electron_capture[2] = (payload_plus_mue->log_electron_capture - payload_minus_mue->log_electron_capture) / mue_denominator;
        result.d_log_neutrino_loss_ec[2] = (payload_plus_mue->log_neutrino_loss_ec - payload_minus_mue->log_neutrino_loss_ec) / mue_denominator;
        result.d_log_positron_capture[2] = (payload_plus_mue->log_positron_capture - payload_minus_mue->log_positron_capture) / mue_denominator;
        result.d_log_antineutrino_loss_bd[2] = (payload_plus_mue->log_antineutrino_loss_bd - payload_minus_mue->log_antineutrino_loss_bd) / mue_denominator;

        return result;
    }


    uint32_t WeakRateInterpolator::pack_isotope_id(const uint16_t A, const uint8_t Z) {
        return (static_cast<uint32_t>(A) << 8) | static_cast<uint32_t>(Z);
    }

}
