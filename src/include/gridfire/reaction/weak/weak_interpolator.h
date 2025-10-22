#pragma once

#include "gridfire/reaction/weak/weak_types.h"
#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/logging/logging.h"

#include <unordered_map>
#include <cstdint>
#include <vector>
#include <expected>
#include <array>



namespace gridfire::rates::weak {
    /**
     * @class WeakRateInterpolator
     * @brief 3D table interpolator for tabulated weak reaction data by isotope.
     *
     * Builds per-isotope 3D grids over (T9, log10(rho*Ye), mu_e) and provides:
     *  - Trilinear interpolation of the tabulated log10(rate) and neutrino-loss fields
     *    into a WeakRatePayload via get_rates().
     *  - Finite-difference estimates of partial derivatives via get_rate_derivatives().
     *
     * Implementation summary (constructor): rows are grouped by (A,Z), then each group's unique
     * axis values are collected and sorted to form the three axes; the 3D payload array is
     * populated at each lattice point with the 6 log10() fields from the raw table.
     */
    class WeakRateInterpolator {
    public:
        /**
         * @brief Raw weak-rate table type expected by the constructor.
         *
         * The size must match the number of rows compiled into the weak-rate library.
         */
        using RowDataTable = std::array<RateDataRow, 77400>; // Total number of entries in the weak rate table NOTE: THIS MUST EQUAL THE VALUE IN weak_rate_library.h

        /**
         * @brief Construct the interpolator from raw weak-rate rows.
         *
         * Groups rows by isotope (A,Z), extracts unique sorted axes for T9, log10(rho*Ye), and mu_e,
         * and fills an internal regular grid with the log10(rate) and neutrino-loss payloads at each node.
         * No interpolation occurs at construction time.
         */
        explicit WeakRateInterpolator(const RowDataTable& raw_data);

        /**
         * @brief List isotopes for which tables are available.
         * @return Vector of available Species (A,Z) derived from internal tables.
         * @throws std::runtime_error If any packed (A,Z) cannot be converted to Species.
         * @par Example
         * @code
         * WeakRateInterpolator interp(rows);
         * auto isotopes = interp.available_isotopes();
         * @endcode
         */
        [[nodiscard]] std::vector<fourdst::atomic::Species> available_isotopes() const;

        /**
         * @brief Trilinear interpolation of weak-rate payload at a state.
         *
         * Interpolates the 6 log10() fields (rates and neutrino losses) at the given state
         * for the requested isotope. If the isotope is unknown or the state lies outside
         * the tabulated ranges, returns an error via std::expected with detailed bounds info.
         *
         * @param A Mass number of the isotope.
         * @param Z Proton number of the isotope.
         * @param t9 Temperature in GK (10^9 K).
         * @param log_rhoYe Log10 of rho*Ye (cgs density times electron fraction).
         * @return expected<WeakRatePayload, InterpolationError>: payload on success;
         *         InterpolationError::UNKNOWN_SPECIES_ERROR if (A,Z) not present; or
         *         InterpolationError::BOUNDS_ERROR if any coordinate is outside the table
         *         (with per-axis bounds included).
         * @par Example
         * @code
         * if (auto res = interp.get_rates(52, 26, 3.0, 6.0, 2.0); res) {
         *   const WeakRatePayload& p = *res;
         * } else {
         *   // inspect res.error().type and optional bounds info
         * }
         * @endcode
         */
        [[nodiscard]] std::expected<WeakRatePayload, InterpolationError> get_rates(
            uint16_t A,
            uint8_t Z,
            double t9,
            double log_rhoYe
        ) const;

        /**
         * @brief Finite-difference partial derivatives of the log10() fields.
         *
         * Uses central differences with small fixed (1e-6) perturbations in each variable
         * (T9, log10(rho*Ye), mu_e) and returns arrays of d(log10(field))/d(var) for all fields.
         * If any perturbed state falls outside the table, returns a BOUNDS_ERROR with per-axis
         * bounds; if the isotope is unknown, returns UNKNOWN_SPECIES_ERROR.
         *
         * @param A Mass number of the isotope.
         * @param Z Proton number of the isotope.
         * @param t9 Temperature in GK (10^9 K).
         * @param log_rhoYe Log10 of rho*Ye (cgs density times electron fraction).
         * @return expected<WeakRateDerivatives, InterpolationError>: derivative payload on success;
         *         otherwise an InterpolationError as described above.
         * @par Example
         * @code
         * if (auto d = interp.get_rate_derivatives(52, 26, 3.0, 6.0, 2.0); d) {
         *   // use d->d_log_beta_minus[0..2], etc.
         * }
         * @endcode
         */
        [[nodiscard]] std::expected<WeakRateDerivatives, InterpolationError> get_rate_derivatives(
            uint16_t A,
            uint8_t Z,
            double t9,
            double log_rhoYe
        ) const;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        /**
         * @brief Per-isotope grids over (T9, log10(rho*Ye), mu_e) with payloads at lattice nodes.
         */
        std::unordered_map<uint32_t, IsotopeGrid> m_rate_table;
    };


}