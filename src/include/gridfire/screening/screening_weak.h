#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/types/types.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

#include "cppad/cppad.hpp"

namespace gridfire::screening {
    /**
     * @class WeakScreeningModel
     * @brief Implements the weak screening model based on the Debye-Hückel approximation.
     *
     * This class provides a concrete implementation of the `ScreeningModel`
     * interface for the weak screening regime, following the formulation of
     * Salpeter (1954). This approach applies the Debye-Hückel theory to model the
     * electrostatic shielding of nuclei in a plasma. It is applicable to
     * non-degenerate, non-relativistic plasmas where thermal energy dominates
     * the electrostatic potential energy.
     *
     * @implements ScreeningModel
     */
    class WeakScreeningModel final : public ScreeningModel {
    public:
        /**
         * @brief Calculates weak screening factors for a set of reactions.
         *
         * This method computes the screening enhancement factor for each reaction
         * based on the Salpeter (1954) formula.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species involved in the network.
         * @param Y The composition object giving the current molar abundances.
         * @param T9 The temperature in units of 10^9 K.
         * @param rho The plasma density in g/cm^3.
         * @return A vector of screening factors (dimensionless), one for each reaction.
         *
         * @b Usage
         * @code
         * WeakScreeningModel weak_model;
         * // ... (initialize reactions, species, Y, T9, rho)
         * std::vector<double> factors = weak_model.calculateScreeningFactors(
         *     reactions, species, Y, T9, rho
         * );
         * @endcode
         */
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Calculates weak screening factors using CppAD types.
         *
         * This is the automatic differentiation-compatible version of the method.
         * It allows the derivatives of the screening factors to be computed with
         * respect to plasma conditions.
         *
         * @param reactions The set of logical reactions in the network.
         * @param species A vector of all atomic species involved in the network.
         * @param Y The composition object giving the current molar abundances.
         * @param T9 The temperature as an AD type.
         * @param rho The plasma density as an AD type.
         * @return A vector of screening factors as AD types.
         */
        [[nodiscard]] std::vector<CppAD::AD<double>> calculateScreeningFactors(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<CppAD::AD<double>> &Y,
            CppAD::AD<double> T9,
            CppAD::AD<double> rho
        ) const override;
        [[nodiscard]] double calculateScreeningFactor(
            const reaction::Reaction& reaction,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /**
         * @brief Analytic d(ln f_r)/dY_k for the weak screening factors.
         *
         * ln f_r = H_r, with H_r proportional to sqrt(zeta) and zeta = sum_k (Z_k^2 + Z_k) Y_k, so
         * d(ln f_r)/dY_k = H_r (Z_k^2 + Z_k) / (2 zeta). It is zero where H_r is capped, below the low-temperature
         * cutoff, and when zeta <= 0 (the factor is then pinned at 1).
         */
        [[nodiscard]] std::optional<std::vector<double>> calculateScreeningLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        /// ln f = H_12 with H_12 proportional to rho^(1/2) T^(-3/2), so d ln f/d ln T = -3/2 H and d ln f/d ln rho = H/2
        /// (zero where the factor is capped or the reaction is unscreened).
        [[nodiscard]] std::optional<std::vector<StateLogDerivatives>> calculateScreeningStateLogDerivatives(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;
    private:
        /// @brief Salpeter prefactor 0.188 zeta^(1/2) rho^(1/2) T_6^(-3/2), times 0 below the low-temperature cutoff.
        template <IsArithmeticOrAD T>
        [[nodiscard]] static T prefactor_impl(const std::vector<fourdst::atomic::Species>& species, const std::vector<T>& Y, T T9, T rho);

        /// @brief Uncapped H_12 for one reaction given the prefactor (0 for one-body reactions).
        template <IsArithmeticOrAD T>
        [[nodiscard]] static T H_impl(const reaction::Reaction& reaction, T prefactor);

        /// @brief exp(min(H_12, H_cap)) for one reaction given the prefactor.
        template <IsArithmeticOrAD T>
        [[nodiscard]] static T factor_impl(const reaction::Reaction& reaction, T prefactor);

        /// @brief Largest H_12 applied (factor e^2), to avoid numerical issues outside the weak regime.
        static constexpr double H_CAP = 2.0;

        /// @brief Logger instance for recording trace and debug information.
        [[maybe_unused]] quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");

    private:
        /**
         * @brief Template implementation for calculating weak screening factors.
         *
         * This private helper function contains the core logic for calculating
         * weak screening factors. It is templated to handle both `double` and
         * `CppAD::AD<double>` numeric types, avoiding code duplication.
         *
         * @tparam T The numeric type, either `double` or `CppAD::AD<double>`.
         * @param reactions The set of reactions.
         * @param species A vector of all species in the network.
         * @param Y The composition object with current molar abundances.
         * @param T9 The temperature in 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A vector of screening factors of type `T`.
         */
        template <IsArithmeticOrAD T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::ReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            T T9,
            T rho
        ) const;
    };

    /**
     * @brief Core implementation of the weak screening calculation (Debye-Hückel model).
     *
     * This function calculates the screening factor `exp(H_12)` for each reaction,
     * based on the Debye-Hückel approximation as formulated by Salpeter (1954).
     *
     * @tparam T The numeric type (`double` or `CppAD::AD<double>`).
     * @param reactions The set of reactions to be screened.
     * @param species The list of all species in the network.
     * @param Y The composition object providing current molar abundances.
     * @param T9 The temperature in 10^9 K.
     * @param rho The density in g/cm^3.
     * @return A vector of screening factors, one for each reaction.
     *
     * @b Algorithm
     * 1.  **Low-Temperature Cutoff**: If T9 is below a small threshold (1e-9),
     *     screening is effectively turned off to prevent numerical instability.
     * 2.  **Zeta Factor (ζ)**: A composition-dependent term is calculated:
     *     `ζ = ∑(Z_i² + Z_i) * Y_i`, where Z_i is the charge and Y_i is the
     *     molar abundance of species i.
     * 3.  **Prefactor**: A key prefactor is computed:
     *     `prefactor = 0.188 * sqrt(ρ / T₇³) * sqrt(ζ)`,
     *     where T₇ is the temperature in units of 10^7 K.
     * 4.  **Screening Term (H_12)**: For each reaction, the term H_12 is calculated:
     *     - For a two-body reaction (reactants Z₁ and Z₂): `H_12 = prefactor * Z₁ * Z₂`.
     *     - For the triple-alpha reaction (3 * He4): `H_12 = 3 * (prefactor * Z_α * Z_α)`.
     *     - For one-body reactions (decays), H_12 is 0, so the factor is 1.
     * 5.  **Capping**: The value of H_12 is capped at 2.0 to prevent excessively large
     *     and unphysical screening factors (exp(2) ≈ 7.4).
     * 6.  **Final Factor**: The screening factor for the reaction is `exp(H_12)`.
     */
    template <IsArithmeticOrAD T>
    T WeakScreeningModel::prefactor_impl(
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<T>& Y,
        const T T9,
        const T rho
    ) {
        // --- CppAD safe low temperature check ---
        const T zero(0.0);
        const T one(1.0);
        const T low_temp_threshold(1e-9);
        const T low_T_flag = CppAD::CondExpLt(T9, low_temp_threshold, zero, one);

        // --- Composition-dependent term ---
        // zeta = sum_i (Z_i^2 + Z_i) X_i / A_i = sum_i (Z_i^2 + Z_i) Y_i, with Y_i the molar abundance (mol/g)
        T zeta(0.0);
        for (size_t i = 0; i < species.size(); ++i) {
            const T Z = species[i].m_z;
            zeta += (Z * Z + Z) * Y[i];
        }

        // Salpeter (1954): H_12 = 0.188 Z_1 Z_2 zeta^(1/2) rho^(1/2) T_6^(-3/2) (Clayton 1968, eq. 4-221). The coefficient
        // is for T in units of 10^6 K; evaluating it with T_7 overstated H by 10^1.5 (screening factors ~6 instead of
        // ~1.06 for p + p in the solar core).
        const T T6 = T9 * static_cast<T>(1000.0);
        const T T6_safe = CppAD::CondExpLe(T6, low_temp_threshold, low_temp_threshold, T6);
        return low_T_flag * static_cast<T>(0.188) * CppAD::sqrt(rho / (T6_safe * T6_safe * T6_safe)) * CppAD::sqrt(zeta);
    }

    template <IsArithmeticOrAD T>
    T WeakScreeningModel::H_impl(const reaction::Reaction& reaction, const T prefactor) {
        const auto& reactants = reaction.reactants();
        const bool isTripleAlpha = (
            reactants.size() == 3 &&
            reactants[0].m_z == 2 &&
            reactants[1].m_z == 2 &&
            reactants[2].m_z == 2 &&
            reactants[0] == reactants[1] &&
            reactants[1] == reactants[2]
        );
        if (reactants.size() == 2) {
            return prefactor * static_cast<T>(reactants[0].m_z) * static_cast<T>(reactants[1].m_z);
        }
        if (isTripleAlpha) {
            // alpha + alpha, then alpha + 8Be: H = H(2,2) + H(2,4) = 3 H(2,2)
            return static_cast<T>(3.0) * prefactor * static_cast<T>(4.0);
        }
        // One-body reactions are not screened
        return static_cast<T>(0.0);
    }

    template <IsArithmeticOrAD T>
    T WeakScreeningModel::factor_impl(const reaction::Reaction& reaction, const T prefactor) {
        const T H_12 = H_impl<T>(reaction, prefactor);
        const T cap = static_cast<T>(H_CAP);
        return CppAD::exp(CppAD::CondExpGe(H_12, cap, cap, H_12));
    }

    template <IsArithmeticOrAD T>
    std::vector<T> WeakScreeningModel::calculateFactors_impl(
        const reaction::ReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<T>& Y,
        const T T9,
        const T rho
    ) const {
        LOG_TRACE_L3(m_logger, "Calculating weak screening factors for {} reactions...", reactions.size());
        const T prefactor = prefactor_impl<T>(species, Y, T9, rho);
        std::vector<T> factors;
        factors.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            factors.push_back(factor_impl<T>(*reaction, prefactor));
        }
        return factors;
    }

}