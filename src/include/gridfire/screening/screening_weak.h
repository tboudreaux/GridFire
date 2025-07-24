#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"

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
         * @param Y A vector of the molar abundances (mol/g) for each species.
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
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
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
         * @param Y A vector of the molar abundances as AD types.
         * @param T9 The temperature as an AD type.
         * @param rho The plasma density as an AD type.
         * @return A vector of screening factors as AD types.
         */
        [[nodiscard]] std::vector<CppAD::AD<double>> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<CppAD::AD<double>>& Y,
            const CppAD::AD<double> T9,
            const CppAD::AD<double> rho
        ) const override;
    private:
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
         * @param Y A vector of molar abundances.
         * @param T9 The temperature in 10^9 K.
         * @param rho The density in g/cm^3.
         * @return A vector of screening factors of type `T`.
         */
        template <typename T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            const T T9,
            const T rho
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
     * @param Y The molar abundances of the species.
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
    template <typename T>
    std::vector<T> WeakScreeningModel::calculateFactors_impl(
        const reaction::LogicalReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<T>& Y,
        const T T9,
        const T rho
    ) const {
        LOG_TRACE_L3(
            m_logger,
            "Calculating weak screening factors for {} reactions...",
            reactions.size()
        );
        // --- CppAD Safe low temp checking ---
        const T zero(0.0);
        const T one(1.0);
        const T low_temp_threshold(1e-9);

        const T low_T_flag = CppAD::CondExpLt(T9, low_temp_threshold, zero, one);

        // --- Calculate composition-dependent terms ---
        // ζ = ∑(Z_i^2 + Z_i) * X_i / A_i
        // This simplifies somewhat to ζ = ∑ (Z_i^2 + Z_i) * Y_i
        // Where Y_i is the molar abundance (mol/g)
        T zeta(0.0);
        for (size_t i = 0; i < species.size(); ++i) {
            const T Z = species[i].m_z;
            zeta += (Z * Z + Z) * Y[i];
        }

        // --- Constant prefactors ---
        const T T7 = T9 * static_cast<T>(100.00);
        const T T7_safe = CppAD::CondExpLe(T7, low_temp_threshold, low_temp_threshold, T7);
        const T prefactor = static_cast<T>(0.188) * CppAD::sqrt(rho / (T7_safe * T7_safe * T7_safe)) * CppAD::sqrt(zeta);

        // --- Loop through reactions and calculate screening factors for each ---
        std::vector<T> factors;
        factors.reserve(reactions.size());
        for (const auto& reaction : reactions) {
            T H_12(0.0); // screening abundance term
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
                LOG_TRACE_L3(m_logger, "Calculating screening factor for reaction: {}", reaction.peName());
                const T Z1 = static_cast<T>(reactants[0].m_z);
                const T Z2 = static_cast<T>(reactants[1].m_z);
                H_12 = prefactor * Z1 * Z2;
            }
            else if (isTripleAlpha) {
                LOG_TRACE_L3(m_logger, "Special case for triple alpha process in reaction: {}", reaction.peName());
                // Special case for triple alpha process
                const T Z_alpha = static_cast<T>(2.0);
                const T H_alpha_alpha = prefactor * Z_alpha * Z_alpha;
                H_12 = static_cast<T>(3.0) * H_alpha_alpha; // Triple alpha process
            }
            // For 1 body reactions H_12 remains 0 so e^H_12 will be 1.0 (screening does not apply)
            // Aside from triple alpha, all other astrophysically relevant reactions are 2-body in the weak screening regime

            H_12 *= low_T_flag; // Apply low temperature flag to screening factor
            H_12 = CppAD::CondExpGe(H_12, static_cast<T>(2.0), static_cast<T>(2.0), H_12); // Caps the screening factor at 10 to avoid numerical issues
            factors.push_back(CppAD::exp(H_12));
            // std::cout << "Screening factor: " << reaction.peName() << " : " << factors.back() << "(" << H_12 << ")" << std::endl;
        }
        return factors;
    }

}