#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

#include "cppad/cppad.hpp"

namespace gridfire::screening {
    class WeakScreeningModel final : public ScreeningModel {
    public:
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        [[nodiscard]] std::vector<CppAD::AD<double>> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<CppAD::AD<double>>& Y,
            const CppAD::AD<double> T9,
            const CppAD::AD<double> rho
        ) const override;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");

    private:
        template <typename T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            const T T9,
            const T rho
        ) const;
    };

    template <typename T>
    std::vector<T> WeakScreeningModel::calculateFactors_impl(
        const reaction::LogicalReactionSet& reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<T>& Y,
        const T T9,
        const T rho
    ) const {
        LOG_TRACE_L1(
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