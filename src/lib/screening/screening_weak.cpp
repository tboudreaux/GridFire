#include "gridfire/screening/screening_weak.h"

#include "fourdst/atomic/atomicSpecies.h"

#include "cppad/cppad.hpp"

#include <optional>
#include <vector>


namespace gridfire::screening {
    using ADDouble = CppAD::AD<double>;
    std::vector<ADDouble> WeakScreeningModel::calculateScreeningFactors(
        const reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<screening::ADDouble> &Y,
        const ADDouble T9,
        const ADDouble rho
    ) const {
        return calculateFactors_impl<ADDouble>(reactions, species, Y, T9, rho);
    }

    std::vector<double> WeakScreeningModel::calculateScreeningFactors(
        const reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return calculateFactors_impl<double>(reactions, species, Y, T9, rho);
    }

    double WeakScreeningModel::calculateScreeningFactor(
        const reaction::Reaction &reaction,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return factor_impl<double>(reaction, prefactor_impl<double>(species, Y, T9, rho));
    }

    std::optional<std::vector<double>> WeakScreeningModel::calculateScreeningLogDerivatives(
        const reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        std::vector<double> dlnf_dY(reactions.size() * species.size(), 0.0);
        double zeta = 0.0;
        for (size_t k = 0; k < species.size(); ++k) {
            const double Z = species[k].m_z;
            zeta += (Z * Z + Z) * Y[k];
        }
        if (!(zeta > 0.0)) {
            return dlnf_dY; // sqrt(zeta) has no finite derivative here, and every factor is pinned at 1
        }
        const double prefactor = prefactor_impl<double>(species, Y, T9, rho); // 0 below the low-temperature cutoff
        size_t r = 0;
        for (const auto& reaction : reactions) {
            // Same H as calculateScreeningFactors; the capped branch is constant in Y
            if (const double H_12 = H_impl<double>(*reaction, prefactor); H_12 > 0.0 && H_12 < H_CAP) {
                double* row = dlnf_dY.data() + r * species.size();
                for (size_t k = 0; k < species.size(); ++k) {
                    const double Z = species[k].m_z;
                    row[k] = H_12 * (Z * Z + Z) / (2.0 * zeta);
                }
            }
            ++r;
        }
        return dlnf_dY;
    }
}
