#include "gridfire/screening/screening_weak.h"

#include "fourdst/composition/atomicSpecies.h"

#include "cppad/cppad.hpp"

#include <vector>


namespace gridfire::screening {
    using ADDouble = CppAD::AD<double>;
    std::vector<ADDouble> WeakScreeningModel::calculateScreeningFactors(
        const reaction::LogicalReactionSet &reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<ADDouble> &Y,
        const ADDouble T9,
        const ADDouble rho
    ) const {
        return calculateFactors_impl<ADDouble>(reactions, species, Y, T9, rho);
    }

    std::vector<double> WeakScreeningModel::calculateScreeningFactors(
        const reaction::LogicalReactionSet &reactions,
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double> &Y,
        const double T9,
        const double rho
    ) const {
        return calculateFactors_impl<double>(reactions, species, Y, T9, rho);
    }
}
