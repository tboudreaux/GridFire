#pragma once

#include "gridfire/reaction/reaction.h"

#include "fourdst/composition/atomicSpecies.h"

#include "cppad/cppad.hpp"

#include <vector>

namespace gridfire::screening {
    class ScreeningModel {
    public:
        using ADDouble = CppAD::AD<double>;
        virtual ~ScreeningModel() = default;

        virtual std::vector<double> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
            ) const = 0;

        virtual std::vector<ADDouble> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble>& Y,
            const ADDouble T9,
            const ADDouble rho
            ) const = 0;
    };
}