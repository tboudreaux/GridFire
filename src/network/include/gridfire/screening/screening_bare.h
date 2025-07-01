#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "cppad/cppad.hpp"

namespace gridfire::screening {
    class BareScreeningModel final : public ScreeningModel {
    using ADDouble = CppAD::AD<double>;
    public:
        [[nodiscard]] std::vector<double> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const override;

        [[nodiscard]] std::vector<ADDouble> calculateScreeningFactors(
            const reaction::LogicalReactionSet& reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<ADDouble>& Y,
            const ADDouble T9,
            const ADDouble rho
        ) const override;
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

    template<typename T>
    std::vector<T> BareScreeningModel::calculateFactors_impl(
        const reaction::LogicalReactionSet &reactions,
        const std::vector<fourdst::atomic::Species> &species,
        const std::vector<T> &Y,
        const T T9,
        const T rho
    ) const {
        return std::vector<T>(reactions.size(), T(1.0)); // Bare screening returns 1.0 for all reactions
    }
}
