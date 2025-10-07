#pragma once

#include "gridfire/screening/screening_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "cppad/cppad.hpp"

namespace gridfire::screening {
    class IntermediateScreeningModel final : public ScreeningModel {
    public:
        std::vector<double> calculateScreeningFactors(
            const reaction::ReactionSet &reactions,
            const std::vector<fourdst::atomic::Species> &species,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        std::vector<ADDouble> calculateScreeningFactors(
            const reaction::ReactionSet &reactions,
            const std::vector<fourdst::atomic::Species> &species,
            const std::vector<ADDouble> &Y,
            ADDouble T9,
            ADDouble rho
        ) const override;

    private:
        template<typename T>
        [[nodiscard]] std::vector<T> calculateFactors_impl(
            const reaction::ReactionSet &reactions,
            const std::vector<fourdst::atomic::Species>& species,
            const std::vector<T>& Y,
            T T9,
            T rho
        ) const;
    };

    template<typename T>
    std::vector<T> IntermediateScreeningModel::calculateFactors_impl(
        const reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species> &species,
        const std::vector<T>& Y,
        const T T9,
        const T rho
    ) const {
        // TODO: Implement the intermediate screening model logic here. Follow the prescription from Graboske et al. (1973)
        return std::vector<T>(species.size(), 0);
    }
}
