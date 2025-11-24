#pragma once

#include "gridfire/screening/screening.h"

#include <vector>

#include "cppad/cppad.hpp"

class PyScreening final : public gridfire::screening::ScreeningModel {
    [[nodiscard]] std::vector<double> calculateScreeningFactors(
        const gridfire::reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species> &species,
        const std::vector<double> &Y,
        double T9,
        double rho
    ) const override;

    [[nodiscard]] std::vector<ADDouble> calculateScreeningFactors(
        const gridfire::reaction::ReactionSet &reactions,
        const std::vector<fourdst::atomic::Species> &species,
        const std::vector<ADDouble> &Y,
        ADDouble T9,
        ADDouble rho
    ) const override;
};