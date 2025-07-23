#pragma once

#include "gridfire/screening/screening.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <vector>

#include "cppad/cppad.hpp"

class PyScreening final : public gridfire::screening::ScreeningModel {
    std::vector<double> calculateScreeningFactors(const gridfire::reaction::LogicalReactionSet &reactions, const std::vector<fourdst::atomic::Species> &species, const std::vector<double> &Y, const double T9, const double rho) const override;
    std::vector<ADDouble> calculateScreeningFactors(const gridfire::reaction::LogicalReactionSet &reactions, const std::vector<fourdst::atomic::Species> &species, const std::vector<ADDouble> &Y, const ADDouble T9, const ADDouble rho) const override;
};