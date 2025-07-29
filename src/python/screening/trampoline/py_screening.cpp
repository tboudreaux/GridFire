#include "gridfire/screening/screening.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <vector>

#include "py_screening.h"

#include "cppad/cppad.hpp"


namespace py = pybind11;

std::vector<double> PyScreening::calculateScreeningFactors(const gridfire::reaction::LogicalReactionSet &reactions, const std::vector<fourdst::atomic::Species> &species, const std::vector<double> &Y, const double T9, const double rho) const {
    PYBIND11_OVERLOAD_PURE(
        std::vector<double>,  // Return type
        gridfire::screening::ScreeningModel,
        calculateScreeningFactors // Method name
    );
}

using ADDouble = gridfire::screening::ScreeningModel::ADDouble;
std::vector<ADDouble> PyScreening::calculateScreeningFactors(const gridfire::reaction::LogicalReactionSet &reactions, const std::vector<fourdst::atomic::Species> &species, const std::vector<ADDouble> &Y, const ADDouble T9, const ADDouble rho) const {
    PYBIND11_OVERLOAD_PURE(
        std::vector<ADDouble>,  // Return type
        gridfire::screening::ScreeningModel,
        calculateScreeningFactors // Method name
    );
}
