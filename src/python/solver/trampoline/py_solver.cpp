#include "gridfire/solver/solver.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <vector>
#include <tuple>
#include <string>
#include <any>

#include "py_solver.h"


namespace py = pybind11;

gridfire::NetOut PyDynamicNetworkSolverStrategy::evaluate(const gridfire::NetIn &netIn) {
    PYBIND11_OVERRIDE_PURE(
        gridfire::NetOut,  // Return type
        gridfire::solver::DynamicNetworkSolverStrategy,  // Base class
        evaluate,  // Method name
        netIn  // Arguments
    );
}

void PyDynamicNetworkSolverStrategy::set_callback(const std::any &callback) {
    PYBIND11_OVERRIDE_PURE(
        void,
        gridfire::solver::DynamicNetworkSolverStrategy,  // Base class
        set_callback,  // Method name
        callback  // Arguments
    );
}

std::vector<std::tuple<std::string, std::string>> PyDynamicNetworkSolverStrategy::describe_callback_context() const {
    using DescriptionVector = std::vector<std::tuple<std::string, std::string>>;
    PYBIND11_OVERRIDE_PURE(
        DescriptionVector,  // Return type
        gridfire::solver::DynamicNetworkSolverStrategy,  // Base class
        describe_callback_context  // Method name
    );
}
