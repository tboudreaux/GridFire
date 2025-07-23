#include "gridfire/solver/solver.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h> // Needed for std::function

#include <vector>

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
