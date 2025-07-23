#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include "bindings.h"

#include "gridfire/solver/solver.h"
#include "trampoline/py_solver.h"

namespace py = pybind11;


void register_solver_bindings(py::module &m) {
    auto py_dynamic_network_solving_strategy = py::class_<gridfire::solver::DynamicNetworkSolverStrategy, PyDynamicNetworkSolverStrategy>(m, "DynamicNetworkSolverStrategy");
    auto py_direct_network_solver = py::class_<gridfire::solver::DirectNetworkSolver, gridfire::solver::DynamicNetworkSolverStrategy>(m, "DirectNetworkSolver");

    py_direct_network_solver.def(py::init<gridfire::DynamicEngine&>(),
        py::arg("engine"),
        "Constructor for the DirectNetworkSolver. Takes a DynamicEngine instance to use for evaluating the network.");

    py_direct_network_solver.def("evaluate",
        &gridfire::solver::DirectNetworkSolver::evaluate,
        py::arg("netIn"),
        "Evaluate the network for a given timestep. Returns the output conditions after the timestep.");
}

