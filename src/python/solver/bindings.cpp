#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly
#include <pybind11/numpy.h>

#include <boost/numeric/ublas/vector.hpp>

#include "bindings.h"

#include "gridfire/solver/solver.h"
#include "trampoline/py_solver.h"

namespace py = pybind11;


void register_solver_bindings(const py::module &m) {
    auto py_dynamic_network_solving_strategy = py::class_<gridfire::solver::DynamicNetworkSolverStrategy, PyDynamicNetworkSolverStrategy>(m, "DynamicNetworkSolverStrategy");
    auto py_direct_network_solver = py::class_<gridfire::solver::DirectNetworkSolver, gridfire::solver::DynamicNetworkSolverStrategy>(m, "DirectNetworkSolver");

    py_direct_network_solver.def(py::init<gridfire::DynamicEngine&>(),
        py::arg("engine"),
        "Constructor for the DirectNetworkSolver. Takes a DynamicEngine instance to use for evaluating the network."
    );

    py_direct_network_solver.def("evaluate",
        &gridfire::solver::DirectNetworkSolver::evaluate,
        py::arg("netIn"),
        "Evaluate the network for a given timestep. Returns the output conditions after the timestep."
    );

    py_direct_network_solver.def("set_callback",
        [](gridfire::solver::DirectNetworkSolver &self, const gridfire::solver::DirectNetworkSolver::TimestepCallback& cb) {
            self.set_callback(cb);
        },
        py::arg("callback"),
        "Sets a callback function to be called at each timestep."
    );

    py::class_<gridfire::solver::DirectNetworkSolver::TimestepContext>(py_direct_network_solver, "TimestepContext")
        .def_readonly("t", &gridfire::solver::DirectNetworkSolver::TimestepContext::t, "Current time in the simulation.")
        .def_property_readonly(
            "state", [](const gridfire::solver::DirectNetworkSolver::TimestepContext& ctx) {
                std::vector<double> state(ctx.state.size());
                std::ranges::copy(ctx.state, state.begin());
                return py::array_t<double>(static_cast<ssize_t>(state.size()), state.data());
            })
        .def_readonly("dt", &gridfire::solver::DirectNetworkSolver::TimestepContext::dt, "Current timestep size.")
        .def_readonly("cached_time", &gridfire::solver::DirectNetworkSolver::TimestepContext::cached_time, "Cached time for the last computed result.")
        .def_readonly("last_observed_time", &gridfire::solver::DirectNetworkSolver::TimestepContext::last_observed_time, "Last time the state was observed.")
        .def_readonly("last_step_time", &gridfire::solver::DirectNetworkSolver::TimestepContext::last_step_time, "Last step time taken for the integration.")
        .def_readonly("T9", &gridfire::solver::DirectNetworkSolver::TimestepContext::T9, "Temperature in units of 10^9 K.")
        .def_readonly("rho", &gridfire::solver::DirectNetworkSolver::TimestepContext::rho, "Temperature in units of 10^9 K.")
        .def_property_readonly("cached_result", [](const gridfire::solver::DirectNetworkSolver::TimestepContext& ctx) -> py::object {
            if (ctx.cached_result.has_value()) {
                const auto&[dydt, nuclearEnergyGenerationRate] = ctx.cached_result.value();
                return py::make_tuple(
                    py::array_t<double>(static_cast<ssize_t>(dydt.size()), dydt.data()),
                    nuclearEnergyGenerationRate
                );
            }
            return py::none();
        }, "Cached result of the step derivatives.")
        .def_readonly("num_steps", &gridfire::solver::DirectNetworkSolver::TimestepContext::num_steps, "Total number of steps taken in the simulation.")
        .def_property_readonly("engine", [](const gridfire::solver::DirectNetworkSolver::TimestepContext &ctx) -> const gridfire::DynamicEngine & {
            return ctx.engine;
        }, py::return_value_policy::reference)

        .def_property_readonly("network_species", [](const gridfire::solver::DirectNetworkSolver::TimestepContext &ctx) -> const std::vector<fourdst::atomic::Species> & {
            return ctx.networkSpecies;
        }, py::return_value_policy::reference);
}

