#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <functional>

#include "bindings.h"

#include "gridfire/solver/strategies/CVODE_solver_strategy.h"
#include "trampoline/py_solver.h"

namespace py = pybind11;


void register_solver_bindings(const py::module &m) {
    auto py_solver_context_base = py::class_<gridfire::solver::SolverContextBase>(m, "SolverContextBase");

    auto py_cvode_timestep_context = py::class_<gridfire::solver::CVODESolverStrategy::TimestepContext, gridfire::solver::SolverContextBase>(m, "CVODETimestepContext");
    py_cvode_timestep_context.def_readonly("t", &gridfire::solver::CVODESolverStrategy::TimestepContext::t);
    py_cvode_timestep_context.def_property_readonly(
        "state",
        [](const gridfire::solver::CVODESolverStrategy::TimestepContext& self) -> std::vector<double> {
            const sunrealtype* nvec_data = N_VGetArrayPointer(self.state);
            const sunindextype length = N_VGetLength(self.state);
            return std::vector<double>(nvec_data, nvec_data + length);
        }
    );
    py_cvode_timestep_context.def_readonly("dt", &gridfire::solver::CVODESolverStrategy::TimestepContext::dt);
    py_cvode_timestep_context.def_readonly("last_step_time", &gridfire::solver::CVODESolverStrategy::TimestepContext::last_step_time);
    py_cvode_timestep_context.def_readonly("T9", &gridfire::solver::CVODESolverStrategy::TimestepContext::T9);
    py_cvode_timestep_context.def_readonly("rho", &gridfire::solver::CVODESolverStrategy::TimestepContext::rho);
    py_cvode_timestep_context.def_readonly("num_steps", &gridfire::solver::CVODESolverStrategy::TimestepContext::num_steps);
    py_cvode_timestep_context.def_readonly("currentConvergenceFailures", &gridfire::solver::CVODESolverStrategy::TimestepContext::currentConvergenceFailures);
    py_cvode_timestep_context.def_readonly("currentNonlinearIterations", &gridfire::solver::CVODESolverStrategy::TimestepContext::currentNonlinearIterations);
    py_cvode_timestep_context.def_property_readonly(
        "engine",
        [](const gridfire::solver::CVODESolverStrategy::TimestepContext& self) -> const gridfire::engine::DynamicEngine& {
            return self.engine;
        }
    );
    py_cvode_timestep_context.def_property_readonly(
        "networkSpecies",
        [](const gridfire::solver::CVODESolverStrategy::TimestepContext& self) -> std::vector<fourdst::atomic::Species> {
            return self.networkSpecies;
        }
    );

    auto py_dynamic_network_solver_strategy = py::class_<gridfire::solver::DynamicNetworkSolverStrategy, PyDynamicNetworkSolverStrategy>(m, "DynamicNetworkSolverStrategy");
    py_dynamic_network_solver_strategy.def(
        "evaluate",
        &gridfire::solver::DynamicNetworkSolverStrategy::evaluate,
        py::arg("netIn"),
        "evaluate the dynamic engine using the dynamic engine class"
    );


    py_dynamic_network_solver_strategy.def(
        "describe_callback_context",
        &gridfire::solver::DynamicNetworkSolverStrategy::describe_callback_context,
        "Get a structure representing what data is in the callback context in a human readable format"
    );

    auto py_cvode_solver_strategy = py::class_<gridfire::solver::CVODESolverStrategy, gridfire::solver::DynamicNetworkSolverStrategy>(m, "CVODESolverStrategy");

    py_cvode_solver_strategy.def(
        py::init<gridfire::engine::DynamicEngine&>(),
        py::arg("engine"),
        "Initialize the CVODESolverStrategy object."
    );

    py_cvode_solver_strategy.def(
        "evaluate",
        py::overload_cast<const gridfire::NetIn&, bool>(&gridfire::solver::CVODESolverStrategy::evaluate),
        py::arg("netIn"),
        py::arg("display_trigger") = false,
        "evaluate the dynamic engine using the dynamic engine class"
    );

    py_cvode_solver_strategy.def(
        "get_stdout_logging_enabled",
        &gridfire::solver::CVODESolverStrategy::get_stdout_logging_enabled,
        "Check if solver logging to standard output is enabled."
    );

    py_cvode_solver_strategy.def(
        "set_stdout_logging_enabled",
        &gridfire::solver::CVODESolverStrategy::set_stdout_logging_enabled,
        py::arg("logging_enabled"),
        "Enable logging to standard output."
    );

    py_cvode_solver_strategy.def(
        "set_absTol",
        &gridfire::solver::CVODESolverStrategy::set_absTol,
        py::arg("absTol"),
        "Set the absolute tolerance for the CVODE solver."
    );

    py_cvode_solver_strategy.def(
        "set_relTol",
        &gridfire::solver::CVODESolverStrategy::set_relTol,
        py::arg("relTol"),
        "Set the relative tolerance for the CVODE solver."
    );

    py_cvode_solver_strategy.def(
        "get_absTol",
        &gridfire::solver::CVODESolverStrategy::get_absTol,
        "Get the absolute tolerance for the CVODE solver."
    );

    py_cvode_solver_strategy.def(
        "get_relTol",
        &gridfire::solver::CVODESolverStrategy::get_relTol,
        "Get the relative tolerance for the CVODE solver."
    );

    py_cvode_solver_strategy.def(
    "set_callback",
    [](
        gridfire::solver::CVODESolverStrategy& self,
        std::function<void(const gridfire::solver::CVODESolverStrategy::TimestepContext&)> cb
    ) {
        self.set_callback(std::any(cb));
    },
    py::arg("cb"),
    "Set a callback function which will run at the end of every successful timestep"
    );
}

