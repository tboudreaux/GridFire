#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <functional>

#include "bindings.h"

#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/engine/scratchpads/blob.h"
#include "trampoline/py_solver.h"

namespace py = pybind11;


void register_solver_bindings(const py::module &m) {
    auto py_cvode_timestep_context = py::class_<gridfire::solver::PointSolverTimestepContext>(m, "PointSolverTimestepContext");
    py_cvode_timestep_context.def_readonly("t", &gridfire::solver::PointSolverTimestepContext::t);
    py_cvode_timestep_context.def_readonly("dt", &gridfire::solver::PointSolverTimestepContext::dt);
    py_cvode_timestep_context.def_readonly("last_step_time", &gridfire::solver::PointSolverTimestepContext::last_step_time);
    py_cvode_timestep_context.def_readonly("T9", &gridfire::solver::PointSolverTimestepContext::T9);
    py_cvode_timestep_context.def_readonly("rho", &gridfire::solver::PointSolverTimestepContext::rho);
    py_cvode_timestep_context.def_readonly("num_steps", &gridfire::solver::PointSolverTimestepContext::num_steps);
    py_cvode_timestep_context.def_readonly("currentConvergenceFailures", &gridfire::solver::PointSolverTimestepContext::currentConvergenceFailures);
    py_cvode_timestep_context.def_readonly("currentNonlinearIterations", &gridfire::solver::PointSolverTimestepContext::currentNonlinearIterations);
    py_cvode_timestep_context.def_property_readonly(
        "engine",
        [](const gridfire::solver::PointSolverTimestepContext& self) -> const gridfire::engine::DynamicEngine& {
            return self.engine;
        }
    );
    py_cvode_timestep_context.def_property_readonly(
        "networkSpecies",
        [](const gridfire::solver::PointSolverTimestepContext& self) -> std::vector<fourdst::atomic::Species> {
            return self.networkSpecies;
        }
    );
    py_cvode_timestep_context.def_property_readonly(
        "state_ctx",
        [](const gridfire::solver::PointSolverTimestepContext& self) {
            return &(self.state_ctx);
        },
        py::return_value_policy::reference_internal
    );
    py_cvode_timestep_context.def_property_readonly(
        "composition",
        [](const gridfire::solver::PointSolverTimestepContext& self) -> fourdst::composition::Composition {
            return self.getPhysicalComposition();
        }
    );
    py_cvode_timestep_context.def_property_readonly(
        "rawState",
        [](const gridfire::solver::PointSolverTimestepContext& self) -> std::vector<double> {
            const std::span<const double> s = self.rawState();
            return std::vector<double>(s.begin(), s.end());
        }
    );
    py_cvode_timestep_context.def("abundance",
        py::overload_cast<size_t>(&gridfire::solver::PointSolverTimestepContext::abundance, py::const_),
        py::arg("species_index"));
    py_cvode_timestep_context.def("abundance",
        py::overload_cast<const fourdst::atomic::Species&>(&gridfire::solver::PointSolverTimestepContext::abundance, py::const_),
        py::arg("species"));
    py_cvode_timestep_context.def_property_readonly("accumulatedSpecificEnergy", &gridfire::solver::PointSolverTimestepContext::accumulatedSpecificEnergy);



    auto py_solver_context_base = py::class_<gridfire::solver::SolverContextBase>(m, "SolverContextBase");
    auto py_point_solver_context = py::class_<gridfire::solver::PointSolverContext, gridfire::solver::SolverContextBase>(m, "PointSolverContext");

    py_point_solver_context
        .def_readonly(
            "sun_ctx", &gridfire::solver::PointSolverContext::sun_ctx
        )
        .def_readonly(
            "cvode_mem", &gridfire::solver::PointSolverContext::cvode_mem
        )
        .def_readonly(
            "Y", &gridfire::solver::PointSolverContext::Y
        )
        .def_readonly(
            "YErr",  &gridfire::solver::PointSolverContext::YErr
        )
        .def_readonly(
            "J", &gridfire::solver::PointSolverContext::J
        )
        .def_readonly(
            "LS", &gridfire::solver::PointSolverContext::LS
        )
        .def_property_readonly(
            "engine_ctx",
            [](const gridfire::solver::PointSolverContext& self) -> gridfire::engine::scratch::StateBlob& {
                return *(self.engine_ctx);
            },
            py::return_value_policy::reference
        )
        .def_readonly(
            "num_steps", &gridfire::solver::PointSolverContext::num_steps
        )
        .def_property(
            "abs_tol",
            [](const gridfire::solver::PointSolverContext& self) -> double {
                return self.abs_tol.value();
            },
            [](gridfire::solver::PointSolverContext& self, double abs_tol) -> void {
                self.abs_tol = abs_tol;
            }
        )
        .def_property(
            "rel_tol",
            [](const gridfire::solver::PointSolverContext& self) -> double {
                return self.rel_tol.value();
            },
            [](gridfire::solver::PointSolverContext& self, double rel_tol) -> void {
                self.rel_tol = rel_tol;
            }
        )
        .def_property(
            "stdout_logging",
            [](const gridfire::solver::PointSolverContext& self) -> bool {
                return self.stdout_logging;
            },
            [](gridfire::solver::PointSolverContext& self, const bool enable) -> void {
                self.stdout_logging = enable;
            }
        )
        .def_property(
            "detailed_logging",
            [](const gridfire::solver::PointSolverContext& self) -> bool {
                return self.detailed_step_logging;
            },
            [](gridfire::solver::PointSolverContext& self, const bool enable) -> void {
                self.detailed_step_logging = enable;
            }
        )
        .def_property(
            "callback",
            [](const gridfire::solver::PointSolverContext& self) -> std::optional<std::function<void(const gridfire::solver::PointSolverTimestepContext&)>> {
                return self.callback;
            },
            [](gridfire::solver::PointSolverContext& self, const std::optional<std::function<void(const gridfire::solver::PointSolverTimestepContext&)>>& cb) {
                self.callback = cb;
            }
        )
        .def("reset_all", &gridfire::solver::PointSolverContext::reset_all)
        .def("reset_user", &gridfire::solver::PointSolverContext::reset_user)
        .def("reset_cvode", &gridfire::solver::PointSolverContext::reset_cvode)
        .def("clear_context", &gridfire::solver::PointSolverContext::clear_context)
        .def("init_context", &gridfire::solver::PointSolverContext::init_context)
        .def("has_context", &gridfire::solver::PointSolverContext::has_context)
        .def("init", &gridfire::solver::PointSolverContext::init)
        .def(py::init<const gridfire::engine::scratch::StateBlob&>(), py::arg("engine_ctx"));



    auto py_single_zone_dynamic_network_solver = py::class_<gridfire::solver::SingleZoneDynamicNetworkSolver, PySingleZoneDynamicNetworkSolver>(m, "SingleZoneDynamicNetworkSolver");
    py_single_zone_dynamic_network_solver.def(
        "evaluate",
        &gridfire::solver::SingleZoneDynamicNetworkSolver::evaluate,
        py::arg("solver_ctx"),
        py::arg("netIn"),
        "evaluate the dynamic engine using the dynamic engine class for a single zone"
    );
    auto py_multi_zone_dynamic_network_solver = py::class_<gridfire::solver::MultiZoneDynamicNetworkSolver, PyMultiZoneDynamicNetworkSolver>(m, "MultiZoneDynamicNetworkSolver");
    py_multi_zone_dynamic_network_solver.def(
        "evaluate",
        &gridfire::solver::MultiZoneDynamicNetworkSolver::evaluate,
        py::arg("solver_ctx"),
        py::arg("netIns"),
        "evaluate the dynamic engine using the dynamic engine class for multiple zones (using openmp if available)"
    );

    auto py_point_solver = py::class_<gridfire::solver::PointSolver, gridfire::solver::SingleZoneDynamicNetworkSolver>(m, "PointSolver");

    py_point_solver.def(
        py::init<gridfire::engine::DynamicEngine&>(),
        py::arg("engine"),
        "Initialize the PointSolver object."
    );

    py_point_solver.def(
        py::init<gridfire::engine::DynamicEngine&, gridfire::config::GridFireConfig&>(),
        py::arg("engine"),
        py::arg("config"),
        "Initialize the PointSolver object with a configuration set."
    );

    py_point_solver.def(
        "getConfig",
        &gridfire::solver::PointSolver::getConfig,
        "Get a copy of the config object"
    );


    py_point_solver.def(
        "evaluate",
        py::overload_cast<gridfire::solver::SolverContextBase&, const gridfire::NetIn&, bool, bool>(&gridfire::solver::PointSolver::evaluate, py::const_),
        py::arg("solver_ctx"),
        py::arg("netIn"),
        py::arg("display_trigger") = false,
        py::arg("force_reinitialization") = false,
        "evaluate the dynamic engine using the dynamic engine class"
    );

    auto py_grid_solver_context = py::class_<gridfire::solver::GridSolverContext, gridfire::solver::SolverContextBase>(m, "GridSolverContext");
    py_grid_solver_context.def(py::init<const gridfire::engine::scratch::StateBlob&>(), py::arg("ctx_template"));
    py_grid_solver_context.def("init", &gridfire::solver::GridSolverContext::init);
    py_grid_solver_context.def("reset", &gridfire::solver::GridSolverContext::reset);
    py_grid_solver_context.def("set_callback", py::overload_cast<const std::function<void(const gridfire::solver::TimestepContextBase&)>&>(&gridfire::solver::GridSolverContext::set_callback) , py::arg("callback"));
    py_grid_solver_context.def("set_callback", py::overload_cast<const std::function<void(const gridfire::solver::TimestepContextBase&)>&, size_t>(&gridfire::solver::GridSolverContext::set_callback) , py::arg("callback"), py::arg("zone_idx"));
    py_grid_solver_context.def("clear_callback", py::overload_cast<>(&gridfire::solver::GridSolverContext::clear_callback));
    py_grid_solver_context.def("clear_callback", py::overload_cast<size_t>(&gridfire::solver::GridSolverContext::clear_callback), py::arg("zone_idx"));
    py_grid_solver_context.def_property(
        "stdout_logging",
        [](const gridfire::solver::GridSolverContext& self) -> bool {
            return self.zone_stdout_logging;
        },
        [](gridfire::solver::GridSolverContext& self, const bool enable) -> void {
            self.zone_stdout_logging = enable;
        }
    )
    .def_property(
        "detailed_logging",
        [](const gridfire::solver::GridSolverContext& self) -> bool {
            return self.zone_detailed_logging;
        },
        [](gridfire::solver::GridSolverContext& self, const bool enable) -> void {
            self.zone_detailed_logging = enable;
        }
    )
    .def_property(
        "zone_completion_logging",
        [](const gridfire::solver::GridSolverContext& self) -> bool {
            return self.zone_completion_logging;
        },
        [](gridfire::solver::GridSolverContext& self, const bool enable) -> void {
            self.zone_completion_logging = enable;
        }
    );

    auto py_grid_solver = py::class_<gridfire::solver::GridSolver, gridfire::solver::MultiZoneDynamicNetworkSolver>(m, "GridSolver");
    py_grid_solver.def(
        py::init<const gridfire::engine::DynamicEngine&, const gridfire::solver::SingleZoneDynamicNetworkSolver&>(),
        py::arg("engine"),
        py::arg("solver"),
        "Initialize the GridSolver object."
    );

    py_grid_solver.def(
        "evaluate",
        &gridfire::solver::GridSolver::evaluate,
        py::arg("solver_ctx"),
        py::arg("netIns"),
        "evaluate the dynamic engine using the dynamic engine class"
    );

}

