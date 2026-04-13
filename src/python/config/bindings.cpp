#include "bindings.h"

#include "gridfire/config/config.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void register_config_bindings(pybind11::module &m) {

    py::class_<gridfire::config::BoundaryFluxConfig>(m, "BoundaryFluxConfig")
        .def(py::init<>())
        .def_readwrite("relativeThreshold", &gridfire::config::BoundaryFluxConfig::relativeThreshold)
        .def_readwrite("absoluteThreshold", &gridfire::config::BoundaryFluxConfig::absoluteThreshold);

    py::class_<gridfire::config::TriggerConfig>(m, "TriggerConfig")
        .def(py::init<>())
        .def_readwrite("offDiagonalThreshold", &gridfire::config::TriggerConfig::offDiagonalThreshold)
        .def_readwrite("timestepCollapseRatio", &gridfire::config::TriggerConfig::timestepCollapseRatio)
        .def_readwrite("maxConvergenceFailures", &gridfire::config::TriggerConfig::maxConvergenceFailures)
        .def_readwrite("boundaryFlux", &gridfire::config::TriggerConfig::boundaryFlux);

    py::class_<gridfire::config::PointSolverConfig>(m, "PointSolverConfig")
        .def(py::init<>())
        .def_readwrite("absTol", &gridfire::config::PointSolverConfig::absTol)
        .def_readwrite("relTol", &gridfire::config::PointSolverConfig::relTol)
        .def_readwrite("trigger", &gridfire::config::PointSolverConfig::trigger);

    py::class_<gridfire::config::SolverConfig>(m, "SolverConfig")
        .def(py::init<>())
        .def_readwrite("pointSolver", &gridfire::config::SolverConfig::pointSolver);

    py::class_<gridfire::config::AdaptiveEngineViewConfig>(m, "AdaptiveEngineViewConfig")
        .def(py::init<>())
        .def_readwrite("relativeCullingThreshold", &gridfire::config::AdaptiveEngineViewConfig::relativeCullingThreshold);

    py::class_<gridfire::config::EngineViewConfig>(m, "EngineViewConfig")
        .def(py::init<>())
        .def_readwrite("adaptiveEngineView", &gridfire::config::EngineViewConfig::adaptiveEngineView);

    py::class_<gridfire::config::EngineConfig>(m, "EngineConfig")
        .def(py::init<>())
        .def_readwrite("views", &gridfire::config::EngineConfig::views);

    py::class_<gridfire::config::GridFireConfig>(m, "GridFireConfig")
        .def(py::init<>())
        .def_readwrite("solver", &gridfire::config::GridFireConfig::solver)
        .def_readwrite("engine", &gridfire::config::GridFireConfig::engine);
}