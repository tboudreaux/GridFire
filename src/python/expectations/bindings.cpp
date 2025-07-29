#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/expectations/expectations.h"

void register_expectation_bindings(py::module &m) {
    py::enum_<gridfire::expectations::EngineErrorTypes>(m, "EngineErrorTypes")
        .value("FAILURE", gridfire::expectations::EngineErrorTypes::FAILURE)
        .value("INDEX", gridfire::expectations::EngineErrorTypes::INDEX)
        .value("STALE", gridfire::expectations::EngineErrorTypes::STALE)
        .export_values();

    py::enum_<gridfire::expectations::StaleEngineErrorTypes>(m, "StaleEngineErrorTypes")
        .value("SYSTEM_RESIZED", gridfire::expectations::StaleEngineErrorTypes::SYSTEM_RESIZED)
        .export_values();

    // Bind the base class
    py::class_<gridfire::expectations::EngineError>(m, "EngineError")
        .def_readonly("message", &gridfire::expectations::EngineError::m_message)
        .def_readonly("type", &gridfire::expectations::EngineError::type)
        .def("__str__", [](const gridfire::expectations::EngineError &e) {return e.m_message;});

    // Bind the EngineIndexError, specifying EngineError as the base
    py::class_<gridfire::expectations::EngineIndexError, gridfire::expectations::EngineError>(m, "EngineIndexError")
        .def(py::init<int>(), py::arg("index"))
        .def_readonly("index", &gridfire::expectations::EngineIndexError::m_index)
        .def("__str__", [](const gridfire::expectations::EngineIndexError &e) {
            return e.m_message + " at index " + std::to_string(e.m_index);
        });

    // Bind the StaleEngineError, specifying EngineError as the base
    py::class_<gridfire::expectations::StaleEngineError, gridfire::expectations::EngineError>(m, "StaleEngineError")
        .def(py::init<gridfire::expectations::StaleEngineErrorTypes>(), py::arg("stale_type"))
        .def_readonly("stale_type", &gridfire::expectations::StaleEngineError::staleType)
        .def("__str__", [](const gridfire::expectations::StaleEngineError &e) {
            return static_cast<std::string>(e);
        });
}
