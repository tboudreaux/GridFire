#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <iostream>

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/exceptions/exceptions.h"

void register_exception_bindings(py::module &m) {
    py::register_exception<gridfire::exceptions::EngineError>(m, "GridFireEngineError");

    // TODO: Make it so that we can grab the stale state in python
    // m.attr("StaleEngineTrigger") = py::register_exception<gridfire::exceptions::StaleEngineTrigger>(m, "StaleEngineTrigger", m.attr("GridFireEngineError"));
    m.attr("StaleEngineError") = py::register_exception<gridfire::exceptions::StaleEngineError>(m, "StaleEngineError", m.attr("GridFireEngineError"));
    m.attr("FailedToPartitionEngineError") = py::register_exception<gridfire::exceptions::FailedToPartitionEngineError>(m, "FailedToPartitionEngineError", m.attr("GridFireEngineError"));
    m.attr("NetworkResizedError") = py::register_exception<gridfire::exceptions::NetworkResizedError>(m, "NetworkResizedError", m.attr("GridFireEngineError"));
    m.attr("UnableToSetNetworkReactionsError") = py::register_exception<gridfire::exceptions::UnableToSetNetworkReactionsError>(m, "UnableToSetNetworkReactionsError", m.attr("GridFireEngineError"));

    py::class_<gridfire::exceptions::StaleEngineTrigger::state>(m, "StaleEngineState")
        .def(py::init<>())
        .def_readwrite("T9", &gridfire::exceptions::StaleEngineTrigger::state::m_T9)
        .def_readwrite("rho", &gridfire::exceptions::StaleEngineTrigger::state::m_rho)
        .def_readwrite("Y", &gridfire::exceptions::StaleEngineTrigger::state::m_Y)
        .def_readwrite("t", &gridfire::exceptions::StaleEngineTrigger::state::m_t)
        .def_readwrite("total_steps", &gridfire::exceptions::StaleEngineTrigger::state::m_total_steps)
        .def_readwrite("eps_nuc", &gridfire::exceptions::StaleEngineTrigger::state::m_eps_nuc);

    py::class_<gridfire::exceptions::StaleEngineTrigger>(m, "StaleEngineTrigger")
        .def(py::init<const gridfire::exceptions::StaleEngineTrigger::state &>())
        .def("getState", &gridfire::exceptions::StaleEngineTrigger::getState)
        .def("numSpecies", &gridfire::exceptions::StaleEngineTrigger::numSpecies)
        .def("totalSteps", &gridfire::exceptions::StaleEngineTrigger::totalSteps)
        .def("energy", &gridfire::exceptions::StaleEngineTrigger::energy)
        .def("getMolarAbundance", &gridfire::exceptions::StaleEngineTrigger::getMolarAbundance)
        .def("temperature", &gridfire::exceptions::StaleEngineTrigger::temperature)
        .def("density", &gridfire::exceptions::StaleEngineTrigger::density)
        .def("__repr__", [&](const gridfire::exceptions::StaleEngineTrigger& self) {
           return self.what();
        });

}
