#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <string_view>
#include <vector>

#include "bindings.h"

#include "gridfire/screening/screening.h"
#include "trampoline/py_screening.h"

namespace py = pybind11;

void register_screening_bindings(py::module &m) {
    py::class_<gridfire::screening::ScreeningModel, PyScreening>(m, "ScreeningModel");

    py::enum_<gridfire::screening::ScreeningType>(m, "ScreeningType")
        .value("BARE", gridfire::screening::ScreeningType::BARE)
        .value("WEAK", gridfire::screening::ScreeningType::WEAK)
        .export_values();

    m.def("selectScreeningModel", &gridfire::screening::selectScreeningModel,
          py::arg("type"),
          "Select a screening model based on the specified type. Returns a pointer to the selected model.");

    py::class_<gridfire::screening::BareScreeningModel>(m, "BareScreeningModel")
        .def(py::init<>())
        .def("calculateScreeningFactors",
             py::overload_cast<const gridfire::reaction::LogicalReactionSet&, const std::vector<fourdst::atomic::Species>&, const std::vector<double>&, double, double>(&gridfire::screening::BareScreeningModel::calculateScreeningFactors, py::const_),
             py::arg("reactions"), py::arg("species"), py::arg("Y"), py::arg("T9"), py::arg("rho"),
             "Calculate the bare plasma screening factors. This always returns 1.0 (bare)"
            );

    py::class_<gridfire::screening::WeakScreeningModel>(m, "WeakScreeningModel")
        .def(py::init<>())
        .def("calculateScreeningFactors",
             py::overload_cast<const gridfire::reaction::LogicalReactionSet&, const std::vector<fourdst::atomic::Species>&, const std::vector<double>&, double, double>(&gridfire::screening::WeakScreeningModel::calculateScreeningFactors, py::const_),
             py::arg("reactions"), py::arg("species"), py::arg("Y"), py::arg("T9"), py::arg("rho"),
             "Calculate the weak plasma screening factors using the Salpeter (1954) model."
            );
}