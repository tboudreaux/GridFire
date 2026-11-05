#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include <vector>

#include "bindings.h"

#include "gridfire/screening/screening.h"
#include "trampoline/py_screening.h"

namespace py = pybind11;

void register_screening_bindings(py::module &m) {
    auto screening_model = py::class_<gridfire::screening::ScreeningModel, PyScreening>(m, "ScreeningModel");

    py::enum_<gridfire::screening::ScreeningType>(m, "ScreeningType")
        .value("BARE", gridfire::screening::ScreeningType::BARE)
        .value("WEAK", gridfire::screening::ScreeningType::WEAK)
        .value("CD09", gridfire::screening::ScreeningType::CD09)
        .value("BLENDED", gridfire::screening::ScreeningType::BLENDED)
        .export_values();

    py::class_<gridfire::screening::ScreeningBlendWindow>(m, "ScreeningBlendWindow")
        .def(py::init<>())
        .def_readwrite("weakExponent", &gridfire::screening::ScreeningBlendWindow::weakExponent,
            "Below this weak-screening exponent H the weak prescription is used alone")
        .def_readwrite("strongExponent", &gridfire::screening::ScreeningBlendWindow::strongExponent,
            "Above this weak-screening exponent H the strong prescription is used alone (at most 2)");

    m.def(
        "selectScreeningModel",
        &gridfire::screening::selectScreeningModel,
        py::arg("type"),
        py::arg("window") = gridfire::screening::ScreeningBlendWindow{},
        "Select a screening model based on the specified type. Returns a pointer to the selected model."
    );

    py::class_<gridfire::screening::BareScreeningModel>(m, "BareScreeningModel")
        .def(py::init<>())
        .def("calculateScreeningFactors",
             py::overload_cast<
                const gridfire::reaction::ReactionSet&,
                const std::vector<fourdst::atomic::Species>&,
                const std::vector<double>&,
                double,
                double
            >(&gridfire::screening::BareScreeningModel::calculateScreeningFactors, py::const_),
            py::arg("reactions"),
            py::arg("species"),
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate the bare plasma screening factors. This always returns 1.0 (bare)"
        );

    py::class_<gridfire::screening::WeakScreeningModel>(m, "WeakScreeningModel")
        .def(py::init<>())
        .def("calculateScreeningFactors",
             py::overload_cast<
                const gridfire::reaction::ReactionSet&,
                const std::vector<fourdst::atomic::Species>&,
                const std::vector<double>&,
                double,
                double
            >(&gridfire::screening::WeakScreeningModel::calculateScreeningFactors, py::const_),
            py::arg("reactions"),
            py::arg("species"),
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate the weak plasma screening factors using the Salpeter (1954) model."
        );

    py::class_<gridfire::screening::CD09ScreeningModel>(m, "CD09ScreeningModel")
        .def(py::init<>())
        .def("calculateScreeningFactors",
             py::overload_cast<
                const gridfire::reaction::ReactionSet&,
                const std::vector<fourdst::atomic::Species>&,
                const std::vector<double>&,
                double,
                double
            >(&gridfire::screening::CD09ScreeningModel::calculateScreeningFactors, py::const_),
            py::arg("reactions"),
            py::arg("species"),
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate plasma screening factors with Chugunov & DeWitt (2009), valid from weak to strong coupling."
        );

    py::class_<gridfire::screening::BlendedScreeningModel>(m, "BlendedScreeningModel")
        .def(py::init([](const gridfire::screening::ScreeningBlendWindow& window) {
            return gridfire::screening::BlendedScreeningModel(
                std::make_unique<gridfire::screening::WeakScreeningModel>(),
                std::make_unique<gridfire::screening::CD09ScreeningModel>(),
                window);
        }), py::arg("window") = gridfire::screening::ScreeningBlendWindow{},
            "Salpeter weak screening handing over to Chugunov & DeWitt (2009) across the blend window")
        .def("weight", &gridfire::screening::BlendedScreeningModel::weight, py::arg("H"),
            "Blend weight for a weak-screening exponent H (0: weak only, 1: CD09 only)")
        .def("calculateScreeningFactors",
             py::overload_cast<
                const gridfire::reaction::ReactionSet&,
                const std::vector<fourdst::atomic::Species>&,
                const std::vector<double>&,
                double,
                double
            >(&gridfire::screening::BlendedScreeningModel::calculateScreeningFactors, py::const_),
            py::arg("reactions"),
            py::arg("species"),
            py::arg("Y"),
            py::arg("T9"),
            py::arg("rho"),
            "Calculate plasma screening factors, blending Salpeter and Chugunov & DeWitt (2009)."
        );
}
