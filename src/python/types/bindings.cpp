#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/network.h"

void register_type_bindings(pybind11::module &m) {
    py::class_<gridfire::NetIn>(m, "NetIn")
    .def(py::init<>())
    .def_readwrite("composition", &gridfire::NetIn::composition)
    .def_readwrite("tMax", &gridfire::NetIn::tMax)
    .def_readwrite("dt0", &gridfire::NetIn::dt0)
    .def_readwrite("temperature", &gridfire::NetIn::temperature)
    .def_readwrite("density", &gridfire::NetIn::density)
    .def_readwrite("energy", &gridfire::NetIn::energy)
    .def("__repr__", [](const gridfire::NetIn &netIn) {
        std::stringstream ss;
        ss << "NetIn(composition=" << netIn.composition
           << ", tMax=" << netIn.tMax
           << ", dt0=" << netIn.dt0
           << ", temperature=" << netIn.temperature
           << ", density=" << netIn.density
           << ", energy=" << netIn.energy << ")";
        return ss.str();
    });

    py::class_<gridfire::NetOut>(m, "NetOut")
        .def_readonly("composition", &gridfire::NetOut::composition)
        .def_readonly("num_steps", &gridfire::NetOut::num_steps)
        .def_readonly("energy", &gridfire::NetOut::energy)
        .def("__repr__", [](const gridfire::NetOut &netOut) {
            std::stringstream ss;
            ss << "NetOut(composition=" << netOut.composition
               << ", num_steps=" << netOut.num_steps
               << ", energy=" << netOut.energy << ")";
            return ss.str();
        });

}
