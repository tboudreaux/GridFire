#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/utils/logging.h"

void register_utils_bindings(py::module &m) {
    m.def("formatNuclearTimescaleLogString",
          &gridfire::utils::formatNuclearTimescaleLogString,
          py::arg("engine"),
          py::arg("Y"),
          py::arg("T9"),
          py::arg("rho"),
          "Format a string for logging nuclear timescales based on temperature, density, and energy generation rate."
    );
}
