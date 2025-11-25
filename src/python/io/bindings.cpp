#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include <vector>

#include "bindings.h"

#include "gridfire/io/io.h"
#include "trampoline/py_io.h"

namespace py = pybind11;

auto register_io_bindings(const py::module &m) -> void {
    auto ParsedNetworkData = py::class_<gridfire::io::ParsedNetworkData>(m, "ParsedNetworkData");

    auto NetworkFileParser = py::class_<gridfire::io::NetworkFileParser, PyNetworkFileParser>(m, "NetworkFileParser");

    auto SimpleReactionListFileParser = py::class_<gridfire::io::SimpleReactionListFileParser, gridfire::io::NetworkFileParser>(m, "SimpleReactionListFileParser")
        .def("parse", &gridfire::io::SimpleReactionListFileParser::parse,
            py::arg("filename"),
            "Parse a simple reaction list file and return a ParsedNetworkData object.");
}