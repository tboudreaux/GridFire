#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc if needed directly

#include <string_view>
#include <vector>

#include "bindings.h"

#include "gridfire/io/io.h"
#include "trampoline/py_io.h"

namespace py = pybind11;

void register_io_bindings(py::module &m) {
    py::class_<gridfire::io::ParsedNetworkData>(m, "ParsedNetworkData");

    py::class_<gridfire::io::NetworkFileParser, PyNetworkFileParser>(m, "NetworkFileParser");

    py::class_<gridfire::io::SimpleReactionListFileParser, gridfire::io::NetworkFileParser>(m, "SimpleReactionListFileParser")
        .def("parse", &gridfire::io::SimpleReactionListFileParser::parse,
            py::arg("filename"),
            "Parse a simple reaction list file and return a ParsedNetworkData object.");

    // py::class_<gridfire::io::MESANetworkFileParser, gridfire::io::NetworkFileParser>(m, "MESANetworkFileParser")
    //     .def("parse", &gridfire::io::MESANetworkFileParser::parse,
    //         py::arg("filename"),
    //         "Parse a MESA network file and return a ParsedNetworkData object.");
}