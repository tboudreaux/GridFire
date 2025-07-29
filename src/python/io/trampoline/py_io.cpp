#include "gridfire/io/io.h"
#include "py_io.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

gridfire::io::ParsedNetworkData PyNetworkFileParser::parse(const std::string &filename) const {
    PYBIND11_OVERLOAD_PURE(
        gridfire::io::ParsedNetworkData,
        gridfire::io::NetworkFileParser,
        parse  // Method name
    );
}