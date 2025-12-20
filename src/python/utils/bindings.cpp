#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly

#include "bindings.h"

#include "gridfire/utils/hashing.h"

namespace py = pybind11;

#include "gridfire/utils/logging.h"

void register_utils_bindings(py::module &m) {
    m.def("formatNuclearTimescaleLogString",
          &gridfire::utils::formatNuclearTimescaleLogString,
          py::arg("ctx"),
          py::arg("engine"),
          py::arg("Y"),
          py::arg("T9"),
          py::arg("rho"),
          "Format a string for logging nuclear timescales based on temperature, density, and energy generation rate."
    );

    m.def(
        "hash_atomic",
        &gridfire::utils::hash_atomic,
        py::arg("a"),
        py::arg("z")
    );

    auto hashing_module = m.def_submodule("hashing", "module for gridfire hashing functions");
    auto reaction_hashing_module = hashing_module.def_submodule("reaction", "utility module for hashing gridfire reaction functions");

    reaction_hashing_module.def(
        "splitmix64",
        &gridfire::utils::hashing::reaction::splitmix64,
        py::arg("x")
    );

    reaction_hashing_module.def(
        "mix_species",
        &gridfire::utils::hashing::reaction::mix_species,
        py::arg("a"),
        py::arg("z")
    );

    reaction_hashing_module.def(
        "multiset_combine",
        &gridfire::utils::hashing::reaction::multiset_combine,
        py::arg("acc"),
        py::arg("x")
    );

    m.def(
        "hash_reaction",
        &gridfire::utils::hash_reaction,
        py::arg("reaction")
    );

}
