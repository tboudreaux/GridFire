#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

#include "types/bindings.h"
#include "partition/bindings.h"
#include "expectations/bindings.h"
#include "engine/bindings.h"
#include "exceptions/bindings.h"
#include "io/bindings.h"
#include "reaction/bindings.h"
#include "screening/bindings.h"
#include "solver/bindings.h"
#include "utils/bindings.h"

PYBIND11_MODULE(gridfire, m) {
    m.doc() = "Python bindings for the fourdst utility modules which are a part of the 4D-STAR project.";

    pybind11::module::import("fourdst.constants");
    pybind11::module::import("fourdst.composition");
    pybind11::module::import("fourdst.config");
    pybind11::module::import("fourdst.atomic");

    auto typeMod  = m.def_submodule("type", "GridFire type bindings");
    register_type_bindings(typeMod);

    auto partitionMod  = m.def_submodule("partition", "GridFire partition function bindings");
    register_partition_bindings(partitionMod);

    auto expectationMod  = m.def_submodule("expectations", "GridFire expectations bindings");
    register_expectation_bindings(expectationMod);

    auto reactionMod  = m.def_submodule("reaction", "GridFire reaction bindings");
    register_reaction_bindings(reactionMod);

    auto screeningMod  = m.def_submodule("screening", "GridFire plasma screening bindings");
    register_screening_bindings(screeningMod);

    auto ioMod  = m.def_submodule("io", "GridFire io bindings");
    register_io_bindings(ioMod);

    auto exceptionMod  = m.def_submodule("exceptions", "GridFire exceptions bindings");
    register_exception_bindings(exceptionMod);

    auto engineMod  = m.def_submodule("engine", "Engine and Engine View bindings");
    register_engine_bindings(engineMod);

    auto solverMod  = m.def_submodule("solver", "GridFire numerical solver bindings");
    register_solver_bindings(solverMod);

    auto utilsMod  = m.def_submodule("utils", "GridFire utility method bindings");
    register_utils_bindings(utilsMod);
}