#pragma once

#include <pybind11/pybind11.h>

void register_scratchpad_types_bindings(pybind11::module_& m);
void register_scratchpad_bindings(pybind11::module_& m);
void register_state_blob_bindings(pybind11::module_& m);