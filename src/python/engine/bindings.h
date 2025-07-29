#pragma once

#include <pybind11/pybind11.h>

void register_engine_bindings(pybind11::module &m);

void register_base_engine_bindings(pybind11::module &m);

void register_engine_view_bindings(pybind11::module &m);

void abs_stype_register_engine_bindings(pybind11::module &m);
void abs_stype_register_dynamic_engine_bindings(pybind11::module &m);

void con_stype_register_graph_engine_bindings(pybind11::module &m);


