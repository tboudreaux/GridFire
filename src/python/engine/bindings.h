#pragma once

#include <pybind11/pybind11.h>

void register_engine_bindings(pybind11::module &m);

void register_base_engine_bindings(const pybind11::module &m);

void register_engine_view_bindings(const pybind11::module &m);

void abs_stype_register_engine_bindings(const pybind11::module &m);
void abs_stype_register_dynamic_engine_bindings(const pybind11::module &m);

void con_stype_register_graph_engine_bindings(const pybind11::module &m);

void register_engine_diagnostic_bindings(pybind11::module &m);
void register_engine_procedural_bindings(pybind11::module &m);

void register_engine_construction_bindings(pybind11::module &m);
void register_engine_priming_bindings(pybind11::module &m);

void register_engine_type_bindings(pybind11::module &m);
void register_engine_building_type_bindings(pybind11::module &m);
void register_engine_reporting_type_bindings(pybind11::module &m);


