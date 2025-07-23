#pragma once

#include <pybind11/pybind11.h>

void register_partition_bindings(pybind11::module &m);

void register_partition_types_bindings(pybind11::module &m);

void register_ground_state_partition_bindings(pybind11::module &m);

void register_rauscher_thielemann_partition_data_record_bindings(pybind11::module &m);

void register_rauscher_thielemann_partition_bindings(pybind11::module &m);

void register_composite_partition_bindings(pybind11::module &m);

