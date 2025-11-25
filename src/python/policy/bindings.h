#pragma once

#include <pybind11/pybind11.h>

void register_policy_bindings(pybind11::module& m);
void register_reaction_chain_policy_bindings(pybind11::module& m);
void register_network_policy_bindings(pybind11::module& m);