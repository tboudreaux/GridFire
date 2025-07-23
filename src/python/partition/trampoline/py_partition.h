#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gridfire/partition/partition.h"


class PyPartitionFunction final : public gridfire::partition::PartitionFunction {
    double evaluate(int z, int a, double T9) const override;
    double evaluateDerivative(int z, int a, double T9) const override;
    bool supports(int z, int a) const override;
    std::string type() const override;
    std::unique_ptr<gridfire::partition::PartitionFunction> clone() const override;
};