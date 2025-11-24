#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gridfire/partition/partition.h"


class PyPartitionFunction final : public gridfire::partition::PartitionFunction {
    [[nodiscard]] double evaluate(int z, int a, double T9) const override;
    [[nodiscard]] double evaluateDerivative(int z, int a, double T9) const override;
    [[nodiscard]] bool supports(int z, int a) const override;
    [[nodiscard]] std::string type() const override;
    [[nodiscard]] std::unique_ptr<PartitionFunction> clone() const override;
};