#include "py_partition.h"

#include "gridfire/partition/partition.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <string>
#include <memory>


namespace py = pybind11;


double PyPartitionFunction::evaluate(int z, int a, double T9) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::partition::PartitionFunction,
        evaluate,
        z, a, T9
    );
}

double PyPartitionFunction::evaluateDerivative(int z, int a, double T9) const {
    PYBIND11_OVERRIDE_PURE(
        double,
        gridfire::partition::PartitionFunction,
        evaluateDerivative,
        z, a, T9
    );
}

bool PyPartitionFunction::supports(int z, int a) const {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::partition::PartitionFunction,
        supports,
        z, a
    );
}

std::string PyPartitionFunction::type() const {
    PYBIND11_OVERRIDE_PURE(
        std::string,
        gridfire::partition::PartitionFunction,
        type
    );
}

std::unique_ptr<gridfire::partition::PartitionFunction> PyPartitionFunction::clone() const {
    PYBIND11_OVERRIDE_PURE(
        std::unique_ptr<gridfire::partition::PartitionFunction>,
        gridfire::partition::PartitionFunction,
        clone
    );
}

