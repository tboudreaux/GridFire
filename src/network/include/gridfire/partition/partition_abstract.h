#pragma once

#include <string>

namespace gridfire::partition {
    class PartitionFunction {
    public:
        virtual ~PartitionFunction() = default;
        virtual double evaluate(int z, int a, double T9) const = 0;
        virtual double evaluateDerivative(int z, int a, double T9) const = 0;
        virtual bool supports(int z, int a) const = 0;
        virtual std::string type() const = 0;
    };
}