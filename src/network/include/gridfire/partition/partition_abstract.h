#pragma once

#include <string>
#include <memory>

namespace gridfire::partition {
    class PartitionFunction {
    public:
        virtual ~PartitionFunction() = default;
        [[nodiscard]] virtual double evaluate(int z, int a, double T9) const = 0;
        [[nodiscard]] virtual double evaluateDerivative(int z, int a, double T9) const = 0;
        [[nodiscard]] virtual bool supports(int z, int a) const = 0;
        [[nodiscard]] virtual std::string type() const = 0;
        [[nodiscard]] virtual std::unique_ptr<PartitionFunction> clone() const = 0;
    };
}