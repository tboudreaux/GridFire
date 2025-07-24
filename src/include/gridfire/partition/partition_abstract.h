#pragma once

#include <string>
#include <memory>

namespace gridfire::partition {

    /**
     * @class PartitionFunction
     * @brief Abstract interface for evaluating nuclear partition functions.
     *
     * Provides methods to compute the partition function and its temperature derivative
     * for a given isotope, to query if the function supports that isotope, and to
     * clone the function object. Concrete implementations must provide temperature-
     * dependent statistical models.
     */
    class PartitionFunction {
    public:
        /**
         * @brief Virtual destructor.
         *
         * Ensures proper cleanup in derived classes.
         */
        virtual ~PartitionFunction() = default;

        /**
         * @brief Evaluate the partition function for a given isotope.
         *
         * @param z Proton number (atomic number) of the isotope; must be >= 1.
         * @param a Mass number of the isotope; must be >= z.
         * @param T9 Temperature in units of 10^9 K; must be > 0.
         * @return Partition function value (dimensionless) at the specified temperature.
         * @pre Derived implementation supports (z, a) and T9 > 0.
         * @post No side effects; pure function.
         */
        [[nodiscard]] virtual double evaluate(int z, int a, double T9) const = 0;

        /**
         * @brief Evaluate the temperature derivative of the partition function.
         *
         * Computes d/dT (partition function) at the given parameters.
         *
         * @param z Proton number (atomic number) of the isotope; must be >= 1.
         * @param a Mass number of the isotope; must be >= z.
         * @param T9 Temperature in units of 10^9 K; must be > 0.
         * @return Temperature derivative of the partition function.
         * @pre Derived implementation supports (z, a) and T9 > 0.
         * @post No side effects; pure function.
         */
        [[nodiscard]] virtual double evaluateDerivative(int z, int a, double T9) const = 0;

        /**
         * @brief Check if this partition function supports an isotope.
         *
         * @param z Proton number of the isotope.
         * @param a Mass number of the isotope.
         * @return true if evaluate and evaluateDerivative can be called for this isotope; false otherwise.
         * @post No side effects.
         */
        [[nodiscard]] virtual bool supports(int z, int a) const = 0;

        /**
         * @brief Get the human-readable type of this partition function.
         *
         * @return String identifier for the partition function implementation.
         * @post No side effects.
         */
        [[nodiscard]] virtual std::string type() const = 0;

        /**
         * @brief Create a deep copy of this PartitionFunction.
         *
         * @return Unique pointer to a new PartitionFunction instance with identical state.
         * @post The caller owns the returned object and must manage its lifetime.
         */
        [[nodiscard]] virtual std::unique_ptr<PartitionFunction> clone() const = 0;
    };
}