#pragma once

#include "gridfire/partition/partition_abstract.h"
#include "gridfire/partition/partition_types.h"

#include "fourdst/logging/logging.h"

#include <string>
#include <unordered_map>
#include <vector>

#include <memory>


namespace gridfire::partition {
    /**
     * @class CompositePartitionFunction
     * @brief Combines multiple PartitionFunction instances into a single composite strategy.
     *
     * Maintains an ordered list of sub-functions and delegates evaluation and derivative calls
     * to the first function that supports the requested isotope.
     *
     * See partition_composite.cpp for details on sub-function selection and error logging.
     *
     * @throws std::runtime_error If no sub-function supports a given (z,a,T9) in evaluate or evaluateDerivative.
     */
     class CompositePartitionFunction final : public PartitionFunction {
     public:
        /**
         * @brief Construct a composite function from specified types.
         *
         * Instantiates sub-functions according to the order of types provided.
         * @param partitionFunctions List of BasePartitionType identifiers for sub-functions.
         * @pre partitionFunctions must not be empty.
         * @post m_partitionFunctions contains instances matching each type.
+         */
         explicit CompositePartitionFunction(const std::vector<BasePartitionType>& partitionFunctions);
        /**
         * @brief Copy constructor.
         *
         * Creates deep clones of the sub-functions in another composite.
         * @param other Existing composite to copy from.
         * @post m_partitionFunctions contains clones of other's sub-functions.
+         */
         CompositePartitionFunction(const CompositePartitionFunction& other);
        /**
         * @brief Evaluate the composite partition function.
         *
         * Calls evaluate on the first sub-function supporting the isotope.
         * @param z Atomic number (>=1).
         * @param a Mass number (>=z).
         * @param T9 Temperature in 10^9 K.
         * @return Partition function value from supporting sub-function.
         * @throws std::runtime_error If no sub-function supports (z,a,T9).
+         */
         [[nodiscard]] double evaluate(int z, int a, double T9) const override;
        /**
         * @brief Evaluate temperature derivative of the composite function.
         *
         * Delegates to the first supporting sub-function's derivative.
         * @param z Atomic number.
         * @param a Mass number.
         * @param T9 Temperature in 10^9 K.
         * @return d/dT9 of the partition function.
         * @throws std::runtime_error If no sub-function supports (z,a,T9).
+         */
         [[nodiscard]] double evaluateDerivative(int z, int a, double T9) const override;
        /**
         * @brief Check support across all sub-functions.
         *
         * @param z Atomic number.
         * @param a Mass number.
         * @return true if any sub-function supports (z,a); false otherwise.
+         */
         [[nodiscard]] bool supports(int z, int a) const override;
        /**
         * @brief Get composite type identifier.
         *
         * Concatenates the type() strings of all sub-functions.
         * @return A string like "CompositePartitionFunction(func1, func2, ...)".
+         */
         [[nodiscard]] std::string type() const override;
        /**
         * @brief Clone this composite partition function.
         *
         * @return Unique pointer to a deep copy of this object.
+         */
         [[nodiscard]] std::unique_ptr<PartitionFunction> clone() const override {
             return std::make_unique<CompositePartitionFunction>(*this);
         }
     private:
         quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
         std::vector<std::unique_ptr<PartitionFunction>> m_partitionFunctions; ///< Set of partition functions to use in the composite partition function.
     private:
        /**
         * @brief Instantiate a sub-function by its type.
         *
         * @param type Enumeration value selecting the desired function implementation.
         * @return Unique pointer to a new PartitionFunction instance of the given type.
         * @throws std::runtime_error If the given type is not recognized.
+         */
         std::unique_ptr<PartitionFunction> selectPartitionFunction(const BasePartitionType type) const;
     };
 }