#pragma once

#include "gridfire/partition/partition_abstract.h"
#include "gridfire/partition/partition_types.h"

#include "fourdst/logging/logging.h"

#include <string>
#include <unordered_map>
#include <vector>

#include <memory>


namespace gridfire::partition {
    class CompositePartitionFunction final : public PartitionFunction {
    public:
        explicit CompositePartitionFunction(const std::vector<BasePartitionType>& partitionFunctions);
        CompositePartitionFunction(const CompositePartitionFunction& other);
        [[nodiscard]] double evaluate(int z, int a, double T9) const override;
        [[nodiscard]] double evaluateDerivative(int z, int a, double T9) const override;
        [[nodiscard]] bool supports(int z, int a) const override;
        [[nodiscard]] std::string type() const override;
        [[nodiscard]] std::unique_ptr<PartitionFunction> clone() const override {
            return std::make_unique<CompositePartitionFunction>(*this);
        }
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<std::unique_ptr<PartitionFunction>> m_partitionFunctions; ///< Set of partition functions to use in the composite partition function.
    private:
        std::unique_ptr<PartitionFunction> selectPartitionFunction(const BasePartitionType type) const;
    };
}