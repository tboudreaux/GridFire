#pragma once

#include "gridfire/partition/partition_abstract.h"

#include "fourdst/logging/logging.h"

#include <string>
#include <unordered_map>
#include <vector>

#include <memory>


namespace gridfire::partition {
    enum BasePartitionType {
        RauscherThielemann, ///< Rauscher-Thielemann partition function
        GroundState,        ///< Ground state partition function
    };

    inline std::unordered_map<BasePartitionType, std::string> basePartitionTypeToString = {
        {RauscherThielemann, "RauscherThielemann"},
        {GroundState, "GroundState"}
    };

    inline std::unordered_map<std::string, BasePartitionType> stringToBasePartitionType = {
        {"RauscherThielemann", RauscherThielemann},
        {"GroundState", GroundState}
    };

    class CompositePartitionFunction final : public PartitionFunction {
    public:
        explicit CompositePartitionFunction(const std::vector<BasePartitionType>& partitionFunctions);
        double evaluate(int z, int a, double T9) const override;
        double evaluateDerivative(int z, int a, double T9) const override;
        bool supports(int z, int a) const override;
        std::string type() const override;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::vector<std::unique_ptr<PartitionFunction>> m_partitionFunctions; ///< Set of partition functions to use in the composite partition function.
    private:
        std::unique_ptr<PartitionFunction> selectPartitionFunction(const BasePartitionType type) const;
    };
}