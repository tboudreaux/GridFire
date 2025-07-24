#include "gridfire/partition/composite/partition_composite.h"

#include <vector>
#include <set>

#include "gridfire/partition/partition_ground.h"
#include "gridfire/partition/partition_rauscher_thielemann.h"
#include "quill/LogMacros.h"

namespace gridfire::partition {
    CompositePartitionFunction::CompositePartitionFunction(
        const std::vector<BasePartitionType>& partitionFunctions
    ) {
        for (const auto& type : partitionFunctions) {
            LOG_TRACE_L2(m_logger, "Adding partition function of type: {}", basePartitionTypeToString[type]);
            m_partitionFunctions.push_back(selectPartitionFunction(type));
        }
    }

    CompositePartitionFunction::CompositePartitionFunction(const CompositePartitionFunction &other) {
        m_partitionFunctions.reserve(other.m_partitionFunctions.size());
        for (const auto& pf : other.m_partitionFunctions) {
            m_partitionFunctions.push_back(pf->clone());
        }
    }

    double CompositePartitionFunction::evaluate(int z, int a, double T9) const {
        LOG_TRACE_L3(m_logger, "Evaluating partition function for Z={} A={} T9={}", z, a, T9);
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
                LOG_TRACE_L3(m_logger, "Partition function of type {} supports Z={} A={}", partitionFunction->type(), z, a);
                return partitionFunction->evaluate(z, a, T9);
            } else {
                LOG_TRACE_L3(m_logger, "Partition function of type {} does not support Z={} A={}", partitionFunction->type(), z, a);
            }
        }
        LOG_ERROR(
            m_logger,
            "No partition function supports Z={} A={} T9={}. Tried: {}",
            z,
            a,
            T9,
            type()
        );
        throw std::runtime_error("No partition function supports the given Z, A, and T9 values.");
    }

    double CompositePartitionFunction::evaluateDerivative(int z, int a, double T9) const {
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
                LOG_TRACE_L3(m_logger, "Evaluating derivative of partition function for Z={} A={} T9={}", z, a, T9);
                return partitionFunction->evaluateDerivative(z, a, T9);
            }
        }
        LOG_ERROR(
            m_logger,
            "No partition function supports Z={} A={} T9={}. Tried: {}",
            z,
            a,
            T9,
            type()
        );
        throw std::runtime_error("No partition function supports the given Z, A, and T9 values.");
    }

    bool CompositePartitionFunction::supports(int z, int a) const {
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
                LOG_TRACE_L2(m_logger, "Partition function supports Z={} A={}", z, a);
                return true;
            }
        }
        return false;
    }

    std::string CompositePartitionFunction::type() const {
        std::stringstream ss;
        ss << "CompositePartitionFunction(";
        size_t count = 0;
        for (const auto& partitionFunction : m_partitionFunctions) {
            ss << partitionFunction->type();
            if (count < m_partitionFunctions.size() - 1) {
                ss << ", ";
            }
            count++;
        }
        ss << ")";
        std::string types = ss.str();
        return types;
    }

    std::unique_ptr<PartitionFunction> CompositePartitionFunction::selectPartitionFunction(
        const BasePartitionType type
    ) const {
        switch (type) {
            case RauscherThielemann: {
                return std::make_unique<RauscherThielemannPartitionFunction>();
            }
            case GroundState: {
                return std::make_unique<GroundStatePartitionFunction>();
            }
            default: {
                LOG_ERROR(m_logger, "Unknown partition function type");
                throw std::runtime_error("Unknown partition function type");
            }
        }
    }
}
