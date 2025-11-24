#include "gridfire/partition/composite/partition_composite.h"

#include <vector>
#include <set>

#include "gridfire/partition/partition_ground.h"
#include "gridfire/partition/partition_rauscher_thielemann.h"
#include "gridfire/utils/hashing.h"
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
        const uint_fast32_t hash = utils::hash_atomic(a, z);
        if (m_supportCache.contains(hash)) {
            return m_supportCache.at(hash).evaluate(z, a, T9);
        }
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
                m_supportCache.emplace(hash, *partitionFunction);
                return partitionFunction->evaluate(z, a, T9);
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
        const uint_fast32_t hash = utils::hash_atomic(a, z);
        if (m_supportCache.contains(hash)) {
            return m_supportCache.at(hash).evaluateDerivative(z, a, T9);
        }
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
                m_supportCache.emplace(hash, *partitionFunction);
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
        const uint_fast32_t hash = utils::hash_atomic(a, z);
        if (m_supportCache.contains(hash)) {
            return true;
        }
        for (const auto& partitionFunction : m_partitionFunctions) {
            if (partitionFunction->supports(z, a)) {
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
