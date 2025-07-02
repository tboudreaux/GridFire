#pragma once

#include "gridfire/partition/partition_abstract.h"

#include "fourdst/logging/logging.h"

#include <unordered_map>

#include "quill/Logger.h"

namespace gridfire::partition {
    class GroundStatePartitionFunction final : public PartitionFunction {
    public:
        GroundStatePartitionFunction();
        double evaluate(
            const int z,
            const int a,
            const double T9
        ) const override;
        double evaluateDerivative(
            const int z,
            const int a,
            const double T9
        ) const override;
        bool supports(
            const int z,
            const int a
        ) const override;
        std::string type() const override { return "GroundState"; }
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::unordered_map<int, double> m_ground_state_spin;
        static constexpr int make_key(
            const int z,
        const int a);
    };

}
