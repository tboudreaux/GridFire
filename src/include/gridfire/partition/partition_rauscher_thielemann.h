#pragma once

#include "gridfire/partition/partition_abstract.h"

#include "fourdst/logging/logging.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>
#include <memory>

namespace gridfire::partition {
    class RauscherThielemannPartitionFunction final : public PartitionFunction {
    public:
        RauscherThielemannPartitionFunction();
        double evaluate(int z, int a, double T9) const override;
        double evaluateDerivative(int z, int a, double T9) const override;
        bool supports(int z, int a) const override;
        std::string type() const override { return "RauscherThielemann"; }
    private:
        enum Bounds {
            FRONT,
            BACK,
            MIDDLE
        };
    private:
        struct IsotopeData {
            double ground_state_spin;
            std::array<double, 24> normalized_g_values;
        };
        struct InterpolationPoints {
            double T9_high;
            double G_norm_high;
            double T9_low;
            double G_norm_low;
        };
        struct IdentifiedIsotope {
            Bounds bound;
            const IsotopeData& data;
            size_t upperIndex;
            size_t lowerIndex;
        };
        std::unique_ptr<PartitionFunction> clone() const override {
            return std::make_unique<RauscherThielemannPartitionFunction>(*this);
        }
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        std::unordered_map<int, IsotopeData> m_partitionData;
    private:
        static InterpolationPoints get_interpolation_points(
            const size_t upper_index,
            const size_t lower_index,
            const std::array<double, 24>& normalized_g_values
        );
        IdentifiedIsotope find(int z, int a, double T9) const;
        static constexpr int make_key(int z, int a);
    };
}
