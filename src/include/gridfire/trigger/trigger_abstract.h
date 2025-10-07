#pragma once

#include "gridfire/trigger/trigger_result.h"

#include <string>

namespace gridfire::trigger {
    template <typename TriggerContextStruct>
    class Trigger {
    public:
        virtual ~Trigger() = default;

        virtual bool check(const TriggerContextStruct& ctx) const = 0;

        virtual void update(const TriggerContextStruct& ctx) = 0;
        virtual void reset() = 0;

        [[nodiscard]] virtual std::string name() const = 0;
        [[nodiscard]] virtual std::string describe() const = 0;
        [[nodiscard]] virtual TriggerResult why(const TriggerContextStruct& ctx) const = 0;

        [[nodiscard]] virtual size_t numTriggers() const = 0;
        [[nodiscard]] virtual size_t numMisses() const = 0;
    };
}