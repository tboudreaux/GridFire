#pragma once

#include "gridfire/trigger/trigger_result.h"

#include <string>

namespace gridfire::trigger {
    /**
     * @brief Generic trigger interface for signaling events/conditions during integration.
     *
     * A Trigger encapsulates a condition evaluated against a user-defined context struct
     * (TriggerContextStruct). The interface supports:
     *  - check(): evaluate the current condition without mutating external state
     *  - update(): feed the context to update internal state (counters, sliding windows, etc.)
     *  - reset(): clear internal state and counters
     *  - name()/describe(): human-readable identification and textual description
     *  - why(): structured explanation (tree of TriggerResult) of the last evaluation
     *  - numTriggers()/numMisses(): simple counters for diagnostics
     *
     * Logical compositions (AND/OR/NOT/EveryNth) are provided in trigger_logical.h and implement
     * this interface for any TriggerContextStruct.
     *
     * @tparam TriggerContextStruct User-defined context passed to triggers (e.g., timestep info).
     */
    template <typename TriggerContextStruct>
    class Trigger {
    public:
        /**
         * @brief Virtual destructor for polymorphic use.
         */
        virtual ~Trigger() = default;

        /**
         * @brief Evaluate the trigger condition against the provided context.
         * @param ctx Context snapshot (immutable view) used to evaluate the condition.
         * @return true if the condition is satisfied; false otherwise.
         */
        virtual bool check(const TriggerContextStruct& ctx) const = 0;

        /**
         * @brief Update any internal state with the given context (e.g., counters, windows).
         * @param ctx Context snapshot used to update state.
         */
        virtual void update(const TriggerContextStruct& ctx) = 0;
        /**
         * @brief Reset internal state and diagnostics counters.
         */
        virtual void reset() = 0;

        /**
         * @brief Short, stable name for this trigger (suitable for logs/UI).
         */
        [[nodiscard]] virtual std::string name() const = 0;
        /**
         * @brief Human-readable description of this trigger's logic.
         */
        [[nodiscard]] virtual std::string describe() const = 0;
        /**
         * @brief Explain why the last evaluation would be true/false in a structured way.
         * @param ctx Context snapshot for the explanation.
         * @return A TriggerResult tree with a boolean value and per-cause details.
         */
        [[nodiscard]] virtual TriggerResult why(const TriggerContextStruct& ctx) const = 0;

        /**
         * @brief Total number of times this trigger evaluated to true since last reset.
         */
        [[nodiscard]] virtual size_t numTriggers() const = 0;
        /**
         * @brief Total number of times this trigger evaluated to false since last reset.
         */
        [[nodiscard]] virtual size_t numMisses() const = 0;
    };
}