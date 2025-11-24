#pragma once

#include <vector>
#include <string>

namespace gridfire::trigger {
    /**
     * @file trigger_result.h
     * @brief Structured explanation node for trigger evaluations.
     *
     * TriggerResult represents a tree describing why a trigger evaluated to true/false.
     * Each node contains a boolean value, a short name, a human-readable description,
     * and optional nested causes for composite triggers (e.g., AND/OR/NOT).
     *
     * @par Example
     * @code
     * // Produce a result and pretty-print it
     * gridfire::trigger::TriggerResult r{"A>5", "Threshold passed", true, {}};
     * // See procedures/trigger_pprint.h for printWhy()
     * // gridfire::trigger::printWhy(r);
     * @endcode
     */
    struct TriggerResult {
        std::string name;        ///< Short identifier for the condition (e.g., "Temperature Rise").
        std::string description; ///< Human-readable reason summarizing the outcome at this node.
        bool value;              ///< Evaluation result for this node (true/false).
        std::vector<TriggerResult> causes; ///< Sub-reasons for composite triggers.
    };
}