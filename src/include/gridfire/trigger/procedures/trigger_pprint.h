#pragma once

#include "gridfire/trigger/trigger_result.h"

#include <iostream>

namespace gridfire::trigger {
    /**
     * @brief Pretty-print a TriggerResult explanation tree to std::cout.
     *
     * Prints one line per node prefixed with Unicode bullets and indentation to visualize
     * the explanation hierarchy. Each line shows [TRUE|FALSE], the node name, and description.
     *
     * @param result Root TriggerResult to print.
     * @param indent Current indentation level (number of two-space indents); callers typically
     *               omit this parameter and let recursion handle it.
     *
     * @par Example
     * @code
     * using gridfire::trigger::TriggerResult;
     * TriggerResult leaf{"A>5", "Threshold passed", true, {}};
     * TriggerResult root{"AND", "Both conditions true", true, {leaf}};
     * gridfire::trigger::printWhy(root);
     * @endcode
     */
    inline void printWhy(const TriggerResult& result, const int indent = 0) { // NOLINT(*-no-recursion)
        const std::string prefix(indent * 2, ' ');
        std::cout << prefix << "• [" << (result.value ? "TRUE" : "FALSE")
                  << "] " << result.name << ": " << result.description << std::endl;

        for (const auto& cause : result.causes) {
            printWhy(cause, indent + 1);
        }
    }
}