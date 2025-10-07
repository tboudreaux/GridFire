#pragma once

#include "gridfire/trigger/trigger_result.h"

#include <iostream>

namespace gridfire::trigger {
    inline void printWhy(const TriggerResult& result, const int indent = 0) { // NOLINT(*-no-recursion)
        const std::string prefix(indent * 2, ' ');
        std::cout << prefix << "• [" << (result.value ? "TRUE" : "FALSE")
                  << "] " << result.name << ": " << result.description << std::endl;

        for (const auto& cause : result.causes) {
            printWhy(cause, indent + 1);
        }
    }
}