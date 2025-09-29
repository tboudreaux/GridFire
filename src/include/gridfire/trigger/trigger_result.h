#pragma once

#include <vector>
#include <string>

namespace gridfire::trigger {
    struct TriggerResult {
        std::string name;
        std::string description;
        bool value;
        std::vector<TriggerResult> causes;
    };
}