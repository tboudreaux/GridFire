#pragma once

#include "gridfire/engine/engine_abstract.h"

#include <string>
#include <vector>

namespace gridfire::utils {
    std::string formatNuclearTimescaleLogString(
        const DynamicEngine& engine,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    );
}