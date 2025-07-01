#pragma once

#include "gridfire/screening/screening_abstract.h"

#include <memory>

namespace gridfire::screening {
    enum class ScreeningType {
        BARE,    ///< No screening applied
        WEAK,   ///< Weak screening model
    };

    std::unique_ptr<ScreeningModel> selectScreeningModel(ScreeningType type);
}