#pragma once

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    class UtilityError : public GridFireError {
        using GridFireError::GridFireError;
    };

    class HashingError final : public UtilityError {
        using UtilityError::UtilityError;
    };
}
