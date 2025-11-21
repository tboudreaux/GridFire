#pragma once

#include <exception>
#include <string>
#include <utility>

#include "gridfire/exceptions/gridfire_exception.h"

namespace gridfire::exceptions {
    class UtilityError : public GridFireError {
        using GridFireError::GridFireError;
    };

    class HashingError final : public UtilityError {
        using UtilityError::UtilityError;
    };
}
