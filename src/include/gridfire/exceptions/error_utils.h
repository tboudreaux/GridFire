#pragma once

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    /**
     * @brief Base class for utility-related errors in GridFire.
     *
     * This class serves as a base for all exceptions related to utility functions
     */
    class UtilityError : public GridFireError {
        using GridFireError::GridFireError;
    };

    /**
     * @brief Exception class for hashing-related errors in GridFire.
     *
     * This class is used to represent errors that occur during hashing operations.
     */
    class HashingError final : public UtilityError {
        using UtilityError::UtilityError;
    };
}
