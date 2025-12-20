#pragma once

#include <string>
#include <stdexcept>

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    class ScratchPadError : public GridFireError {
    public:
        explicit ScratchPadError(const std::string& msg)
            : GridFireError("ScratchPadError: " + msg) {}
    };
}