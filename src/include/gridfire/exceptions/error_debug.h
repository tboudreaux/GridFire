#pragma once

#include <source_location>
#include <string>

#include "gridfire/exceptions/error_gridfire.h"

namespace gridfire::exceptions {
    class DebugException final : public GridFireError{
    public:
    #ifdef NDEBUG
        #if defined(__clang__)
            __attribute__((unavailable("DebugExceptions may not be used in release builds")))

        #elif defined(__GNUC__)
            __attribute__((error("DebugExceptions may not be used in release builds")))
        #endif
    #endif

    explicit DebugException(const std::string_view message,
                            const std::source_location loc = std::source_location::current())
                : GridFireError(format_error(message, loc)) {
            }
        private:
            static std::string format_error(std::string_view message, const std::source_location loc) {
                return std::format("[DEBUG HALT] {}:{}: {}",
                    loc.file_name(),
                    loc.line(),
                    message);
            }
        };
}