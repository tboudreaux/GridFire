#pragma once

#include "gridfire/reaction/reaction.h"

namespace gridfire::reaclib {

    [[maybe_unused]] static bool s_initialized = false;

    /**
     * @brief Provides global access to the fully initialized REACLIB reaction set.
     *
     * This function uses a singleton pattern to ensure the reaction data is loaded
     * from the embedded binary only once, the very first time it is requested.
     *
     * @return A constant reference to the application-wide reaction set.
     */
    const reaction::ReactionSet &get_all_reaclib_reactions();

    /// Name of the REACLIB snapshot compiled into this build (the file utils/reaclib/format.py converted).
    inline constexpr const char* DATA_SOURCE = "REACLIB (pynucastro snapshot reaclib_pynucastro_latest_Z26)";

    /// XXH64 of the REACLIB data compiled into this build: identifies the rate coefficients exactly.
    uint64_t data_hash();
} // namespace gridfire::reaclib
