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
    const reaction::LogicalReactionSet& get_all_reactions();

} // namespace gridfire::reaclib