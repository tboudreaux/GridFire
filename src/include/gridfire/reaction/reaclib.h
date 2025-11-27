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

    // Simple heuristic to check if a reaclib reaction is a strong or weak reaction
    /*  A weak reaction is defined here as one where:
        - The number of reactants is equal to the number of products
        - There is only one reactant and one product
        - The mass number (A) of the reactant is equal to the mass number (A) of the product
    */
    bool reaction_is_weak(const reaction::Reaction& reaction);

} // namespace gridfire::reaclib
