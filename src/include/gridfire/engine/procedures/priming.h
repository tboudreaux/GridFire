#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/types/types.h"

#include "fourdst/atomic/atomicSpecies.h"

#include "gridfire/engine/scratchpads/blob.h"

#include <optional>


namespace gridfire::engine {

    /**
     * @brief Primes absent species in the network to their equilibrium abundances.
     *
     * Executes a network priming algorithm that iteratively rebuilds the reaction network,
     * calculates equilibrium mass fractions for species with zero initial abundance,
     * and applies mass transfers based on reaction flows.
     *
     * Refer to priming.cpp for implementation details on logging, algorithmic steps, and error handling.
     *
     * @param ctx
     * @param netIn Input network data containing initial composition, temperature, and density.
     * @param engine DynamicEngine used to build and evaluate the reaction network.
     * @param ignoredReactionTypes Types of reactions to ignore during priming (e.g., weak reactions).
     * @pre netIn.composition defines species and their mass fractions; engine is constructed with a valid network.
     * @post engine.networkReactions restored to its initial state; returned report contains primedComposition,
     *       massFractionChanges for each species, success flag, and status code.
     * @return PrimingReport encapsulating the results of the priming operation.
     */
     PrimingReport primeNetwork(
         scratch::StateBlob &ctx,
         const NetIn& netIn,
         const GraphEngine& engine, const std::optional<std::vector<reaction::ReactionType>>& ignoredReactionTypes
     );
 }