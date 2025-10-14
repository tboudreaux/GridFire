#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/composition/atomicSpecies.h"


namespace gridfire {

    /**
     * @brief Primes absent species in the network to their equilibrium abundances.
     *
     * Executes a network priming algorithm that iteratively rebuilds the reaction network,
     * calculates equilibrium mass fractions for species with zero initial abundance,
     * and applies mass transfers based on reaction flows.
     *
     * Refer to priming.cpp for implementation details on logging, algorithmic steps, and error handling.
     *
     * @param netIn Input network data containing initial composition, temperature, and density.
     * @param engine DynamicEngine used to build and evaluate the reaction network.
     * @param ignoredReactionTypes Types of reactions to ignore during priming (e.g., weak reactions).
     * @pre netIn.composition defines species and their mass fractions; engine is constructed with a valid network.
     * @post engine.networkReactions restored to its initial state; returned report contains primedComposition,
     *       massFractionChanges for each species, success flag, and status code.
     * @return PrimingReport encapsulating the results of the priming operation.
     */
     PrimingReport primeNetwork(
         const NetIn& netIn,
         DynamicEngine& engine,
         const std::optional<std::vector<reaction::ReactionType>>& ignoredReactionTypes
     );

    /**
     * @brief Computes the destruction rate constant for a specific species.
     *
     * Calculates the sum of molar reaction flows for all reactions where the species
     * is a reactant (negative stoichiometry) after scaling its abundance to unity.
     *
     * @param engine Engine providing the current set of network reactions and flow calculations.
     * @param species The atomic species whose destruction rate is computed.
     * @param composition Current composition providing abundances for all species.
     * @param T9 Temperature in units of 10^9 K.
     * @param rho Density of the medium.
     * @param reactionTypesToIgnore types of reactions to ignore during calculation.
     * @pre Y.size() matches engine.getNetworkReactions().size() mapping species order.
     * @post Returned rate constant is non-negative.
     * @return Sum of absolute stoichiometry-weighted destruction flows for the species.
     */
     double calculateDestructionRateConstant(
         const DynamicEngine& engine,
         const fourdst::atomic::Species& species,
         const fourdst::composition::Composition& composition,
         double T9,
         double rho, const std::optional<std::vector<reaction::ReactionType>> &
         reactionTypesToIgnore
     );

    /**
     * @brief Computes the creation rate for a specific species.
     *
     * Sums molar reaction flows for all reactions where the species
     * appears as a product (positive stoichiometry).
     *
     * @param engine Engine providing the current set of network reactions and flow calculations.
     * @param species The atomic species whose creation rate is computed.
     * @param composition Composition object containing current abundances.
     * @param T9 Temperature in units of 10^9 K.
     * @param rho Density of the medium.
     * @param reactionTypesToIgnore types of reactions to ignore during calculation.
     * @pre Y.size() matches engine.getNetworkReactions().size() mapping species order.
     * @post Returned creation rate is non-negative.
     * @return Sum of stoichiometry-weighted creation flows for the species.
     */
     double calculateCreationRate(
         const DynamicEngine& engine,
         const fourdst::atomic::Species& species,
         const fourdst::composition::Composition& composition,
         double T9,
         double rho, const std::optional<std::vector<reaction::ReactionType>> &reactionTypesToIgnore
     );
 }