#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"

#include <map>
#include <ranges>
#include <sstream>


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
     * @pre netIn.composition defines species and their mass fractions; engine is constructed with a valid network.
     * @post engine.networkReactions restored to its initial state; returned report contains primedComposition,
     *       massFractionChanges for each species, success flag, and status code.
     * @return PrimingReport encapsulating the results of the priming operation.
     */
     PrimingReport primeNetwork(
         const NetIn&,
         DynamicEngine& engine
     );

    /**
     * @brief Computes the destruction rate constant for a specific species.
     *
     * Calculates the sum of molar reaction flows for all reactions where the species
     * is a reactant (negative stoichiometry) after scaling its abundance to unity.
     *
     * @param engine Engine providing the current set of network reactions and flow calculations.
     * @param species The atomic species whose destruction rate is computed.
     * @param Y Vector of molar abundances for all species in the engine.
     * @param T9 Temperature in units of 10^9 K.
     * @param rho Density of the medium.
     * @pre Y.size() matches engine.getNetworkReactions().size() mapping species order.
     * @post Returned rate constant is non-negative.
     * @return Sum of absolute stoichiometry-weighted destruction flows for the species.
     */
     double calculateDestructionRateConstant(
         const DynamicEngine& engine,
         const fourdst::atomic::Species& species,
         const std::vector<double>& Y,
         double T9,
         double rho
     );

    /**
     * @brief Computes the creation rate for a specific species.
     *
     * Sums molar reaction flows for all reactions where the species
     * appears as a product (positive stoichiometry).
     *
     * @param engine Engine providing the current set of network reactions and flow calculations.
     * @param species The atomic species whose creation rate is computed.
     * @param Y Vector of molar abundances for all species in the engine.
     * @param T9 Temperature in units of 10^9 K.
     * @param rho Density of the medium.
     * @pre Y.size() matches engine.getNetworkReactions().size() mapping species order.
     * @post Returned creation rate is non-negative.
     * @return Sum of stoichiometry-weighted creation flows for the species.
     */
     double calculateCreationRate(
         const DynamicEngine& engine,
         const fourdst::atomic::Species& species,
         const std::vector<double>& Y,
         double T9,
         double rho
     );
 }