#include "gridfire/engine/procedures/priming.h"
#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/engine/procedures/construction.h"
#include "gridfire/solver/solver.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

namespace gridfire {
    using fourdst::composition::Composition;
    using fourdst::atomic::Species;

    const reaction::Reaction* findDominantCreationChannel (
        const DynamicEngine& engine,
        const Species& species,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) {
        const reaction::Reaction* dominateReaction = nullptr;
        double maxFlow = -1.0;
        for (const auto& reaction : engine.getNetworkReactions()) {
            if (reaction->contains(species) && reaction->stoichiometry(species) > 0) {
                const double flow = engine.calculateMolarReactionFlow(*reaction, comp, T9, rho);
                if (flow > maxFlow) {
                    maxFlow = flow;
                    dominateReaction = reaction.get();
                }
            }
        }
        return dominateReaction;
    }

    /**
     * @brief Primes absent species in the network to their equilibrium abundances using a robust, two-stage approach.
     *
     * @details This function implements a robust network priming algorithm that avoids the pitfalls of
     * sequential, one-by-one priming. The previous, brittle method could allow an early priming
     * reaction to consume all of a shared reactant, starving later reactions. This new, two-stage
     * method ensures that all priming reactions are considered collectively, competing for the
     * same limited pool of initial reactants in a physically consistent manner.
     *
     * The algorithm proceeds in three main stages:
     * 1.  **Calculation Stage:** It first loops through all species that need priming. For each one,
     * it calculates its theoretical equilibrium mass fraction and identifies the dominant
     * creation channel. Crucially, it *does not* modify any abundances at this stage. Instead,
     * it stores these calculations as a list of "mass transfer requests".
     *
     * 2.  **Collective Scaling Stage:** It then processes the full list of requests to determine the
     * total "debit" required from each reactant. By comparing these total debits to the
     * initially available mass of each reactant, it calculates a single, global `scalingFactor`.
     * If any reactant is overdrawn, this factor will be less than 1.0, ensuring that no
     * reactant's abundance can go negative.
     *
     * 3.  **Application Stage:** Finally, it loops through the requests again, applying the mass
     * transfers. Each calculated equilibrium mass fraction and corresponding reactant debit is
     * multiplied by the global `scalingFactor` before being applied to the final composition.
     * This ensures that if resources are limited, all primed species are scaled down proportionally.
     *
     * @param netIn Input network data containing initial composition, temperature, and density.
     * @param engine DynamicEngine used to build and evaluate the reaction network.
     * @return PrimingReport encapsulating the results of the priming operation, including the new
     * robustly primed composition.
     */
    PrimingReport primeNetwork(const NetIn& netIn, DynamicEngine& engine) {
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        // --- Initial Setup ---
        // Identify all species with zero initial mass fraction that need to be primed.
        std::vector<Species> speciesToPrime;
        for (const auto &entry: netIn.composition | std::views::values) {
            if (entry.mass_fraction() == 0.0) {
                speciesToPrime.push_back(entry.isotope());
            }
        }
        LOG_DEBUG(logger, "Priming {} species in the network.", speciesToPrime.size());

        // If no species need priming, return immediately.
        PrimingReport report;
        if (speciesToPrime.empty()) {
            report.primedComposition = netIn.composition;
            report.success = true;
            report.status = PrimingReportStatus::NO_SPECIES_TO_PRIME;
            return report;
        }

        const double T9 = netIn.temperature / 1e9;
        const double rho = netIn.density;
        const auto initialReactionSet = engine.getNetworkReactions();
        report.status = PrimingReportStatus::FULL_SUCCESS;
        report.success = true;

        // Create a mutable map of the mass fractions that we will modify.
        std::unordered_map<Species, double> currentMassFractions;
        for (const auto& entry : netIn.composition | std::views::values) {
            currentMassFractions[entry.isotope()] = entry.mass_fraction();
        }
        // Ensure all species to be primed exist in the map, initialized to zero.
        for (const auto& entry : speciesToPrime) {
            currentMassFractions[entry] = 0.0;
        }

        // Rebuild the engine with the full network to ensure all possible creation channels are available.
        engine.rebuild(netIn.composition, NetworkBuildDepth::Full);

        // --- STAGE 1: Calculation and Bookkeeping (No Modifications) ---
        // In this stage, we calculate all required mass transfers but do not apply them yet.

        // A struct to hold the result of each individual priming calculation.
        struct MassTransferRequest {
            Species species_to_prime;
            double equilibrium_mass_fraction;
            std::vector<Species> reactants;
        };
        std::vector<MassTransferRequest> requests;

        for (const auto& primingSpecies : speciesToPrime) {
            // Create a temporary composition reflecting the current state for rate calculations.
            Composition tempComp;
            for(const auto& [sp, mf] : currentMassFractions) {
                tempComp.registerSymbol(std::string(sp.name()));
                tempComp.setMassFraction(sp, std::max(0.0, mf));
            }
            tempComp.finalize(true);

            NetworkPrimingEngineView primer(primingSpecies, engine);
            if (primer.getNetworkReactions().size() == 0) {
                LOG_ERROR(logger, "No priming reactions found for species {}.", primingSpecies.name());
                report.success = false;
                report.status = PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS;
                continue;
            }

            const double destructionRateConstant = calculateDestructionRateConstant(primer, primingSpecies, tempComp, T9, rho);

            if (destructionRateConstant > 1e-99) {
                const double creationRate = calculateCreationRate(primer, primingSpecies, tempComp, T9, rho);
                double equilibriumMassFraction = (creationRate / destructionRateConstant) * primingSpecies.mass();
                if (std::isnan(equilibriumMassFraction)) equilibriumMassFraction = 0.0;

                if (const reaction::Reaction* dominantChannel = findDominantCreationChannel(primer, primingSpecies, tempComp, T9, rho)) {
                    // Store the request instead of applying it immediately.
                    requests.push_back({primingSpecies, equilibriumMassFraction, dominantChannel->reactants()});
                } else {
                    LOG_ERROR(logger, "Failed to find dominant creation channel for {}.", primingSpecies.name());
                    report.status = PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL;
                }
            } else {
                LOG_WARNING(logger, "No destruction channel found for {}. Using fallback abundance.", primingSpecies.name());
                // For species with no destruction, we can't calculate an equilibrium.
                // We add a request with a tiny fallback abundance to ensure it exists in the network.
                requests.push_back({primingSpecies, 1e-40, {}});
            }
        }

        // --- STAGE 2: Collective Scaling Based on Reactant Availability ---
        // Now, we determine the total demand for each reactant and find a global scaling factor.

        std::unordered_map<Species, double> total_mass_debits;
        for (const auto& req : requests) {
            if (req.reactants.empty()) continue; // Skip fallbacks which don't consume reactants.

            double totalReactantMass = 0.0;
            for (const auto& reactant : req.reactants) {
                totalReactantMass += reactant.mass();
            }
            if (totalReactantMass == 0.0) continue;

            for (const auto& reactant : req.reactants) {
                const double massToSubtract = req.equilibrium_mass_fraction * (reactant.mass() / totalReactantMass);
                total_mass_debits[reactant] += massToSubtract;
            }
        }

        double globalScalingFactor = 1.0;
        for (const auto& [reactant, total_debit] : total_mass_debits) {
            double availableMass;
            if (currentMassFractions.contains(reactant)) {
                availableMass = currentMassFractions.at(reactant);
            } else {
                availableMass = 0.0;
            }
            if (total_debit > availableMass && availableMass > 0) {
                globalScalingFactor = std::min(globalScalingFactor, availableMass / total_debit);
            }
        }

        if (globalScalingFactor < 1.0) {
            LOG_WARNING(logger, "Priming was limited by reactant availability. All transfers will be scaled by {:.4e}", globalScalingFactor);
        }

        // --- STAGE 3: Application of Scaled Mass Transfers ---
        // Finally, apply all the transfers, scaled by our global factor.

        std::unordered_map<Species, double> totalMassFractionChanges;
        for (const auto&[species_to_prime, equilibrium_mass_fraction, reactants] : requests) {
            const double scaled_equilibrium_mf = equilibrium_mass_fraction * globalScalingFactor;

            // Add the scaled mass to the primed species.
            currentMassFractions.at(species_to_prime) += scaled_equilibrium_mf;
            totalMassFractionChanges[species_to_prime] += scaled_equilibrium_mf;

            // Subtract the scaled mass from the reactants.
            if (!reactants.empty()) {
                double totalReactantMass = 0.0;
                for (const auto& reactant : reactants) {
                    totalReactantMass += reactant.mass();
                }
                if (totalReactantMass == 0.0) continue;

                for (const auto& reactant : reactants) {
                    const double massToSubtract = scaled_equilibrium_mf * (reactant.mass() / totalReactantMass);
                    if (massToSubtract != 0) {
                        currentMassFractions.at(reactant) -= massToSubtract;
                        totalMassFractionChanges[reactant] -= massToSubtract;
                    }
                }
            }
        }

        // --- Final Composition Construction ---
        std::vector<std::string> final_symbols;
        std::vector<double> final_mass_fractions;
        for(const auto& [species, mass_fraction] : currentMassFractions) {
            final_symbols.emplace_back(species.name());
            final_mass_fractions.push_back(std::max(0.0, mass_fraction)); // Ensure no negative mass fractions.
        }

        Composition primedComposition(final_symbols, final_mass_fractions, true);

        report.primedComposition = primedComposition;
        for (const auto& [species, change] : totalMassFractionChanges) {
            report.massFractionChanges.emplace_back(species, change);
        }

        // Restore the engine to its original, smaller network state.
        engine.setNetworkReactions(initialReactionSet);
        return report;
    }



    // PrimingReport primeNetwork(const NetIn& netIn, DynamicEngine& engine) {
    //     auto logger = LogManager::getInstance().getLogger("log");
    //
    //     std::vector<Species> speciesToPrime;
    //     for (const auto &entry: netIn.composition | std::views::values) {
    //         std::cout << "Checking species: " << entry.isotope().name() << " with mass fraction: " << entry.mass_fraction() << std::endl;
    //         if (entry.mass_fraction() == 0.0) {
    //             speciesToPrime.push_back(entry.isotope());
    //         }
    //     }
    //     LOG_DEBUG(logger, "Priming {} species in the network.", speciesToPrime.size());
    //
    //     PrimingReport report;
    //     if (speciesToPrime.empty()) {
    //         report.primedComposition = netIn.composition;
    //         report.success = true;
    //         report.status = PrimingReportStatus::NO_SPECIES_TO_PRIME;
    //         return report;
    //     }
    //
    //     const double T9 = netIn.temperature / 1e9;
    //     const double rho = netIn.density;
    //     const auto initialReactionSet = engine.getNetworkReactions();
    //
    //     report.status = PrimingReportStatus::FULL_SUCCESS;
    //     report.success = true;
    //
    //     // --- 1: pack composition into internal map ---
    //     std::unordered_map<Species, double> currentMassFractions;
    //     for (const auto& entry : netIn.composition | std::views::values) {
    //         currentMassFractions[entry.isotope()] = entry.mass_fraction();
    //     }
    //     for (const auto& entry : speciesToPrime) {
    //         currentMassFractions[entry] = 0.0; // Initialize priming species with 0 mass fraction
    //     }
    //
    //     std::unordered_map<Species, double> totalMassFractionChanges;
    //
    //     engine.rebuild(netIn.composition, NetworkBuildDepth::Full);
    //
    //     for (const auto& primingSpecies : speciesToPrime) {
    //         LOG_TRACE_L3(logger, "Priming species: {}", primingSpecies.name());
    //
    //         // Create a temporary composition from the current internal state for the primer
    //         Composition tempComp;
    //         for(const auto& [sp, mf] : currentMassFractions) {
    //             tempComp.registerSymbol(std::string(sp.name()));
    //             if (mf < 0.0 && std::abs(mf) < 1e-16) {
    //                 tempComp.setMassFraction(sp, 0.0); // Avoid negative mass fractions
    //             } else {
    //                 tempComp.setMassFraction(sp, mf);
    //             }
    //         }
    //         tempComp.finalize(true); // Finalize with normalization
    //
    //         NetIn tempNetIn = netIn;
    //         tempNetIn.composition = tempComp;
    //
    //         NetworkPrimingEngineView primer(primingSpecies, engine);
    //
    //         if (primer.getNetworkReactions().size() == 0) {
    //             LOG_ERROR(logger, "No priming reactions found for species {}.", primingSpecies.name());
    //             report.success = false;
    //             report.status = PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS;
    //             continue;
    //         }
    //
    //         const auto Y = primer.mapNetInToMolarAbundanceVector(tempNetIn);
    //         const double destructionRateConstant = calculateDestructionRateConstant(primer, primingSpecies, Y, T9, rho);
    //
    //         if (destructionRateConstant > 1e-99) {
    //             double equilibriumMassFraction = 0.0;
    //             const double creationRate = calculateCreationRate(primer, primingSpecies, Y, T9, rho);
    //             equilibriumMassFraction = (creationRate / destructionRateConstant) * primingSpecies.mass();
    //             if (std::isnan(equilibriumMassFraction)) {
    //                 LOG_WARNING(logger, "Equilibrium mass fraction for {} is NaN. Setting to 0.0. This is likely not an issue. It probably originates from all reactions leading to creation and destruction being frozen out. In that case 0.0 should be a good approximation. Hint: This happens often when the network temperature is very the low. ", primingSpecies.name());
    //                 equilibriumMassFraction = 0.0;
    //             }
    //             LOG_TRACE_L3(logger, "Found equilibrium for {}: X_eq = {:.4e}", primingSpecies.name(), equilibriumMassFraction);
    //
    //             if (const reaction::Reaction* dominantChannel = findDominantCreationChannel(primer, primingSpecies, Y, T9, rho)) {
    //                 LOG_TRACE_L3(logger, "Dominant creation channel for {}: {}", primingSpecies.name(), dominantChannel->id());
    //
    //                 double totalReactantMass = 0.0;
    //                 for (const auto& reactant : dominantChannel->reactants()) {
    //                     totalReactantMass += reactant.mass();
    //                 }
    //
    //                 double scalingFactor = 1.0;
    //                 for (const auto& reactant : dominantChannel->reactants()) {
    //                     const double massToSubtract = equilibriumMassFraction * (reactant.mass() / totalReactantMass);
    //                     double availableMass = 0.0;
    //                     if (currentMassFractions.contains(reactant)) {
    //                         availableMass = currentMassFractions.at(reactant);
    //                     }
    //                     if (massToSubtract > availableMass && availableMass > 0) {
    //                         scalingFactor = std::min(scalingFactor, availableMass / massToSubtract);
    //                     }
    //                 }
    //
    //                 if (scalingFactor < 1.0) {
    //                     LOG_WARNING(logger, "Priming for {} was limited by reactant availability. Scaling transfer by {:.4e}", primingSpecies.name(), scalingFactor);
    //                     equilibriumMassFraction *= scalingFactor;
    //                 }
    //
    //                 // Update the internal mass fraction map and accumulate total changes
    //                 totalMassFractionChanges[primingSpecies] += equilibriumMassFraction;
    //                 currentMassFractions.at(primingSpecies) += equilibriumMassFraction;
    //
    //                 for (const auto& reactant : dominantChannel->reactants()) {
    //                     const double massToSubtract = equilibriumMassFraction * (reactant.mass() / totalReactantMass);
    //                     std::cout << "[Priming: " << primingSpecies.name() << ", Channel: " << dominantChannel->id() <<  "] Subtracting " << massToSubtract << " from reactant " << reactant.name() << std::endl;
    //                     totalMassFractionChanges[reactant] -= massToSubtract;
    //                     currentMassFractions[reactant] -= massToSubtract;
    //                 }
    //             } else {
    //                  LOG_ERROR(logger, "Failed to find dominant creation channel for {}.", primingSpecies.name());
    //                  report.status = PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL;
    //                  totalMassFractionChanges[primingSpecies] += 1e-40;
    //                  currentMassFractions.at(primingSpecies) += 1e-40;
    //             }
    //         } else {
    //             LOG_WARNING(logger, "No destruction channel found for {}. Using fallback abundance.", primingSpecies.name());
    //             totalMassFractionChanges.at(primingSpecies) += 1e-40;
    //             currentMassFractions.at(primingSpecies) += 1e-40;
    //             report.status = PrimingReportStatus::BASE_NETWORK_TOO_SHALLOW;
    //         }
    //     }
    //
    //     // --- Final Composition Construction ---
    //     std::vector<std::string> final_symbols;
    //     std::vector<double> final_mass_fractions;
    //     for(const auto& [species, mass_fraction] : currentMassFractions) {
    //         final_symbols.emplace_back(species.name());
    //         if (mass_fraction < 0.0 && std::abs(mass_fraction) < 1e-16) {
    //             final_mass_fractions.push_back(0.0); // Avoid negative mass fractions
    //         } else {
    //             final_mass_fractions.push_back(mass_fraction);
    //         }
    //     }
    //
    //     // Create the final composition object from the pre-normalized mass fractions
    //     Composition primedComposition(final_symbols, final_mass_fractions, true);
    //
    //     report.primedComposition = primedComposition;
    //     for (const auto& [species, change] : totalMassFractionChanges) {
    //         report.massFractionChanges.emplace_back(species, change);
    //     }
    //
    //     engine.setNetworkReactions(initialReactionSet);
    //     return report;
    // }

  double calculateDestructionRateConstant(
        const DynamicEngine& engine,
        const Species& species,
        const Composition& comp,
        const double T9,
        const double rho
    ) {
        //TODO: previously (when using raw vectors) I set y[speciesIndex] = 1.0 to let there be enough so that just the destruction rate could be found (without bottlenecks from abundance) we will need to do a similar thing here.
        double destructionRateConstant = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            if (reaction->contains_reactant(species)) {
                const int stoichiometry = reaction->stoichiometry(species);
                destructionRateConstant += std::abs(stoichiometry) * engine.calculateMolarReactionFlow(*reaction, comp, T9, rho);
            }
        }
        return destructionRateConstant;
    }

    double calculateCreationRate(
        const DynamicEngine& engine,
        const Species& species,
        const Composition& comp,
        const double T9,
        const double rho
    ) {
        double creationRate = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            const int stoichiometry = reaction->stoichiometry(species);
            if (stoichiometry > 0) {
                if (engine.calculateMolarReactionFlow(*reaction, comp, T9, rho) > 0.0) {
                }
                creationRate += stoichiometry * engine.calculateMolarReactionFlow(*reaction, comp, T9, rho);
            }
        }
        return creationRate;
    }
}