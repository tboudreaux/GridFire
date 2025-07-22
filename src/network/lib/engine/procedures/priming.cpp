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
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) {
        const reaction::Reaction* dominateReaction = nullptr;
        double maxFlow = -1.0;
        for (const auto& reaction : engine.getNetworkReactions()) {
            if (reaction.contains(species) && reaction.stoichiometry(species) > 0) {
                const double flow = engine.calculateMolarReactionFlow(reaction, Y, T9, rho);
                if (flow > maxFlow) {
                    maxFlow = flow;
                    dominateReaction = &reaction;
                }
            }
        }
        return dominateReaction;
    }


    PrimingReport primeNetwork(const NetIn& netIn, DynamicEngine& engine) {
        auto logger = LogManager::getInstance().getLogger("log");

        std::vector<Species> speciesToPrime;
        for (const auto &entry: netIn.composition | std::views::values) {
            if (entry.mass_fraction() == 0.0) {
                speciesToPrime.push_back(entry.isotope());
            }
        }
        LOG_DEBUG(logger, "Priming {} species in the network.", speciesToPrime.size());

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

        // --- 1: pack composition into internal map ---
        std::unordered_map<Species, double> currentMassFractions;
        for (const auto& entry : netIn.composition | std::views::values) {
            currentMassFractions[entry.isotope()] = entry.mass_fraction();
        }
        for (const auto& entry : speciesToPrime) {
            currentMassFractions[entry] = 0.0; // Initialize priming species with 0 mass fraction
        }

        std::unordered_map<Species, double> totalMassFractionChanges;

        engine.rebuild(netIn.composition, NetworkBuildDepth::Full);

        for (const auto& primingSpecies : speciesToPrime) {
            LOG_TRACE_L3(logger, "Priming species: {}", primingSpecies.name());

            // Create a temporary composition from the current internal state for the primer
            Composition tempComp;
            for(const auto& [sp, mf] : currentMassFractions) {
                tempComp.registerSymbol(std::string(sp.name()));
                if (mf < 0.0 && std::abs(mf) < 1e-16) {
                    tempComp.setMassFraction(sp, 0.0); // Avoid negative mass fractions
                } else {
                    tempComp.setMassFraction(sp, mf);
                }
            }
            tempComp.finalize(true); // Finalize with normalization

            NetIn tempNetIn = netIn;
            tempNetIn.composition = tempComp;

            NetworkPrimingEngineView primer(primingSpecies, engine);

            if (primer.getNetworkReactions().size() == 0) {
                LOG_ERROR(logger, "No priming reactions found for species {}.", primingSpecies.name());
                report.success = false;
                report.status = PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS;
                continue;
            }

            const auto Y = primer.mapNetInToMolarAbundanceVector(tempNetIn);
            const double destructionRateConstant = calculateDestructionRateConstant(primer, primingSpecies, Y, T9, rho);

            double equilibriumMassFraction = 0.0;

            if (destructionRateConstant > 1e-99) {
                const double creationRate = calculateCreationRate(primer, primingSpecies, Y, T9, rho);
                equilibriumMassFraction = (creationRate / destructionRateConstant) * primingSpecies.mass();
                if (std::isnan(equilibriumMassFraction)) {
                    LOG_WARNING(logger, "Equilibrium mass fraction for {} is NaN. Setting to 0.0. This is likely not an issue. It probably originates from all reactions leading to creation and destruction being frozen out. In that case 0.0 should be a good approximation. Hint: This happens often when the network temperature is very the low. ", primingSpecies.name());
                    equilibriumMassFraction = 0.0;
                }
                LOG_TRACE_L3(logger, "Found equilibrium for {}: X_eq = {:.4e}", primingSpecies.name(), equilibriumMassFraction);

                const reaction::Reaction* dominantChannel = findDominantCreationChannel(primer, primingSpecies, Y, T9, rho);

                if (dominantChannel) {
                    LOG_TRACE_L3(logger, "Dominant creation channel for {}: {}", primingSpecies.name(), dominantChannel->peName());

                    double totalReactantMass = 0.0;
                    for (const auto& reactant : dominantChannel->reactants()) {
                        totalReactantMass += reactant.mass();
                    }

                    double scalingFactor = 1.0;
                    for (const auto& reactant : dominantChannel->reactants()) {
                        const double massToSubtract = equilibriumMassFraction * (reactant.mass() / totalReactantMass);
                        double availableMass = 0.0;
                        if (currentMassFractions.contains(reactant)) {
                            availableMass = currentMassFractions.at(reactant);
                        }
                        if (massToSubtract > availableMass && availableMass > 0) {
                            scalingFactor = std::min(scalingFactor, availableMass / massToSubtract);
                        }
                    }

                    if (scalingFactor < 1.0) {
                        LOG_WARNING(logger, "Priming for {} was limited by reactant availability. Scaling transfer by {:.4e}", primingSpecies.name(), scalingFactor);
                        equilibriumMassFraction *= scalingFactor;
                    }

                    // Update the internal mass fraction map and accumulate total changes
                    totalMassFractionChanges[primingSpecies] += equilibriumMassFraction;
                    currentMassFractions[primingSpecies] += equilibriumMassFraction;

                    for (const auto& reactant : dominantChannel->reactants()) {
                        const double massToSubtract = equilibriumMassFraction * (reactant.mass() / totalReactantMass);
                        totalMassFractionChanges[reactant] -= massToSubtract;
                        currentMassFractions[reactant] -= massToSubtract;
                    }
                } else {
                     LOG_ERROR(logger, "Failed to find dominant creation channel for {}.", primingSpecies.name());
                     report.status = PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL;
                     totalMassFractionChanges[primingSpecies] += 1e-40;
                     currentMassFractions[primingSpecies] += 1e-40;
                }
            } else {
                LOG_WARNING(logger, "No destruction channel found for {}. Using fallback abundance.", primingSpecies.name());
                totalMassFractionChanges[primingSpecies] += 1e-40;
                currentMassFractions[primingSpecies] += 1e-40;
                report.status = PrimingReportStatus::BASE_NETWORK_TOO_SHALLOW;
            }
        }

        // --- Final Composition Construction ---
        std::vector<std::string> final_symbols;
        std::vector<double> final_mass_fractions;
        for(const auto& [species, mass_fraction] : currentMassFractions) {
            final_symbols.push_back(std::string(species.name()));
            if (mass_fraction < 0.0 && std::abs(mass_fraction) < 1e-16) {
                final_mass_fractions.push_back(0.0); // Avoid negative mass fractions
            } else {
                final_mass_fractions.push_back(mass_fraction);
            }
        }

        // Create the final composition object from the pre-normalized mass fractions
        Composition primedComposition(final_symbols, final_mass_fractions, true);

        report.primedComposition = primedComposition;
        for (const auto& [species, change] : totalMassFractionChanges) {
            report.massFractionChanges.emplace_back(species, change);
        }

        engine.setNetworkReactions(initialReactionSet);
        return report;
    }

  double calculateDestructionRateConstant(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) {
        const int speciesIndex = engine.getSpeciesIndex(species);
        std::vector<double> Y_scaled(Y.begin(), Y.end());
        Y_scaled[speciesIndex] = 1.0; // Set the abundance of the species to 1.0 for rate constant calculation
        double destructionRateConstant = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            if (reaction.contains_reactant(species)) {
                const int stoichiometry = reaction.stoichiometry(species);
                destructionRateConstant += std::abs(stoichiometry) * engine.calculateMolarReactionFlow(reaction, Y_scaled, T9, rho);
            }
        }
        return destructionRateConstant;
    }

    double calculateCreationRate(
        const DynamicEngine& engine,
        const fourdst::atomic::Species& species,
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) {
        double creationRate = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            const int stoichiometry = reaction.stoichiometry(species);
            if (stoichiometry > 0) {
                creationRate += stoichiometry * engine.calculateMolarReactionFlow(reaction, Y, T9, rho);
            }
        }
        return creationRate;
    }
}