#include "gridfire/engine/procedures/priming.h"
#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/engine/procedures/construction.h"
#include "gridfire/solver/solver.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

namespace {
    // Create a dummy wrapper Composition to measure the unrestricted flow rate of species
    class UnrestrictedComposition final : public fourdst::composition::Composition {
    private:
        const fourdst::atomic::Species& m_unrestrictedSpecies;
    public:
        explicit UnrestrictedComposition(const Composition& baseComposition, const fourdst::atomic::Species& unrestrictedSpecies) :
            Composition(baseComposition),
            m_unrestrictedSpecies(unrestrictedSpecies) {}

        double getMolarAbundance(const fourdst::atomic::Species &species) const override {
            if (species == m_unrestrictedSpecies) {
                return 1.0; // Set to a high value to simulate unrestricted abundance
            }
            return Composition::getMolarAbundance(species);
        }

        double getMolarAbundance(const std::string &symbol) const override {
            if (symbol == m_unrestrictedSpecies.name()) {
                return 1.0; // Set to a high value to simulate unrestricted abundance
            }
            return Composition::getMolarAbundance(symbol);
        }
    };

    bool isReactionIgnorable(
        const gridfire::reaction::Reaction& reaction,
        const std::optional<std::vector<gridfire::reaction::ReactionType>>& reactionsTypesToIgnore
    ) {
        if (reactionsTypesToIgnore.has_value()) {
            for (const auto& type : reactionsTypesToIgnore.value()) {
                if (reaction.type() == type) {
                    return true;
                }
            }
        }
        return false;
    }
}

namespace gridfire {
    using fourdst::composition::Composition;
    using fourdst::atomic::Species;

    const reaction::Reaction* findDominantCreationChannel (
        const DynamicEngine& engine,
        const Species& species,
        const Composition &comp,
        const double T9,
        const double rho,
        const std::optional<std::vector<reaction::ReactionType>> &reactionsTypesToIgnore
    ) {
        const reaction::Reaction* dominateReaction = nullptr;
        double maxFlow = -1.0;
        for (const auto& reaction : engine.getNetworkReactions()) {
            if (isReactionIgnorable(*reaction, reactionsTypesToIgnore)) continue;

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
     * @param ignoredReactionTypes Types of reactions to ignore during priming (e.g., weak reactions).
     * @return PrimingReport encapsulating the results of the priming operation, including the new
     * robustly primed composition.
     */
    PrimingReport primeNetwork(
        const NetIn& netIn,
        DynamicEngine& engine,
        const std::optional<std::vector<reaction::ReactionType>>& ignoredReactionTypes
    ) {
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        // --- Initial Setup ---
        // Identify all species with zero initial mass fraction that need to be primed.
        std::vector<Species> speciesToPrime;
        for (const auto &entry: netIn.composition | std::views::values) {
            if (entry.mass_fraction() == 0.0) {
                speciesToPrime.push_back(entry.isotope());
            }
        }

        // sort primingSpecies by mass number, lightest to heaviest. This ensures we prime in a physically sensible order.
        std::ranges::sort(speciesToPrime, [](const Species& a, const Species& b) {
            return a.mass() < b.mass();
        });

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
            bool didFinalize = tempComp.finalize(true);
            if (!didFinalize) {
                LOG_ERROR(logger, "Failed to finalize temporary composition during priming.");
                report.success = false;
                report.status = PrimingReportStatus::FAILED_TO_FINALIZE_COMPOSITION;
                continue;
            }

            NetworkPrimingEngineView primer(primingSpecies, engine);
            if (primer.getNetworkReactions().size() == 0) {
                LOG_ERROR(logger, "No priming reactions found for species {}.", primingSpecies.name());
                report.success = false;
                report.status = PrimingReportStatus::FAILED_TO_FIND_PRIMING_REACTIONS;
                continue;
            }

            const double destructionRateConstant = calculateDestructionRateConstant(
                primer,
                primingSpecies,
                tempComp,
                T9,
                rho,
                ignoredReactionTypes
            );

            if (destructionRateConstant > 1e-99) {
                const double creationRate = calculateCreationRate(
                    primer,
                    primingSpecies,
                    tempComp,
                    T9,
                    rho,
                    ignoredReactionTypes
                );
                double equilibriumMassFraction = (creationRate / destructionRateConstant) * primingSpecies.mass();

                // ReSharper disable once CppDFAUnusedValue
                if (std::isnan(equilibriumMassFraction)) equilibriumMassFraction = 0.0;

                if (const reaction::Reaction* dominantChannel = findDominantCreationChannel(
                    primer,
                    primingSpecies,
                    tempComp,
                    T9,
                    rho,
                    ignoredReactionTypes)
                ) {
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

  double calculateDestructionRateConstant(
      const DynamicEngine& engine,
      const Species& species,
      const Composition& composition,
      const double T9,
      const double rho,
      const std::optional<std::vector<reaction::ReactionType>> &reactionTypesToIgnore
  ) {
        const UnrestrictedComposition unrestrictedComp(composition, species); // Create a composition that simulates an enormous source abundance of the target species (getMolarAbundance(species) always returns 1.0)
        double destructionRateConstant = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            if (isReactionIgnorable(*reaction, reactionTypesToIgnore)) continue;

            const int stoichiometry = reaction->stoichiometry(species);
            if (stoichiometry < 0) {
                destructionRateConstant += std::abs(stoichiometry) * engine.calculateMolarReactionFlow(*reaction, unrestrictedComp, T9, rho);
            }
        }
        return destructionRateConstant;
    }

    double calculateCreationRate(
        const DynamicEngine& engine,
        const Species& species,
        const Composition& composition,
        const double T9,
        const double rho,
        const std::optional<std::vector<reaction::ReactionType>> &reactionTypesToIgnore
    ) {
        double creationRate = 0.0;
        for (const auto& reaction: engine.getNetworkReactions()) {
            if (isReactionIgnorable(*reaction, reactionTypesToIgnore)) continue;

            const int stoichiometry = reaction->stoichiometry(species);
            if (stoichiometry > 0) {
                creationRate += stoichiometry * engine.calculateMolarReactionFlow(*reaction, composition, T9, rho);
            }
        }
        return creationRate;
    }
}