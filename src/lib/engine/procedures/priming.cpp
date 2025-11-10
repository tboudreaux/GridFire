#include "gridfire/engine/procedures/priming.h"

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"
#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/engine/procedures/construction.h"
#include "gridfire/solver/solver.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "quill/Logger.h"
#include "quill/LogMacros.h"

namespace {
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
        GraphEngine& engine,
        const std::optional<std::vector<reaction::ReactionType>>& ignoredReactionTypes
    ) {
        auto logger = LogManager::getInstance().getLogger("log");

        // --- Initial Setup ---
        // Identify all species with zero initial abundance that need to be primed.
        std::vector<Species> speciesToPrime;
        for (const auto &[sp, y]: netIn.composition) {
            if (y == 0.0) {
                speciesToPrime.push_back(sp);
            }
        }

        // Sort priming species by mass number, lightest to heaviest.
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

        // Build full set of species
        std::set<Species> allSpecies;
        for (const auto &sp: netIn.composition | std::views::keys) {
            allSpecies.insert(sp);
        }
        for (const auto& sp : speciesToPrime) {
            allSpecies.insert(sp);
        }

        // Rebuild the engine with the full network to ensure all possible creation channels are available.
        engine.rebuild(netIn.composition, NetworkBuildDepth::Full);

        // Initialize mutable molar abundances for all species
        std::map<Species, double> molarAbundances;
        for (const auto& sp : allSpecies) {
            molarAbundances[sp] = netIn.composition.contains(sp) ? netIn.composition.getMolarAbundance(sp) : 0.0;
        }

        // --- Prime Each Species ---
        // Since molar abundances are independent, we can directly calculate and apply changes
        std::unordered_map<Species, double> totalMolarAbundanceChanges;

        for (const auto& primingSpecies : speciesToPrime) {
            // Create a temporary composition reflecting the current state for rate calculations.
            Composition tempComp(allSpecies);
            for (const auto& [sp, abundance] : molarAbundances) {
                tempComp.setMolarAbundance(sp, abundance);
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

            double equilibriumMolarAbundance = 0.0;
            std::vector<Species> reactants;

            if (destructionRateConstant > 1e-99) {
                const double creationRate = calculateCreationRate(
                    primer,
                    primingSpecies,
                    tempComp,
                    T9,
                    rho,
                    ignoredReactionTypes
                );

                // Equilibrium: creation rate = destruction rate constant * molar abundance
                equilibriumMolarAbundance = creationRate / destructionRateConstant;

                if (std::isnan(equilibriumMolarAbundance)) {
                    equilibriumMolarAbundance = 0.0;
                }

                if (const reaction::Reaction* dominantChannel = findDominantCreationChannel(
                    primer,
                    primingSpecies,
                    tempComp,
                    T9,
                    rho,
                    ignoredReactionTypes)
                ) {
                    reactants = dominantChannel->reactants();
                } else {
                    LOG_TRACE_L1(logger, "Failed to find dominant creation channel for {}.", primingSpecies.name());
                    report.status = PrimingReportStatus::FAILED_TO_FIND_CREATION_CHANNEL;
                    reactants.clear(); // Use fallback
                }
            } else {
                LOG_TRACE_L2(logger, "No destruction channel found for {}. Using fallback abundance.", primingSpecies.name());
                // For species with no destruction, use a tiny fallback abundance
                equilibriumMolarAbundance = 1e-40;
            }

            // Add the equilibrium molar abundance to the primed species
            molarAbundances.at(primingSpecies) += equilibriumMolarAbundance;
            totalMolarAbundanceChanges[primingSpecies] += equilibriumMolarAbundance;

            // Subtract from reactants proportionally to their stoichiometry
            if (!reactants.empty()) {
                const double totalStoichiometry = static_cast<double>(reactants.size());
                const double abundancePerReactant = equilibriumMolarAbundance / totalStoichiometry;

                for (const auto& reactant : reactants) {
                    LOG_TRACE_L1(logger, "Transferring {:.4e} molar abundance from {} to {} during priming.",
                        abundancePerReactant, reactant.name(), primingSpecies.name());

                    if (!molarAbundances.contains(reactant)) {
                        continue;
                    }
                    molarAbundances.at(reactant) -= abundancePerReactant;
                    totalMolarAbundanceChanges[reactant] -= abundancePerReactant;

                    // Ensure non-negative abundances
                    if (molarAbundances.at(reactant) < 0.0) {
                        LOG_WARNING(logger, "Species {} went negative during priming. Clamping to zero.", reactant.name());
                        totalMolarAbundanceChanges[reactant] += molarAbundances.at(reactant); // Adjust change to reflect clamp
                        molarAbundances.at(reactant) = 0.0;
                    }
                }
            }
        }

        // --- Final Composition Construction ---
        std::vector<Species> final_species;
        std::vector<double> final_molar_abundances;

        for (const auto& [species, abundance] : molarAbundances) {
            final_species.emplace_back(species);
            final_molar_abundances.push_back(std::max(0.0, abundance));

            LOG_TRACE_L1(logger, "After priming, species {} has molar abundance {:.4e} (had {:0.4e} prior to priming).",
                species.name(),
                std::max(0.0, abundance),
                netIn.composition.getMolarAbundance(species));
        }

        Composition primedComposition(final_species, final_molar_abundances);

        report.primedComposition = primedComposition;

        // Convert molar abundance changes to mass fraction changes for reporting
        for (const auto& [species, molarChange] : totalMolarAbundanceChanges) {
            double massFractionChange = molarChange * species.mass();
            report.massFractionChanges.emplace_back(species, massFractionChange);
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
        Composition unrestrictedComp(composition);
        unrestrictedComp.setMolarAbundance(species, 1.0);
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