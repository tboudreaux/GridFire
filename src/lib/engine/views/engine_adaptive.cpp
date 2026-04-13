#include "gridfire/engine/views/engine_adaptive.h"

#include <ranges>
#include <queue>
#include <algorithm>


#include "gridfire/types/types.h"
#include "gridfire/exceptions/error_engine.h"
#include "gridfire/utils/hashing.h"

#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/utils.h"
#include "gridfire/engine/scratchpads/engine_adaptive_scratchpad.h"

#include "quill/LogMacros.h"
#include "quill/Logger.h"

namespace gridfire::engine {
    using fourdst::atomic::Species;
    AdaptiveEngineView::AdaptiveEngineView(
        DynamicEngine &baseEngine
    ) :
    m_baseEngine(baseEngine) {}

    fourdst::composition::Composition AdaptiveEngineView::project(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        auto *state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        state->active_reactions.clear();
        state->active_species.clear();

        fourdst::composition::Composition baseUpdatedComposition = m_baseEngine.project(ctx, netIn);
        NetIn updatedNetIn = netIn;

        updatedNetIn.composition = baseUpdatedComposition;

        LOG_TRACE_L1(m_logger, "Updating AdaptiveEngineView with new network input...");

        auto [allFlows, composition] = calculateAllReactionFlows(ctx, updatedNetIn);

        double maxFlow = 0.0;

        for (const auto&[reactionPtr, flowRate]: allFlows) {
            if (flowRate > maxFlow) {
                maxFlow = flowRate;
            }
        }
        LOG_DEBUG(m_logger, "Maximum reaction flow rate in adaptive engine view: {:0.3E} [mol/s]", maxFlow);

        const std::unordered_set<Species> reachableSpecies = findReachableSpecies(ctx, updatedNetIn);
        LOG_DEBUG(m_logger, "Found {} reachable species in adaptive engine view.", reachableSpecies.size());

        const std::vector<const reaction::Reaction*> finalReactions = cullReactionsByFlow(ctx, allFlows, reachableSpecies, composition, maxFlow);

        finalizeActiveSet(ctx, finalReactions);

        auto [rescuedReactions, rescuedSpecies] = rescueEdgeSpeciesDestructionChannel(
            ctx,
            composition,
            netIn.temperature/1e9,
            netIn.density
        );

        for (const auto& reactionPtr : rescuedReactions) {
            state->active_reactions.add_reaction(*reactionPtr);
        }

        for (const auto& species : rescuedSpecies) {
            if (!std::ranges::contains(state->active_species, species) && m_baseEngine.getSpeciesStatus(ctx, species) == SpeciesStatus::ACTIVE) {
                state->active_species.push_back(species);
            }
        }

        LOG_INFO(m_logger, "AdaptiveEngineView updated successfully with {} active species and {} active reactions.", state->active_species.size(), state->active_reactions.size());

        return updatedNetIn.composition;
    }

    const std::vector<Species> & AdaptiveEngineView::getNetworkSpecies(scratch::StateBlob& ctx) const {
        return scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx)->active_species;
    }

    std::expected<StepDerivatives<double>, EngineStatus> AdaptiveEngineView::calculateRHSAndEnergy(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho, bool trust
    ) const {
        LOG_TRACE_L2(m_logger, "Calculating RHS and Energy in AdaptiveEngineView at T9 = {}, rho = {}.", T9, rho);

        const fourdst::composition::Composition collectedComp = collectComposition(ctx, comp, T9, rho);

        auto result = m_baseEngine.calculateRHSAndEnergy(ctx, collectedComp, T9, rho, true);
        LOG_TRACE_L2(m_logger, "Base engine calculation of RHS and Energy complete.");

        if (!result) {
            LOG_TRACE_L2(m_logger, "Base engine returned stale error during RHS and Energy calculation.");
            return std::unexpected{result.error()};
        }

        return result.value();
    }

    EnergyDerivatives AdaptiveEngineView::calculateEpsDerivatives(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateEpsDerivatives(ctx, comp, T9, rho);
    }

    NetworkJacobian AdaptiveEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const auto *state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        return generateJacobianMatrix(ctx, comp, T9, rho, state->active_species);
    }

    NetworkJacobian AdaptiveEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::vector<Species> &activeSpecies
    ) const {
        return m_baseEngine.generateJacobianMatrix(ctx, comp, T9, rho, activeSpecies);

    }

    NetworkJacobian AdaptiveEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const SparsityPattern &sparsityPattern
    ) const {
        return m_baseEngine.generateJacobianMatrix(ctx, comp, T9, rho, sparsityPattern);
    }

    double AdaptiveEngineView::calculateMolarReactionFlow(
        scratch::StateBlob& ctx,
        const reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const auto *state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        if (!state->active_reactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the adaptive engine view.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }

        return m_baseEngine.calculateMolarReactionFlow(ctx, reaction, comp, T9, rho);
    }

    const reaction::ReactionSet & AdaptiveEngineView::getNetworkReactions(
        scratch::StateBlob& ctx
    ) const {
        return scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx) -> active_reactions;
    }

    reaction::ReactionSet AdaptiveEngineView::getInactiveNetworkReactions(scratch::StateBlob &ctx) const {
        const reaction::ReactionSet& baseEngineReactions = m_baseEngine.getNetworkReactions(ctx);
        const reaction::ReactionSet baseEngineInactiveReactions = m_baseEngine.getInactiveNetworkReactions(ctx);


        reaction::ReactionSet inactiveReactions = baseEngineInactiveReactions;
        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        const reaction::ReactionSet& activeReactions = state->active_reactions;

        for (const auto& active_reaction : baseEngineReactions) {
            if (!inactiveReactions.contains(*active_reaction) && !activeReactions.contains(*active_reaction)) {
                inactiveReactions.add_reaction(*active_reaction);
            }
        }

        return inactiveReactions;

    }

    double AdaptiveEngineView::getInactiveReactionMolarReactionFlow(
        scratch::StateBlob &ctx,
        const reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        return m_baseEngine.calculateMolarReactionFlow(ctx, reaction, comp, T9, rho);
    }

    std::expected<std::unordered_map<Species, double>, EngineStatus> AdaptiveEngineView::getSpeciesTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {

        const auto result = m_baseEngine.getSpeciesTimescales(ctx, comp, T9, rho);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        const std::unordered_map<Species, double>& fullTimescales = result.value();


        std::unordered_map<Species, double> culledTimescales;
        culledTimescales.reserve(state->active_species.size());
        for (const auto& active_species : state->active_species) {
            if (fullTimescales.contains(active_species)) {
                culledTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return culledTimescales;

    }

    std::expected<std::unordered_map<Species, double>, EngineStatus> AdaptiveEngineView::getSpeciesDestructionTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {

        const auto result = m_baseEngine.getSpeciesDestructionTimescales(ctx, comp, T9, rho);
        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        const std::unordered_map<Species, double>& destructionTimescales = result.value();

        std::unordered_map<Species, double> culledTimescales;
        culledTimescales.reserve(state->active_species.size());
        for (const auto& active_species : state->active_species) {
            if (destructionTimescales.contains(active_species)) {
                culledTimescales[active_species] = destructionTimescales.at(active_species);
            }
        }
        return culledTimescales;
    }

    screening::ScreeningType AdaptiveEngineView::getScreeningModel(
        scratch::StateBlob& ctx
    ) const {
        return m_baseEngine.getScreeningModel(ctx);
    }

    PrimingReport AdaptiveEngineView::primeEngine(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        return m_baseEngine.primeEngine(ctx, netIn);
    }

    fourdst::composition::Composition AdaptiveEngineView::collectComposition(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        fourdst::composition::Composition result = m_baseEngine.collectComposition(ctx, comp, T9, rho);

        for (const auto& species : state->active_species) {
            if (!result.contains(species)) {
                result.registerSpecies(species);
            }
        }

        return result;
    }

    SpeciesStatus AdaptiveEngineView::getSpeciesStatus(
        scratch::StateBlob& ctx,
        const Species &species
    ) const {
        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        const SpeciesStatus status = m_baseEngine.getSpeciesStatus(ctx, species);
        if (status == SpeciesStatus::ACTIVE && std::ranges::find(state->active_species, species) == state->active_species.end()) {
            return SpeciesStatus::INACTIVE_FLOW;
        }
        return status;
    }

    std::optional<StepDerivatives<double>> AdaptiveEngineView::getMostRecentRHSCalculation(
        scratch::StateBlob &ctx
    ) const {
        return m_baseEngine.getMostRecentRHSCalculation(ctx);
    }

    std::unique_ptr<scratch::StateBlob> AdaptiveEngineView::constructStateBlob(const scratch::StateBlob *blob) const {
        std::unique_ptr<scratch::StateBlob> i_blob;
        if (blob) {
            i_blob = blob->clone_structure();
        } else {
            i_blob = std::make_unique<scratch::StateBlob>();
        }
        i_blob->enroll<scratch::AdaptiveEngineViewScratchPad>();
        auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, false>(*i_blob);
        state->initialize(*this);
        return i_blob;
    }

    size_t AdaptiveEngineView::getSpeciesIndex(
        scratch::StateBlob& ctx,
        const Species &species
    ) const {
        const auto *state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        const auto it = std::ranges::find(state->active_species, species);
        if (it != state->active_species.end()) {
            return static_cast<int>(std::distance(state->active_species.begin(), it));
        } else {
            LOG_ERROR(m_logger, "Species '{}' not found in active species list.", species.name());
            m_logger->flush_log();
            throw std::runtime_error("Species not found in active species list: " + std::string(species.name()));
        }
    }

    std::pair<std::vector<AdaptiveEngineView::ReactionFlow>, fourdst::composition::Composition> AdaptiveEngineView::calculateAllReactionFlows(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        const auto& fullSpeciesList = m_baseEngine.getNetworkSpecies(ctx);
        fourdst::composition::Composition composition = netIn.composition;

        for (const auto& species: fullSpeciesList) {
            if (!netIn.composition.contains(species)) {
                LOG_TRACE_L2(m_logger, "Species '{}' not found in composition. Registering", species.name());
                composition.registerSpecies(species);
            }
        }

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        std::vector<ReactionFlow> reactionFlows;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions(ctx);
        reactionFlows.reserve(fullReactionSet.size());
        for (const auto& reaction : fullReactionSet) {
            const double flow = m_baseEngine.calculateMolarReactionFlow(ctx, *reaction, composition, T9, rho);
            reactionFlows.push_back({reaction.get(), flow});
            LOG_TRACE_L3(m_logger, "Reaction '{}' has flow rate: {:0.3E} [mol/s/g]", reaction->id(), flow);
        }
        return {reactionFlows, composition};
    }

    std::unordered_set<Species> AdaptiveEngineView::findReachableSpecies(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        std::unordered_set<Species> reachable;
        std::queue<Species> to_vist;

        constexpr double ABUNDANCE_FLOOR = 1e-12; // Abundance floor for a species to be considered part of the initial fuel
        for (const auto& species: m_baseEngine.getNetworkSpecies(ctx)) {
            if (netIn.composition.contains(species) && netIn.composition.getMassFraction(std::string(species.name())) > ABUNDANCE_FLOOR) {
                if (!reachable.contains(species)) {
                    to_vist.push(species);
                    reachable.insert(species);
                    LOG_TRACE_L2(m_logger, "Network Connectivity Analysis: Species {:5} is part of the initial fuel", species.name());
                }
            }
        }

        bool new_species_found_in_pass = true;
        while (new_species_found_in_pass) {
            new_species_found_in_pass = false;
            for (const auto& reaction: m_baseEngine.getNetworkReactions(ctx)) {
                bool all_reactants_reachable = true;
                for (const auto& reactant: reaction->reactants()) {
                    if (!reachable.contains(reactant)) {
                        all_reactants_reachable = false;
                        break;
                    }
                }
                if (all_reactants_reachable) {
                    for (const auto& product: reaction->products()) {
                        if (!reachable.contains(product)) {
                            reachable.insert(product);
                            new_species_found_in_pass = true;
                            LOG_TRACE_L2(m_logger, "Network Connectivity Analysis: Species {:5} is reachable via reaction {:20}", product.name(), reaction->id());
                        }
                    }
                }
            }
        }

        return reachable;
    }

    std::vector<const reaction::Reaction *> AdaptiveEngineView::cullReactionsByFlow(
        scratch::StateBlob& ctx,
        const std::vector<ReactionFlow> &allFlows,
        const std::unordered_set<fourdst::atomic::Species> &reachableSpecies,
        const fourdst::composition::Composition &comp,
        const double maxFlow
    ) const {
        LOG_TRACE_L1(m_logger, "Culling reactions based on flow rates...");

        const auto relative_culling_threshold = m_config->engine.views.adaptiveEngineView.relativeCullingThreshold;

        double absoluteCullingThreshold = relative_culling_threshold * maxFlow;
        LOG_DEBUG(m_logger, "Relative culling threshold: {:7.3E} ({:7.3E})", relative_culling_threshold, absoluteCullingThreshold);
        std::vector<const reaction::Reaction*> culledReactions;
        for (const auto& [reactionPtr, flowRate]: allFlows) {
            bool keepReaction = false;
            if (flowRate > absoluteCullingThreshold) {
                LOG_TRACE_L2(m_logger, "Maintaining reaction '{:20}' with relative (abs) flow rate: {:7.3E} ({:7.3E} [mol/s])", reactionPtr->id(), flowRate/maxFlow, flowRate);
                keepReaction = true;
            } else {
                bool zero_flow_due_to_reachable_reactants = false;
                if (flowRate < 1e-99 && flowRate > 0.0) {
                    for (const auto& reactant: reactionPtr->reactants()) {
                        if (comp.getMolarAbundance(reactant) < 1e-99 && reachableSpecies.contains(reactant)) {
                            LOG_TRACE_L1(m_logger, "Maintaining reaction {:20} with low flow ({:7.3E} [mol/s/g]) due to reachable reactant '{:6}'.", reactionPtr->id(), flowRate, reactant.name());
                            zero_flow_due_to_reachable_reactants = true;
                            break;
                        }
                    }
                }
                if (zero_flow_due_to_reachable_reactants) {
                    keepReaction = true;
                }
            }
            if (keepReaction) {
                culledReactions.push_back(reactionPtr);
            } else {
                LOG_TRACE_L1(m_logger, "Culling     reaction '{:20}' due to low flow rate or lack of connectivity.", reactionPtr->id());
            }
        }
        LOG_DEBUG(m_logger, "Selected {:5} (total: {:5}, culled: {:5}) reactions based on flow rates.", culledReactions.size(), allFlows.size(), allFlows.size() - culledReactions.size());
        return culledReactions;
    }

    AdaptiveEngineView::RescueSet AdaptiveEngineView::rescueEdgeSpeciesDestructionChannel(
        scratch::StateBlob& ctx,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho
    ) const {
        const auto result = m_baseEngine.getSpeciesTimescales(ctx, comp, T9, rho);
        if (!result) {
            LOG_CRITICAL(m_logger, "Failed to get species timescales due to base engine failure");
            m_logger->flush_log();
            throw exceptions::EngineError("Failed to get species timescales due base engine failure");
        }

        const auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        std::unordered_map<Species, double> timescales = result.value();
        std::set<Species> onlyProducedSpecies;
        for (const auto& reaction : state->active_reactions) {
            const std::vector<Species>& products = reaction->products();
            onlyProducedSpecies.insert(products.begin(), products.end());
        }

        // Remove species that are consumed by any one of the active reactions.
        std::erase_if(
            onlyProducedSpecies,
            [&](const Species &species) {
                for (const auto& reaction : state->active_reactions) {
                    if (reaction->contains_reactant(species)) {
                        return true; // If any active reaction consumes the species then erase it from the set.
                    }
                }
                return false;
            }
        );

        // Remove species that have a non-zero timescale (these are expected to be +inf as they should be the equilibrium species if using with a MultiscalePartitioningEngineView)
        std::erase_if(
            onlyProducedSpecies,
            [&](const Species &species) {
                return std::isinf(timescales.at(species));
            }
        );

        std::unordered_map<Species, const reaction::Reaction*> reactionsToRescue;
        for (const auto& species : onlyProducedSpecies) {
            double maxSpeciesConsumptionRate = 0.0;
            for (const auto& reaction : m_baseEngine.getNetworkReactions(ctx)) {
                const bool speciesToCheckIsConsumed = reaction->contains_reactant(species);
                if (!speciesToCheckIsConsumed) {
                    continue; // If the species is not consumed by this reaction, skip it.
                }
                bool allOtherReactantsAreAvailable = true;
                for (const auto& reactant : reaction->reactants()) {
                    const bool reactantIsAvailable = std::ranges::contains(state->active_species, reactant);
                    if (!reactantIsAvailable && reactant != species) {
                        allOtherReactantsAreAvailable = false;
                    }
                }
                if (allOtherReactantsAreAvailable) {
                    std::vector<double> Y = comp.getMolarAbundanceVector();

                    const double Ye = comp.getElectronAbundance();

                    std::unordered_map<Species, double> speciesMassMap;
                    for (const auto &sp: comp | std::views::keys) {
                        speciesMassMap[sp] = sp.mass();
                    }
                    std::unordered_map<size_t, Species> speciesIndexMap;
                    for (const auto& sp: comp | std::views::keys) {
                        size_t distance = std::distance(speciesMassMap.begin(), speciesMassMap.find(sp));
                        speciesIndexMap.emplace(distance, sp);
                    }
                    double rate = reaction->calculate_rate(T9, rho, Ye, 0.0, Y, speciesIndexMap);
                    if (rate > maxSpeciesConsumptionRate) {
                        maxSpeciesConsumptionRate = rate;
                        reactionsToRescue[species] = reaction.get();
                    }
                }
            }
        }
        LOG_TRACE_L1(
            m_logger,
            "Rescuing {} {}reactions",
            reactionsToRescue.size(),
            [&]() -> std::string {
                std::ostringstream ss;
                if (reactionsToRescue.empty()) {
                    return "";
                }
                int count = 0;
                ss << "(";
                for (const auto &reaction : reactionsToRescue | std::views::values) {
                    ss << reaction->id();
                    if (count < reactionsToRescue.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }
                ss << ") ";
                return ss.str();
            }()
        );

        LOG_TRACE_L1(
            m_logger,
            "Timescale adjustments due to reaction rescue: {}",
            [&]() -> std::string {
                std::stringstream ss;
                if (reactionsToRescue.empty()) {
                    return "No reactions rescued...";
                }
                int count = 0;
                for (const auto& [species, reaction] : reactionsToRescue) {
                    ss << "(Species: " << species.name() << " started with a timescale of " << timescales.at(species);
                    ss << ", rescued by reaction: " << reaction->id();
                    ss << " whose product timescales are --- [";
                    int icount = 0;
                    for (const auto& product : reaction->products()) {
                        ss << product.name() << ": " << timescales.at(product);
                        if (icount < reaction->products().size() - 1) {
                            ss << ", ";
                        }
                        icount++;
                    }
                    ss << "])";

                    if (count < reactionsToRescue.size() - 1) {
                        ss << ", ";
                    }
                    count++;
                }

                return ss.str();
            }()
        );

        std::unordered_set<const reaction::Reaction*> newReactions;
        std::unordered_set<Species> newSpecies;

        for (const auto &reactionPtr: reactionsToRescue | std::views::values) {
            newReactions.insert(reactionPtr);
            for (const auto& product : reactionPtr->products()) {
                newSpecies.insert(product);
            }
        }
        return {std::move(newReactions), std::move(newSpecies)};
    }

    void AdaptiveEngineView::finalizeActiveSet(
        scratch::StateBlob& ctx,
        const std::vector<const reaction::Reaction *> &finalReactions
    ) const {
        auto* state = scratch::get_state<scratch::AdaptiveEngineViewScratchPad, true>(ctx);
        std::unordered_set<Species>finalSpeciesSet;
        state->active_reactions.clear();
        for (const auto* reactionPtr: finalReactions) {
            state->active_reactions.add_reaction(*reactionPtr);
            for (const auto& reactant : reactionPtr->reactants()) {
                const SpeciesStatus reactantStatus = m_baseEngine.getSpeciesStatus(ctx, reactant);
                if (!finalSpeciesSet.contains(reactant) && (reactantStatus == SpeciesStatus::ACTIVE || reactantStatus == SpeciesStatus::EQUILIBRIUM)) {
                    LOG_TRACE_L3(m_logger, "Adding reactant '{}' to active species set through reaction {}.", reactant.name(), reactionPtr->id());
                    finalSpeciesSet.insert(reactant);
                }
            }
            for (const auto& product : reactionPtr->products()) {
                const SpeciesStatus productStatus = m_baseEngine.getSpeciesStatus(ctx, product);
                if (!finalSpeciesSet.contains(product) && (productStatus == SpeciesStatus::ACTIVE || productStatus == SpeciesStatus::EQUILIBRIUM)) {
                    LOG_TRACE_L3(m_logger, "Adding product '{}' to active species set through reaction {}.", product.name(), reactionPtr->id());
                    finalSpeciesSet.insert(product);
                }
            }
        }
        state->active_species.clear();
        state->active_species = std::vector<Species>(finalSpeciesSet.begin(), finalSpeciesSet.end());
        std::ranges::sort(
            state->active_species,
            [](const Species &a, const Species &b) { return a.mass() < b.mass(); }
        );
    }
}

