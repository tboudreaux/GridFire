#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/engine/engine_graph.h"

#include "fourdst/atomic/species.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/composition/decorators/composition_masked.h"

#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/engine_defined_scratchpad.h"
#include "gridfire/engine/scratchpads/utils.h"

#include "quill/LogMacros.h"

#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace gridfire::engine {
    using fourdst::atomic::Species;

    DefinedEngineView::DefinedEngineView(
        const std::vector<std::string>& peNames,
        GraphEngine& baseEngine
    ) :
    m_baseEngine(baseEngine) {
        // collect(peNames);
    }

    const DynamicEngine & DefinedEngineView::getBaseEngine() const {
        return m_baseEngine;
    }

    const std::vector<Species> & DefinedEngineView::getNetworkSpecies(
         scratch::StateBlob& ctx
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        if (state->active_species_vector_cache.has_value()) {
            return state->active_species_vector_cache.value();
        }
        state->active_species_vector_cache = std::vector<Species>(state->active_species.begin(), state->active_species.end());
        return state->active_species_vector_cache.value();
    }

    std::expected<StepDerivatives<double>, EngineStatus> DefinedEngineView::calculateRHSAndEnergy(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho, bool trust
    ) const {
        auto *state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);

        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());
        const auto result = m_baseEngine.calculateRHSAndEnergy(ctx, masked, T9, rho, state->active_reactions);

        if (!result) {
            return std::unexpected{result.error()};
        }

        return result.value();
    }

    EnergyDerivatives DefinedEngineView::calculateEpsDerivatives(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());

        return m_baseEngine.calculateEpsDerivatives(ctx, masked, T9, rho, state->active_reactions);
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);

        if (!state->active_species_vector_cache.has_value()) {
            state->active_species_vector_cache = std::vector<Species>(state->active_species.begin(), state->active_species.end());
        }
        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());
        return m_baseEngine.generateJacobianMatrix(ctx, masked, T9, rho, state->active_species_vector_cache.value());
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const std::vector<Species> &activeSpecies
    ) const {

        const std::set<fourdst::atomic::Species> activeSpeciesSet(
            activeSpecies.begin(),
            activeSpecies.end()
        );

        const fourdst::composition::MaskedComposition masked(comp, activeSpeciesSet | std::ranges::to<std::vector>());
        return m_baseEngine.generateJacobianMatrix(ctx, masked, T9, rho, activeSpecies);
    }

    NetworkJacobian DefinedEngineView::generateJacobianMatrix(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho,
        const SparsityPattern &sparsityPattern
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());
        return m_baseEngine.generateJacobianMatrix(ctx, masked, T9, rho, sparsityPattern);
    }

    double DefinedEngineView::calculateMolarReactionFlow(
        scratch::StateBlob& ctx,
        const reaction::Reaction &reaction,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);

        if (!state->active_reactions.contains(reaction)) {
            LOG_ERROR(m_logger, "Reaction '{}' is not part of the active reactions in the DefinedEngineView.", reaction.id());
            m_logger -> flush_log();
            throw std::runtime_error("Reaction not found in active reactions: " + std::string(reaction.id()));
        }

        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());
        return m_baseEngine.calculateMolarReactionFlow(ctx, reaction, masked, T9, rho);
    }

    const reaction::ReactionSet & DefinedEngineView::getNetworkReactions(
        scratch::StateBlob& ctx
    ) const {

        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        return state->active_reactions;
    }

    std::expected<std::unordered_map<Species, double>, EngineStatus> DefinedEngineView::getSpeciesTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        const fourdst::composition::MaskedComposition masked(comp, state->active_species | std::ranges::to<std::vector>());

        const auto result = m_baseEngine.getSpeciesTimescales(ctx, masked, T9, rho, state->active_reactions);
        if (!result) {
            return std::unexpected{result.error()};
        }
        const auto& fullTimescales = result.value();

        std::unordered_map<Species, double> definedTimescales;
        for (const auto& active_species : state->active_species) {
            if (fullTimescales.contains(active_species)) {
                definedTimescales[active_species] = fullTimescales.at(active_species);
            }
        }
        return definedTimescales;
    }

    std::expected<std::unordered_map<Species, double>, EngineStatus> DefinedEngineView::getSpeciesDestructionTimescales(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        const fourdst::composition::MaskedComposition masked(comp,  state->active_species| std::ranges::to<std::vector>());

        const auto result = m_baseEngine.getSpeciesDestructionTimescales(ctx, masked, T9, rho, state->active_reactions);

        if (!result) {
            return std::unexpected{result.error()};
        }

        const auto& destructionTimescales = result.value();

        std::unordered_map<Species, double> definedTimescales;
        for (const auto& active_species : state->active_species){
            if (destructionTimescales.contains(active_species)) {
                definedTimescales[active_species] = destructionTimescales.at(active_species);
            }
        }
        return definedTimescales;
    }

    fourdst::composition::Composition DefinedEngineView::project(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        return m_baseEngine.project(ctx, netIn);
    }

    screening::ScreeningType DefinedEngineView::getScreeningModel(
        scratch::StateBlob& ctx
    ) const {
        return m_baseEngine.getScreeningModel(ctx);
    }

    size_t DefinedEngineView::getSpeciesIndex(
        scratch::StateBlob& ctx,
        const Species &species
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);

        const auto it = std::ranges::find(state->active_species, species);
        if (it != state->active_species.end()) {
            return static_cast<int>(std::distance(state->active_species.begin(), it));
        } else {
            LOG_ERROR(m_logger, "Species '{}' not found in active species list.", species.name());
            m_logger->flush_log();
            throw std::runtime_error("Species not found in active species list: " + std::string(species.name()));
        }
    }

    PrimingReport DefinedEngineView::primeEngine(
        scratch::StateBlob& ctx,
        const NetIn &netIn
    ) const {
        return m_baseEngine.primeEngine(ctx, netIn);
    }

    fourdst::composition::Composition DefinedEngineView::collectComposition(
        scratch::StateBlob& ctx,
        const fourdst::composition::CompositionAbstract &comp,
        const double T9,
        const double rho
    ) const {
        fourdst::composition::Composition result = m_baseEngine.collectComposition(ctx, comp, T9, rho);
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);

        for (const auto& species : state->active_species) {
            if (!result.contains(species)) {
                result.registerSpecies(species);
            }
        }
        return result;
    }

    SpeciesStatus DefinedEngineView::getSpeciesStatus(
        scratch::StateBlob& ctx,
        const Species &species
    ) const {
        const auto *state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        const SpeciesStatus status = m_baseEngine.getSpeciesStatus(ctx, species);
        if (status == SpeciesStatus::ACTIVE && !state->active_species.contains(species)) {
            return SpeciesStatus::INACTIVE_FLOW;
        }
        return status;
    }

    std::optional<StepDerivatives<double>> DefinedEngineView::getMostRecentRHSCalculation(
        scratch::StateBlob &ctx
    ) const {
        return m_baseEngine.getMostRecentRHSCalculation(ctx);
    }

    std::vector<size_t> DefinedEngineView::constructSpeciesIndexMap(
        scratch::StateBlob& ctx
    ) const {
        LOG_TRACE_L3(m_logger, "Constructing species index map for DefinedEngineView...");
        std::unordered_map<Species, size_t> fullSpeciesReverseMap;
        const auto& fullSpeciesList = m_baseEngine.getNetworkSpecies(ctx);

        fullSpeciesReverseMap.reserve(fullSpeciesList.size());

        for (size_t i = 0; i < fullSpeciesList.size(); ++i) {
            fullSpeciesReverseMap[fullSpeciesList[i]] = i;
        }

        std::vector<size_t> speciesIndexMap;
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        speciesIndexMap.reserve(state->active_species.size());

        for (const auto& active_species : state->active_species) {
            auto it = fullSpeciesReverseMap.find(active_species);
            if (it != fullSpeciesReverseMap.end()) {
                speciesIndexMap.push_back(it->second);
            } else {
                LOG_ERROR(m_logger, "Species '{}' not found in full species map.", active_species.name());
                m_logger -> flush_log();
                throw std::runtime_error("Species not found in full species map: " + std::string(active_species.name()));
            }
        }
        LOG_TRACE_L3(m_logger, "Species index map constructed with {} entries.", speciesIndexMap.size());
        return speciesIndexMap;

    }

    std::vector<size_t> DefinedEngineView::constructReactionIndexMap(
        scratch::StateBlob& ctx
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        LOG_TRACE_L3(m_logger, "Constructing reaction index map for DefinedEngineView...");

        // --- Step 1: Create a reverse map using the reaction's unique ID as the key. ---
        std::unordered_map<std::string_view, size_t> fullReactionReverseMap;
        const auto& fullReactionSet = m_baseEngine.getNetworkReactions(ctx);
        fullReactionReverseMap.reserve(fullReactionSet.size());

        for (size_t i_full = 0; i_full < fullReactionSet.size(); ++i_full) {
            fullReactionReverseMap[fullReactionSet[i_full].id()] = i_full;
        }

        // --- Step 2: Build the final index map using the active reaction set. ---
        std::vector<size_t> reactionIndexMap;
        reactionIndexMap.reserve(state->active_reactions.size());

        for (const auto& active_reaction_ptr : state->active_reactions) {
            auto it = fullReactionReverseMap.find(active_reaction_ptr->id());

            if (it != fullReactionReverseMap.end()) {
                reactionIndexMap.push_back(it->second);
            } else {
                LOG_ERROR(m_logger, "Active reaction '{}' not found in base engine during reaction index map construction.", active_reaction_ptr->id());
                m_logger->flush_log();
                throw std::runtime_error("Mismatch between active reactions and base engine.");
            }
        }

        LOG_TRACE_L3(m_logger, "Reaction index map constructed with {} entries.", reactionIndexMap.size());
        return reactionIndexMap;
    }

    std::vector<double> DefinedEngineView::mapViewToFull(
        scratch::StateBlob& ctx,
        const std::vector<double>& culled
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        std::vector<double> full(m_baseEngine.getNetworkSpecies(ctx).size(), 0.0);
        for (size_t i_culled = 0; i_culled < culled.size(); ++i_culled) {
            const size_t i_full = state->species_index_map[i_culled];
            full[i_full] += culled[i_culled];
        }
        return full;
    }

    std::vector<double> DefinedEngineView::mapFullToView(
        scratch::StateBlob& ctx,
        const std::vector<double>& full
    ) {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        std::vector<double> culled(state->active_species.size(), 0.0);
        for (size_t i_culled = 0; i_culled < state->active_species.size(); ++i_culled) {
            const size_t i_full = state->species_index_map[i_culled];
            culled[i_culled] = full[i_full];
        }
        return culled;
    }

    size_t DefinedEngineView::mapViewToFullSpeciesIndex(
        scratch::StateBlob& ctx,
        size_t culledSpeciesIndex
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        if (culledSpeciesIndex >= state->species_index_map.size()) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for species index map of size {}.", culledSpeciesIndex, state->species_index_map.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledSpeciesIndex) + " is out of bounds for species index map of size " + std::to_string(state->species_index_map.size()) + ".");
        }
        return state->species_index_map[culledSpeciesIndex];
    }

    size_t DefinedEngineView::mapViewToFullReactionIndex(
        scratch::StateBlob& ctx,
        size_t culledReactionIndex
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        if (culledReactionIndex >= state->reaction_index_map.size()) {
            LOG_ERROR(m_logger, "Defined index {} is out of bounds for reaction index map of size {}.", culledReactionIndex, state->reaction_index_map.size());
            m_logger->flush_log();
            throw std::out_of_range("Defined index " + std::to_string(culledReactionIndex) + " is out of bounds for reaction index map of size " + std::to_string(state->reaction_index_map.size()) + ".");
        }
        return state->reaction_index_map[culledReactionIndex];
    }

    void DefinedEngineView::collect(
        scratch::StateBlob& ctx,
        const std::vector<std::string> &peNames
    ) const {
        auto* state = scratch::get_state<scratch::DefinedEngineViewScratchPad, true>(ctx);
        std::unordered_set<Species> seenSpecies;

        const auto& fullNetworkReactionSet = m_baseEngine.getNetworkReactions(ctx);
        for (const auto& peName : peNames) {
            if (!fullNetworkReactionSet.contains(peName)) {
                LOG_ERROR(m_logger, "Reaction with name '{}' not found in the base engine's network reactions. Aborting...", peName);
                m_logger->flush_log();
                throw std::runtime_error("Reaction with name '" + std::string(peName) + "' not found in the base engine's network reactions.");
            }
            const reaction::Reaction* reaction = &fullNetworkReactionSet[peName];
            for (const auto& reactant : reaction->reactants()) {
                if (!seenSpecies.contains(reactant)) {
                    seenSpecies.insert(reactant);
                    state->active_species.emplace(reactant);
                }
            }
            for (const auto& product : reaction->products()) {
                if (!seenSpecies.contains(product)) {
                    seenSpecies.insert(product);
                    state->active_species.emplace(product);
                }
            }
            state->active_reactions.add_reaction(*reaction);
        }
        LOG_TRACE_L3(m_logger, "DefinedEngineView built with {} active species and {} active reactions.", m_activeSpecies.size(), m_activeReactions.size());
        LOG_TRACE_L3(m_logger, "Active species: {}", [this]() -> std::string {
            std::string result;
            for (const auto& species : m_activeSpecies) {
                result += std::string(species.name()) + ", ";
            }
            if (!result.empty()) {
                result.pop_back(); // Remove last space
                result.pop_back(); // Remove last comma
            }
            return result;
        }());
        LOG_TRACE_L3(m_logger, "Active reactions: {}", [this]() -> std::string {
            std::string result;
            for (const auto& reaction : m_activeReactions) {
                result += std::string(reaction->id()) + ", ";
            }
            if (!result.empty()) {
                result.pop_back(); // Remove last space
                result.pop_back(); // Remove last comma
            }
            return result;
        }());
        state->species_index_map = constructSpeciesIndexMap(ctx);
        state->reaction_index_map = constructReactionIndexMap(ctx);
    }


    ////////////////////////////////////////////
    /// FileDefinedEngineView Implementation ///
    /////////////////////////////////////////////

    FileDefinedEngineView::FileDefinedEngineView(
        GraphEngine &baseEngine,
        const std::string &fileName,
        const io::NetworkFileParser &parser
    ):
    DefinedEngineView(parser.parse(fileName), baseEngine),
    m_fileName(fileName),
    m_parser(parser) {}
}
