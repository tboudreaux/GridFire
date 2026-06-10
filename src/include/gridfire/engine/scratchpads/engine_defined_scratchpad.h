/**
 * @file engine_defined_scratchpad.h
 * @brief Scratchpad implementation for the DefinedEngineView.
 *
 * This header defines the DefinedEngineViewScratchPad, a concrete implementation
 * of AbstractScratchPad designed for use with engines that have a statically
 * defined reaction network. It provides storage for the active species, reactions,
 * and index mappings that enable efficient lookups during computations.
 *
 * @par Purpose
 * The DefinedEngineViewScratchPad stores:
 * - The set of active species in the defined network
 * - The set of active reactions in the defined network
 * - Index mappings for efficient species and reaction lookups
 * - A cached vector representation of active species for performance
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/engine_defined_scratchpad.h"
 *
 * // Create a scratchpad for a defined engine
 * gridfire::engine::scratch::DefinedEngineViewScratchPad scratch;
 *
 * // Populate active species
 * scratch.active_species.insert(species1);
 * scratch.active_species.insert(species2);
 *
 * // Set up index mappings for efficient lookups
 * scratch.species_index_map = {0, 1, 2, 3};
 * scratch.reaction_index_map = {0, 1};
 *
 * if (scratch.is_initialized()) {
 *     // Use the scratchpad for computations
 *     for (size_t idx : scratch.species_index_map) {
 *         // Process species by index
 *     }
 * }
 *
 * // Clone for parallel execution
 * auto worker_scratch = scratch.clone();
 * @endcode
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread should have its own instance
 * of DefinedEngineViewScratchPad. Use clone() to create independent copies
 * for parallel workers.
 *
 * @see AbstractScratchPad
 * @see DefinedEngineView
 */

#pragma once

#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "fourdst/atomic/atomicSpecies.h"
#include <vector>
#include <memory>
#include <set>
#include <optional>

namespace gridfire::engine::scratch {

/**
 * @brief Scratchpad for storing working memory used by defined reaction network engines.
 *
 * DefinedEngineViewScratchPad provides storage for species and reaction data
 * in engines with statically defined reaction networks. It includes index mappings
 * for efficient lookups and an optional cache for the active species vector.
 *
 * Unlike adaptive scratchpads, the defined scratchpad is considered initialized
 * by default (has_initialized = true), as the network structure is known at
 * construction time.
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread should operate on its own
 * independent instance. Use clone() to create copies for parallel execution.
 */
struct DefinedEngineViewScratchPad final : AbstractScratchPad {
    /// @brief Unique identifier for this scratchpad type.
    constexpr static auto ID = ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD;

    /// @brief Flag indicating whether the scratchpad is initialized (default: true).
    bool has_initialized = true;

    /// @brief Set of species active in the defined network.
    std::set<fourdst::atomic::Species> active_species;

    /// @brief Set of reactions active in the defined network.
    reaction::ReactionSet active_reactions;

    /// @brief Mapping from local indices to global species indices.
    std::vector<size_t> species_index_map;

    /// @brief Mapping from local indices to global reaction indices.
    std::vector<size_t> reaction_index_map;

    /// @brief Cached vector of active species for performance optimization.
    /// @details This optional cache avoids repeated conversion from set to vector.
    std::optional<std::vector<fourdst::atomic::Species>> active_species_vector_cache = std::nullopt;

    /**
     * @brief Check whether the scratchpad has been initialized.
     * @return true if initialized (always true by default for defined networks).
     */
    bool is_initialized() const override { return has_initialized; }

    /**
     * @brief Create a deep copy of this scratchpad.
     *
     * Creates an independent copy of all internal state, including
     * active species, reactions, index mappings, and the species vector cache.
     * The clone can be modified without affecting the original.
     *
     * @return A unique pointer to the cloned scratchpad.
     *
     * @par Examples
     * @code{.cpp}
     * DefinedEngineViewScratchPad scratch;
     * scratch.active_species.insert(species);
     * scratch.species_index_map = {0, 1, 2};
     *
     * // Create independent copy for a worker thread
     * auto worker_copy = scratch.clone();
     * @endcode
     */
    std::unique_ptr<AbstractScratchPad> clone() const override {
        auto pad = std::make_unique<DefinedEngineViewScratchPad>();
        pad->has_initialized = this->has_initialized;
        pad->active_species = this->active_species;
        pad->active_reactions = this->active_reactions;
        pad->species_index_map = this->species_index_map;
        pad->reaction_index_map = this->reaction_index_map;
        pad->active_species_vector_cache = this->active_species_vector_cache;
        return pad;
    }
};

} // namespace gridfire::engine::scratch
