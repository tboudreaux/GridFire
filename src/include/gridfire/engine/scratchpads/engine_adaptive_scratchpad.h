/**
 * @file engine_adaptive_scratchpad.h
 * @brief Scratchpad implementation for the AdaptiveEngineView.
 *
 * This header defines the AdaptiveEngineViewScratchPad, a concrete implementation
 * of AbstractScratchPad designed for use with the AdaptiveEngineView. It provides
 * thread-local storage for active species and reactions that are dynamically
 * determined during adaptive network computations.
 *
 * @par Purpose
 * The AdaptiveEngineViewScratchPad stores:
 * - The set of currently active species in the adaptive network
 * - The set of currently active reactions based on network topology
 * - Initialization state to track whether the scratchpad is ready for use
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/engine_adaptive_scratchpad.h"
 * #include "gridfire/engine/views/engine_adaptive.h"
 *
 * // Create and initialize the scratchpad from an AdaptiveEngineView
 * gridfire::engine::scratch::AdaptiveEngineViewScratchPad scratch;
 * AdaptiveEngineView engine = create_adaptive_engine();
 * scratch.initialize(engine);
 *
 * if (scratch.is_initialized()) {
 *     // Access active species for computation
 *     for (const auto& species : scratch.active_species) {
 *         // Process each active species
 *     }
 * }
 *
 * // Clone for parallel execution
 * auto worker_scratch = scratch.clone();
 * @endcode
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread should have its own instance
 * of AdaptiveEngineViewScratchPad. Use clone() to create independent copies
 * for parallel workers.
 *
 * @see AbstractScratchPad
 * @see AdaptiveEngineView
 */

#pragma once

#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/types.h"

#include "gridfire/engine/views/engine_adaptive.h"

namespace gridfire::engine::scratch {

/**
 * @brief Scratchpad for storing working memory used by AdaptiveEngineView computations.
 *
 * AdaptiveEngineViewScratchPad provides temporary storage for the active species
 * and reactions determined by the adaptive network algorithm. This allows the
 * engine to avoid recalculating network topology on every evaluation.
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread should operate on its own
 * independent instance. Use clone() to create copies for parallel execution.
 */
struct AdaptiveEngineViewScratchPad final : AbstractScratchPad {
    /// @brief Unique identifier for this scratchpad type.
    static constexpr auto ID = ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD;

    /// @brief Flag indicating whether the scratchpad has been initialized.
    bool has_initialized = false;


    /// @brief Vector of species currently active in the adaptive network.
    std::vector<fourdst::atomic::Species> active_species;

    /// @brief Set of reactions currently active in the adaptive network.
    reaction::ReactionSet active_reactions;

    /**
     * @brief Check whether the scratchpad has been initialized.
     * @return true if initialized, false otherwise.
     */
    [[nodiscard]] bool is_initialized() const override { return has_initialized; }

    /**
     * @brief Initialize the scratchpad from an AdaptiveEngineView.
     *
     * Clears any existing state and prepares the scratchpad for use.
     * This method is idempotent; calling it multiple times has no effect
     * after the first successful initialization.
     *
     * @param engine The AdaptiveEngineView to initialize from.
     *
     * @par Examples
     * @code{.cpp}
     * AdaptiveEngineViewScratchPad scratch;
     * AdaptiveEngineView engine = create_engine();
     * scratch.initialize(engine);
     * @endcode
     */
    void initialize(const AdaptiveEngineView& engine) {
        if (has_initialized) return;
        active_species.clear();
        active_reactions.clear();
        has_initialized = true;
    }

    /**
     * @brief Create a deep copy of this scratchpad.
     *
     * Creates an independent copy of all internal state, including
     * active species and reactions. The clone can be modified without
     * affecting the original.
     *
     * @return A unique pointer to the cloned scratchpad.
     *
     * @par Examples
     * @code{.cpp}
     * auto original = std::make_unique<AdaptiveEngineViewScratchPad>();
     * original->initialize(engine);
     *
     * // Create independent copy for a worker thread
     * auto worker_copy = original->clone();
     * @endcode
     */
    std::unique_ptr<AbstractScratchPad> clone() const override {
        auto ptr = std::make_unique<AdaptiveEngineViewScratchPad>();
        ptr->has_initialized = has_initialized;
        ptr->active_species = active_species;
        ptr->active_reactions = active_reactions;
        return ptr;
    }
};

} // namespace gridfire::engine::scratch
