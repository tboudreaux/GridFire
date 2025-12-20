/**
 * @file engine_multiscale_scratchpad.h
 * @brief Scratchpad implementation for the MultiscalePartitioningEngineView.
 *
 * This header defines the MultiscalePartitioningEngineViewScratchPad, a concrete
 * implementation of AbstractScratchPad designed for use with multiscale partitioning
 * algorithms. It provides thread-local storage for QSE (Quasi-Static Equilibrium)
 * groups, solvers, species classifications, and SUNDIALS contexts required for
 * solving stiff subsystems.
 *
 * @par Purpose
 * The MultiscalePartitioningEngineViewScratchPad stores:
 * - QSE groups representing clusters of species in quasi-static equilibrium
 * - QSE solvers for each group (managed via unique_ptr)
 * - Classification of species into dynamic vs. algebraic categories
 * - A composition cache for efficient lookup of computed compositions
 * - A SUNDIALS context for numerical solver operations
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/engine_multiscale_scratchpad.h"
 *
 * // Create and initialize the scratchpad
 * gridfire::engine::scratch::MultiscalePartitioningEngineViewScratchPad scratch;
 * scratch.initialize();
 *
 * if (scratch.is_initialized()) {
 *     // Add QSE groups and species classifications
 *     scratch.dynamic_species.push_back(hydrogen);
 *     scratch.dynamic_species.push_back(helium);
 *     scratch.algebraic_species.push_back(carbon);
 *
 *     // Use the SUNDIALS context for solver operations
 *     SUNContext ctx = scratch.sun_ctx;
 * }
 *
 * // Clone for parallel execution (note: solvers not cloned)
 * auto worker_scratch = scratch.clone();
 * worker_scratch->initialize();  // Must re-initialize for new SUNContext
 * @endcode
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread must have its own instance
 * because SUNDIALS contexts and QSE solvers maintain internal state. When
 * cloning, the QSE solvers are not copied and the SUNContext is not cloned;
 * the new instance must call initialize() to create its own context.
 *
 * @see AbstractScratchPad
 * @see MultiscalePartitioningEngineView
 * @see SUNContext
 */

#pragma once

#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/types.h"

#include "fourdst/atomic/atomicSpecies.h"

#include <vector>
#include <memory>
#include <unordered_map>

#include "sundials/sundials_context.h"

namespace gridfire::engine::scratch {

/**
 * @brief Scratchpad for multiscale partitioning engine computations with QSE groups.
 *
 * MultiscalePartitioningEngineViewScratchPad provides thread-local storage for
 * the multiscale partitioning algorithm, which separates species into fast
 * (algebraic/QSE) and slow (dynamic) timescale groups. This enables efficient
 * integration of stiff reaction networks by solving QSE subsystems separately.
 *
 * @par Thread Safety
 * This class is **not thread-safe**. SUNDIALS contexts and QSE solver objects
 * maintain internal state that cannot be shared across threads. Each thread
 * must have its own scratchpad instance with its own SUNContext. When using
 * clone(), the new instance starts uninitialized and must call initialize()
 * to create a new SUNContext.
 *
 * @note The destructor properly cleans up SUNDIALS resources by freeing the
 *       SUNContext after clearing all solvers that depend on it.
 *
 * @warning QSE solvers are **not cloned** - they must be re-created in the
 *          cloned scratchpad as needed.
 */
struct MultiscalePartitioningEngineViewScratchPad final : AbstractScratchPad {
    /// @brief Type alias for QSE group from the multiscale engine view.
    using QSEGroup = MultiscalePartitioningEngineView::QSEGroup;

    /// @brief Type alias for QSE solver from the multiscale engine view.
    using QSESolver = MultiscalePartitioningEngineView::QSESolver;

    /// @brief Type alias for atomic species.
    using Species = fourdst::atomic::Species;

    /// @brief Unique identifier for this scratchpad type.
    static constexpr auto ID = ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD;

    /// @brief Flag indicating whether the scratchpad has been initialized.
    bool has_initialized = false;

    /// @brief Vector of QSE groups representing equilibrium clusters.
    std::vector<QSEGroup> qse_groups;

    /// @brief Vector of QSE solvers, one per QSE group.
    /// @note These are not cloned; new instances must create their own solvers.
    std::vector<std::unique_ptr<QSESolver>> qse_solvers;

    /// @brief Species that evolve on the dynamic (slow) timescale.
    std::vector<Species> dynamic_species;

    /// @brief Species that are solved algebraically (fast timescale/QSE).
    std::vector<Species> algebraic_species;

    /// @brief Cache of computed compositions indexed by a hash key.
    std::unordered_map<uint64_t, fourdst::composition::Composition> composition_cache;

    /// @brief SUNDIALS context for solver operations.
    /// @note Must be freed in destructor after clearing solvers.
    SUNContext sun_ctx = nullptr;

    /**
     * @brief Check whether the scratchpad has been initialized.
     * @return true if initialized with a valid SUNContext, false otherwise.
     */
    [[nodiscard]] bool is_initialized() const override { return has_initialized; }

    /**
     * @brief Initialize the scratchpad by creating a SUNDIALS context.
     *
     * Creates a new SUNContext for use with SUNDIALS solvers. This method
     * is idempotent; calling it multiple times has no effect after the
     * first successful initialization.
     *
     * @throws std::runtime_error if SUNContext creation fails.
     *
     * @par Examples
     * @code{.cpp}
     * MultiscalePartitioningEngineViewScratchPad scratch;
     * scratch.initialize();
     *
     * // Now safe to use SUNContext for solver creation
     * SUNContext ctx = scratch.sun_ctx;
     * @endcode
     */
    void initialize() {
        if (has_initialized) return;

        const int flag = SUNContext_Create(SUN_COMM_NULL, &sun_ctx);
        if (flag != 0) {
            throw std::runtime_error("Failed to create SUNContext in MultiscalePartitioningEngineViewScratchPad.");
        }

        has_initialized = true;
    }

    /**
     * @brief Destructor that properly releases SUNDIALS resources.
     *
     * Clears all QSE solvers before freeing the SUNContext to ensure
     * proper cleanup order and avoid dangling references.
     */
    ~MultiscalePartitioningEngineViewScratchPad() override {
        qse_solvers.clear();
        if (sun_ctx != nullptr) {
            SUNContext_Free(&sun_ctx);
            sun_ctx = nullptr;
        }
    }

    /**
     * @brief Create a partial copy of this scratchpad.
     *
     * Creates a copy with the QSE groups, species classifications, and
     * composition cache.
     *
     * @return A unique pointer to the cloned scratchpad.
     *
     * @note The new instance will automatically initialize a new SUNContext and clone the QSE solvers.
     *
     * @par Examples
     * @code{.cpp}
     * MultiscalePartitioningEngineViewScratchPad scratch;
     * scratch.initialize();
     * scratch.dynamic_species = {hydrogen, helium};
     *
     * // Clone for a worker thread
     * auto worker = scratch.clone();
     *
     * @endcode
     */
    [[nodiscard]] std::unique_ptr<AbstractScratchPad> clone() const override {
        auto clone_pad = std::make_unique<MultiscalePartitioningEngineViewScratchPad>();
        clone_pad->qse_groups = this->qse_groups;
        clone_pad->dynamic_species = this->dynamic_species;
        clone_pad->algebraic_species = this->algebraic_species;
        clone_pad->composition_cache = this->composition_cache;

        clone_pad->initialize();
        clone_pad->qse_solvers.reserve(this->qse_solvers.size());
        for (const auto& solver : qse_solvers) {
            clone_pad->qse_solvers.push_back(solver->clone(clone_pad->sun_ctx)); // Must rebind context to new SUNContext
        }
        return clone_pad;
    }
};

} // namespace gridfire::engine::scratch
