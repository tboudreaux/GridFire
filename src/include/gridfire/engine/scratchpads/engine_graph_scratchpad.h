/**
 * @file engine_graph_scratchpad.h
 * @brief Scratchpad implementation for the GraphEngine using CppAD automatic differentiation.
 *
 * This header defines the GraphEngineScratchPad, a concrete implementation of
 * AbstractScratchPad designed for use with the GraphEngine. It provides thread-local
 * storage for CppAD automatic differentiation functions, Jacobian computation work
 * structures, and cached derivatives used during ODE integration.
 *
 * @par Purpose
 * The GraphEngineScratchPad stores:
 * - A local copy of the CppAD ADFun for RHS evaluation
 * - Work structures for sparse Jacobian calculations
 * - Cached abundance values for efficient reuse
 * - Cached step derivatives and Jacobian subsets by timestep
 * - The most recent RHS calculation for warm-starting
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
 * #include "gridfire/engine/engine_graph.h"
 *
 * // Create and initialize the scratchpad from a GraphEngine
 * gridfire::engine::scratch::GraphEngineScratchPad scratch;
 * GraphEngine engine = create_graph_engine();
 * scratch.initialize(engine);
 *
 * if (scratch.is_initialized()) {
 *     // Access the local ADFun for thread-safe evaluation
 *     auto& adfun = scratch.rhsADFun.value();
 *
 *     // Use cached Jacobian work for efficient sparse computations
 *     auto& jac_work = scratch.jac_work;
 * }
 *
 * // Clone for parallel execution
 * auto worker_scratch = scratch.clone();
 * @endcode
 *
 * @par Thread Safety
 * This class is **not thread-safe**. Each thread must have its own instance
 * of GraphEngineScratchPad because CppAD ADFun objects maintain internal state
 * that is modified during evaluation. Use clone() to create independent copies
 * for parallel workers, ensuring each thread has its own ADFun instance.
 *
 * @see AbstractScratchPad
 * @see GraphEngine
 * @see CppAD::ADFun
 */

#pragma once
#include <vector>

#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_abstract.h"

#include "cppad/cppad.hpp"

#include <optional>


namespace gridfire::engine::scratch {

/**
 * @brief Scratchpad for storing CppAD automatic differentiation state for GraphEngine.
 *
 * GraphEngineScratchPad provides thread-local storage for all CppAD-related
 * objects needed during ODE integration with the GraphEngine. This includes
 * the ADFun object for evaluating the right-hand side of the ODE and computing
 * Jacobians, as well as various caches to improve performance.
 *
 * @par Thread Safety
 * This class is **not thread-safe**. CppAD ADFun objects maintain internal
 * state that is modified during Forward and Reverse mode operations. Each
 * thread must have its own scratchpad instance. Use clone() to create
 * independent copies for parallel execution.
 *
 * @note When cloning, if the rhsADFun has not been initialized, the clone
 *       will also be uninitialized and has_initialized will be false.
 */
struct GraphEngineScratchPad final : AbstractScratchPad {
    /**
     * @brief Result codes for ADFun registration operations.
     */
    enum class ADFunRegistrationResult : uint8_t {
        SUCCESS,           ///< Registration completed successfully.
        ALREADY_REGISTERED ///< ADFun was already registered; no action taken.
    };

    /// @brief CppAD function object for evaluating the ODE right-hand side.
    /// @details Contains the computational graph for automatic differentiation.
    std::optional<CppAD::ADFun<double>> rhsADFun;

    /// @brief Work structure for sparse Jacobian calculations.
    /// @details Reused across Jacobian evaluations to avoid reallocation.
    CppAD::sparse_jac_work jac_work;

    /// @brief Local cache of abundance values for efficient RHS evaluation.
    std::vector<double> local_abundance_cache;

    /// @brief Cache of step derivatives indexed by timestep identifier.
    std::unordered_map<size_t, StepDerivatives<double>> stepDerivativesCache;

    /// @brief Cache of sparse Jacobian subsets indexed by timestep identifier.
    std::unordered_map<size_t, CppAD::sparse_rcv<std::vector<size_t>, std::vector<double>>> jacobianSubsetCache;

    /// @brief Cache of Jacobian work structures indexed by timestep identifier.
    std::unordered_map<size_t, CppAD::sparse_jac_work> jacWorkCache;

    /// @brief The most recent RHS calculation result for warm-starting.
    std::optional<StepDerivatives<double>> most_recent_rhs_calculation;

    /// @brief Flag indicating whether the scratchpad has been initialized.
    bool has_initialized = false;

    /// @brief Unique identifier for this scratchpad type.
    static constexpr auto ID = ScratchPadType::GRAPH_ENGINE_SCRATCHPAD;

    /**
     * @brief Check whether the scratchpad has been initialized.
     * @return true if initialized with a valid ADFun, false otherwise.
     */
    [[nodiscard]] bool is_initialized() const override { return has_initialized; }

    /**
     * @brief Initialize the scratchpad from a GraphEngine.
     *
     * Copies the authoritative ADFun from the engine and clears all caches.
     * This method is idempotent; calling it multiple times has no effect
     * after the first successful initialization.
     *
     * @param engine The GraphEngine to initialize from.
     *
     * @par Examples
     * @code{.cpp}
     * GraphEngineScratchPad scratch;
     * GraphEngine engine = create_engine();
     * scratch.initialize(engine);
     *
     * // Now safe to use for thread-local computations
     * auto& adfun = scratch.rhsADFun.value();
     * @endcode
     */
    void initialize(const GraphEngine& engine) {
        if (has_initialized) return;

        const auto& sourceTape = engine.getAuthoritativeADFun();
        rhsADFun.emplace();
        *rhsADFun = sourceTape;
        jac_work.clear();
        local_abundance_cache.clear();
        stepDerivativesCache.clear();
        jacobianSubsetCache.clear();
        jacWorkCache.clear();
        most_recent_rhs_calculation = std::nullopt;
        has_initialized = true;
    }

    /**
     * @brief Create a deep copy of this scratchpad.
     *
     * Creates an independent copy of all internal state, including the
     * CppAD ADFun object and all caches. The clone can be safely used
     * in a separate thread without affecting the original.
     *
     * @return A unique pointer to the cloned scratchpad.
     *
     * @note If rhsADFun is not initialized, the clone will also be
     *       uninitialized (has_initialized = false).
     *
     * @par Examples
     * @code{.cpp}
     * GraphEngineScratchPad scratch;
     * scratch.initialize(engine);
     *
     * // Create independent copies for parallel workers
     * std::vector<std::unique_ptr<AbstractScratchPad>> worker_pads;
     * for (int i = 0; i < num_threads; ++i) {
     *     worker_pads.push_back(scratch.clone());
     * }
     * @endcode
     */
    [[nodiscard]] std::unique_ptr<AbstractScratchPad> clone() const override {
        auto ptr = std::make_unique<GraphEngineScratchPad>();
        if (!rhsADFun.has_value()) {
            ptr->rhsADFun = std::nullopt;
            ptr->has_initialized = false;
        } else {
            ptr->rhsADFun.emplace();
            *ptr->rhsADFun = rhsADFun.value();
            ptr->has_initialized = true;
        }
        ptr->jac_work = jac_work;
        ptr->local_abundance_cache = local_abundance_cache;
        ptr->stepDerivativesCache = stepDerivativesCache;
        ptr->jacobianSubsetCache = jacobianSubsetCache;
        ptr->jacWorkCache = jacWorkCache;
        ptr->most_recent_rhs_calculation = most_recent_rhs_calculation;
        return ptr;
    }
};

} // namespace gridfire::engine::scratch
