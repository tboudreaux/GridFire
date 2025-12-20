/**
 * @file scratchpads.h
 * @brief Unified header for the scratchpad memory management system.
 *
 * This is the main include file for the scratchpad subsystem. It provides
 * a single include point for all scratchpad-related functionality, including
 * the abstract base class, concrete implementations, type definitions,
 * the StateBlob container, utility functions, and formatters for debugging.
 *
 * @par What are Scratchpads?
 * Scratchpads are temporary working memory containers used by computational
 * engines during ODE integration and reaction network calculations. They serve
 * several critical purposes:
 *
 * - **Performance**: Pre-allocate memory once and reuse across iterations,
 *   avoiding repeated heap allocations during time-critical computations.
 * - **Caching**: Store intermediate results (Jacobians, derivatives, etc.)
 *   that can be reused across solver steps.
 * - **Thread Safety**: Provide thread-local storage for parallel execution,
 *   where each thread operates on its own independent scratchpad instance.
 * - **State Management**: Encapsulate engine-specific working state separate
 *   from the engine's persistent configuration.
 *
 * @par Architecture Overview
 * The scratchpad system consists of:
 *
 * - **AbstractScratchPad**: Base interface defining `is_initialized()` and `clone()`
 * - **Concrete Scratchpads**: Engine-specific implementations
 *   - `GraphEngineScratchPad`: CppAD ADFun and Jacobian caches
 *   - `AdaptiveEngineViewScratchPad`: Active species/reactions for adaptive networks
 *   - `DefinedEngineViewScratchPad`: Species/reactions for static networks
 *   - `MultiscalePartitioningEngineViewScratchPad`: QSE groups and SUNDIALS context
 * - **StateBlob**: Container managing multiple scratchpads with type-safe access
 * - **Utilities**: Helper functions for exception-based retrieval
 * - **Formatters**: std::format specializations for debugging output
 *
 * @par Why Use Scratchpads?
 * During numerical integration of stiff reaction networks, engines must:
 * 1. Evaluate right-hand side functions (species derivatives)
 * 2. Compute sparse Jacobian matrices
 * 3. Solve linear systems within Newton iterations
 *
 * These operations require substantial temporary memory. Without scratchpads,
 * each evaluation would allocate and deallocate working buffers, causing:
 * - Memory fragmentation
 * - Cache thrashing
 * - Unnecessary allocation overhead
 *
 * Scratchpads solve this by providing persistent, reusable working memory
 * that lives for the duration of an integration step (or longer).
 *
 * @par Thread Safety Model
 * The scratchpad system is designed for thread-local usage:
 *
 * - Scratchpads are **not thread-safe** by design
 * - Each thread must have its own scratchpad instances
 * - Use `StateBlob::clone_structure()` to create independent copies for workers
 * - The original scratchpad/blob remains usable by the main thread
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/scratchpads.h"
 *
 * using namespace gridfire::engine::scratch;
 *
 * // === Basic Usage ===
 * // Create a StateBlob and enroll scratchpads
 * StateBlob blob;
 * blob.enroll<GraphEngineScratchPad>();
 * blob.enroll<AdaptiveEngineViewScratchPad>();
 *
 * // Initialize scratchpads
 * auto* graph_scratch = get_state<GraphEngineScratchPad>(blob);
 * graph_scratch->initialize(engine);
 *
 * // Use initialized scratchpad
 * auto* scratch = get_state<GraphEngineScratchPad, true>(blob);  // Throws if not initialized
 * auto& adfun = scratch->rhsADFun.value();
 *
 * // === Parallel Execution ===
 * // Clone for worker threads
 * std::vector<std::unique_ptr<StateBlob>> worker_blobs;
 * for (int i = 0; i < num_threads; ++i) {
 *     worker_blobs.push_back(blob.clone_structure());
 * }
 *
 * // Each worker uses its own blob
 * #pragma omp parallel for
 * for (int i = 0; i < work_items; ++i) {
 *     int tid = omp_get_thread_num();
 *     StateBlob& my_blob = *worker_blobs[tid];
 *     auto* my_scratch = get_state<GraphEngineScratchPad>(my_blob);
 *     // Thread-safe: each thread has its own scratchpad
 *     compute_with_scratchpad(*my_scratch);
 * }
 *
 * // === Debugging ===
 * // Use formatters for logging
 * std::cout << std::format("Blob state: {}\n", blob);
 * // Output: StateBlob(Enrolled: GraphEngineScratchPad(...), AdaptiveEngineViewScratchPad(...))
 *
 * // Check status
 * auto status = blob.get_status<GraphEngineScratchPad>();
 * if (status == StateBlob::ScratchPadStatus::ENROLLED_NOT_INITIALIZED) {
 *     // Need to initialize before use
 * }
 * @endcode
 *
 * @par Included Headers
 * This header includes:
 * - scratchpad_abstract.h - AbstractScratchPad base class
 * - engine_graph_scratchpad.h - GraphEngineScratchPad
 * - engine_adaptive_scratchpad.h - AdaptiveEngineViewScratchPad
 * - engine_multiscale_scratchpad.h - MultiscalePartitioningEngineViewScratchPad
 * - engine_defined_scratchpad.h - DefinedEngineViewScratchPad
 * - types.h - ScratchPadType enumeration
 * - blob.h - StateBlob container
 * - utils.h - get_state() helper functions
 * - formatters.h - std::format specializations
 *
 * @see AbstractScratchPad
 * @see StateBlob
 * @see ScratchPadType
 */

/**
 * @namespace gridfire::engine::scratch
 * @brief Scratchpad memory management for computational engines.
 *
 * The scratch namespace contains all components related to temporary working
 * memory management for GridFire's computational engines. This includes the
 * abstract scratchpad interface, concrete implementations for each engine type,
 * the StateBlob container for managing multiple scratchpads, and utilities
 * for convenient access.
 *
 * @par Key Components
 * - AbstractScratchPad: Interface for all scratchpad types
 * - GraphEngineScratchPad: Working memory for CppAD-based graph engines
 * - AdaptiveEngineViewScratchPad: Storage for adaptive network computations
 * - DefinedEngineViewScratchPad: Storage for static network computations
 * - MultiscalePartitioningEngineViewScratchPad: QSE solver state management
 * - StateBlob: Type-safe container for multiple scratchpads
 * - ScratchPadType: Enumeration of registered scratchpad types
 * - IsScratchPad: Concept constraining valid scratchpad types
 * - get_state(): Utility functions for exception-based retrieval
 *
 * @par Design Philosophy
 * The scratchpad system follows these design principles:
 *
 * 1. **Separation of Concerns**: Working memory is separate from engine configuration
 * 2. **Type Safety**: Compile-time verification of scratchpad types via concepts
 * 3. **Performance**: O(1) access via enum-indexed arrays
 * 4. **Thread Locality**: Each thread owns its scratchpad instances
 * 5. **Clonability**: Deep copying enables parallel execution patterns
 */

#pragma once

#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_adaptive_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_multiscale_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_defined_scratchpad.h"
#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/utils.h"
#include "gridfire/engine/scratchpads/formatters.h"