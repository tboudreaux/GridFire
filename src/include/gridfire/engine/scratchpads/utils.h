/**
 * @file utils.h
 * @brief Utility functions for convenient scratchpad retrieval with exception handling.
 *
 * This header provides helper functions that wrap StateBlob's get() methods,
 * converting error codes into exceptions for simpler error handling. These
 * utilities eliminate the need to manually check std::expected results and
 * switch on error codes at every call site.
 *
 * @par Purpose
 * The utility functions provide:
 * - Exception-based error handling instead of std::expected
 * - Consistent error messages across the codebase
 * - Both mutable and const-correct overloads
 * - Optional initialization checking via template parameter
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/utils.h"
 * #include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
 *
 * using namespace gridfire::engine::scratch;
 *
 * void compute(StateBlob& blob) {
 *     // Simple retrieval - throws if not found
 *     GraphEngineScratchPad* scratch = get_state<GraphEngineScratchPad>(blob);
 *     scratch->initialize(engine);
 *
 *     // Retrieval with initialization check - throws if not initialized
 *     auto* initialized_scratch = get_state<GraphEngineScratchPad, true>(blob);
 *     // Safe to use - guaranteed to be initialized
 * }
 *
 * void read_only_access(const StateBlob& blob) {
 *     // Const overload for read-only access
 *     const GraphEngineScratchPad* scratch = get_state<GraphEngineScratchPad>(blob);
 *     bool ready = scratch->is_initialized();
 * }
 * @endcode
 *
 * @par Thread Safety
 * These functions inherit the thread safety characteristics of StateBlob.
 * They are **not thread-safe** - each thread should operate on its own
 * StateBlob instance.
 *
 * @see StateBlob
 * @see IsScratchPad
 * @see exceptions::ScratchPadError
 */

#pragma once
#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/exceptions/error_scratchpad.h"

namespace gridfire::engine::scratch {

/**
 * @brief Retrieve a scratchpad from a StateBlob, throwing on error.
 *
 * Convenience wrapper around StateBlob::get() that converts error codes
 * into exceptions. Use this when you expect the scratchpad to exist and
 * want exception-based error handling.
 *
 * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
 *
 * @param ctx The StateBlob to retrieve the scratchpad from.
 *
 * @return Pointer to the requested scratchpad.
 *
 * @throws exceptions::ScratchPadError if the scratchpad is not found,
 *         cannot be cast to the requested type, or any other error occurs.
 *
 * @par Examples
 * @code{.cpp}
 * StateBlob blob;
 * blob.enroll<GraphEngineScratchPad>();
 *
 * // Retrieve the scratchpad - throws if not enrolled
 * GraphEngineScratchPad* scratch = get_state<GraphEngineScratchPad>(blob);
 * scratch->initialize(engine);
 * @endcode
 */
template <IsScratchPad CTX>
CTX* get_state(StateBlob& ctx) {
    auto result = ctx.get<CTX>();
    if (!result.has_value()) {
        switch (result.error()) {
            case StateBlob::Error::SCRATCHPAD_NOT_FOUND:
                throw exceptions::ScratchPadError("Requested scratchpad not found in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_BAD_CAST:
                throw exceptions::ScratchPadError("Failed to cast scratchpad to the requested type.");
            case StateBlob::Error::SCRATCHPAD_TYPE_COLLISION:
                throw exceptions::ScratchPadError("Scratchpad type collision detected in StateBlob.");
            default:
                throw exceptions::ScratchPadError("Unknown error occurred while retrieving scratchpad from StateBlob.");
        }
    }
    return result.value();
}

/**
 * @brief Retrieve a const scratchpad from a const StateBlob, throwing on error.
 *
 * Const-correct overload of get_state() for read-only access to scratchpads.
 * Use this when you have a const reference to the StateBlob and only need
 * to read from the scratchpad.
 *
 * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
 *
 * @param ctx The const StateBlob to retrieve the scratchpad from.
 *
 * @return Const pointer to the requested scratchpad.
 *
 * @throws exceptions::ScratchPadError if the scratchpad is not found,
 *         cannot be cast to the requested type, or any other error occurs.
 *
 * @par Examples
 * @code{.cpp}
 * void inspect(const StateBlob& blob) {
 *     const GraphEngineScratchPad* scratch = get_state<GraphEngineScratchPad>(blob);
 *     if (scratch->is_initialized()) {
 *         // Read from scratch...
 *     }
 * }
 * @endcode
 */
template <IsScratchPad CTX>
const CTX* get_state(const StateBlob& ctx) {
    auto result = ctx.get<CTX>();
    if (!result.has_value()) {
        switch (result.error()) {
            case StateBlob::Error::SCRATCHPAD_NOT_FOUND:
                throw exceptions::ScratchPadError("Requested scratchpad not found in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_BAD_CAST:
                throw exceptions::ScratchPadError("Failed to cast scratchpad to the requested type.");
            case StateBlob::Error::SCRATCHPAD_TYPE_COLLISION:
                throw exceptions::ScratchPadError("Scratchpad type collision detected in StateBlob.");
            default:
                throw exceptions::ScratchPadError("Unknown error occurred while retrieving scratchpad from StateBlob.");
        }
    }
    return result.value();
}

/**
 * @brief Retrieve a scratchpad with optional initialization check, throwing on error.
 *
 * Extended version of get_state() that can optionally verify the scratchpad
 * is initialized before returning it. When MUST_BE_INITIALIZED is true, an
 * exception is thrown if the scratchpad exists but hasn't been initialized.
 *
 * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
 * @tparam MUST_BE_INITIALIZED If true, throws when the scratchpad is not initialized.
 *
 * @param ctx The StateBlob to retrieve the scratchpad from.
 *
 * @return Pointer to the requested scratchpad (guaranteed initialized if MUST_BE_INITIALIZED is true).
 *
 * @throws exceptions::ScratchPadError if the scratchpad is not found,
 *         cannot be cast, is not initialized (when required), or any other error.
 *
 * @par Examples
 * @code{.cpp}
 * // Ensure scratchpad is initialized before use
 * try {
 *     auto* scratch = get_state<GraphEngineScratchPad, true>(blob);
 *     // Guaranteed to be initialized here
 *     use_scratchpad(*scratch);
 * } catch (const exceptions::ScratchPadError& e) {
 *     // Handle missing or uninitialized scratchpad
 * }
 * @endcode
 */
template <IsScratchPad CTX, bool MUST_BE_INITIALIZED>
CTX* get_state(StateBlob& ctx) {
    auto result = ctx.get<CTX, MUST_BE_INITIALIZED>();
    if (!result.has_value()) {
        switch (result.error()) {
            case StateBlob::Error::SCRATCHPAD_NOT_FOUND:
                throw exceptions::ScratchPadError("Requested scratchpad not found in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_BAD_CAST:
                throw exceptions::ScratchPadError("Failed to cast scratchpad to the requested type.");
            case StateBlob::Error::SCRATCHPAD_TYPE_COLLISION:
                throw exceptions::ScratchPadError("Scratchpad type collision detected in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_NOT_INITIALIZED:
                throw exceptions::ScratchPadError("Requested scratchpad not initialized in StateBlob. If this is acceptable behavior, use scratch::get_state<>() without the MUST_BE_INITIALIZED template parameter.");
            default:
                throw exceptions::ScratchPadError("Unknown error occurred while retrieving scratchpad from StateBlob.");
        }
    }
    return result.value();
}

/**
 * @brief Retrieve a const scratchpad with optional initialization check, throwing on error.
 *
 * Const-correct overload of the initialization-checking get_state(). Combines
 * read-only access with optional initialization verification.
 *
 * @tparam CTX The scratchpad type to retrieve (must satisfy IsScratchPad).
 * @tparam MUST_BE_INITIALIZED If true, throws when the scratchpad is not initialized.
 *
 * @param ctx The const StateBlob to retrieve the scratchpad from.
 *
 * @return Const pointer to the requested scratchpad (guaranteed initialized if MUST_BE_INITIALIZED is true).
 *
 * @throws exceptions::ScratchPadError if the scratchpad is not found,
 *         cannot be cast, is not initialized (when required), or any other error.
 *
 * @par Examples
 * @code{.cpp}
 * void validate(const StateBlob& blob) {
 *     // Get const access, ensuring initialization
 *     const auto* scratch = get_state<GraphEngineScratchPad, true>(blob);
 *     // Safe to read - guaranteed initialized
 *     const auto& adfun = scratch->rhsADFun.value();
 * }
 * @endcode
 */
template <IsScratchPad CTX, bool MUST_BE_INITIALIZED>
const CTX* get_state(const StateBlob& ctx) {
    auto result = ctx.get<CTX, MUST_BE_INITIALIZED>();
    if (!result.has_value()) {
        switch (result.error()) {
            case StateBlob::Error::SCRATCHPAD_NOT_FOUND:
                throw exceptions::ScratchPadError("Requested scratchpad not found in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_BAD_CAST:
                throw exceptions::ScratchPadError("Failed to cast scratchpad to the requested type.");
            case StateBlob::Error::SCRATCHPAD_TYPE_COLLISION:
                throw exceptions::ScratchPadError("Scratchpad type collision detected in StateBlob.");
            case StateBlob::Error::SCRATCHPAD_NOT_INITIALIZED:
                throw exceptions::ScratchPadError("Requested scratchpad is not initialized. If this is acceptable behavior, use scratch::get_state<>() without the MUST_BE_INITIALIZED template parameter.");
            default:
                throw exceptions::ScratchPadError("Unknown error occurred while retrieving scratchpad from StateBlob.");
        }
    }
    return result.value();
}

} // namespace gridfire::engine::scratch
