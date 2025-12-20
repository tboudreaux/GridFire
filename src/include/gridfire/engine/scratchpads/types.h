/**
 * @file types.h
 * @brief Type definitions and utilities for the scratchpad system.
 *
 * This header defines the ScratchPadType enumeration which identifies all
 * registered scratchpad types in the system, along with utility functions
 * for querying scratchpad type information at compile-time and runtime.
 *
 * @par Purpose
 * The types header provides:
 * - A centralized enumeration of all scratchpad types
 * - Compile-time constant for the maximum number of scratchpad types
 * - Runtime conversion of scratchpad types to human-readable names
 *
 * @par Adding New Scratchpad Types
 * To add a new scratchpad type:
 * 1. Add a new enumerator before `_COUNT` in ScratchPadType
 * 2. Add a corresponding case in get_scratchpad_type_name()
 * 3. Create the concrete scratchpad class with a static `ID` member set to the new type
 *
 * @par Examples
 * @code{.cpp}
 * #include "gridfire/engine/scratchpads/types.h"
 *
 * using namespace gridfire::engine::scratch;
 *
 * // Get the maximum number of scratchpad types at compile time
 * constexpr size_t max_types = get_max_scratchpad_types();
 * std::array<bool, max_types> enrolled_flags{};
 *
 * // Get a human-readable name for a scratchpad type
 * ScratchPadType type = ScratchPadType::GRAPH_ENGINE_SCRATCHPAD;
 * std::string_view name = get_scratchpad_type_name(type);
 * // name == "GraphEngineScratchPad"
 *
 * // Iterate over all scratchpad types
 * for (size_t i = 0; i < get_max_scratchpad_types(); ++i) {
 *     auto type = static_cast<ScratchPadType>(i);
 *     std::cout << get_scratchpad_type_name(type) << "\n";
 * }
 * @endcode
 *
 * @see AbstractScratchPad
 * @see StateBlob
 */

#pragma once
#include <cstdint>
#include <string_view>

namespace gridfire::engine::scratch {

/**
 * @brief Enumeration of all registered scratchpad types.
 *
 * Each scratchpad implementation must have a unique type identifier in this
 * enumeration. The concrete scratchpad class should define a static `ID`
 * member initialized to its corresponding ScratchPadType value.
 *
 * @note The `_COUNT` enumerator is a sentinel value used to determine the
 *       total number of scratchpad types. It must always be the last entry.
 *       Do not use `_COUNT` as an actual scratchpad type.
 */
enum class ScratchPadType : uint8_t {
    GRAPH_ENGINE_SCRATCHPAD,                          ///< GraphEngineScratchPad for CppAD-based engines.
    MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD,   ///< MultiscalePartitioningEngineViewScratchPad for QSE partitioning.
    ADAPTIVE_ENGINE_VIEW_SCRATCHPAD,                  ///< AdaptiveEngineViewScratchPad for adaptive networks.
    DEFINED_ENGINE_VIEW_SCRATCHPAD,                   ///< DefinedEngineViewScratchPad for static networks.
    PRIMING_ENGINE_VIEW_SCRATCHPAD,                   ///< PrimingEngineViewScratchPad for engine priming.

    _COUNT  ///< Sentinel value representing the total number of scratchpad types. Do not use as a type.
};

/**
 * @brief Get the maximum number of scratchpad types at compile time.
 *
 * Returns the total count of registered scratchpad types, derived from
 * the ScratchPadType::_COUNT sentinel value. This is useful for sizing
 * fixed-size arrays that need a slot for each scratchpad type.
 *
 * @return The number of valid scratchpad types (excluding _COUNT).
 *
 * @par Examples
 * @code{.cpp}
 * // Use at compile time for array sizing
 * constexpr size_t NUM_TYPES = get_max_scratchpad_types();
 * std::array<std::unique_ptr<AbstractScratchPad>, NUM_TYPES> scratchpads;
 *
 * // Use in static_assert
 * static_assert(get_max_scratchpad_types() > 0, "No scratchpad types defined");
 * @endcode
 */
consteval size_t get_max_scratchpad_types() {
    return static_cast<size_t>(ScratchPadType::_COUNT);
}

/**
 * @brief Convert a ScratchPadType to a human-readable name.
 *
 * Returns a string view containing the class name associated with
 * the given scratchpad type. Useful for logging, debugging, and
 * error messages.
 *
 * @param scratchpad_type The scratchpad type to convert.
 *
 * @return A string view containing the scratchpad class name, or
 *         "UnknownScratchPadType" for unrecognized values.
 *
 * @par Examples
 * @code{.cpp}
 * ScratchPadType type = ScratchPadType::GRAPH_ENGINE_SCRATCHPAD;
 * std::cout << "Using: " << get_scratchpad_type_name(type) << "\n";
 * // Output: "Using: GraphEngineScratchPad"
 *
 * // Use in error messages
 * throw std::runtime_error(
 *     std::format("Failed to initialize {}", get_scratchpad_type_name(type))
 * );
 * @endcode
 */
constexpr std::string_view get_scratchpad_type_name(const ScratchPadType scratchpad_type) {
    if constexpr (get_max_scratchpad_types() == 0) {
        return {""};
    }
    switch (scratchpad_type) {
        case ScratchPadType::GRAPH_ENGINE_SCRATCHPAD:
            return "GraphEngineScratchPad";
        case ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD:
            return "MultiscalePartitioningEngineViewScratchPad";
        case ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD:
            return "AdaptiveEngineViewScratchPad";
        case ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD:
            return "DefinedEngineViewScratchPad";
        case ScratchPadType::PRIMING_ENGINE_VIEW_SCRATCHPAD:
            return "PrimingEngineViewScratchPad";
        default:
            return "UnknownScratchPadType";
    }
}

} // namespace gridfire::engine::scratch
