#pragma once

#include "gridfire/engine/scratchpads/blob.h"
#include "gridfire/engine/scratchpads/scratchpad_abstract.h"
#include "gridfire/engine/scratchpads/engine_graph_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_adaptive_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_multiscale_scratchpad.h"
#include "gridfire/engine/scratchpads/engine_defined_scratchpad.h"
#include "gridfire/engine/scratchpads/types.h"
#include "gridfire/utils/logging.h"

#include <format>
#include <string>
#include <string_view>

// 1. ScratchPadType: Inherit from string_view formatter for efficiency
template <>
struct std::formatter<gridfire::engine::scratch::ScratchPadType> : std::formatter<std::string_view> {
    // Note: NOT static, marked const
    auto format(gridfire::engine::scratch::ScratchPadType type, auto& ctx) const {
        // Convert to string_view
        std::string_view name = gridfire::engine::scratch::get_scratchpad_type_name(type);
        // Delegate to the base class to handle width/fill/alignment
        return std::formatter<std::string_view>::format(name, ctx);
    }
};

// 2. AbstractScratchPad
template <>
struct std::formatter<gridfire::engine::scratch::AbstractScratchPad> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::AbstractScratchPad& pad, auto& ctx) const {
        std::string str = std::format("AbstractScratchPad(Initialized: {})",
                           pad.is_initialized());
        return std::formatter<std::string>::format(str, ctx);
    }
};

// 3. GraphEngineScratchPad
template <>
struct std::formatter<gridfire::engine::scratch::GraphEngineScratchPad> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::GraphEngineScratchPad& pad, auto& ctx) const {
        std::string str = std::format("GraphEngineScratchPad(Initialized: {}, HasADFun: {}, CachedStepDerivatives: {}, CachedJacobians: {})",
                           pad.has_initialized,
                           pad.rhsADFun.has_value(),
                           pad.stepDerivativesCache.size(),
                           pad.jacobianSubsetCache.size());
        return std::formatter<std::string>::format(str, ctx);
    }
};

// 4. AdaptiveEngineViewScratchPad
template<>
struct std::formatter<gridfire::engine::scratch::AdaptiveEngineViewScratchPad> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::AdaptiveEngineViewScratchPad& pad, auto& ctx) const {
        std::string str = std::format("AdaptiveEngineViewScratchPad(Initialized: {}, Active Species: {}, Active Reactions: {})",
                           pad.has_initialized,
                           pad.active_species.size(),
                           pad.active_reactions.size());
        return std::formatter<std::string>::format(str, ctx);
    }
};

// 5. MultiscalePartitioningEngineViewScratchPad
template <>
struct std::formatter<gridfire::engine::scratch::MultiscalePartitioningEngineViewScratchPad> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::MultiscalePartitioningEngineViewScratchPad& pad, auto& ctx) const {
        std::string str = std::format("MultiscalePartitioningEngineViewScratchPad(Initialized: {}, QSE Groups: {}, Dynamic Species: {}, Algebraic Species: {}, Cached Compositions: {})",
                           pad.has_initialized,
                           pad.qse_groups.size(),
                           pad.dynamic_species.size(),
                           pad.algebraic_species.size(),
                           pad.composition_cache.size());
        return std::formatter<std::string>::format(str, ctx);
    }
};

// 6. DefinedEngineViewScratchPad
template <>
struct std::formatter<gridfire::engine::scratch::DefinedEngineViewScratchPad> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::DefinedEngineViewScratchPad& pad, auto& ctx) const {
        std::string str = std::format("DefinedEngineViewScratchPad(Initialized: {}, Active Species: {}, Active Reactions: {})",
                           pad.has_initialized,
                           pad.active_species.size(),
                           pad.active_reactions.size());
        return std::formatter<std::string>::format(str, ctx);
    }
};

// 7. StateBlob
template <>
struct std::formatter<gridfire::engine::scratch::StateBlob> : std::formatter<std::string> {
    auto format(const gridfire::engine::scratch::StateBlob& blob, auto& ctx) const {
        // Construct the full string representation
        std::string str = std::format("StateBlob(Enrolled: {})",
            gridfire::utils::iterable_to_delimited_string(
                blob.get_registered_scratchpads(),
                ", ",
                [&blob](const auto& type) {
                    auto result = blob.get(type);
                    if (!result.has_value()) {
                        return std::format("{}(Error: {})",
                            gridfire::engine::scratch::get_scratchpad_type_name(type),
                            gridfire::engine::scratch::StateBlob::error_to_string(result.error()));
                    }

                    gridfire::engine::scratch::AbstractScratchPad* scratchpad = result.value();

                    // We can reuse the formatters we defined above by dereferencing the cast pointers!
                    switch (type) {
                        case gridfire::engine::scratch::ScratchPadType::GRAPH_ENGINE_SCRATCHPAD : {
                            auto* cast_pad = dynamic_cast<gridfire::engine::scratch::GraphEngineScratchPad*>(scratchpad);
                            // This works because we defined a formatter for GraphEngineScratchPad above
                            return std::format("{}", *cast_pad);
                        }
                        case gridfire::engine::scratch::ScratchPadType::MULTISCALE_PARTITIONING_ENGINE_VIEW_SCRATCHPAD : {
                            auto* cast_pad = dynamic_cast<gridfire::engine::scratch::MultiscalePartitioningEngineViewScratchPad*>(scratchpad);
                            return std::format("{}", *cast_pad);
                        }
                        case gridfire::engine::scratch::ScratchPadType::ADAPTIVE_ENGINE_VIEW_SCRATCHPAD : {
                            auto* cast_pad = dynamic_cast<gridfire::engine::scratch::AdaptiveEngineViewScratchPad*>(scratchpad);
                            return std::format("{}", *cast_pad);
                        }
                        case gridfire::engine::scratch::ScratchPadType::DEFINED_ENGINE_VIEW_SCRATCHPAD : {
                            auto* cast_pad = dynamic_cast<gridfire::engine::scratch::DefinedEngineViewScratchPad*>(scratchpad);
                            return std::format("{}", *cast_pad);
                        }
                        default: {
                            return std::format("{}(Unknown ScratchPad Type)", gridfire::engine::scratch::get_scratchpad_type_name(type));
                        }
                    }
                }
            )
        );

        // Delegate to base class to write the string to the context
        return std::formatter<std::string>::format(str, ctx);
    }
};