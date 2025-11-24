#pragma once

namespace gridfire::engine {
    /**
     * @enum EngineTypes
     * @brief Enumeration of different engine types available in GridFire.
     *
     * Values:
     *   - GRAPH_ENGINE: The standard graph-based engine.
     *   - ADAPTIVE_ENGINE_VIEW: An engine that adapts based on certain criteria.
     *   - MULTISCALE_PARTITIONING_ENGINE_VIEW: An engine that partitions the system at multiple scales.
     *   - PRIMING_ENGINE_VIEW: An engine that uses a priming strategy for simulations.
     *   - DEFINED_ENGINE_VIEW: An engine defined by user specifications.
     *   - FILE_DEFINED_ENGINE_VIEW: An engine defined through external files.
     */
    enum class EngineTypes {
        GRAPH_ENGINE,
        ADAPTIVE_ENGINE_VIEW,
        MULTISCALE_PARTITIONING_ENGINE_VIEW,
        PRIMING_ENGINE_VIEW,
        DEFINED_ENGINE_VIEW,
        FILE_DEFINED_ENGINE_VIEW
    };

    /**
     * @brief Converts an EngineTypes enum value to its corresponding string representation.
     *
     * @param type The EngineTypes enum value to convert.
     * @return A string_view representing the name of the engine type.
     */
    constexpr std::string_view engine_type_to_string(const EngineTypes type) {
        switch (type) {
            case EngineTypes::GRAPH_ENGINE:
                return "GraphEngine";
            case EngineTypes::ADAPTIVE_ENGINE_VIEW:
                return "AdaptiveEngineView";
            case EngineTypes::MULTISCALE_PARTITIONING_ENGINE_VIEW:
                return "MultiscalePartitioningEngineView";
            case EngineTypes::PRIMING_ENGINE_VIEW:
                return "PrimingEngineView";
            case EngineTypes::DEFINED_ENGINE_VIEW:
                return "DefinedEngineView";
            case EngineTypes::FILE_DEFINED_ENGINE_VIEW:
                return "FileDefinedEngineView";
        }
        return "UnknownEngineType";
    }
}