#pragma once

#include "gridfire/engine/types/building.h"
#include "gridfire/engine/types/reporting.h"

namespace gridfire {
    enum class EngineTypes {
        GRAPH_ENGINE,
        ADAPTIVE_ENGINE_VIEW,
        MULTISCALE_PARTITIONING_ENGINE_VIEW,
        PRIMING_ENGINE_VIEW,
        DEFINED_ENGINE_VIEW,
        FILE_DEFINED_ENGINE_VIEW
    };

    inline std::string engine_type_to_string(const EngineTypes type) {
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