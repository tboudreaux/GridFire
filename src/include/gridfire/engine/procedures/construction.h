#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/types/building.h"

#include "fourdst/composition/composition.h"

#include <variant>

namespace gridfire {

    reaction::LogicalReactionSet build_reaclib_nuclear_network(
        const fourdst::composition::Composition &composition,
        BuildDepthType maxLayers = NetworkBuildDepth::Full,
        bool reverse = false
    );
}