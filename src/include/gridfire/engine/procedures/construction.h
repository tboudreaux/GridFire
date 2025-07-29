#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/types/building.h"

#include "fourdst/composition/composition.h"

#include <variant>

namespace gridfire {

    /**
     * @brief Builds a nuclear reaction network from the Reaclib library based on an initial composition.
     *
     * Constructs a layered reaction network by collecting reactions up to the specified depth
     * from the Reaclib dataset. Starting species are those with non-zero mass fractions in the input
     * composition. Layers expand by including products of collected reactions until the depth limit.
     * Optionally selects reverse reactions instead of forward.
     *
     * See implementation in construction.cpp for details on the layering algorithm, logging, and performance.
     *
     * @param composition Mapping of isotopic species to their mass fractions; species with positive
     *        mass fraction seed the network.
     * @param maxLayers Variant specifying either a predefined NetworkBuildDepth or a custom integer depth;
     *        negative depth (Full) collects all reactions, zero is invalid.
     * @param reverse If true, collects reverse reactions (decays or back-reactions); if false, uses forward reactions.
     * @pre composition must have at least one species with positive mass fraction.
     * @pre Resolved integer depth from maxLayers must not be zero.
     * @post Returned network includes only reactions satisfying the depth and reverse criteria.
     * @return A LogicalReactionSet encapsulating the collected reactions for graph-based engines.
     * @throws std::logic_error If the resolved network depth is zero (no reactions can be collected).
     */
    reaction::LogicalReactionSet build_reaclib_nuclear_network(
        const fourdst::composition::Composition &composition,
        BuildDepthType maxLayers = NetworkBuildDepth::Full,
        bool reverse = false
    );
}