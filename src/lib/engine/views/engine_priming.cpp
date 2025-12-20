#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/engine/scratchpads/blob.h"

#include "fourdst/atomic/species.h"

#include "quill/LogMacros.h"
#include "quill/Logger.h"

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

namespace gridfire::engine {
    using fourdst::atomic::species;

    NetworkPrimingEngineView::NetworkPrimingEngineView(
        scratch::StateBlob& ctx,
        const std::string &primingSymbol,
        GraphEngine &baseEngine
    ) :
    DefinedEngineView(
        constructPrimingReactionSet(
            ctx,
            species.at(primingSymbol),
            baseEngine
        ),
        baseEngine
    ),
    m_primingSpecies(species.at(primingSymbol)) {}

    NetworkPrimingEngineView::NetworkPrimingEngineView(
        scratch::StateBlob& ctx,
        const fourdst::atomic::Species &primingSpecies,
        GraphEngine &baseEngine
    ) :
    DefinedEngineView(
        constructPrimingReactionSet(
            ctx,
            primingSpecies,
            baseEngine
        ),
        baseEngine
    ),
    m_primingSpecies(primingSpecies) {}

    std::vector<std::string> NetworkPrimingEngineView::constructPrimingReactionSet(
        scratch::StateBlob& ctx,
        const fourdst::atomic::Species &primingSpecies,
        const GraphEngine &baseEngine
    ) const {
        std::unordered_set<std::string> primeReactions;
        for (const auto &reaction : baseEngine.getNetworkReactions(ctx)) {
            if (reaction->contains(primingSpecies)) {
                primeReactions.insert(std::string(reaction->id()));
            }
        }
        if (primeReactions.empty()) {
            LOG_INFO(m_logger, "No priming reactions found for species '{}', returning empty peName set.", primingSpecies.name());
            return std::vector<std::string>{};
            // m_logger->flush_log();
            // throw std::runtime_error("No priming reactions found for species '" + std::string(primingSpecies.name()) + "'.");
        }
        std::vector<std::string> primingReactionSet(primeReactions.begin(), primeReactions.end());
        // LOG_INFO(m_logger, "Constructed priming reaction set with {} reactions for species '{}'.", primingReactionSet.size(), primingSpecies.name());
        return primingReactionSet;
    }



}
