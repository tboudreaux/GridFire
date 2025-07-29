#include "gridfire/engine/views/engine_priming.h"
#include "gridfire/solver/solver.h"

#include "fourdst/composition/species.h"
#include "fourdst/logging/logging.h"

#include "quill/LogMacros.h"
#include "quill/Logger.h"

#include <vector>
#include <string>
#include <unordered_set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <ranges>
#include <cmath>


namespace gridfire {
    using fourdst::atomic::species;

    NetworkPrimingEngineView::NetworkPrimingEngineView(
        const std::string &primingSymbol,
        DynamicEngine &baseEngine
    ) :
    DefinedEngineView(
        constructPrimingReactionSet(
            species.at(primingSymbol),
            baseEngine
        ),
        baseEngine
    ),
    m_primingSpecies(species.at(primingSymbol)) {}

    NetworkPrimingEngineView::NetworkPrimingEngineView(
        const fourdst::atomic::Species &primingSpecies,
        DynamicEngine &baseEngine
    ) :
    DefinedEngineView(
        constructPrimingReactionSet(
            primingSpecies,
            baseEngine
        ),
        baseEngine
    ),
    m_primingSpecies(primingSpecies) {
    }


    std::vector<std::string> NetworkPrimingEngineView::constructPrimingReactionSet(
        const fourdst::atomic::Species &primingSpecies,
        const DynamicEngine &baseEngine
    ) const {
        std::unordered_set<std::string> primeReactions;
        for (const auto &reaction : baseEngine.getNetworkReactions()) {
            if (reaction.contains(primingSpecies)) {
                primeReactions.insert(std::string(reaction.peName()));
            }
        }
        if (primeReactions.empty()) {
            LOG_ERROR(m_logger, "No priming reactions found for species '{}'.", primingSpecies.name());
            m_logger->flush_log();
            throw std::runtime_error("No priming reactions found for species '" + std::string(primingSpecies.name()) + "'.");
        }
        std::vector<std::string> primingReactionSet(primeReactions.begin(), primeReactions.end());
        // LOG_INFO(m_logger, "Constructed priming reaction set with {} reactions for species '{}'.", primingReactionSet.size(), primingSpecies.name());
        return primingReactionSet;
    }



}
