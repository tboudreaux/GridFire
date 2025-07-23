#include "gridfire/engine/procedures/construction.h"

#include <ranges>
#include <stdexcept>

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/reaclib.h"

#include "fourdst/composition/composition.h"

#include "fourdst/logging/logging.h"

#include "quill/Logger.h"
#include "quill/LogMacros.h"

namespace gridfire {
    using reaction::LogicalReactionSet;
    using reaction::ReactionSet;
    using reaction::Reaction;
    using fourdst::composition::Composition;
    using fourdst::atomic::Species;


    LogicalReactionSet build_reaclib_nuclear_network(
        const Composition &composition,
        BuildDepthType maxLayers,
        bool reverse
    ) {
        int depth;
        if (std::holds_alternative<NetworkBuildDepth>(maxLayers)) {
            depth = static_cast<int>(std::get<NetworkBuildDepth>(maxLayers));
        } else {
            depth = std::get<int>(maxLayers);
        }
        auto logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        if (depth == 0) {
            LOG_ERROR(logger, "Network build depth is set to 0. No reactions will be collected.");
            throw std::logic_error("Network build depth is set to 0. No reactions will be collected.");
        }

        const auto allReactions = reaclib::get_all_reactions();
        std::vector<Reaction> remainingReactions;
        for (const auto& reaction : allReactions) {
            if (reaction.is_reverse() == reverse) {
                remainingReactions.push_back(reaction);
            }
        }

        if (depth == static_cast<int>(NetworkBuildDepth::Full)) {
            LOG_INFO(logger, "Building full nuclear network with a total of {} reactions.", allReactions.size());
            const ReactionSet reactionSet(remainingReactions);
            return reaction::packReactionSetToLogicalReactionSet(reactionSet);
        }

        std::unordered_set<Species> availableSpecies;
        for (const auto &entry: composition | std::views::values) {
            if (entry.mass_fraction() > 0.0) {
                availableSpecies.insert(entry.isotope());
            }
        }


        std::vector<Reaction> collectedReactions;

        LOG_INFO(logger, "Starting network construction with {} available species.", availableSpecies.size());
        for (int layer = 0; layer < depth && !remainingReactions.empty(); ++layer) {
            LOG_TRACE_L1(logger, "Collecting reactions for layer {} with {} remaining reactions. Currently there are {} available species", layer, remainingReactions.size(), availableSpecies.size());
            std::vector<Reaction> reactionsForNextPass;
            std::unordered_set<Species> newProductsThisLayer;
            bool newReactionsAdded = false;

            reactionsForNextPass.reserve(remainingReactions.size());

            for (const auto &reaction : remainingReactions) {
                bool allReactantsAvailable = true;
                for (const auto& reactant : reaction.reactants()) {
                    if (!availableSpecies.contains(reactant)) {
                        allReactantsAvailable = false;
                        break;
                    }
                }

                if (allReactantsAvailable) {
                    collectedReactions.push_back(reaction);
                    newReactionsAdded = true;

                    for (const auto& product : reaction.products()) {
                        newProductsThisLayer.insert(product);
                    }
                } else {
                    reactionsForNextPass.push_back(reaction);
                }
            }

            if (!newReactionsAdded) {
                LOG_INFO(logger, "No new reactions added in layer {}. Stopping network construction with {} reactions collected.", layer, collectedReactions.size());
                break;
            }

            LOG_TRACE_L1(logger, "Layer {}: Collected {} reactions. New products this layer: {}", layer, collectedReactions.size(), newProductsThisLayer.size());
            availableSpecies.insert(newProductsThisLayer.begin(), newProductsThisLayer.end());

            remainingReactions = std::move(reactionsForNextPass);
        }

        LOG_INFO(logger, "Network construction completed with {} reactions collected.", collectedReactions.size());
        const ReactionSet reactionSet(collectedReactions);
        return reaction::packReactionSetToLogicalReactionSet(reactionSet);
    }



}