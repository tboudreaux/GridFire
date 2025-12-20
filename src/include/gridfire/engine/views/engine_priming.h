#pragma once

#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/engine/scratchpads/blob.h"


#include "fourdst/logging/logging.h"
#include "fourdst/atomic/atomicSpecies.h"


#include "quill/Logger.h"

#include <vector>
#include <string>

namespace gridfire::engine {

    /**
     * @class NetworkPrimingEngineView
     * @brief Provides a view of a DynamicEngine filtered to reactions involving a specified priming species.
     *
     * This view constructs a subset of the network reactions from the base engine that
     * contain the given priming species and delegates all engine operations to the underlying engine.
     *
     * See implementation in engine_priming.cpp for details on reaction set construction.
     *
     * @note Throws std::runtime_error if no priming reactions are found for the species.
     */
    class NetworkPrimingEngineView final : public DefinedEngineView {
    public:
        /**
         * @brief Constructs the view by looking up the priming species by symbol.
         *
         * @param primingSymbol Symbol string of the species to prime.
         * @param baseEngine Reference to the base DynamicEngine to wrap.
         * @pre primingSymbol must correspond to a valid species in atomic::species registry.
         * @post The view will contain only reactions that involve the priming species.
         * @throws std::out_of_range If primingSymbol is not found in the species registry.
         * @throws std::runtime_error If no reactions contain the priming species.
         */
        NetworkPrimingEngineView(
            scratch::StateBlob& ctx,
            const std::string& primingSymbol,
            GraphEngine& baseEngine
        );
        /**
         * @brief Constructs the view using an existing Species object.
         *
         * @param primingSpecies The species object to prime.
         * @param baseEngine Reference to the base DynamicEngine to wrap.
         * @pre primingSpecies must be valid and present in the network of baseEngine.
         * @post The view will contain only reactions that involve the priming species.
         * @throws std::runtime_error If no reactions contain the priming species.
         */
        NetworkPrimingEngineView(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species& primingSpecies,
            GraphEngine& baseEngine
        );


    private:
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        fourdst::atomic::Species m_primingSpecies; ///< The priming species, if specified.
    private:
        /**
         * @brief Constructs the set of reaction names that involve the priming species.
         *
         * @param primingSpecies Species for which to collect priming reactions.
         * @param baseEngine Base engine containing the full network of reactions.
         * @pre baseEngine.getNetworkReactions() returns a valid iterable set of reactions.
         * @post Returns a vector of unique reaction name strings containing the priming species.
         * @return Vector of reaction name strings containing the priming species.
         * @throws std::runtime_error If no reactions involve the priming species.
         */
        [[nodiscard]] std::vector<std::string> constructPrimingReactionSet(
            scratch::StateBlob& ctx,
            const fourdst::atomic::Species& primingSpecies,
            const GraphEngine& baseEngine
        ) const;
    };

}