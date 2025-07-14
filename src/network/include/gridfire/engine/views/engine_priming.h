#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_defined.h"

#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"

#include "quill/Logger.h"

#include <vector>
#include <string>

namespace gridfire {

    class NetworkPrimingEngineView final : public DefinedEngineView {
    public:
        NetworkPrimingEngineView(const std::string& primingSymbol, DynamicEngine& baseEngine);
        NetworkPrimingEngineView(const fourdst::atomic::Species& primingSpecies, DynamicEngine& baseEngine);


    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        fourdst::atomic::Species m_primingSpecies; ///< The priming species, if specified.
    private:
        std::vector<std::string> constructPrimingReactionSet(
            const fourdst::atomic::Species& primingSpecies,
            const DynamicEngine& baseEngine
        ) const;
    };

}