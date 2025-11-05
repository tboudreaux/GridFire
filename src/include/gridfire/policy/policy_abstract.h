#pragma once

#include "fourdst/composition/atomicSpecies.h"
#include "gridfire/engine/types/building.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"

#include <string>
#include <set>


namespace gridfire::policy {
    enum class NetworkPolicyStatus {
        UNINITIALIZED,
        INITIALIZED_UNVERIFIED,
        MISSING_KEY_REACTION,
        MISSING_KEY_SPECIES,
        INITIALIZED_VERIFIED
    };

    class NetworkPolicy {
    public:
        virtual ~NetworkPolicy() = default;
        virtual std::string name() const = 0;

        virtual const std::set<fourdst::atomic::Species> get_seed_species() const = 0;
        virtual const reaction::ReactionSet& get_seed_reactions() const = 0;

        virtual DynamicEngine& construct() = 0;

        virtual NetworkPolicyStatus getStatus() const = 0;
    };

    class ReactionChainPolicy {
    public:
        virtual ~ReactionChainPolicy() = default;
        virtual const reaction::ReactionSet& get_reactions() const = 0;
    };

    class MultiReactionChainPolicy : public ReactionChainPolicy {
    public:
        virtual const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const = 0;
    };

}
