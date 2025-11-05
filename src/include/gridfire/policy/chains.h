#pragma once

#include "gridfire/policy/policy_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "gridfire/reaction/reaclib.h"

#include "gridfire/exceptions/error_policy.h"

#include <memory>

namespace gridfire::policy {

    class ProtonProtonChainPolicy final: public ReactionChainPolicy {
    public:
        ProtonProtonChainPolicy();

        const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
    private:
        std::vector<std::string> m_reactionIDs = {
            "p(p,e+)d",
            "d(p,g)he3",
            "he3(he3,2p)he4",
            "he4(he3,g)be7",
            "be7(e-,)li7",
            "li7(p,a)he4",
            "be7(p,g)b8",
            "b8(,e+)be8",
            "be8(,a)he4"
        };
        reaction::ReactionSet m_reactions;
    };

    class CNOChainPolicy final: public ReactionChainPolicy {
    public:
        CNOChainPolicy();
        const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
    private:
        std::set<std::string> m_reactionIDs = {
            "c12(p,g)n13",
            "n13(,e+)c13",
            "c13(p,g)n14",
            "n14(p,g)o15",
            "o15(,e+)n15",
            "n15(p,a)c12",

            "n15(p,g)o16",
            "o16(p,g)f17",
            "f17(,e+)o17",
            "o17(p,a)n14",
            "n14(p,g)o15",
            "o15(,e+)n15",

            "o17(p,g)f18",
            "f18(,e+)o18",
            "o18(p,a)n15",
            "n15(p,g)o16",
            "o16(p,g)f17",
            "f17(,e+)o17",

            "o18(p,g)f19",
            "f19(p,a)o16",
            "o16(p,g)f17",
            "f17(,e+)o17",
            "o17(p,g)f18",
            "f18(,e+)o18"
        };
        reaction::ReactionSet m_reactions;
    };

    class HotCNOChainPolicy final : public ReactionChainPolicy {
    public:
        HotCNOChainPolicy();
        const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
    private:
        std::set<std::string> m_reactionIDs = {
            "c12(p,g)n13",
            "n13(p,g)o14",
            "o14(,e+)n14",
            "n14(p,g)o15",
            "o15(,e+)n15",
            "n15(p,a)c12",

            "n15(p,g)o16",
            "o16(p,g)f17",
            "f17(p,g)ne18",
            "ne18(,e+)f18",
            "f18(p,a)o15",
            "o15(,e+)n15",

            "f18(p,g)ne19",
            "ne19(,e+)f19",
            "f19(p,a)o16",
            "o16(p,g)f17",
            "f17(p,g)ne18",
            "ne18(,e+)f18"
        };

        reaction::ReactionSet m_reactions;
    };

    class LowMassMainSequenceReactionChainPolicy final : public MultiReactionChainPolicy {
    public:
        LowMassMainSequenceReactionChainPolicy();

        const reaction::ReactionSet & get_reactions() const override;

        const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const override;

    private:
        std::vector<std::unique_ptr<ReactionChainPolicy>> m_chain_policies;
        reaction::ReactionSet m_reactions;
    };

    inline ProtonProtonChainPolicy::ProtonProtonChainPolicy() {
        const auto& all_reaclib_reactions = reaclib::get_all_reaclib_reactions();

        for (const auto& reactionID : m_reactionIDs) {
            auto reaction = all_reaclib_reactions.get(reactionID);
            if (!reaction) {
                throw exceptions::MissingBaseReactionError("The Underlying REACLIB reaction set is missing the reaction " + std::string(reactionID) + " needed for the proton-proton chain. This indicates that there is an issue with the GridFire binary you are using. Please try to recompile and if that fails please report this issue to the developers.");
            }
            m_reactions.add_reaction(reaction.value()->clone());
        }
    }

    inline CNOChainPolicy::CNOChainPolicy() {
        const auto& all_reaclib_reactions = reaclib::get_all_reaclib_reactions();
        for (const auto& reactionID : m_reactionIDs) {
            auto reaction = all_reaclib_reactions.get(reactionID);
            if (!reaction) {
                throw exceptions::MissingBaseReactionError("The Underlying REACLIB reaction set is missing the reaction " + std::string(reactionID) + " needed for the CNO cycle. This indicates that there is an issue with the GridFire binary you are using. Please try to recompile and if that fails please report this issue to the developers.");
            }
            m_reactions.add_reaction(reaction.value()->clone());
        }
    }

    inline HotCNOChainPolicy::HotCNOChainPolicy() {
        const auto& all_reaclib_reactions = reaclib::get_all_reaclib_reactions();
        for (const auto& reactionID : m_reactionIDs) {
            auto reaction = all_reaclib_reactions.get(reactionID);
            if (!reaction) {
                throw exceptions::MissingBaseReactionError("The Underlying REACLIB reaction set is missing the reaction " + std::string(reactionID) + " needed for the Hot CNO cycle. This indicates that there is an issue with the GridFire binary you are using. Please try to recompile and if that fails please report this issue to the developers.");
            }
            m_reactions.add_reaction(reaction.value()->clone());
        }
    }


    inline LowMassMainSequenceReactionChainPolicy::LowMassMainSequenceReactionChainPolicy() {
        m_chain_policies.emplace_back(std::make_unique<ProtonProtonChainPolicy>());
        m_chain_policies.emplace_back(std::make_unique<CNOChainPolicy>());
        for (const auto& policy_ptr : m_chain_policies) {
            m_reactions.extend(policy_ptr->get_reactions());
        }
    }

    inline const reaction::ReactionSet & LowMassMainSequenceReactionChainPolicy::get_reactions() const {
        return m_reactions;
    }

    inline const std::vector<std::unique_ptr<ReactionChainPolicy>>& LowMassMainSequenceReactionChainPolicy::get_chain_policies() const {
        return m_chain_policies;
    }



}

