/**
 * @file chains.h
 * @brief Concrete implementations of ReactionChainPolicy for key stellar reaction chains.
 *
 * This file provides concrete policies for fundamental nuclear reaction chains, such as the
 * Proton-Proton chain and the CNO cycle. These policies inherit from ReactionChainPolicy (see
 * `policy_abstract.h`) and provide a pre-defined set of reactions.
 *
 * They are typically used by higher-level NetworkPolicy implementations (e.g., `LowMassMainSequencePolicy`
 * in `stellar_policy.h`) to compose a complete set of required reactions for a particular
 * stellar environment.
 *
 * @note Constructors for these policies may throw `gridfire::exceptions::MissingBaseReactionError`
 *       if a required reaction is not found in the underlying REACLIB database. This usually
 *       indicates an issue with the GridFire installation or the bundled reaction data.
 */
#pragma once

#include "gridfire/policy/policy_abstract.h"
#include "gridfire/reaction/reaction.h"

#include "gridfire/reaction/reaclib.h"

#include "gridfire/exceptions/error_policy.h"

#include <memory>

namespace gridfire::policy {

    /**
     * @class ProtonProtonChainPolicy
     * @brief A ReactionChainPolicy for the Proton-Proton (PP) chain.
     *
     * Encapsulates the set of reactions that constitute the three branches of the PP chain,
     * which is the primary energy generation mechanism in stars like the Sun.
     *
     * @throws gridfire::exceptions::MissingBaseReactionError if a required reaction for the PP chain
     *         is not found in the REACLIB database during construction.
     */
    class ProtonProtonChainPolicy final: public ReactionChainPolicy {
    public:
        /**
         * @brief Constructs the policy and initializes its reaction set from REACLIB.
         */
        ProtonProtonChainPolicy();

        /**
         * @brief Returns the set of reactions in the PP chain.
         * @return const reaction::ReactionSet&
         *
         * @par Example
         * @code
         * ProtonProtonChainPolicy pp_policy;
         * const auto& reactions = pp_policy.get_reactions();
         * std::cout << "PP chain has " << reactions.size() << " reactions." << std::endl;
         * @endcode
         */
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
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

    /**
     * @class CNOChainPolicy
     * @brief A ReactionChainPolicy for the Carbon-Nitrogen-Oxygen (CNO) cycle.
     *
     * Encapsulates the reactions of the CNO cycle, a catalytic cycle that is the dominant
     * source of energy in massive stars.
     *
     * @throws gridfire::exceptions::MissingBaseReactionError if a required reaction for the CNO cycle
     *         is not found in the REACLIB database during construction.
     */
    class CNOChainPolicy final: public ReactionChainPolicy {
    public:
        /**
         * @brief Constructs the policy and initializes its reaction set from REACLIB.
         */
        CNOChainPolicy();
        /**
         * @brief Returns the set of reactions in the CNO cycle.
         * @return const reaction::ReactionSet&
         *
         * @par Example
         * @code
         * CNOChainPolicy cno_policy;
         * const auto& reactions = cno_policy.get_reactions();
         * assert(reactions.contains("c12(p,g)n13"));
         * @endcode
         */
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
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

    /**
     * @class HotCNOChainPolicy
     * @brief A ReactionChainPolicy for the Hot CNO (HCNO) cycle.
     *
     * Encapsulates the reactions of the HCNO cycle, which becomes significant at higher
     * temperatures and densities than the standard CNO cycle, often in explosive scenarios.
     *
     * @throws gridfire::exceptions::MissingBaseReactionError if a required reaction for the HCNO cycle
     *         is not found in the REACLIB database during construction.
     */
    class HotCNOChainPolicy final : public ReactionChainPolicy {
    public:
        /**
         * @brief Constructs the policy and initializes its reaction set from REACLIB.
         */
        HotCNOChainPolicy();
        /**
         * @brief Returns the set of reactions in the HCNO cycle.
         * @return const reaction::ReactionSet&
         */
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override { return m_reactions; }
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

    /**
     * @class LowMassMainSequenceReactionChainPolicy
     * @brief A MultiReactionChainPolicy for low-mass main-sequence stars.
     *
     * This policy composes the `ProtonProtonChainPolicy` and `CNOChainPolicy` to represent the
     * key energy-generating reaction chains active in low-mass stars like the Sun.
     */
    class LowMassMainSequenceReactionChainPolicy final : public MultiReactionChainPolicy {
    public:
        /**
         * @brief Constructs the policy and initializes its child policies.
         */
        LowMassMainSequenceReactionChainPolicy();

        /**
         * @brief Returns the combined set of reactions from all child policies (PP and CNO).
         * @return const reaction::ReactionSet&
         */
        [[nodiscard]] const reaction::ReactionSet & get_reactions() const override;

        /**
         * @brief Returns the vector of child policies.
         * @return const std::vector<std::unique_ptr<ReactionChainPolicy>>&
         *
         * @par Example
         * @code
         * LowMassMainSequenceReactionChainPolicy lmms_policy;
         * const auto& child_policies = lmms_policy.get_chain_policies();
         * std::cout << "Low-mass policy has " << child_policies.size() << " child policies." << std::endl;
         * @endcode
         */
        [[nodiscard]] const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const override;

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
