#pragma once

#include "gridfire/policy/policy_abstract.h"

#include <vector>
#include <memory>

namespace gridfire::policy {
    /**
     * @class MultiReactionChainPolicy
     * @brief A ReactionChainPolicy composed of multiple child ReactionChainPolicy instances.
     *
     * Useful for policies that represent a union of several reaction chains (for example the
     * `LowMassMainSequenceReactionChainPolicy` composes the proton-proton and CNO chains).
     *
     * @par Example
     * @code
     * LowMassMainSequenceReactionChainPolicy multi;
     * const auto &chains = multi.get_chain_policies();
     * for (const auto &ch : chains) { std::cout << ch->get_reactions().size() << " reactions\n"; }
     * @endcode
     */
    class MultiReactionChainPolicy : public ReactionChainPolicy {
    public:
        explicit MultiReactionChainPolicy(std::vector<std::unique_ptr<ReactionChainPolicy>>&& chain_policies);

        [[nodiscard]] const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const;

        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override;

        [[nodiscard]] bool contains(const std::string &id) const override;
        [[nodiscard]] bool contains(const reaction::Reaction &reaction) const override;

        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;
        [[nodiscard]] std::string name() const override;
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;

        [[nodiscard]] bool operator==(const ReactionChainPolicy& other) const override;
        [[nodiscard]] bool operator!=(const ReactionChainPolicy& other) const override;

        [[nodiscard]] size_t size() const;

        auto begin() { return m_chain_policies.begin(); }
        [[nodiscard]] auto begin() const { return m_chain_policies.cbegin(); }
        auto end() { return m_chain_policies.end(); }
        [[nodiscard]] auto end() const { return m_chain_policies.cend(); }
    protected:
        std::vector<std::unique_ptr<ReactionChainPolicy>> m_chain_policies;
        reaction::ReactionSet m_reactions;
    };
}