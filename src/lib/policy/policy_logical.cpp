#include "gridfire/policy/policy_abstract.h"
#include "gridfire/policy/policy_logical.h"

#include "xxhash64.h"
#include "gridfire/utils/hashing.h"

namespace gridfire::policy {
    MultiReactionChainPolicy::MultiReactionChainPolicy(
        std::vector<std::unique_ptr<ReactionChainPolicy>> &&chain_policies
    ) : m_chain_policies(std::move(chain_policies)) {
        for (const auto &ch : m_chain_policies) {
            m_reactions.extend(ch->get_reactions());
        }
    }

    const std::vector<std::unique_ptr<ReactionChainPolicy>> & MultiReactionChainPolicy::get_chain_policies() const {
        return m_chain_policies;
    }

    const reaction::ReactionSet & MultiReactionChainPolicy::get_reactions() const {
        return m_reactions;
    }

    bool MultiReactionChainPolicy::contains(const std::string &id) const {
        return m_reactions.contains(id);
    }

    bool MultiReactionChainPolicy::contains(const reaction::Reaction &reaction) const {
        return m_reactions.contains(reaction);
    }

    std::unique_ptr<ReactionChainPolicy> MultiReactionChainPolicy::clone() const {
        return std::make_unique<MultiReactionChainPolicy>(
            [this]() {
                std::vector<std::unique_ptr<ReactionChainPolicy>> chain_policies;
                for (const auto &ch : m_chain_policies) {
                    chain_policies.push_back(ch->clone());
                }
                return chain_policies;
            }()
        );
    }

    std::string MultiReactionChainPolicy::name() const {
        return "MultiReactionChainPolicy";
    }

    uint64_t MultiReactionChainPolicy::hash(const uint64_t seed) const {
        std::vector<size_t> reaction_hashes;
        for (const auto& reaction : m_reactions) {
            reaction_hashes.push_back(reaction->hash(seed));
        }

        return XXHash64::hash(reaction_hashes.data(), reaction_hashes.size(), seed);
    }

    bool MultiReactionChainPolicy::operator==(const ReactionChainPolicy &other) const {
        return this->hash(0) == other.hash(0);
    }

    bool MultiReactionChainPolicy::operator!=(const ReactionChainPolicy &other) const {
        return this->hash(0) != other.hash(0);
    }


    size_t MultiReactionChainPolicy::size() const {
        return m_chain_policies.size();
    }
}
