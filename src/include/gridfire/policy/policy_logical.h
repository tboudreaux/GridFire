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
        /**
         * @brief Constructs a MultiReactionChainPolicy from a vector of ReactionChainPolicy instances.
         *
         * The provided chain policies are moved into the new MultiReactionChainPolicy instance.
         *
         * @param chain_policies vector of unique_ptr to ReactionChainPolicy instances.
         *
         * @par Example
         * @code
         * std::vector<std::unique_ptr<ReactionChainPolicy>> chains;
         * chains.push_back(std::make_unique<ProtonProtonChainPolicy>());
         * chains.push_back(std::make_unique<CNOChainPolicy>());
         * MultiReactionChainPolicy multi(std::move(chains));
         * @endcode
         */
        explicit MultiReactionChainPolicy(std::vector<std::unique_ptr<ReactionChainPolicy>>&& chain_policies);

        /**
         * @brief Returns the vector of child ReactionChainPolicy instances.
         *
         * @return const std::vector<std::unique_ptr<ReactionChainPolicy>>& reference to the child chain policies.
         *
         * @par Example
         * @code
         * const auto &chains = multi.get_chain_policies();
         * for (const auto &ch : chains) { std::cout << ch->get_reactions().size() << " reactions\n"; }
         * @endcode
         */
        [[nodiscard]] const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const;

        /**
         * @brief Returns the combined ReactionSet of all child chain policies.
         *
         * @return const reaction::ReactionSet& reference to the combined reactions.
         *
         * @par Example
         * @code
         * const auto &reactions = multi.get_reactions();
         * std::cout << "Multi chain contains " << reactions.size() << " reactions\n";
         * @endcode
         */
        [[nodiscard]] const reaction::ReactionSet& get_reactions() const override;

        /**
         * @brief Checks if the MultiReactionChainPolicy contains a reaction by ID
         *
         * @param id the reaction ID to check for.
         * @return true if the reaction ID is present in the combined ReactionSet, false otherwise.
         *
         * @par Example
         * @code
         * bool has_pp = multi.contains("p(p,e+)d");
         * @endcode
         */
        [[nodiscard]] bool contains(const std::string &id) const override;

        /**
         * @brief Checks if the MultiReactionChainPolicy contains a specific reaction.
         *
         * @param reaction the Reaction to check for.
         * @return true if the reaction is present in the combined ReactionSet, false otherwise.
         *
         * @par Example
         * @code
         * reaction::Reaction r = ...; // obtain a Reaction instance
         * bool has_reaction = multi.contains(r);
         * @endcode
         */
        [[nodiscard]] bool contains(const reaction::Reaction &reaction) const override;

        /**
         * @brief Creates a deep copy of the MultiReactionChainPolicy.
         *
         * @return std::unique_ptr<ReactionChainPolicy> unique pointer to the cloned instance.
         *
         * @par Example
         * @code
         * std::unique_ptr<ReactionChainPolicy> clone = multi.clone();
         * @endcode
         */
        [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

        /**
         * @brief Returns the name of the MultiReactionChainPolicy.
         *
         * @return std::string the name of the policy.
         *
         * @par Example
         * @code
         * std::string n = multi.name();
         * std::cout << "Using policy: " << n << std::endl;
         * @endcode
         */
        [[nodiscard]] std::string name() const override;

        /**
         * @brief Computes a hash value for the MultiReactionChainPolicy.
         *
         * @param seed the seed value for the hash computation.
         * @return uint64_t the computed hash value.
         *
         * @par Example
         * @code
         * uint64_t h = multi.hash(0);
         * std::cout << "Policy hash: " << h << std::endl;
         * @endcode
         */
        [[nodiscard]] uint64_t hash(uint64_t seed) const override;

        /**
         * @brief Equality comparison operator.
         *
         * @param other the other ReactionChainPolicy to compare against.
         * @return true if the policies are equal, false otherwise.
         *
         * @par Example
         * @code
         * bool equal = (multi1 == multi2);
         * @endcode
         */
        [[nodiscard]] bool operator==(const ReactionChainPolicy& other) const override;

        /**
         * @brief Inequality comparison operator.
         *
         * @param other the other ReactionChainPolicy to compare against.
         * @return true if the policies are not equal, false otherwise.
         *
         * @par Example
         * @code
         * bool not_equal = (multi1 != multi2);
         * @endcode
         */
        [[nodiscard]] bool operator!=(const ReactionChainPolicy& other) const override;

        /**
         * @brief Returns the number of child ReactionChainPolicy instances.
         *
         * @return size_t the number of child chain policies.
         *
         * @par Example
         * @code
         * size_t n = multi.size();
         * std::cout << "Multi chain has " << n << " child policies\n";
         * @endcode
         */
        [[nodiscard]] size_t size() const;

        /**
         * @brief Returns iterator to the beginning of the child ReactionChainPolicy instances.
         *
         * @return iterator to the beginning.
         *
         * @par Example
         * @code
         * for (auto it = multi.begin(); it != multi.end(); ++it) {
         *     std::cout << (*it)->name() << std::endl;
         * }
         * @endcode
         */
        auto begin() { return m_chain_policies.begin(); }

        /**
         * @brief Returns const iterator to the beginning of the child ReactionChainPolicy instances.
         *
         * @return const iterator to the beginning.
         *
         * @par Example
         * @code
         * for (auto it = multi.begin(); it != multi.end(); ++it) {
         *     std::cout << (*it)->name() << std::endl;
         * }
         * @endcode
         */
        [[nodiscard]] auto begin() const { return m_chain_policies.cbegin(); }

        /**
         * @brief Returns iterator to the end of the child ReactionChainPolicy instances.
         *
         * @return iterator to the end.
         *
         * @par Example
         * @code
         * for (auto it = multi.begin(); it != multi.end(); ++it) {
         *     std::cout << (*it)->name() << std::endl;
         * }
         * @endcode
         */
        auto end() { return m_chain_policies.end(); }

        /**
         * @brief Returns const iterator to the end of the child ReactionChainPolicy instances.
         *
         * @return const iterator to the end.
         *
         * @par Example
         * @code
         * for (auto it = multi.begin(); it != multi.end(); ++it) {
         *     std::cout << (*it)->name() << std::endl;
         * }
         * @endcode
         */
        [[nodiscard]] auto end() const { return m_chain_policies.cend(); }
    protected:
        std::vector<std::unique_ptr<ReactionChainPolicy>> m_chain_policies{}; ///< Child chain policies
        reaction::ReactionSet m_reactions; ///< Combined reactions from all child policies
    };
}