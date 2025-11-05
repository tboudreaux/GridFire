/**
 * @file policy_abstract.h
 * @brief Abstract policy interfaces used to construct reaction networks (DynamicEngine) from seed compositions.
 *
 * This header declares the base interfaces for high-level "policies" that drive automatic network
 * construction and verification. A NetworkPolicy encapsulates the information required to build
 * a production-ready reaction network from an initial composition (seed species and seed reactions)
 * and to return a constructed DynamicEngine ready for use by the solvers.
 *
 * Concrete policy implementations live in this directory and provide the real behaviour. See for
 * example:
 *  - gridfire/policy/stellar_policy.h  (concrete stellar network policy used in many examples)
 *  - gridfire/policy/chains.h          (reaction-chain helper policies such as proton-proton, CNO)
 *
 * An example of using a concrete policy to construct and run an engine is available at:
 *  - tests/graphnet_sandbox/main.cpp
 *
 * @note Doxygen comments on public methods include @par Example usage blocks. Methods that may
 *       throw in concrete implementations include @throws tags.
 */

#pragma once

#include "fourdst/composition/atomicSpecies.h"
#include "gridfire/engine/types/building.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine_abstract.h"

#include <string>
#include <set>


namespace gridfire::policy {
    /**
     * @enum NetworkPolicyStatus
     * @brief Lifecycle/verification status for a NetworkPolicy instance.
     *
     * Used by concrete NetworkPolicy implementations to report the state of the policy and the
     * constructed network during/after construction. Tests and callers can inspect the status
     * to determine whether the resulting network meets the policy's requirements.
     */
    enum class NetworkPolicyStatus {
        UNINITIALIZED,
        INITIALIZED_UNVERIFIED,
        MISSING_KEY_REACTION,
        MISSING_KEY_SPECIES,
        INITIALIZED_VERIFIED
    };

    /**
     * @class NetworkPolicy
     * @brief Abstract interface for policies that construct DynamicEngine networks from a seed composition.
     *
     * A NetworkPolicy provides three main pieces of information:
     *  - A name identifying the policy.
     *  - A set of seed species and seed reactions required to initialize the network.
     *  - A constructor method that returns a fully constructed DynamicEngine (or view stack) built
     *    to satisfy the policy.
     *
     * Concrete implementations include `LowMassMainSequencePolicy` (see `stellar_policy.h`) and may
     * throw policy-specific exceptions during construction (for example when required reactions or
     * species are missing).
     *
     * @par Example
     * @code
     * // Example usage (see tests/graphnet_sandbox/main.cpp for a complete example):
     * // fourdst::composition::Composition comp = prepared composition
     * // gridfire::policy::LowMassMainSequencePolicy policy(comp);
     * // // construct() returns a reference to a DynamicEngine (could be a view stack)
     * // DynamicEngine &engine = policy.construct();
     * // // engine can now be passed to a solver (e.g. CVODESolverStrategy)
     * @endcode
     */
    class NetworkPolicy {
    public:
        virtual ~NetworkPolicy() = default;

        /**
         * @brief Human-readable name for the policy.
         *
         * @return a std::string identifying the policy implementation (e.g. "LowMassMainSequencePolicy").
         *
         * @par Example
         * @code
         * std::string n = policy.name();
         * std::cout << "Using policy: " << n << std::endl;
         * @endcode
         */
        [[nodiscard]] virtual std::string name() const = 0;

        /**
         * @brief Returns the seed species the policy requires to initialize the network.
         *
         * The returned set contains atomic species identifiers (fourdst::atomic::Species) which the
         * policy expects to be present in the initial composition used to build the network.
         *
         * Implementations should return a copy or an immutable reference to their internal set of
         * required seed species.
         *
         * @par Example
         * @code
         * const auto seeds = policy.get_seed_species();
         * for (const auto &s : seeds) { std::cout << s.name() << std::endl; }
         * @endcode
         */
        [[nodiscard]] virtual const std::set<fourdst::atomic::Species> get_seed_species() const = 0;

        /**
         * @brief Returns the set of seed reactions the policy requires.
         *
         * The ReactionSet describes reactions that must be present in the constructed network for the
         * policy to be considered satisfied. Concrete policies often implement their reaction
         * requirements by composing one or more ReactionChainPolicy instances (see `chains.h`).
         *
         * @par Example
         * @code
         * const reaction::ReactionSet &reacs = policy.get_seed_reactions();
         * // inspect reaction IDs or count
         * std::cout << "Policy requires " << reacs.size() << " reactions" << std::endl;
         * @endcode
         */
        [[nodiscard]] virtual const reaction::ReactionSet& get_seed_reactions() const = 0;

        /**
         * @brief Construct and return a DynamicEngine instance (or engine view stack) satisfying the policy.
         *
         * Implementations typically build one or more engine layers (GraphEngine, MultiscalePartitioningEngineView,
         * AdaptiveEngineView, etc.) and return a reference to the top-most DynamicEngine. The storage lifetime of
         * the returned reference is implementation-defined (usually owned by the policy instance).
         *
         * @return DynamicEngine& reference to a running/constructed engine ready for solver consumption.
         *
         * @throws gridfire::exceptions::MissingKeyReactionError if required reactions are not present in the
         *         constructed network (see `gridfire/exceptions/error_policy.h`).
         * @throws gridfire::exceptions::MissingSeedSpeciesError if required seed species are missing from the
         *         initializing composition.
         * @throws gridfire::exceptions::PolicyError for other construction/verification failures.
         *
         * @par Example
         * @code
         * DynamicEngine &engine = policy.construct();
         * solver::CVODESolverStrategy solver(engine);
         * NetOut out = solver.evaluate(netIn, true);
         * @endcode
         */
        [[nodiscard]] virtual DynamicEngine& construct() = 0;

        /**
         * @brief Returns the current verification/construction status of the policy.
         *
         * The status reports whether the policy has been initialized and whether the constructed
         * network satisfies the policy's key requirements.
         *
         * @par Example
         * @code
         * NetworkPolicyStatus s = policy.getStatus();
         * if (s != NetworkPolicyStatus::INITIALIZED_VERIFIED) { // handle error }
         * @endcode
         */
        [[nodiscard]] virtual NetworkPolicyStatus getStatus() const = 0;
    };

    /**
     * @class ReactionChainPolicy
     * @brief Abstract interface encapsulating a set of reactions representing a single chain or pathway.
     *
     * ReactionChainPolicy implementations (see `chains.h`) supply a ReactionSet describing the reactions
     * that comprise a nuclear reaction chain (for example the proton-proton chain, CNO cycle, etc.).
     *
     * @par Example
     * @code
     * ProtonProtonChainPolicy pp;
     * const auto &reacs = pp.get_reactions();
     * @endcode
     *
     * @note Concrete implementations may throw exceptions on construction if the underlying reaction
     *       database (e.g. REACLIB) does not include requested reactions. See `chains.h` for details.
     */
    class ReactionChainPolicy {
    public:
        virtual ~ReactionChainPolicy() = default;

        /**
         * @brief Returns the ReactionSet describing this chain.
         *
         * @return const reaction::ReactionSet& reference to the chain's reactions.
         *
         * @par Example
         * @code
         * const reaction::ReactionSet &set = chainPolicy.get_reactions();
         * std::cout << "Chain contains " << set.size() << " reactions\n";
         * @endcode
         *
         * @throws gridfire::exceptions::MissingBaseReactionError may be thrown by concrete implementations
         *         at construction time if the required reactions cannot be found in the base reaction set.
         */
        [[nodiscard]] virtual const reaction::ReactionSet& get_reactions() const = 0;
    };

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
         * @brief Returns the vector of child ReactionChainPolicy instances.
         * @return const std::vector<std::unique_ptr<ReactionChainPolicy>>&
         */
        [[nodiscard]] virtual const std::vector<std::unique_ptr<ReactionChainPolicy>>& get_chain_policies() const = 0;
    };
}
