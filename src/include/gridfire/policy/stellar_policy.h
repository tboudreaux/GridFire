/**
 * @file stellar_policy.h
 * @brief High-level concrete NetworkPolicy for specific stellar environments.
 *
 * This file defines the `MainSequencePolicy`, a concrete implementation of the `NetworkPolicy`
 * interface. This policy is designed to construct a reaction network suitable for simulating
 * nucleosynthesis in low-mass main-sequence stars (like the Sun).
 *
 * It serves as a practical example of how to create a high-level policy by:
 * 1.  Defining a set of required seed species.
 * 2.  Composing lower-level `ReactionChainPolicy` instances (from `chains.h`) to specify
 *     the necessary reactions (e.g., PP-chain and CNO-cycle).
 * 3.  Implementing the `construct()` method to build a sophisticated, multi-layered `DynamicEngine`
 *     stack, ready for use with a solver.
 *
 */
#pragma once

#include <vector>
#include <set>

#include "gridfire/policy/policy_abstract.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/reaction/reaction.h"


#include "fourdst/composition/composition.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "gridfire/partition/composite/partition_composite.h"

#include "gridfire/policy/chains.h"

namespace gridfire::policy {
    /**
     * @class MainSequencePolicy
     * @brief A NetworkPolicy for building reaction networks suitable for low-mass main-sequence stars.
     *
     * This policy ensures that a constructed network contains all necessary species and reactions
     * for modeling the core hydrogen burning phase in low-mass stars, primarily the PP-chain and
     * the CNO-cycle.
     *
     * The `construct()` method builds a stack of engine views (`GraphEngine` ->
     * `MultiscalePartitioningEngineView` -> `AdaptiveEngineView`) to provide a dynamically
     * adjusting, performance-optimized reaction network.
     *
     * This policy requires the following seed species:
     *  - H-1
     *  - He-3
     *  - He-4
     *  - C-12
     *  - N-14
     *  - O-16
     *  - Ne-20
     *  - Mg-24
     *
     *
     * This policy composes the `ProtonProtonChainPolicy` and `CNOChainPolicy` to define
     * the required reactions.
     */
    class MainSequencePolicy final: public NetworkPolicy {
    public:
        /**
         * @brief Constructs the policy from an existing composition object.
         *
         * @param composition The initial composition, which must contain all seed species required by the policy.
         *
         * @throws exceptions::MissingSeedSpeciesError if the provided `composition` is missing a required species.
         *
         * @par Example
         * @code
         * fourdst::composition::Composition comp;
         * // ... populate composition ...
         * LowMassMainSequencePolicy policy(comp);
         * @endcode
         */
        explicit MainSequencePolicy(const fourdst::composition::Composition& composition);
        /**
         * @brief Constructs the policy from a list of species and their mass fractions.
         *
         * @param seed_species A vector of atomic species.
         * @param mass_fractions A vector of corresponding mass fractions.
         *
         * @throws exceptions::MissingSeedSpeciesError if the provided species list is missing a required seed species.
         * @throws fourdst::composition::exceptions::CompositionNotFinalizedError if the internal composition fails to finalize.
         *
         * @par Example
         * @code
         * using namespace fourdst::atomic;
         * std::vector<Species> species = {H_1, He_4, C_12, O_16};
         * std::vector<double> mass_fractions = {0.7, 0.28, 0.01, 0.01};
         * LowMassMainSequencePolicy policy(species, mass_fractions);
         * @endcode
         */
        explicit MainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, const std::vector<double> &mass_fractions);

        /**
         * @brief Returns the name of the policy.
         * @return "MainSequencePolicy"
         */
        [[nodiscard]] std::string name() const override { return "MainSequencePolicy"; }

        /**
         * @brief Returns the set of seed species required by this policy.
         * @return const std::set<fourdst::atomic::Species>&
         */
        [[nodiscard]] const std::set<fourdst::atomic::Species>& get_seed_species() const override { return m_seed_species; }
        /**
         * @brief Returns the set of seed reactions required by this policy (from the PP and CNO chains).
         * @return const reaction::ReactionSet&
         */
        [[nodiscard]] const reaction::ReactionSet& get_seed_reactions() const override { return m_reaction_policy->get_reactions(); }

        /**
         * @brief Constructs and returns the complete, multi-layered dynamic engine.
         *
         * This method builds the full network engine stack:
         * - A base `GraphEngine` is created with the initial composition. This is constructed three layers deep which is sufficient to capture all required reactions. Further a
         *   composite ground state and Rauscher Thielemann partition function is used.
         * - A `MultiscalePartitioningEngineView` is layered on top for performance optimization. This will put some species into equilibrium groups based on their reaction timescales.
         * - An `AdaptiveEngineView` is added as the final layer to dynamically cull reaction pathways based on molar reaction flows.
         *
         * After construction, it verifies that the resulting network meets the policy's requirements.
         *
         * @return DynamicEngine& A reference to the top-level `AdaptiveEngineView`.
         *
         * @throws exceptions::MissingKeyReactionError if the final network is missing a required reaction.
         * @throws exceptions::MissingSeedSpeciesError if the final network is missing a required species.
         * @throws exceptions::PolicyError on other verification failures.
         *
         * @par Example
         * @code
         * LowMassMainSequencePolicy enginePolicy(composition);
         * DynamicEngine& engine = enginePolicy.construct();
         * solver::CVODESolverStrategy solver(engine);
         * // ... run solver ...
         * @endcode
         */
        engine::DynamicEngine& construct() override;

        /**
         * @brief Gets the current status of the policy.
         * @return NetworkPolicyStatus The construction and verification status.
         */
        [[nodiscard]] NetworkPolicyStatus get_status() const override;

        [[nodiscard]] const std::vector<std::unique_ptr<engine::DynamicEngine>> &get_engine_stack() const override;

        [[nodiscard]] std::vector<engine::EngineTypes> get_engine_types_stack() const override;
        [[nodiscard]] const std::unique_ptr<partition::PartitionFunction>& get_partition_function() const override;

    private:
        std::set<fourdst::atomic::Species> m_seed_species; ///< The set of seed species required by this policy. These are H-1, He-3, He-4, C-12, N-14, O-16, Ne-20, Mg-24.

        std::unique_ptr<ReactionChainPolicy> m_reaction_policy = std::make_unique<MainSequenceReactionChainPolicy>(); ///< The composed reaction chain policy (PP-chain + CNO-cycle).
        fourdst::composition::Composition m_initializing_composition; ///< The initial composition used to build the network.
        std::unique_ptr<partition::PartitionFunction> m_partition_function; ///< The partition function used in network construction.
        std::vector<std::unique_ptr<engine::DynamicEngine>> m_network_stack; ///< The stack of dynamic engines constructed by the policy.

        NetworkPolicyStatus m_status = NetworkPolicyStatus::UNINITIALIZED; ///< The current status of the policy.
    private:
        static std::unique_ptr<partition::PartitionFunction> build_partition_function();
        [[nodiscard]] NetworkPolicyStatus check_status() const;

    public:

    };



}