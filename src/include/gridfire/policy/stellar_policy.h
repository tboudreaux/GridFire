/**
 * @file stellar_policy.h
 * @brief High-level concrete NetworkPolicy for specific stellar environments.
 *
 * This file defines the `LowMassMainSequencePolicy`, a concrete implementation of the `NetworkPolicy`
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

#include "gridfire/exceptions/error_policy.h"

#include "fourdst/composition/composition.h"
#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/partition/composite/partition_composite.h"

#include "gridfire/policy/chains.h"

namespace gridfire::policy {
    /**
     * @class LowMassMainSequencePolicy
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
    class LowMassMainSequencePolicy final: public NetworkPolicy {
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
        explicit LowMassMainSequencePolicy(const fourdst::composition::Composition& composition);
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
        explicit LowMassMainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, std::vector<double> mass_fractions);

        /**
         * @brief Returns the name of the policy.
         * @return "LowMassMainSequencePolicy"
         */
        [[nodiscard]] std::string name() const override { return "LowMassMainSequencePolicy"; }

        /**
         * @brief Returns the set of seed species required by this policy.
         * @return const std::set<fourdst::atomic::Species>&
         */
        [[nodiscard]] const std::set<fourdst::atomic::Species> get_seed_species() const override { return m_seed_species; }
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
        DynamicEngine& construct() override;

        /**
         * @brief Gets the current status of the policy.
         * @return NetworkPolicyStatus The construction and verification status.
         */
        [[nodiscard]] NetworkPolicyStatus getStatus() const override;
    private:
        std::set<fourdst::atomic::Species> m_seed_species = {
            fourdst::atomic::H_1,
            fourdst::atomic::He_3,
            fourdst::atomic::He_4,
            fourdst::atomic::C_12,
            fourdst::atomic::N_14,
            fourdst::atomic::O_16,
            fourdst::atomic::Ne_20,
            fourdst::atomic::Mg_24
        };

        std::unique_ptr<ReactionChainPolicy> m_reaction_policy = std::make_unique<LowMassMainSequenceReactionChainPolicy>();
        fourdst::composition::Composition m_initializing_composition;
        std::unique_ptr<partition::PartitionFunction> m_partition_function;
        std::vector<std::unique_ptr<DynamicEngine>> m_network_stack;

        NetworkPolicyStatus m_status = NetworkPolicyStatus::UNINITIALIZED;
    private:
        static std::unique_ptr<partition::PartitionFunction> build_partition_function();
        NetworkPolicyStatus check_status() const;

    };

    inline LowMassMainSequencePolicy::LowMassMainSequencePolicy(const fourdst::composition::Composition& composition) {
        for (const auto& species : m_seed_species) {
            if (!composition.hasSpecies(species)) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize LowMassMainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }
        m_initializing_composition = composition;
        m_partition_function = build_partition_function();
    }

    inline LowMassMainSequencePolicy::LowMassMainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, std::vector<double> mass_fractions) {
        for (const auto& species : m_seed_species) {
            if (std::ranges::find(seed_species, species) == seed_species.end()) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize LowMassMainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }

        for (const auto& [species, x] : std::views::zip(seed_species, mass_fractions)) {
            m_initializing_composition.registerSpecies(species);
            m_initializing_composition.setMassFraction(species, x);
        }

        const bool didFinalize = m_initializing_composition.finalize(true);
        if (!didFinalize) {
            throw fourdst::composition::exceptions::CompositionNotFinalizedError("Failed to finalize initial composition for LowMassMainSequencePolicy.");
        }

        m_partition_function = build_partition_function();
    }

    inline DynamicEngine& LowMassMainSequencePolicy::construct() {
        m_network_stack.clear();

        m_network_stack.emplace_back(
            std::make_unique<GraphEngine>(m_initializing_composition, *m_partition_function, NetworkBuildDepth::ThirdOrder, NetworkConstructionFlags::DEFAULT)
        );
        m_network_stack.emplace_back(
            std::make_unique<MultiscalePartitioningEngineView>(*m_network_stack.back().get())
        );
        m_network_stack.emplace_back(
            std::make_unique<AdaptiveEngineView>(*m_network_stack.back().get())
        );

        m_status = NetworkPolicyStatus::INITIALIZED_UNVERIFIED;
        m_status = check_status();

        switch (m_status) {
            case NetworkPolicyStatus::MISSING_KEY_REACTION:
                throw exceptions::MissingKeyReactionError("LowMassMainSequencePolicy construction failed: The constructed network is missing key reactions required by the policy.");
            case NetworkPolicyStatus::MISSING_KEY_SPECIES:
                throw exceptions::MissingSeedSpeciesError("LowMassMainSequencePolicy construction failed: The constructed network is missing key seed species required by the policy.");
            case NetworkPolicyStatus::UNINITIALIZED:
                throw exceptions::PolicyError("LowMassMainSequencePolicy construction failed: The network policy is uninitialized.");
            case NetworkPolicyStatus::INITIALIZED_UNVERIFIED:
                throw exceptions::PolicyError("LowMassMainSequencePolicy construction failed: The network policy status could not be verified.");
            case NetworkPolicyStatus::INITIALIZED_VERIFIED:
                break;
        }
        return *m_network_stack.back();
    }

    inline std::unique_ptr<partition::PartitionFunction> LowMassMainSequencePolicy::build_partition_function() {
        using partition::BasePartitionType;
        const auto partitionFunction = partition::CompositePartitionFunction({
            BasePartitionType::RauscherThielemann,
            BasePartitionType::GroundState
        });
        return std::make_unique<partition::CompositePartitionFunction>(partitionFunction);
    }

    inline NetworkPolicyStatus LowMassMainSequencePolicy::getStatus() const {
        return m_status;
    }

    inline NetworkPolicyStatus LowMassMainSequencePolicy::check_status() const {
        for (const auto& species : m_seed_species) {
            if (!m_initializing_composition.hasSpecies(species)) {
                return NetworkPolicyStatus::MISSING_KEY_SPECIES;
            }
        }
        const reaction::ReactionSet& baseReactions = m_network_stack.front()->getNetworkReactions();
        for (const auto& reaction : m_reaction_policy->get_reactions()) {
            const bool result = baseReactions.contains(*reaction);
            if (!result) {
                return NetworkPolicyStatus::MISSING_KEY_REACTION;
            }
        }
        return NetworkPolicyStatus::INITIALIZED_VERIFIED;
    }


}