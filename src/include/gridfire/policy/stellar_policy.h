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
    class LowMassMainSequencePolicy final: public NetworkPolicy {
    public:
        explicit LowMassMainSequencePolicy(const fourdst::composition::Composition& composition);
        explicit LowMassMainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, std::vector<double> mass_fractions);

        std::string name() const override { return "LowMassMainSequencePolicy"; }

        const std::set<fourdst::atomic::Species> get_seed_species() const override { return m_seed_species; }
        const reaction::ReactionSet& get_seed_reactions() const override { return m_reaction_policy->get_reactions(); }

        DynamicEngine& construct() override;

        NetworkPolicyStatus getStatus() const override;
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