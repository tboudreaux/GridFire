#include "gridfire/policy/stellar_policy.h"
#include "gridfire/policy/policy_abstract.h"
#include "gridfire/exceptions/error_policy.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_views.h"

#include "fourdst/composition/exceptions/exceptions_composition.h"

namespace gridfire::policy {
    MainSequencePolicy::MainSequencePolicy(const fourdst::composition::Composition& composition) {
        for (const auto& species : m_seed_species) {
            if (!composition.hasSpecies(species)) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize MainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }
        m_initializing_composition = composition;
        m_partition_function = build_partition_function();
    }

    MainSequencePolicy::MainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, std::vector<double> mass_fractions) {
        for (const auto& species : m_seed_species) {
            if (std::ranges::find(seed_species, species) == seed_species.end()) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize MainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }

        for (const auto& [species, x] : std::views::zip(seed_species, mass_fractions)) {
            m_initializing_composition.registerSpecies(species);
            m_initializing_composition.setMassFraction(species, x);
        }

        const bool didFinalize = m_initializing_composition.finalize(true);
        if (!didFinalize) {
            throw fourdst::composition::exceptions::CompositionNotFinalizedError("Failed to finalize initial composition for MainSequencePolicy.");
        }

        m_partition_function = build_partition_function();
    }

    DynamicEngine& MainSequencePolicy::construct() {
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
                throw exceptions::MissingKeyReactionError("MainSequencePolicy construction failed: The constructed network is missing key reactions required by the policy.");
            case NetworkPolicyStatus::MISSING_KEY_SPECIES:
                throw exceptions::MissingSeedSpeciesError("MainSequencePolicy construction failed: The constructed network is missing key seed species required by the policy.");
            case NetworkPolicyStatus::UNINITIALIZED:
                throw exceptions::PolicyError("MainSequencePolicy construction failed: The network policy is uninitialized.");
            case NetworkPolicyStatus::INITIALIZED_UNVERIFIED:
                throw exceptions::PolicyError("MainSequencePolicy construction failed: The network policy status could not be verified.");
            case NetworkPolicyStatus::INITIALIZED_VERIFIED:
                break;
        }
        return *m_network_stack.back();
    }

    inline std::unique_ptr<partition::PartitionFunction> MainSequencePolicy::build_partition_function() {
        using partition::BasePartitionType;
        const auto partitionFunction = partition::CompositePartitionFunction({
            BasePartitionType::RauscherThielemann,
            BasePartitionType::GroundState
        });
        return std::make_unique<partition::CompositePartitionFunction>(partitionFunction);
    }

    inline NetworkPolicyStatus MainSequencePolicy::getStatus() const {
        return m_status;
    }

    inline NetworkPolicyStatus MainSequencePolicy::check_status() const {
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