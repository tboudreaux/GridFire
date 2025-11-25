#include "gridfire/policy/stellar_policy.h"
#include "gridfire/policy/policy_abstract.h"
#include "gridfire/exceptions/error_policy.h"

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_views.h"

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

namespace {
    std::set<fourdst::atomic::Species> initialize_seed_species() {
        return {
            fourdst::atomic::H_1,
            fourdst::atomic::He_3,
            fourdst::atomic::He_4,
            fourdst::atomic::C_12,
            fourdst::atomic::N_14,
            fourdst::atomic::O_16,
            fourdst::atomic::Ne_20,
            fourdst::atomic::Mg_24
        };
    }
}

namespace gridfire::policy {
    MainSequencePolicy::MainSequencePolicy(const fourdst::composition::Composition& composition) : m_seed_species(initialize_seed_species()) {
        for (const auto& species : m_seed_species) {
            if (!composition.contains(species)) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize MainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }
        m_initializing_composition = composition;
        m_partition_function = build_partition_function();
    }

    MainSequencePolicy::MainSequencePolicy(std::vector<fourdst::atomic::Species> seed_species, const std::vector<double> &mass_fractions) {
        for (const auto& species : m_seed_species) {
            if (std::ranges::find(seed_species, species) == seed_species.end()) {
                throw exceptions::MissingSeedSpeciesError("Cannot initialize MainSequencePolicy: Required Seed species " + std::string(species.name()) + " is missing from the provided composition.");
            }
        }

        m_initializing_composition = fourdst::composition::buildCompositionFromMassFractions(seed_species, mass_fractions);

        m_partition_function = build_partition_function();
    }

    engine::DynamicEngine& MainSequencePolicy::construct() {
        m_network_stack.clear();

        m_network_stack.emplace_back(
            std::make_unique<engine::GraphEngine>(m_initializing_composition, *m_partition_function, engine::NetworkBuildDepth::ThirdOrder, engine::NetworkConstructionFlags::DEFAULT)
        );

        auto& graphRepr = dynamic_cast<engine::GraphEngine&>(*m_network_stack.back().get());
        graphRepr.setUseReverseReactions(false);


        m_network_stack.emplace_back(
            std::make_unique<engine::MultiscalePartitioningEngineView>(*m_network_stack.back().get())
        );
        m_network_stack.emplace_back(
            std::make_unique<engine::AdaptiveEngineView>(*m_network_stack.back().get())
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

    inline NetworkPolicyStatus MainSequencePolicy::get_status() const {
        return m_status;
    }

    const std::vector<std::unique_ptr<engine::DynamicEngine>> &MainSequencePolicy::get_engine_stack() const {
        if (m_status != NetworkPolicyStatus::INITIALIZED_VERIFIED) {
            throw exceptions::PolicyError("Cannot get engine stack from MainSequencePolicy: Policy is not initialized and verified. Call construct() first.");
        }
        return m_network_stack;
    }

    std::vector<engine::EngineTypes> MainSequencePolicy::get_engine_types_stack() const {
        return {
            engine::EngineTypes::GRAPH_ENGINE,
            engine::EngineTypes::MULTISCALE_PARTITIONING_ENGINE_VIEW,
            engine::EngineTypes::ADAPTIVE_ENGINE_VIEW
        };
    }

    const std::unique_ptr<partition::PartitionFunction>& MainSequencePolicy::get_partition_function() const {
        return m_partition_function;
    }

    inline NetworkPolicyStatus MainSequencePolicy::check_status() const {
        for (const auto& species : m_seed_species) {
            if (!m_initializing_composition.contains(species)) {
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