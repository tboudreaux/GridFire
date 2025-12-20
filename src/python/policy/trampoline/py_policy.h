#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gridfire/policy/policy.h"

class PyNetworkPolicy final : public gridfire::policy::NetworkPolicy {
public:
    [[nodiscard]] std::string name() const override;

    [[nodiscard]] const std::set<fourdst::atomic::Species>& get_seed_species() const override;

    [[nodiscard]] const gridfire::reaction::ReactionSet& get_seed_reactions() const override;

    [[nodiscard]] gridfire::policy::ConstructionResults construct() override;

    [[nodiscard]] gridfire::policy::NetworkPolicyStatus get_status() const override;

    [[nodiscard]] const std::vector<std::unique_ptr<gridfire::engine::DynamicEngine>> &get_engine_stack() const override;

    [[nodiscard]] std::vector<gridfire::engine::EngineTypes> get_engine_types_stack() const override;

    [[nodiscard]] const std::unique_ptr<gridfire::partition::PartitionFunction>& get_partition_function() const override;

    [[nodiscard]] std::unique_ptr<gridfire::engine::scratch::StateBlob> get_stack_scratch_blob() const override;
};

class PyReactionChainPolicy final : public gridfire::policy::ReactionChainPolicy {
public:
    [[nodiscard]] const gridfire::reaction::ReactionSet & get_reactions() const override;

    [[nodiscard]] bool contains(const std::string &id) const override;

    [[nodiscard]] bool contains(const gridfire::reaction::Reaction &reaction) const override;

    [[nodiscard]] std::unique_ptr<ReactionChainPolicy> clone() const override;

    [[nodiscard]] std::string name() const override;

    [[nodiscard]] uint64_t hash(uint64_t seed) const override;

    [[nodiscard]] bool operator==(const ReactionChainPolicy &other) const override;

    [[nodiscard]] bool operator!=(const ReactionChainPolicy &other) const override;
};