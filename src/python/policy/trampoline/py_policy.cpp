#include "py_policy.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fourdst/atomic/atomicSpecies.h"

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/engine.h"

#include "gridfire/policy/policy.h"

#include <string>
#include <set>

namespace py = pybind11;

std::string PyNetworkPolicy::name() const {
    PYBIND11_OVERRIDE_PURE(
        std::string,
        gridfire::policy::NetworkPolicy,
        name
    );
}

const std::set<fourdst::atomic::Species>& PyNetworkPolicy::get_seed_species() const {
    PYBIND11_OVERRIDE_PURE(
        const std::set<fourdst::atomic::Species>&,
        gridfire::policy::NetworkPolicy,
        get_seed_species
    );
}

const gridfire::reaction::ReactionSet& PyNetworkPolicy::get_seed_reactions() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::reaction::ReactionSet&,
        gridfire::policy::NetworkPolicy,
        get_seed_reactions
    );
}

gridfire::engine::DynamicEngine& PyNetworkPolicy::construct() {
    PYBIND11_OVERRIDE_PURE(
        gridfire::engine::DynamicEngine&,
        gridfire::policy::NetworkPolicy,
        construct
    );
}

gridfire::policy::NetworkPolicyStatus PyNetworkPolicy::get_status() const {
    PYBIND11_OVERRIDE_PURE(
        gridfire::policy::NetworkPolicyStatus,
        gridfire::policy::NetworkPolicy,
        getStatus
    );
}

const std::vector<std::unique_ptr<gridfire::engine::DynamicEngine>> &PyNetworkPolicy::get_engine_stack() const {
    PYBIND11_OVERRIDE_PURE(
        const std::vector<std::unique_ptr<gridfire::engine::DynamicEngine>> &,
        gridfire::policy::NetworkPolicy,
        get_engine_stack
    );
}

std::vector<gridfire::engine::EngineTypes> PyNetworkPolicy::get_engine_types_stack() const {
    PYBIND11_OVERRIDE_PURE(
        std::vector<gridfire::engine::EngineTypes>,
        gridfire::policy::NetworkPolicy,
        get_engine_types_stack
    );
}

const std::unique_ptr<gridfire::partition::PartitionFunction>& PyNetworkPolicy::get_partition_function() const {
    PYBIND11_OVERRIDE_PURE(
        const std::unique_ptr<gridfire::partition::PartitionFunction>&,
        gridfire::policy::NetworkPolicy,
        get_partition_function
    );
}

const gridfire::reaction::ReactionSet &PyReactionChainPolicy::get_reactions() const {
    PYBIND11_OVERRIDE_PURE(
        const gridfire::reaction::ReactionSet &,
        gridfire::policy::ReactionChainPolicy,
        get_reactions
    );
}

bool PyReactionChainPolicy::contains(const std::string &id) const {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::policy::ReactionChainPolicy,
        contains,
        id
    );
}

bool PyReactionChainPolicy::contains(const gridfire::reaction::Reaction &reaction) const {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::policy::ReactionChainPolicy,
        contains,
        reaction
    );
}

std::unique_ptr<gridfire::policy::ReactionChainPolicy> PyReactionChainPolicy::clone() const {
    PYBIND11_OVERRIDE_PURE(
        std::unique_ptr<gridfire::policy::ReactionChainPolicy>,
        gridfire::policy::ReactionChainPolicy,
        clone
    );
}

std::string PyReactionChainPolicy::name() const {
    PYBIND11_OVERRIDE_PURE(
        std::string,
        gridfire::policy::ReactionChainPolicy,
        name
    );
}

uint64_t PyReactionChainPolicy::hash(uint64_t seed) const {
    PYBIND11_OVERRIDE_PURE(
        uint64_t,
        gridfire::policy::ReactionChainPolicy,
        hash,
        seed
    );
}

bool PyReactionChainPolicy::operator==(const ReactionChainPolicy &other) const {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::policy::ReactionChainPolicy,
        operator==,
        other
    );
}

bool PyReactionChainPolicy::operator!=(const ReactionChainPolicy &other) const {
    PYBIND11_OVERRIDE_PURE(
        bool,
        gridfire::policy::ReactionChainPolicy,
        operator!=,
        other
    );
}
