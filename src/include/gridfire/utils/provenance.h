#pragma once

#include "gridfire/config/config.h"
#include "gridfire/engine/engine_abstract.h"

#include <cstdint>
#include <string>
#include <string_view>

namespace gridfire {
    /// GridFire's version, e.g. "1.2.0".
    [[nodiscard]] std::string version_string();

    /// Order-independent XXH64 of the reaction ids in a network: two networks with the same reactions agree.
    [[nodiscard]] uint64_t network_hash(const reaction::ReactionSet& reactions);

    /**
     * @brief A description of the physics a GridFire calculation uses, for a host code to record with its results.
     *
     * One "key: value" line per item: GridFire version, network policy, network size and hash, the rate libraries
     * and the hashes of their compiled data, screening, Jacobian method, passive species, energy channels and solver
     * tolerances. The network lines are omitted when engine is null (no network built yet).
     *
     * Two runs that print the same description used the same code version, rates, network and physics options.
     */
    [[nodiscard]] std::string describe_provenance(
        const config::GridFireConfig& config,
        const engine::DynamicEngine* engine = nullptr,
        engine::scratch::StateBlob* ctx = nullptr,
        std::string_view policy_name = {}
    );
}
