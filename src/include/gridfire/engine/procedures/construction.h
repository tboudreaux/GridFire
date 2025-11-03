#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/types/building.h"

#include "fourdst/composition/composition.h"

#include <variant>

#include "gridfire/reaction/weak/weak_interpolator.h"

namespace gridfire {


    enum class NetworkConstructionFlags : uint32_t {
        NONE = 0,

        STRONG = 1 << 0, // 1

        BETA_MINUS = 1 << 1, // 2
        BETA_PLUS = 1 << 2, // 4
        ELECTRON_CAPTURE = 1 << 3, // 8
        POSITRON_CAPTURE = 1 << 4, // 16

        WEAK = BETA_MINUS | BETA_PLUS | ELECTRON_CAPTURE | POSITRON_CAPTURE,
        DEFAULT = STRONG,
        ALL = STRONG | WEAK
    };

    constexpr auto to_underlying(NetworkConstructionFlags f) noexcept {
        return static_cast<std::underlying_type_t<NetworkConstructionFlags>>(f);
    }

    inline NetworkConstructionFlags operator|(const NetworkConstructionFlags lhs, const NetworkConstructionFlags rhs) {
        return static_cast<NetworkConstructionFlags>(to_underlying(lhs) | to_underlying(rhs));
    }

    inline NetworkConstructionFlags operator&(const NetworkConstructionFlags lhs, const NetworkConstructionFlags rhs) {
        return static_cast<NetworkConstructionFlags>(to_underlying(lhs) & to_underlying(rhs));
    }

    inline bool has_flag(const NetworkConstructionFlags flags, const NetworkConstructionFlags flag_to_check) {
        return (flags & flag_to_check) != NetworkConstructionFlags::NONE;
    }

    inline std::string NetworkConstructionFlagsToString(NetworkConstructionFlags flags) {
        std::stringstream ss;
        constexpr std::array<NetworkConstructionFlags, 5> bases_flags_array = {
            NetworkConstructionFlags::STRONG,
            NetworkConstructionFlags::BETA_MINUS,
            NetworkConstructionFlags::BETA_PLUS,
            NetworkConstructionFlags::ELECTRON_CAPTURE,
            NetworkConstructionFlags::POSITRON_CAPTURE
        };

        const std::unordered_map<NetworkConstructionFlags, std::string> bases_string_map = {
            {NetworkConstructionFlags::STRONG, "Strong"},
            {NetworkConstructionFlags::BETA_MINUS, "BetaMinus"},
            {NetworkConstructionFlags::BETA_PLUS, "BetaPlus"},
            {NetworkConstructionFlags::ELECTRON_CAPTURE, "ElectronCapture"},
            {NetworkConstructionFlags::POSITRON_CAPTURE, "PositronCapture"}
        };

        size_t i = 0;
        for (const auto& flagType : bases_flags_array) {
            if (has_flag(flags, flagType)) {
                ss << bases_string_map.at(flagType);
                if (i < bases_flags_array.size() - 1) {
                    ss << ", ";
                }
            }
            ++i;
        }

        std::string result = ss.str();
        if (result.empty()) {
            return "No reactions";
        }
        return result;

    }


    /**
     * @brief Builds a nuclear reaction network from the Reaclib library based on an initial composition.
     *
     * Constructs a layered reaction network by collecting reactions up to the specified depth
     * from the Reaclib dataset. Starting species are those with non-zero mass fractions in the input
     * composition. Layers expand by including products of collected reactions until the depth limit.
     * Optionally selects reverse reactions instead of forward.
     *
     * See implementation in construction.cpp for details on the layering algorithm, logging, and performance.
     *
     * @param composition Mapping of isotopic species to their mass fractions; species with positive
     *        mass fraction seed the network.
     * @param weakInterpolator Interpolator to build weak rates from. Must be constructed and owned by the caller.
     * @param maxLayers Variant specifying either a predefined NetworkBuildDepth or a custom integer depth;
     *        negative depth (Full) collects all reactions, zero is invalid.
     * @param ReactionTypes
     * @pre composition must have at least one species with positive mass fraction.
     * @pre Resolved integer depth from maxLayers must not be zero.
     * @post Returned network includes only reactions satisfying the depth and reverse criteria.
     * @return A LogicalReactionSet encapsulating the collected reactions for graph-based engines.
     * @throws std::logic_error If the resolved network depth is zero (no reactions can be collected).
     */
    reaction::ReactionSet build_nuclear_network(
        const fourdst::composition::Composition &composition,
        const rates::weak::WeakRateInterpolator &weakInterpolator,
        BuildDepthType maxLayers = NetworkBuildDepth::Full,
        NetworkConstructionFlags ReactionTypes = NetworkConstructionFlags::DEFAULT
    );
}
