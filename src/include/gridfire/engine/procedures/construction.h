/**
 * @file construction.h
 * @brief Functions for constructing nuclear reaction networks.
 */
#pragma once

#include "gridfire/reaction/reaction.h"
#include "gridfire/engine/types/building.h"

#include "fourdst/composition/composition_abstract.h"

#include <variant>

#include "gridfire/reaction/weak/weak_interpolator.h"

namespace gridfire::engine {
    /**
     * @brief Flags to specify which types of nuclear reactions to include when constructing a reaction network.
     *
     * These flags allow fine-grained control over the inclusion of strong and weak nuclear reactions
     * (beta decay, electron/positron capture) from various sources (Reaclib, WRL) during network construction.
     * They can be combined using bitwise operations to create custom reaction sets.
     */
    enum class NetworkConstructionFlags : uint32_t {
        NONE = 0,

        STRONG = 1 << 0, // 1

        BETA_MINUS = 1 << 1, // 2
        BETA_PLUS = 1 << 2, // 4
        ELECTRON_CAPTURE = 1 << 3, // 8
        POSITRON_CAPTURE = 1 << 4, // 16

        REACLIB_WEAK = 1 << 5,

        WRL_WEAK = BETA_MINUS | BETA_PLUS | ELECTRON_CAPTURE | POSITRON_CAPTURE,

        REACLIB = STRONG | REACLIB_WEAK,

        // Currently we default to just reaclib reactions but include both their strong and weak set
        DEFAULT = REACLIB,

        ALL = STRONG | WRL_WEAK
    };

    /** @brief Helper function to convert NetworkConstructionFlags to their underlying integer type.
     *
     * This function facilitates bitwise operations on NetworkConstructionFlags by converting them
     * to their underlying integer representation.
     *
     * @param f The NetworkConstructionFlags value to convert.
     * @return The underlying integer representation of the flag.
     */
    constexpr auto to_underlying(NetworkConstructionFlags f) noexcept {
        return static_cast<std::underlying_type_t<NetworkConstructionFlags>>(f);
    }

    /** @brief Bitwise OR operator for NetworkConstructionFlags.
     *
     * This operator allows combining two NetworkConstructionFlags values using the bitwise OR operation.
     *
     * @param lhs The left-hand side NetworkConstructionFlags value.
     * @param rhs The right-hand side NetworkConstructionFlags value.
     * @return A new NetworkConstructionFlags value representing the combination of the two inputs.
     */
    inline NetworkConstructionFlags operator|(const NetworkConstructionFlags lhs, const NetworkConstructionFlags rhs) {
        return static_cast<NetworkConstructionFlags>(to_underlying(lhs) | to_underlying(rhs));
    }

    /** @brief Bitwise AND operator for NetworkConstructionFlags.
     *
     * This operator allows checking for common flags between two NetworkConstructionFlags values
     * using the bitwise AND operation.
     *
     * @param lhs The left-hand side NetworkConstructionFlags value.
     * @param rhs The right-hand side NetworkConstructionFlags value.
     * @return A new NetworkConstructionFlags value representing the intersection of the two inputs.
     */
    inline NetworkConstructionFlags operator&(const NetworkConstructionFlags lhs, const NetworkConstructionFlags rhs) {
        return static_cast<NetworkConstructionFlags>(to_underlying(lhs) & to_underlying(rhs));
    }

    /** @brief Checks if a specific flag is set within a NetworkConstructionFlags value.
     *
     * This function determines whether a particular flag is present in a given NetworkConstructionFlags value.
     *
     * @param flags The NetworkConstructionFlags value to check.
     * @param flag_to_check The specific flag to look for.
     * @return True if the flag is set; otherwise, false.
     */
    inline bool has_flag(const NetworkConstructionFlags flags, const NetworkConstructionFlags flag_to_check) {
        return (flags & flag_to_check) != NetworkConstructionFlags::NONE;
    }

    /** @brief Converts NetworkConstructionFlags to a human-readable string.
     *
     * This function generates a comma-separated string representation of the set flags
     * within a NetworkConstructionFlags value. If no flags are set, it returns "No reactions".
     *
     * @param flags The NetworkConstructionFlags value to convert.
     * @return A string listing the set flags or "No reactions" if none are set.
     */
    inline std::string NetworkConstructionFlagsToString(NetworkConstructionFlags flags) {
        std::stringstream ss;
        constexpr std::array<NetworkConstructionFlags, 6> bases_flags_array = {
            NetworkConstructionFlags::STRONG,
            NetworkConstructionFlags::BETA_MINUS,
            NetworkConstructionFlags::BETA_PLUS,
            NetworkConstructionFlags::ELECTRON_CAPTURE,
            NetworkConstructionFlags::POSITRON_CAPTURE,
            NetworkConstructionFlags::REACLIB_WEAK
        };

        const std::unordered_map<NetworkConstructionFlags, std::string> bases_string_map = {
            {NetworkConstructionFlags::STRONG, "Strong"},
            {NetworkConstructionFlags::BETA_MINUS, "BetaMinus"},
            {NetworkConstructionFlags::BETA_PLUS, "BetaPlus"},
            {NetworkConstructionFlags::ELECTRON_CAPTURE, "ElectronCapture"},
            {NetworkConstructionFlags::POSITRON_CAPTURE, "PositronCapture"},
            {NetworkConstructionFlags::REACLIB_WEAK, "ReaclibWeak"}
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
        const fourdst::composition::CompositionAbstract &composition,
        const rates::weak::WeakRateInterpolator &weakInterpolator,
        BuildDepthType maxLayers = NetworkBuildDepth::Full,
        NetworkConstructionFlags ReactionTypes = NetworkConstructionFlags::DEFAULT
    );
}
