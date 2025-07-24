#pragma once

#include <unordered_map>
#include <string>

namespace gridfire::partition {

    /**
     * @enum BasePartitionType
     * @brief Enumerates available partition function implementations.
     *
     * RauscherThielemann: Uses tabulated normalized G-values and linear interpolation.
     * GroundState: Uses ground state spin (J) to compute partition function as 2J+1.
     */
    enum BasePartitionType {
        RauscherThielemann, ///< Rauscher-Thielemann partition function
        GroundState,        ///< Ground state partition function
    };

    /**
     * @brief Mapping from BasePartitionType enum to human-readable string.
     *
     * Used for logging, reporting, or serialization. Ensure that all enum values
     * are represented in this map.
     * @pre Contains entries for all values of BasePartitionType.
     * @post Can convert BasePartitionType to corresponding string.
     */
    inline std::unordered_map<BasePartitionType, std::string> basePartitionTypeToString = {
        {RauscherThielemann, "RauscherThielemann"},
        {GroundState, "GroundState"}
    };

    /**
     * @brief Mapping from string to BasePartitionType enum.
     *
     * Used for parsing configuration or user input. Strings must match exactly
     * to one of the defined partition types.
     * @pre Uses keys that exactly match the outputs of basePartitionTypeToString.
     * @post Can convert valid string identifiers back to BasePartitionType.
     * @throws std::out_of_range if accessed with a non-existing key via at().
     */
    inline std::unordered_map<std::string, BasePartitionType> stringToBasePartitionType = {
        {"RauscherThielemann", RauscherThielemann},
        {"GroundState", GroundState}
    };

}