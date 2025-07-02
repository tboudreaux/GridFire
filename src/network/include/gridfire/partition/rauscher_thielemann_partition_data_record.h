#pragma once

#include <cstdint>

namespace gridfire::partition::record {
#pragma pack(push, 1)
    struct RauscherThielemannPartitionDataRecord {
        uint32_t z; ///< Atomic number
        uint32_t a; ///< Mass number
        double ground_state_spin; ///< Ground state spin
        double partition_function; ///< Partition function value
        double normalized_g_values[24]; ///< Normalized g-values for the first 24 energy levels
    };
#pragma pack(pop)
}