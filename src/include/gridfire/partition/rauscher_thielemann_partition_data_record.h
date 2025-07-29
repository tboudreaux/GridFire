#pragma once

#include <cstdint>

namespace gridfire::partition::record {
#pragma pack(push, 1)
/**
 * @struct RauscherThielemannPartitionDataRecord
 * @brief Packed binary record of Rauscher-Thielemann partition function data for an isotope.
 *
 * Each record stores the atomic number (Z), mass number (A), ground state spin J,
 * and an array of 24 normalized G-values corresponding to fixed temperature grid points.
 * This struct is read directly from embedded binary data and must remain tightly packed.
 *
 * @note Alignment is set to 1 byte to match the binary layout.
 */
    struct RauscherThielemannPartitionDataRecord {
        uint32_t z; ///< Atomic number
        uint32_t a; ///< Mass number
        double ground_state_spin; ///< Ground state spin
        double normalized_g_values[24]; ///< Normalized g-values for the first 24 energy levels
    };
#pragma pack(pop)
}