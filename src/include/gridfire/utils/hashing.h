#pragma once

#include <cstdint>

namespace gridfire::utils {
    /**
     * @brief Generate a unique hash for an isotope given its mass number (A) and atomic number (Z).
     * @details This function combines the mass number and atomic number into a single 32-bit integer
     *          by shifting the mass number 8 bits to the left and OR'ing it with the atomic number.
     *          This ensures a unique representation for each isotope within physically possible ranges.
     * @param a The mass number (A) of the isotope.
     * @param z The atomic number (Z) of the isotope.
     * @return A unique 32-bit hash representing the isotope. This is computed as (A << 8) | Z into an uint32_t.
     */
    inline uint_fast32_t hash_atomic(const uint16_t a, const uint8_t z) noexcept {
        return (static_cast<uint_fast32_t>(a) << 8) | static_cast<uint_fast32_t>(z);
    }

}