#pragma once

#include <cstdint>
#include <functional>

#include "fourdst/composition/utils/composition_hash.h"
#include "gridfire/exceptions/exceptions.h"
#include "gridfire/reaction/reaction.h"

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

    namespace hashing::reaction {
        static std::uint64_t splitmix64(std::uint64_t x) noexcept {
            x += 0x9E3779B97F4A7C15ULL;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
            x ^= (x >> 31);
            return x;
        }

        static std::uint64_t mix_species(const unsigned a, const unsigned z) noexcept {
            const std::uint64_t code = (static_cast<std::uint64_t>(a) << 7) | static_cast<std::uint64_t>(z);
            return splitmix64(code);
        }

        static std::uint64_t multiset_combine(std::uint64_t acc, const std::uint64_t x) noexcept {
            acc += x;
            acc ^= (x << 23) | (x >> (64 - 23));
            acc = splitmix64(acc);
            return acc;
        }
    }

    inline std::uint64_t hash_reaction(const reaction::Reaction& reaction) noexcept {
        using namespace hashing::reaction;

        std::uint64_t hR = 0;

        for (const auto& s : reaction.reactants()) {
            hR = multiset_combine(hR, mix_species(static_cast<unsigned>(s.a()),
                                                  static_cast<unsigned>(s.z())));
        }

        std::uint64_t hP = 0;
        for (const auto& s : reaction.products()) {
            hP = multiset_combine(hP, mix_species(static_cast<unsigned>(s.a()),
                                                  static_cast<unsigned>(s.z())));
        }

        std::uint64_t h = splitmix64(hR ^ 0xC3A5C85C97CB3127ULL);
        h ^= splitmix64((hP << 1) | (hP >> 63));
        return splitmix64(h);
    }

    template <typename T>
    std::size_t hash_combine(std::size_t seed, const T& v) {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }

    inline std::size_t fast_mix(std::size_t h) noexcept {
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return h;
    }

    inline std::size_t hash_state(
        const fourdst::composition::CompositionAbstract& comp,
        const double T9,
        const double rho,
        const reaction::ReactionSet& reactions
    ) noexcept {
        std::size_t hash = comp.hash();
        const std::size_t topology_hash = reactions.hash(0);

        hash ^= topology_hash + 0x517cc1b727220a95 + (hash << 6) + (hash >> 2);

        const std::uint64_t t9_bits = std::bit_cast<std::uint64_t>(T9);
        const std::uint64_t rho_bits = std::bit_cast<std::uint64_t>(rho);

        hash ^= fast_mix(t9_bits) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        hash ^= fast_mix(rho_bits) + 0x9e3779b9 + (hash << 6) + (hash >> 2);

        return hash;
    }
}
