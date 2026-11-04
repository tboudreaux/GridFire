/**
 * @file gf_test_support.h
 * @brief Shared fixtures and helpers for the GridFire test suites.
 *
 * Engines are expensive to construct (network building and AD taping), so they are built once per test process and
 * shared. Every test that needs mutable engine state gets its own StateBlob via fresh_context(). Building the shared
 * engine also runs GF_PAR_INIT, so tests may use OpenMP (but not std::thread: CppAD identifies threads via OpenMP).
 *
 * Randomised tests draw from a seeded generator. The seed defaults to a fixed value and can be overridden with the
 * GF_TEST_SEED environment variable; it is reported on failure so any failure can be reproduced.
 */
#pragma once

#include "gridfire/gridfire.h"
#include "gridfire/engine/engine_graph.h"

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/composition/composition.h"

#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace gridfire::test {

    /// @brief Silence GridFire logging and route log files away from the source tree. Safe to call repeatedly.
    void quiet_logging();

    /// @brief Solar-like composition over H-1, He-3, He-4, C-12, N-14, O-16, Ne-20, Mg-24.
    fourdst::composition::Composition solar_composition();

    /// @brief A main sequence network (MainSequencePolicy) built once per process.
    struct SharedEngine {
        const engine::DynamicEngine& engine; ///< Top of the policy's engine stack
        const engine::GraphEngine& graph; ///< The underlying GraphEngine

        /// @brief An independent StateBlob for this engine.
        [[nodiscard]] std::unique_ptr<engine::scratch::StateBlob> fresh_context() const;
    };

    /// @brief The shared main sequence engine.
    const SharedEngine& main_sequence_engine();

    /// @brief A thermodynamic state and composition at which to probe an engine.
    struct EngineState {
        double T9;
        double rho;
        fourdst::composition::Composition composition;

        [[nodiscard]] std::string describe() const;
    };

    /// @brief Seed for randomised tests (GF_TEST_SEED or a fixed default).
    std::uint64_t test_seed();

    /**
     * @brief Draw reproducible random engine states.
     *
     * T9 and rho are log-uniform. Abundances start from the solar composition for seed species and log-uniform trace
     * abundances for every other network species, then are normalised to unit mass fraction.
     */
    class StateSampler {
    public:
        StateSampler(
            const std::vector<fourdst::atomic::Species>& species,
            std::uint64_t seed,
            double T9_min = 0.01,
            double T9_max = 3.0,
            double rho_min = 1.0,
            double rho_max = 1.0e6
        );

        EngineState next();

    private:
        std::vector<fourdst::atomic::Species> m_species;
        std::mt19937_64 m_rng;
        double m_T9_min, m_T9_max, m_rho_min, m_rho_max;
    };

    /// @brief Composition over exactly @p species with the given molar abundances.
    fourdst::composition::Composition make_composition(
        const std::vector<fourdst::atomic::Species>& species,
        const std::vector<double>& Y
    );

    /// @brief Molar abundances of @p comp ordered like @p species.
    std::vector<double> abundances_in_order(
        const fourdst::composition::CompositionAbstract& comp,
        const std::vector<fourdst::atomic::Species>& species
    );

    /**
     * @brief Relative difference with a floor: |a - b| / max(|a|, |b|, floor).
     *
     * The floor keeps quantities that should both be ~0 from producing huge relative errors.
     */
    double rel_diff(double a, double b, double floor = 0.0);

    /// @brief A finite-difference derivative together with an estimate of its round-off error.
    struct FiniteDifference {
        double value; ///< Richardson-extrapolated central difference
        double noise; ///< Round-off allowance: a few ulps of the largest sampled value, divided by the step
    };

    /**
     * @brief Richardson-extrapolated central difference of f at x with absolute step h.
     *
     * Combines central differences at h and h/2, so truncation error is O(h^4). Steps should be moderate (around
     * 1e-3 relative): very small steps make round-off dominate, which the noise estimate accounts for.
     */
    FiniteDifference richardson_derivative(const std::function<double(double)>& f, double x, double h);

    /// @brief Fundamental constants used to build eps (same source as the engine).
    struct PhysicalConstants {
        double u; ///< Atomic mass unit [g]
        double Na; ///< Avogadro's number [1/mol]
        double c; ///< Speed of light [cm/s]
    };
    const PhysicalConstants& constants();
}
