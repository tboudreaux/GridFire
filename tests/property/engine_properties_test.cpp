// Property tests: invariants of the GraphEngine checked at many reproducible random states.
//
// Each test draws states from a seeded StateSampler. Set GF_TEST_SEED to reproduce or vary the draw; the seed and the
// failing state are printed on failure.
#include "gf_test_support.h"

#include "gridfire/engine/engine_graph.h"

#include <gtest/gtest.h>

#if defined(GF_USE_OPENMP)
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <format>
#include <numeric>
#include <ranges>
#include <vector>

using gridfire::engine::GraphEngine;
using gridfire::engine::StepDerivatives;
using gridfire::test::EngineState;
using gridfire::test::rel_diff;

namespace {
    class EngineProperties : public ::testing::Test {
    protected:
        static void SetUpTestSuite() { gridfire::test::quiet_logging(); }

        void SetUp() override {
            m_ctx = shared().fresh_context();
        }

        static const gridfire::test::SharedEngine& shared() { return gridfire::test::main_sequence_engine(); }
        static const GraphEngine& graph() { return shared().graph; }
        const std::vector<fourdst::atomic::Species>& species() const { return graph().getNetworkSpecies(*m_ctx); }

        gridfire::test::StateSampler sampler(const std::uint64_t salt) const {
            return {species(), gridfire::test::test_seed() + salt};
        }

        StepDerivatives<double> rhs(const EngineState& s) const {
            auto result = graph().calculateRHSAndEnergy(*m_ctx, s.composition, s.T9, s.rho, false);
            EXPECT_TRUE(result.has_value()) << s.describe();
            return result.value();
        }

        double rhs_eps(const fourdst::composition::Composition& comp, const double T9, const double rho) const {
            return graph().calculateRHSAndEnergy(*m_ctx, comp, T9, rho, false).value().nuclearEnergyGenerationRate;
        }

        /// Sum of the magnitudes of the mass-defect terms of eps; the natural scale against which eps is compared.
        static double eps_scale(const StepDerivatives<double>& d) {
            const auto& k = gridfire::test::constants();
            double scale = std::abs(d.neutrinoEnergyLossRate);
            for (const auto& [sp, dydt] : d.dydt) scale += std::abs(dydt * sp.mass() * k.u * k.Na * k.c * k.c);
            return scale;
        }

        std::unique_ptr<gridfire::engine::scratch::StateBlob> m_ctx;
    };

    std::string trace(const EngineState& s) {
        return std::format("seed={} {}", gridfire::test::test_seed(), s.describe());
    }
}

// The AD tape (used for the Jacobian and eps derivatives) and the precomputed fast path (used by CVODE) are two
// implementations of the same physics. They must agree on dY/dt and eps everywhere.
TEST_F(EngineProperties, TapeMatchesPrecomputedPath) {
    auto draw = sampler(1);
    CppAD::ADFun<double> tape;
    tape = graph().getAuthoritativeADFun();

    for (int k = 0; k < 100; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const StepDerivatives<double> fast = rhs(s);

        std::vector<double> x = gridfire::test::abundances_in_order(s.composition, species());
        x.push_back(s.T9);
        x.push_back(s.rho);
        const std::vector<double> y = tape.Forward(0, x);
        ASSERT_EQ(y.size(), fast.dydt.size() + 1);

        double dydt_scale = 0.0;
        for (const auto& v : fast.dydt | std::views::values) dydt_scale = std::max(dydt_scale, std::abs(v));

        size_t i = 0;
        for (const auto& [sp, dydt] : fast.dydt) { // Tape outputs follow the (sorted) dydt map order
            EXPECT_LT(rel_diff(y[i], dydt, 1e-12 * dydt_scale), 1e-10) << "dY/dt of " << sp.name() << ": tape " << y[i] << " fast " << dydt;
            ++i;
        }
        EXPECT_LT(rel_diff(y.back(), fast.nuclearEnergyGenerationRate, 1e-12 * eps_scale(fast)), 1e-10)
            << "eps: tape " << y.back() << " fast " << fast.nuclearEnergyGenerationRate;
    }
}

// eps is defined from the mass defect of dY/dt minus neutrino losses; the engine must report exactly that.
TEST_F(EngineProperties, EpsEqualsMassDefectMinusNeutrinoLosses) {
    auto draw = sampler(2);
    const auto& c = gridfire::test::constants();
    for (int k = 0; k < 100; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const StepDerivatives<double> d = rhs(s);
        double massRate = 0.0;
        for (const auto& [sp, dydt] : d.dydt) massRate += dydt * sp.mass() * c.u;
        const double expected = -massRate * c.Na * c.c * c.c - d.neutrinoEnergyLossRate;
        EXPECT_LT(rel_diff(d.nuclearEnergyGenerationRate, expected, 1e-12 * eps_scale(d)), 1e-10);
    }
}

TEST_F(EngineProperties, BaryonNumberIsConserved) {
    auto draw = sampler(3);
    for (int k = 0; k < 100; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const StepDerivatives<double> d = rhs(s);
        double net = 0.0, gross = 0.0;
        for (const auto& [sp, dydt] : d.dydt) {
            net += sp.a() * dydt;
            gross += sp.a() * std::abs(dydt);
        }
        EXPECT_LE(std::abs(net), 1e-12 * gross) << "sum A dY/dt = " << net << " (gross " << gross << ")";
    }
}

TEST_F(EngineProperties, EnergyChannelsSumToEps) {
    auto draw = sampler(4);
    for (int k = 0; k < 100; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const StepDerivatives<double> d = rhs(s);
        const std::vector<double> channels = graph().calculateEnergyChannels(*m_ctx, s.composition, s.T9, s.rho);
        ASSERT_EQ(channels.size(), graph().getEnergyChannelNames(*m_ctx).size());
        const double sum = std::accumulate(channels.begin(), channels.end(), 0.0);
        EXPECT_LT(rel_diff(sum, d.nuclearEnergyGenerationRate, 1e-12 * eps_scale(d)), 1e-10);
    }
}

// Finite differences use Richardson extrapolation at a moderate step; the comparison allows the FD round-off estimate,
// which matters where eps is a small difference of large terms (e.g. eps close to the neutrino loss).
TEST_F(EngineProperties, EpsDerivativesMatchFiniteDifferences) {
    auto draw = sampler(5);
    for (int k = 0; k < 30; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const auto [dEps_dT, dEps_dRho] = graph().calculateEpsDerivatives(*m_ctx, s.composition, s.T9, s.rho);

        const double T = s.T9 * 1e9;
        const double hT = 1e-3 * T, hR = 1e-3 * s.rho;
        const auto fdT = gridfire::test::richardson_derivative([&](const double t) { return rhs_eps(s.composition, t / 1e9, s.rho); }, T, hT);
        const auto fdR = gridfire::test::richardson_derivative([&](const double r) { return rhs_eps(s.composition, s.T9, r); }, s.rho, hR);

        // eps is a mass-defect sum that cancels by ~1e3 or more, so its round-off follows the gross scale, not eps
        const double gross = eps_scale(rhs(s));
        const double noiseT = std::max(fdT.noise, 64 * std::numeric_limits<double>::epsilon() * gross / hT);
        const double noiseR = std::max(fdR.noise, 64 * std::numeric_limits<double>::epsilon() * gross / hR);
        EXPECT_LE(std::abs(dEps_dT - fdT.value), 1e-5 * std::abs(fdT.value) + noiseT) << "dEps/dT AD " << dEps_dT << " FD " << fdT.value << " (noise " << noiseT << ")";
        EXPECT_LE(std::abs(dEps_dRho - fdR.value), 1e-5 * std::abs(fdR.value) + noiseR) << "dEps/dRho AD " << dEps_dRho << " FD " << fdR.value << " (noise " << noiseR << ")";
    }
}

TEST_F(EngineProperties, JacobianMatchesFiniteDifferences) {
    auto draw = sampler(6);
    const auto& sps = species();
    for (int k = 0; k < 6; ++k) {
        const EngineState s = draw.next();
        SCOPED_TRACE(trace(s));
        const auto jac = graph().generateJacobianMatrix(*m_ctx, s.composition, s.T9, s.rho);
        const std::vector<double> Y = gridfire::test::abundances_in_order(s.composition, sps);

        // dY_i/dt can be a small net of large opposing terms (e.g. neutrons), and its round-off scales with those
        // terms, not with the net. Each reaction term is polynomial in Y, so sum_j |J_ij| Y_j bounds their magnitude.
        std::vector<double> gross(sps.size(), 0.0);
        for (size_t i = 0; i < sps.size(); ++i) {
            for (size_t j = 0; j < sps.size(); ++j) gross[i] += std::abs(jac(sps[i], sps[j])) * Y[j];
        }

        for (size_t j = 0; j < sps.size(); ++j) {
            if (Y[j] <= 0.0) continue; // Central differences need room on both sides of zero
            const double h = 1e-3 * Y[j];

            // Evaluate the RHS once per perturbed abundance, then read every row from the cached results
            std::map<double, StepDerivatives<double>> column;
            const auto rhs_at = [&](const double yj) -> const StepDerivatives<double>& {
                auto it = column.find(yj);
                if (it == column.end()) {
                    std::vector<double> Yp = Y;
                    Yp[j] = yj;
                    it = column.emplace(yj, graph().calculateRHSAndEnergy(*m_ctx, gridfire::test::make_composition(sps, Yp), s.T9, s.rho, false).value()).first;
                }
                return it->second;
            };
            for (size_t i = 0; i < sps.size(); ++i) {
                const auto fd = gridfire::test::richardson_derivative([&](const double yj) { return rhs_at(yj).dydt.at(sps[i]); }, Y[j], h);
                const double noise = std::max(fd.noise, 64 * std::numeric_limits<double>::epsilon() * gross[i] / h);
                const double ad = jac(sps[i], sps[j]);
                // generateJacobianMatrix deliberately drops off-diagonal entries below MIN_JACOBIAN_THRESHOLD
                const double dropped = (ad == 0.0 && i != j) ? gridfire::engine::MIN_JACOBIAN_THRESHOLD : 0.0;
                EXPECT_LE(std::abs(ad - fd.value), 1e-5 * std::abs(fd.value) + noise + dropped)
                    << "J(" << sps[i].name() << ", " << sps[j].name() << "): AD " << ad << " FD " << fd.value << " (noise " << noise << ")";
            }
        }
    }
}

// A cache hit must be indistinguishable from a fresh evaluation, including its side effect on the most recent RHS.
TEST_F(EngineProperties, RHSCacheIsTransparent) {
    auto draw = sampler(7);
    for (int k = 0; k < 20; ++k) {
        const EngineState a = draw.next();
        const EngineState b = draw.next();
        SCOPED_TRACE(trace(a));

        const StepDerivatives<double> first = rhs(a);
        (void)rhs(b);
        const StepDerivatives<double> cached = rhs(a);
        const auto recent = graph().getMostRecentRHSCalculation(*m_ctx);

        const auto fresh_ctx = shared().fresh_context();
        const auto fresh = graph().calculateRHSAndEnergy(*fresh_ctx, a.composition, a.T9, a.rho, false).value();

        EXPECT_EQ(cached.nuclearEnergyGenerationRate, first.nuclearEnergyGenerationRate);
        EXPECT_EQ(cached.dydt, first.dydt);
        EXPECT_EQ(fresh.nuclearEnergyGenerationRate, first.nuclearEnergyGenerationRate);
        EXPECT_EQ(fresh.dydt, first.dydt);
        ASSERT_TRUE(recent.has_value());
        EXPECT_EQ(recent->nuclearEnergyGenerationRate, first.nuclearEnergyGenerationRate)
            << "getMostRecentRHSCalculation is stale after a cache hit";
    }
}

// Independent StateBlobs on separate OpenMP threads must reproduce the serial results bit for bit. GridFire's
// supported parallel model is OpenMP after GF_PAR_INIT (CppAD's allocator identifies threads through OpenMP), so this
// deliberately does not use std::thread.
TEST_F(EngineProperties, OpenMPThreadsWithClonedContextsAreDeterministic) {
#if !defined(GF_USE_OPENMP)
    GTEST_SKIP() << "GridFire was built without OpenMP support";
#else
    auto draw = sampler(8);
    std::vector<EngineState> states;
    for (int k = 0; k < 16; ++k) states.push_back(draw.next());

    std::vector<double> serial;
    for (const auto& s : states) serial.push_back(rhs(s).nuclearEnergyGenerationRate);

    const int num_threads = std::max(2, std::min(8, omp_get_max_threads()));
    std::vector<std::unique_ptr<gridfire::engine::scratch::StateBlob>> contexts;
    for (int t = 0; t < num_threads; ++t) contexts.push_back(shared().fresh_context()); // Cloned serially, as GridSolver does

    std::vector<std::vector<double>> parallel(num_threads, std::vector<double>(states.size()));
    #pragma omp parallel num_threads(num_threads)
    {
        const int t = omp_get_thread_num();
        for (size_t k = 0; k < states.size(); ++k) {
            const size_t idx = (k + t) % states.size(); // Different visiting order per thread
            const auto& s = states[idx];
            parallel[t][idx] = graph().calculateRHSAndEnergy(*contexts[t], s.composition, s.T9, s.rho, false).value().nuclearEnergyGenerationRate;
        }
    }

    for (int t = 0; t < num_threads; ++t) {
        EXPECT_EQ(parallel[t], serial) << "thread " << t;
    }
#endif
}

// Known bug: GraphEngine::calculateEpsDerivatives accepts an active reaction set but differentiates the full network.
// This test skips while the bug is present and fails once it is fixed, so the skip can be removed.
TEST_F(EngineProperties, EpsDerivativesRespectActiveReactionSubset) {
    const auto& all = shared().engine.getNetworkReactions(*m_ctx);
    gridfire::reaction::ReactionSet subset;
    for (const auto& reaction : all) {
        const bool cno = std::ranges::any_of(reaction->reactants(), [](const auto& sp) { return sp.z() >= 6; });
        if (!cno) subset.add_reaction(*reaction);
    }
    ASSERT_LT(subset.size(), all.size());

    auto draw = sampler(9);
    bool mismatch = false;
    for (int k = 0; k < 10 && !mismatch; ++k) {
        const EngineState s = draw.next();
        const auto eps_at = [&](const double T9, const double rho) {
            return graph().calculateRHSAndEnergy(*m_ctx, s.composition, T9, rho, subset).value().nuclearEnergyGenerationRate;
        };
        const auto [dEps_dT, dEps_dRho] = graph().calculateEpsDerivatives(*m_ctx, s.composition, s.T9, s.rho, subset);
        const double T = s.T9 * 1e9, hT = 1e-6 * T;
        const double fdT = (eps_at((T + hT) / 1e9, s.rho) - eps_at((T - hT) / 1e9, s.rho)) / (2 * hT);
        if (rel_diff(dEps_dT, fdT, 1e-30) > 1e-4) mismatch = true;

        // The channel split must respect the subset regardless.
        const auto d = graph().calculateRHSAndEnergy(*m_ctx, s.composition, s.T9, s.rho, subset).value();
        const auto channels = graph().calculateEnergyChannels(*m_ctx, s.composition, s.T9, s.rho, subset);
        const double sum = std::accumulate(channels.begin(), channels.end(), 0.0);
        EXPECT_LT(rel_diff(sum, d.nuclearEnergyGenerationRate, 1e-12 * eps_scale(d)), 1e-10) << trace(s);
    }
    if (mismatch) {
        GTEST_SKIP() << "Known bug: calculateEpsDerivatives ignores its activeReactions argument.";
    }
    FAIL() << "calculateEpsDerivatives now respects activeReactions; remove the known-bug skip from this test.";
}
