// Unit tests for GraphEngine rate multipliers: the factor reaches the RHS, both Jacobian methods and both eps
// derivative methods, and the CppAD tape picks up changes made after it was recorded (dynamic parameters).
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/exceptions/error_reaction.h"
#include "gridfire/policy/stellar_policy.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

using gridfire::config::GridFireConfig;
using gridfire::engine::GraphEngine;

namespace {
    struct BuiltEngine {
        std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
        std::unique_ptr<gridfire::policy::ConstructionResults> results;

        [[nodiscard]] GraphEngine& graph() const {
            return dynamic_cast<GraphEngine&>(*policy->get_engine_stack().front());
        }
        [[nodiscard]] gridfire::engine::scratch::StateBlob& ctx() const { return *results->scratch_blob; }
    };

    BuiltEngine build() {
        gridfire::test::quiet_logging();
        BuiltEngine built;
        built.policy = std::make_unique<gridfire::policy::MainSequencePolicy>(gridfire::test::solar_composition(), GridFireConfig{});
        built.results = std::make_unique<gridfire::policy::ConstructionResults>(built.policy->construct());
        return built;
    }

    fourdst::composition::Composition network_composition(const BuiltEngine& e) {
        const auto& sps = e.graph().getNetworkSpecies(e.ctx());
        return gridfire::test::make_composition(sps, gridfire::test::abundances_in_order(gridfire::test::solar_composition(), sps));
    }

    double rel(const double a, const double b) {
        const double scale = std::max(std::abs(a), std::abs(b));
        return scale == 0.0 ? 0.0 : std::abs(a - b) / scale;
    }

    constexpr double T9 = 0.0157;
    constexpr double RHO = 150.0;
}

TEST(RateMultipliers, BothJacobiansAndEpsDerivativesAgreeWithFactors) {
    const BuiltEngine e = build();
    const auto comp = network_composition(e);
    const auto& sps = e.graph().getNetworkSpecies(e.ctx());

    // Record the Jacobians once without factors, so the tape copy in this state blob is in use before the change
    const auto ad_before = e.graph().generateAutoDiffJacobianMatrix(e.ctx(), comp, T9, RHO);

    e.graph().setRateMultiplier("p(p,e+)d", 1.7);
    e.graph().setRateMultiplier("n14(p,g)o15", 0.3);

    const auto analytic = *e.graph().generateAnalyticJacobianMatrix(e.ctx(), comp, T9, RHO);
    const auto autodiff = e.graph().generateAutoDiffJacobianMatrix(e.ctx(), comp, T9, RHO);
    size_t changed = 0;
    for (const auto& row : sps) {
        for (const auto& col : sps) {
            const double a = analytic(row, col), d = autodiff(row, col);
            ASSERT_LT(rel(a, d), 1e-10) << "J(" << row.name() << ", " << col.name() << "): analytic " << a << ", CppAD " << d;
            if (rel(d, ad_before(row, col)) > 1e-12) ++changed;
        }
    }
    EXPECT_GT(changed, 0u) << "the factors did not change the CppAD Jacobian";

    const auto eps_analytic = *e.graph().calculateAnalyticEpsDerivatives(comp, T9, RHO);
    const auto eps_autodiff = e.graph().calculateAutoDiffEpsDerivatives(e.ctx(), comp, T9, RHO);
    EXPECT_LT(rel(eps_analytic.dEps_dT, eps_autodiff.dEps_dT), 1e-10);
    EXPECT_LT(rel(eps_analytic.dEps_dRho, eps_autodiff.dEps_dRho), 1e-10);
}

TEST(RateMultipliers, ScaleTheReactionInTheRhsAndClearExactly) {
    const BuiltEngine e = build();
    const auto comp = network_composition(e);
    const auto base = *e.graph().calculateRHSAndEnergy(e.ctx(), comp, T9, RHO, false);

    e.graph().setRateMultiplier("p(p,e+)d", 2.0);
    EXPECT_EQ(e.graph().getRateMultiplier("p(p,e+)d"), 2.0);
    const auto scaled = *e.graph().calculateRHSAndEnergy(e.ctx(), comp, T9, RHO, false);
    EXPECT_NE(scaled.nuclearEnergyGenerationRate, base.nuclearEnergyGenerationRate)
        << "a cached RHS from before the factor was returned";

    // H-1 is destroyed by pp (2 per reaction) and by much else; the extra destruction is exactly one more pp flow
    const auto& reactions = e.graph().getNetworkReactions(e.ctx());
    const auto pp = reactions.get("p(p,e+)d");
    ASSERT_TRUE(pp.has_value());
    const double pp_flow = e.graph().calculateMolarReactionFlow(e.ctx(), **pp, comp, T9, RHO);
    e.graph().clearRateMultipliers();
    const double pp_flow_base = e.graph().calculateMolarReactionFlow(e.ctx(), **pp, comp, T9, RHO);
    EXPECT_LT(rel(pp_flow, 2.0 * pp_flow_base), 1e-13);

    const fourdst::atomic::Species& H1 = fourdst::atomic::species.at("H-1");
    EXPECT_LT(rel(scaled.dydt.at(H1) - base.dydt.at(H1), -2.0 * pp_flow_base), 1e-9);

    const auto cleared = *e.graph().calculateRHSAndEnergy(e.ctx(), comp, T9, RHO, false);
    EXPECT_EQ(cleared.nuclearEnergyGenerationRate, base.nuclearEnergyGenerationRate);
    EXPECT_TRUE(e.graph().getRateMultipliers().empty());
}

TEST(RateMultipliers, RejectUnknownReactionsAndBadFactors) {
    const BuiltEngine e = build();
    EXPECT_THROW(e.graph().setRateMultiplier("x(y,z)w", 2.0), gridfire::exceptions::ReactionError);
    EXPECT_THROW(e.graph().setRateMultiplier("p(p,e+)d", -0.1), std::invalid_argument);
    EXPECT_THROW(e.graph().setRateMultiplier("p(p,e+)d", std::nan("")), std::invalid_argument);
    EXPECT_THROW((void)e.graph().getRateMultiplier("x(y,z)w"), gridfire::exceptions::ReactionError);

    e.graph().setRateMultiplier("p(p,e+)d", 1.25);
    e.graph().setRateMultiplier("n14(p,g)o15", 1.0); // 1 is no factor at all
    const auto listed = e.graph().getRateMultipliers();
    ASSERT_EQ(listed.size(), 1u);
    EXPECT_EQ(listed[0].first, "p(p,e+)d");
    EXPECT_EQ(listed[0].second, 1.25);
}
