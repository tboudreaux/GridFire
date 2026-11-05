// Unit tests for run-time configuration: a GridFireConfig passed to a policy, engine or view is the one it uses, and
// the constructors without one use config::GridFireConfig{}.
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/policy/stellar_policy.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

using gridfire::config::GridFireConfig;
using gridfire::engine::GraphEngine;
using gridfire::engine::JacobianMethod;
using gridfire::screening::ScreeningType;

namespace {
    /// A main sequence engine built from its own policy, so each test can choose the configuration.
    struct BuiltEngine {
        std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
        std::unique_ptr<gridfire::policy::ConstructionResults> results;

        [[nodiscard]] const GraphEngine& graph() const { return dynamic_cast<const GraphEngine&>(results->engine); }
        [[nodiscard]] gridfire::engine::scratch::StateBlob& ctx() const { return *results->scratch_blob; }
    };

    BuiltEngine build(const std::optional<GridFireConfig>& config) {
        gridfire::test::quiet_logging();
        BuiltEngine built;
        built.policy = config
            ? std::make_unique<gridfire::policy::MainSequencePolicy>(gridfire::test::solar_composition(), *config)
            : std::make_unique<gridfire::policy::MainSequencePolicy>(gridfire::test::solar_composition());
        built.results = std::make_unique<gridfire::policy::ConstructionResults>(built.policy->construct());
        return built;
    }

    /// Solar composition on the engine's network species (species not in the solar mix are absent).
    fourdst::composition::Composition network_composition(const BuiltEngine& e) {
        const auto& sps = e.graph().getNetworkSpecies(e.ctx());
        return gridfire::test::make_composition(sps, gridfire::test::abundances_in_order(gridfire::test::solar_composition(), sps));
    }

    constexpr double T9 = 0.015;
    constexpr double RHO = 150.0;
}

TEST(RuntimeConfig, ConstructorsWithoutConfigUseDefaults) {
    const BuiltEngine e = build(std::nullopt);
    const GridFireConfig defaults{};
    EXPECT_EQ(e.graph().getConfig().engine.screening, defaults.engine.screening);
    EXPECT_EQ(e.graph().getConfig().engine.jacobian, defaults.engine.jacobian);
    EXPECT_EQ(e.graph().getScreeningModel(e.ctx()), ScreeningType::BLENDED);
    EXPECT_EQ(e.policy->getConfig().engine.jacobian, JacobianMethod::ANALYTIC);
}

TEST(RuntimeConfig, PolicyPassesConfigToEngine) {
    GridFireConfig config;
    config.engine.screening = ScreeningType::BARE;
    config.engine.jacobian = JacobianMethod::AUTOMATIC_DIFFERENTIATION;
    config.engine.energyChannels.splitHydrogenPPCNO = false;
    const BuiltEngine e = build(config);

    EXPECT_EQ(e.graph().getScreeningModel(e.ctx()), ScreeningType::BARE);
    EXPECT_EQ(e.graph().getConfig().engine.jacobian, JacobianMethod::AUTOMATIC_DIFFERENTIATION);
    // The energy-channel option is applied while the engine is built, so this also checks the configuration is in
    // place before syncInternalMaps()
    const auto& channels = e.graph().getEnergyChannelNames(e.ctx());
    EXPECT_EQ(std::ranges::find(channels, "H_pp"), channels.end());
    EXPECT_NE(std::ranges::find(channels, "H"), channels.end());
}

TEST(RuntimeConfig, JacobianMethodSelectsJacobianAndEpsDerivatives) {
    for (const JacobianMethod method : {JacobianMethod::ANALYTIC, JacobianMethod::AUTOMATIC_DIFFERENTIATION}) {
        GridFireConfig config;
        config.engine.jacobian = method;
        const BuiltEngine e = build(config);
        const auto comp = network_composition(e);
        const auto& sps = e.graph().getNetworkSpecies(e.ctx());

        const auto selected = e.graph().generateJacobianMatrix(e.ctx(), comp, T9, RHO);
        const auto expected = method == JacobianMethod::ANALYTIC
            ? *e.graph().generateAnalyticJacobianMatrix(e.ctx(), comp, T9, RHO)
            : e.graph().generateAutoDiffJacobianMatrix(e.ctx(), comp, T9, RHO);
        for (const auto& row : sps) {
            for (const auto& col : sps) {
                ASSERT_EQ(selected(row, col), expected(row, col)) << "J(" << row.name() << ", " << col.name() << ")";
            }
        }

        // The same setting selects how d(eps)/dT and d(eps)/d(rho) are computed
        const auto eps_selected = e.graph().calculateEpsDerivatives(e.ctx(), comp, T9, RHO);
        const auto eps_expected = method == JacobianMethod::ANALYTIC
            ? *e.graph().calculateAnalyticEpsDerivatives(comp, T9, RHO)
            : e.graph().calculateAutoDiffEpsDerivatives(e.ctx(), comp, T9, RHO);
        EXPECT_EQ(eps_selected.dEps_dT, eps_expected.dEps_dT);
        EXPECT_EQ(eps_selected.dEps_dRho, eps_expected.dEps_dRho);
    }
}

TEST(RuntimeConfig, ScreeningSelectsTheRates) {
    GridFireConfig bare;
    bare.engine.screening = ScreeningType::BARE;
    GridFireConfig weak;
    weak.engine.screening = ScreeningType::WEAK;
    const BuiltEngine screened = build(weak);
    const BuiltEngine unscreened = build(bare);
    const BuiltEngine by_default = build(std::nullopt);

    const auto eps = [](const BuiltEngine& e) {
        const auto rhs = e.graph().calculateRHSAndEnergy(e.ctx(), network_composition(e), T9, RHO, false);
        EXPECT_TRUE(rhs.has_value());
        return rhs->nuclearEnergyGenerationRate;
    };
    // Weak screening raises the solar-core pp rate by a few percent
    EXPECT_GT(eps(screened) / eps(unscreened), 1.02);
    // The default (BLENDED) is Salpeter wherever weak screening holds; in the solar core only a few heavy-ion
    // reactions reach the blend window, and they barely contribute to eps
    EXPECT_NEAR(eps(by_default), eps(screened), 1e-9 * eps(screened));
}

TEST(RuntimeConfig, GraphEngineFromReactionsTakesConfig) {
    const BuiltEngine source = build(std::nullopt);
    GridFireConfig config;
    config.engine.screening = ScreeningType::BARE;
    const GraphEngine engine(source.graph().getNetworkReactions(source.ctx()), config);
    const auto ctx = engine.constructStateBlob(nullptr);
    EXPECT_EQ(engine.getScreeningModel(*ctx), ScreeningType::BARE);
    EXPECT_EQ(GraphEngine(source.graph().getNetworkReactions(source.ctx())).getConfig().engine.screening, ScreeningType::BLENDED);
}

TEST(RuntimeConfig, AdaptiveViewTakesConfig) {
    BuiltEngine e = build(std::nullopt);
    GridFireConfig config;
    config.engine.views.adaptiveEngineView.relativeCullingThreshold = 1e-30;
    auto& base = const_cast<gridfire::engine::DynamicEngine&>(e.results->engine);
    const gridfire::engine::AdaptiveEngineView configured(base, config);
    const gridfire::engine::AdaptiveEngineView defaulted(base);
    EXPECT_EQ(configured.getConfig().engine.views.adaptiveEngineView.relativeCullingThreshold, 1e-30);
    EXPECT_EQ(defaulted.getConfig().engine.views.adaptiveEngineView.relativeCullingThreshold,
              GridFireConfig{}.engine.views.adaptiveEngineView.relativeCullingThreshold);
}
