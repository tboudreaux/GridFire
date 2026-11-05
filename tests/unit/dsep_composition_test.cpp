/**
 * Per-reaction contributions (DSEP requirement 8) and passive species (requirement 7).
 *
 * Per reaction: the contributions sum to the totals of calculateRHSAndEnergy, and pep is its own, electron-capture
 * scaled, reaction. Passive species: they sit in the network with dY/dt = 0, affect the rates through screening and
 * Y_e, do not grow the network, and are rejected when they take part in a reaction.
 */
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/reaction/reaclib.h"

#include "fourdst/atomic/species.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <string>

namespace {
    using fourdst::atomic::Species;

    double sum(const std::vector<double>& v) {
        return std::accumulate(v.begin(), v.end(), 0.0);
    }

    size_t index_of(const gridfire::reaction::ReactionSet& reactions, const std::string& id) {
        for (size_t i = 0; i < reactions.size(); ++i) {
            if (reactions[i].id() == id) return i;
        }
        return reactions.size();
    }

    /// Solar composition plus a trace of Fe-56 standing in for the metals outside the network
    fourdst::composition::Composition solar_with_iron(const double Y_Fe) {
        auto comp = gridfire::test::solar_composition();
        comp.registerSpecies(fourdst::atomic::Fe_56);
        comp.setMolarAbundance(fourdst::atomic::Fe_56, Y_Fe);
        return comp;
    }

    gridfire::config::GridFireConfig passive_config(std::vector<std::string> names) {
        gridfire::config::GridFireConfig cfg;
        cfg.engine.passiveSpecies = std::move(names);
        return cfg;
    }
}

// ---------------------------------------------------------------- per-reaction contributions (requirement 8)

TEST(ReactionContributions, SumToTheEngineTotals) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    const auto& species = shared.engine.getNetworkSpecies(*ctx);
    gridfire::test::StateSampler sampler(species, gridfire::test::test_seed(), 0.01, 0.5, 1.0, 1.0e5);

    for (int n = 0; n < 20; ++n) {
        const auto state = sampler.next();
        const auto totals = shared.engine.calculateRHSAndEnergy(*ctx, state.composition, state.T9, state.rho, false);
        ASSERT_TRUE(totals.has_value()) << state.describe();
        const auto parts = shared.engine.calculateReactionContributions(*ctx, state.composition, state.T9, state.rho);
        ASSERT_TRUE(parts.has_value()) << "GraphEngine must provide per-reaction contributions";

        const size_t num_reactions = shared.engine.getNetworkReactions(*ctx).size();
        ASSERT_EQ(parts->molar_flow.size(), num_reactions);
        ASSERT_EQ(parts->energy_rate.size(), num_reactions);
        ASSERT_EQ(parts->neutrino_energy_loss_rate.size(), num_reactions);
        ASSERT_EQ(parts->neutrino_rate.size(), num_reactions);

        EXPECT_LT(gridfire::test::rel_diff(sum(parts->energy_rate), totals->nuclearEnergyGenerationRate, 1e-30), 1e-10)
            << state.describe();
        EXPECT_LT(gridfire::test::rel_diff(sum(parts->neutrino_energy_loss_rate), totals->neutrinoEnergyLossRate, 1e-30), 1e-10)
            << state.describe();
        EXPECT_LT(gridfire::test::rel_diff(sum(parts->neutrino_rate), totals->totalNeutrinoFlux, 1e-30), 1e-10)
            << state.describe();
        for (size_t j = 0; j < num_reactions; ++j) {
            EXPECT_GE(parts->neutrino_rate[j] * parts->molar_flow[j], 0.0) << "neutrino rate against the flow, reaction " << j;
        }
    }
}

// pep (p + e- + p -> d + nu) is its own reaction, an electron capture whose rate carries the rho * Y_e factor, and it
// runs at a few tenths of a percent of pp in the solar core (Bahcall: pep/pp ~ 0.2-0.5%).
TEST(ReactionContributions, PepIsSeparateFromPpInTheSun) {
    gridfire::test::quiet_logging();
    const auto pep = gridfire::reaclib::get_all_reaclib_reactions().get("p(e-p,)d");
    ASSERT_TRUE(pep.has_value()) << "pep is not in the REACLIB database under p(e-p,)d";
    EXPECT_TRUE(dynamic_cast<const gridfire::reaction::ReaclibReaction&>(*pep.value()).is_electron_capture());

    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    const auto& reactions = shared.engine.getNetworkReactions(*ctx);
    const size_t i_pp = index_of(reactions, "p(p,e+)d");
    const size_t i_pep = index_of(reactions, "p(e-p,)d");
    ASSERT_LT(i_pp, reactions.size());
    ASSERT_LT(i_pep, reactions.size());

    const auto parts = shared.engine.calculateReactionContributions(*ctx, gridfire::test::solar_composition(), 0.0157, 150.0);
    ASSERT_TRUE(parts.has_value());
    const double ratio = parts->molar_flow[i_pep] / parts->molar_flow[i_pp];
    EXPECT_GT(ratio, 1e-3);
    EXPECT_LT(ratio, 1e-2);
    // pep neutrinos are monoenergetic at 1.445 MeV, pp neutrinos average 0.267 MeV
    const double E_pep = parts->neutrino_energy_loss_rate[i_pep] / parts->neutrino_rate[i_pep];
    const double E_pp = parts->neutrino_energy_loss_rate[i_pp] / parts->neutrino_rate[i_pp];
    EXPECT_NEAR(E_pep / E_pp, 1.445 / 0.2668, 1e-3);
}

// ---------------------------------------------------------------- passive species (requirement 7)

class PassiveSpecies : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        gridfire::test::quiet_logging();
        policy = std::make_unique<gridfire::policy::MainSequencePolicy>(solar_with_iron(1e-5), passive_config({"Fe-56"}));
        auto results = policy->construct();
        engine = &results.engine;
        blob = std::move(results.scratch_blob);
    }
    static void TearDownTestSuite() {
        blob.reset();
        policy.reset();
    }

    static inline std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
    static inline const gridfire::engine::DynamicEngine* engine = nullptr;
    static inline std::unique_ptr<gridfire::engine::scratch::StateBlob> blob;
};

TEST_F(PassiveSpecies, AreNetworkSpeciesWithNoReactions) {
    const auto& species = engine->getNetworkSpecies(*blob);
    EXPECT_NE(std::ranges::find(species, fourdst::atomic::Fe_56), species.end());
    for (const auto& reaction : engine->getNetworkReactions(*blob)) {
        EXPECT_FALSE(reaction->contains(fourdst::atomic::Fe_56)) << reaction->id();
    }
}

TEST_F(PassiveSpecies, DoNotGrowTheNetwork) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    EXPECT_EQ(engine->getNetworkReactions(*blob).size(), shared.engine.getNetworkReactions(*ctx).size());
    EXPECT_EQ(engine->getNetworkSpecies(*blob).size(), shared.engine.getNetworkSpecies(*ctx).size() + 1);
}

TEST_F(PassiveSpecies, NeverChange) {
    const auto rhs = engine->calculateRHSAndEnergy(*blob, solar_with_iron(1e-5), 0.0157, 150.0, false);
    ASSERT_TRUE(rhs.has_value());
    EXPECT_EQ(rhs->dydt.at(fourdst::atomic::Fe_56), 0.0);
    EXPECT_NE(rhs->dydt.at(fourdst::atomic::H_1), 0.0);
}

// A passive species adds charge to the plasma, so it strengthens screening: more iron means faster burning. At solar
// conditions screening is weak (Salpeter), where ln f grows with zeta = sqrt(sum Z^2 Y + Y_e).
TEST_F(PassiveSpecies, ContributeToScreening) {
    const auto eps = [](const double Y_Fe) {
        const auto rhs = engine->calculateRHSAndEnergy(*blob, solar_with_iron(Y_Fe), 0.0157, 150.0, false);
        EXPECT_TRUE(rhs.has_value());
        return rhs->nuclearEnergyGenerationRate;
    };
    const double without = eps(0.0);
    const double with = eps(2e-3); // X_Fe ~ 0.1, an exaggerated metallicity so the effect is clear
    EXPECT_GT(with / without, 1.001);
    EXPECT_LT(with / without, 1.2);
}

TEST(PassiveSpeciesConfig, SpeciesInAReactionIsRejected) {
    gridfire::test::quiet_logging();
    gridfire::policy::MainSequencePolicy policy(gridfire::test::solar_composition(), passive_config({"C-12"}));
    EXPECT_ANY_THROW((void)policy.construct());
}

TEST(PassiveSpeciesConfig, UnknownSpeciesIsRejected) {
    gridfire::test::quiet_logging();
    gridfire::policy::MainSequencePolicy policy(gridfire::test::solar_composition(), passive_config({"Xx-999"}));
    EXPECT_ANY_THROW((void)policy.construct());
}
