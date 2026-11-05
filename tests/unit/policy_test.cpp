// Unit tests for HydrogenHeliumBurningPolicy: a fixed network with helium burning, independent of the starting
// composition, and actually burning helium.
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/exceptions/error_policy.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/utils/gf_omp.h"

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils/utils.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

using gridfire::policy::HydrogenHeliumBurningPolicy;

namespace {
    const fourdst::atomic::Species& sp(const char* name) { return fourdst::atomic::species.at(name); }

    fourdst::composition::Composition from_mass_fractions(const std::vector<std::string>& names, const std::vector<double>& X) {
        return fourdst::composition::buildCompositionFromMassFractions(names, X);
    }

    struct Built {
        std::unique_ptr<HydrogenHeliumBurningPolicy> policy;
        std::unique_ptr<gridfire::policy::ConstructionResults> results;
        [[nodiscard]] const gridfire::engine::DynamicEngine& engine() const { return results->engine; }
        [[nodiscard]] gridfire::engine::scratch::StateBlob& ctx() const { return *results->scratch_blob; }
    };

    Built build(const fourdst::composition::Composition& composition, const gridfire::config::GridFireConfig& config = {}) {
        gridfire::test::quiet_logging();
        GF_PAR_INIT();
        Built b;
        b.policy = std::make_unique<HydrogenHeliumBurningPolicy>(composition, config);
        b.results = std::make_unique<gridfire::policy::ConstructionResults>(b.policy->construct());
        return b;
    }

    bool has_reaction(const Built& b, const std::string& id) {
        return b.engine().getNetworkReactions(b.ctx()).contains(id);
    }
}

TEST(HydrogenHeliumBurningPolicy, ContainsHydrogenAndHeliumBurning) {
    const Built b = build(gridfire::test::solar_composition());
    EXPECT_EQ(b.policy->get_status(), gridfire::policy::NetworkPolicyStatus::INITIALIZED_VERIFIED);
    for (const char* id : {"p(p,e+)d", "he3(he3,2p)he4", "c12(p,g)n13", "n14(p,g)o15", "he4(2a,)c12", "c12(a,g)o16",
                           "o16(a,g)ne20", "ne20(a,g)mg24", "n14(a,g)f18", "f18(,e+)o18", "o18(a,g)ne22", "c13(a,n)o16"}) {
        EXPECT_TRUE(has_reaction(b, id)) << id;
    }
    for (const auto& reaction : b.policy->get_seed_reactions()) {
        EXPECT_TRUE(b.engine().getNetworkReactions(b.ctx()).contains(*reaction)) << reaction->id();
    }
}

// The network is built from fixed seeds, so it is identical (species, order and reactions) whatever the composition
TEST(HydrogenHeliumBurningPolicy, NetworkIsIndependentOfTheStartingComposition) {
    const Built solar = build(gridfire::test::solar_composition());
    const Built pristine = build(from_mass_fractions({"H-1", "He-4"}, {0.75, 0.25}));
    const Built helium = build(from_mass_fractions({"He-4", "C-12", "O-16"}, {0.98, 0.01, 0.01}));
    const auto& species = solar.engine().getNetworkSpecies(solar.ctx());
    EXPECT_EQ(pristine.engine().getNetworkSpecies(pristine.ctx()), species);
    EXPECT_EQ(helium.engine().getNetworkSpecies(helium.ctx()), species);
    EXPECT_EQ(pristine.engine().getNetworkReactions(pristine.ctx()).size(), solar.engine().getNetworkReactions(solar.ctx()).size());
    for (const char* name : {"O-18", "Ne-22", "F-18", "Mg-24"}) {
        EXPECT_NE(std::ranges::find(species, sp(name)), species.end()) << name;
    }
}

// No species is required, unlike MainSequencePolicy (which needs Ne-20 and Mg-24 among others)
TEST(HydrogenHeliumBurningPolicy, NeedsNoParticularSeedSpecies) {
    const auto pristine = from_mass_fractions({"H-1", "He-4"}, {0.75, 0.25});
    EXPECT_NO_THROW(build(pristine));
    EXPECT_THROW(gridfire::policy::MainSequencePolicy{pristine}, gridfire::exceptions::MissingSeedSpeciesError);
}

TEST(HydrogenHeliumBurningPolicy, RejectsSpeciesOutsideTheNetwork) {
    const auto with_iron = from_mass_fractions({"H-1", "He-4", "Fe-56"}, {0.7, 0.29, 0.01});
    HydrogenHeliumBurningPolicy policy(with_iron);
    EXPECT_THROW((void)policy.construct(), gridfire::exceptions::PolicyError);
}

TEST(HydrogenHeliumBurningPolicy, PassesConfigToTheEngine) {
    gridfire::config::GridFireConfig config;
    config.engine.screening = gridfire::screening::ScreeningType::BARE;
    const Built b = build(gridfire::test::solar_composition(), config);
    EXPECT_EQ(b.engine().getScreeningModel(b.ctx()), gridfire::screening::ScreeningType::BARE);
}

// Core helium burning conditions: helium is turned into carbon and oxygen
TEST(HydrogenHeliumBurningPolicy, BurnsHelium) {
    const auto start = from_mass_fractions({"He-4", "C-12", "N-14", "O-16"}, {0.98, 0.003, 0.012, 0.005});
    const Built b = build(start);
    gridfire::solver::PointSolverContext sctx(b.ctx());
    sctx.set_stdout_logging(false);
    const gridfire::solver::PointSolver solver(b.engine());

    gridfire::NetIn in;
    in.composition = start;
    in.temperature = 1.5e8;
    in.density = 1e4;
    in.tMax = 3e13;
    in.dt0 = 1e-6;
    in.energy = 0.0;
    const gridfire::NetOut out = solver.evaluate(sctx, in);
    ASSERT_TRUE(out.converged);
    EXPECT_GT(out.energy, 0.0);
    EXPECT_GT(out.eps_nuc, 0.0);
    EXPECT_LT(out.composition.getMolarAbundance(sp("He-4")), start.getMolarAbundance(sp("He-4")));
    EXPECT_GT(out.composition.getMolarAbundance(sp("C-12")), start.getMolarAbundance(sp("C-12")));
    EXPECT_GT(out.composition.getMolarAbundance(sp("O-16")), start.getMolarAbundance(sp("O-16")));
    // N-14 is converted by alpha captures towards Ne-22
    EXPECT_LT(out.composition.getMolarAbundance(sp("N-14")), start.getMolarAbundance(sp("N-14")));
    EXPECT_GT(out.composition.getMolarAbundance(sp("Ne-22")), 0.0);
    const auto& names = out.energy_channel_names;
    const auto he = std::ranges::find(names, "He");
    ASSERT_NE(he, names.end());
    EXPECT_GT(out.eps_by_channel[he - names.begin()], 0.9 * out.eps_nuc) << "helium burning dominates eps";
}
