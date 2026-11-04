// Unit tests on the reaction data of the main sequence network and on small engine utilities.
#include "gf_test_support.h"

#include "gridfire/reaction/reaction.h"
#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/weak/neutrino_energies.h"
#include "gridfire/utils/hashing.h"

#include <gtest/gtest.h>

#include <bit>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>

using gridfire::reaction::ReactionType;

class MainSequenceNetworkData : public ::testing::Test {
protected:
    static void SetUpTestSuite() { gridfire::test::quiet_logging(); }
};

TEST_F(MainSequenceNetworkData, EveryReactionConservesBaryonNumber) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    for (const auto& reaction : shared.engine.getNetworkReactions(*ctx)) {
        int dA = 0;
        for (const auto& [species, coeff] : reaction->stoichiometry()) dA += coeff * species.a();
        EXPECT_EQ(dA, 0) << "reaction " << reaction->id();
    }
}

// Charge-changing reactions must be typed as weak: that typing decides neutrino losses and suppresses reverse rates.
// (Regression: two-body weak reactions such as p(p,e+)d were once typed as strong.)
TEST_F(MainSequenceNetworkData, OnlyWeakReactionsChangeNuclearCharge) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    for (const auto& reaction : shared.engine.getNetworkReactions(*ctx)) {
        int dZ = 0;
        for (const auto& [species, coeff] : reaction->stoichiometry()) dZ += coeff * species.z();
        const bool weak = reaction->type() == ReactionType::WEAK ||
                          reaction->type() == ReactionType::REACLIB_WEAK ||
                          reaction->type() == ReactionType::LOGICAL_REACLIB_WEAK;
        EXPECT_EQ(weak, dZ != 0) << "reaction " << reaction->id() << " changes Z by " << dZ << " but is typed as " << (weak ? "weak" : "strong");
    }
}

// The neutrino loss model reads the lepton from the REACLIB id, so every weak REACLIB reaction must name exactly one.
TEST(ReaclibDatabase, EveryWeakReactionNamesExactlyOneLepton) {
    gridfire::test::quiet_logging();
    size_t weak_count = 0;
    for (const auto& reaction : gridfire::reaclib::get_all_reaclib_reactions()) {
        const std::string id(reaction->id());
        const auto paren = id.find('(');
        const auto close = id.find(')');
        ASSERT_NE(paren, std::string::npos) << id;
        ASSERT_NE(close, std::string::npos) << id;
        const std::string inner = id.substr(paren + 1, close - paren - 1);
        size_t leptons = 0;
        for (size_t pos = inner.find("e+"); pos != std::string::npos; pos = inner.find("e+", pos + 1)) ++leptons;
        for (size_t pos = inner.find("e-"); pos != std::string::npos; pos = inner.find("e-", pos + 1)) ++leptons;
        if (gridfire::reaction::reaction_is_weak(*reaction)) {
            ++weak_count;
            EXPECT_EQ(leptons, 1u) << "weak reaction " << id;
        } else {
            EXPECT_EQ(leptons, 0u) << "strong reaction " << id << " names a lepton";
        }
    }
    EXPECT_GT(weak_count, 0u);
}

TEST_F(MainSequenceNetworkData, StoichiometryMatchesReactantAndProductLists) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    for (const auto& reaction : shared.engine.getNetworkReactions(*ctx)) {
        for (const auto& [species, coeff] : reaction->stoichiometry()) {
            const int consumed = static_cast<int>(std::ranges::count(reaction->reactants(), species));
            const int produced = static_cast<int>(std::ranges::count(reaction->products(), species));
            EXPECT_EQ(coeff, produced - consumed) << "reaction " << reaction->id() << ", species " << species.name();
        }
    }
}

TEST(StateHashing, AnySingleBitChangeInTemperatureOrDensityChangesTheHash) {
    gridfire::test::quiet_logging();
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    const auto& reactions = shared.engine.getNetworkReactions(*ctx);
    const auto comp = gridfire::test::solar_composition();
    const double T9 = 0.015, rho = 150.0;
    const std::size_t base = gridfire::utils::hash_state(comp, T9, rho, reactions);
    for (int bit = 0; bit < 64; ++bit) {
        const auto flip = [bit](const double x) { return std::bit_cast<double>(std::bit_cast<std::uint64_t>(x) ^ (1ULL << bit)); };
        EXPECT_NE(gridfire::utils::hash_state(comp, flip(T9), rho, reactions), base) << "T9 bit " << bit;
        EXPECT_NE(gridfire::utils::hash_state(comp, T9, flip(rho), reactions), base) << "rho bit " << bit;
    }
}

// Average neutrino energies must come from the tabulated values (Bahcall 2002, as in MESA), not the fraction-of-Q
// fallback, for the reactions that dominate solar neutrino losses.
TEST(NeutrinoEnergies, SolarReactionsUseBahcallValues) {
    gridfire::test::quiet_logging();
    const auto& all = gridfire::reaclib::get_all_reaclib_reactions();
    const std::vector<std::pair<std::string, double>> expected = {
        {"p(p,e+)d", 0.2668}, {"he3(p,e+)he4", 9.628}, {"b8(,e+ a)he4", 6.735}, {"be7(e-,)li7", 0.814},
        {"n13(,e+)c13", 0.706}, {"o15(,e+)n15", 0.9964},
    };
    for (const auto& [id, energy] : expected) {
        const auto reaction = all.get(id);
        ASSERT_TRUE(reaction.has_value()) << id << " is not in the REACLIB database";
        const auto tabulated = gridfire::reaction::weak::tabulated_neutrino_energy_MeV(*reaction.value());
        ASSERT_TRUE(tabulated.has_value()) << id << " has no tabulated neutrino energy";
        EXPECT_DOUBLE_EQ(*tabulated, energy) << id;
    }
}

TEST(NeutrinoEnergies, PepIsDistinctFromPp) {
    gridfire::test::quiet_logging();
    const auto& all = gridfire::reaclib::get_all_reaclib_reactions();
    for (const auto& reaction : all) {
        const std::string id(reaction->id());
        if (id.starts_with("p(") && id.ends_with(")d") && id.find("e-") != std::string::npos) {
            const auto tabulated = gridfire::reaction::weak::tabulated_neutrino_energy_MeV(*reaction);
            ASSERT_TRUE(tabulated.has_value()) << id;
            EXPECT_DOUBLE_EQ(*tabulated, 1.445) << id << " (pep)";
            return;
        }
    }
    GTEST_SKIP() << "No pep reaction in this REACLIB snapshot";
}

TEST(NeutrinoEnergies, StrongReactionsAreNotTabulated) {
    gridfire::test::quiet_logging();
    const auto reaction = gridfire::reaclib::get_all_reaclib_reactions().get("c12(p,g)n13");
    ASSERT_TRUE(reaction.has_value());
    EXPECT_FALSE(gridfire::reaction::weak::tabulated_neutrino_energy_MeV(*reaction.value()).has_value());
}

// Every pp-chain and (hot) CNO weak reaction in the network must use a tabulated neutrino energy; the fraction-of-Q
// fallback is reserved for reactions that are energetically negligible in hydrogen burning (they are reported).
TEST_F(MainSequenceNetworkData, HydrogenBurningWeakReactionsHaveTabulatedNeutrinoEnergies) {
    const std::set<std::string> h_burning_parents = {"H-1", "He-3", "Be-7", "B-8", "C-11", "N-13", "O-14", "O-15",
                                                     "F-17", "F-18", "Ne-18", "Ne-19"};
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    std::string fallback;
    size_t checked = 0;
    for (const auto& reaction : shared.engine.getNetworkReactions(*ctx)) {
        if (!gridfire::reaction::reaction_is_weak(*reaction)) continue;
        const bool tabulated = gridfire::reaction::weak::tabulated_neutrino_energy_MeV(*reaction).has_value();
        const auto& parent = std::ranges::max(reaction->reactants(), {}, [](const auto& sp) { return sp.a(); });
        const bool endothermic = reaction->qValue() < 0.0;
        if (h_burning_parents.contains(std::string(parent.name())) && !endothermic) {
            ++checked;
            EXPECT_TRUE(tabulated) << reaction->id() << " uses the fraction-of-Q neutrino fallback";
        }
        if (!tabulated) fallback += " " + std::string(reaction->id());
    }
    EXPECT_GT(checked, 5u);
    if (!fallback.empty()) {
        std::cout << "Weak reactions using the fraction-of-Q neutrino fallback:" << fallback << "\n";
    }
}

// REACLIB electron-capture rates are tabulated per unit electron density; the physical rate is the fit times rho * Ye.
// (Regression: without the factor 7Be electron capture was ~100x too slow in the solar core, inflating the 8B branch.)
TEST(ReaclibRates, ElectronCaptureRatesScaleWithElectronDensity) {
    gridfire::test::quiet_logging();
    const auto& all = gridfire::reaclib::get_all_reaclib_reactions();
    const auto be7_ec = all.get("be7(e-,)li7");
    const auto pp = all.get("p(p,e+)d");
    ASSERT_TRUE(be7_ec.has_value());
    ASSERT_TRUE(pp.has_value());
    const auto& ec = dynamic_cast<const gridfire::reaction::ReaclibReaction&>(*be7_ec.value());
    const auto& beta = dynamic_cast<const gridfire::reaction::ReaclibReaction&>(*pp.value());
    EXPECT_TRUE(ec.is_electron_capture());
    EXPECT_FALSE(beta.is_electron_capture());

    const double T9 = 0.015;
    const auto rate = [T9](const auto& r, const double rho, const double Ye) { return r.calculate_rate(T9, rho, Ye, 0.0, {}, {}); };
    EXPECT_NEAR(rate(ec, 150.0, 0.85) / rate(ec, 1.0, 1.0), 150.0 * 0.85, 1e-10);
    EXPECT_NEAR(rate(ec, 300.0, 0.5) / rate(ec, 150.0, 0.85), (300.0 * 0.5) / (150.0 * 0.85), 1e-12);
    EXPECT_DOUBLE_EQ(rate(beta, 150.0, 0.85), rate(beta, 1.0, 1.0));
}
