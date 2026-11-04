// Unit tests for the energy channel classifier (gridfire::engine::classify_energy_channels).
#include "gf_test_support.h"

#include "gridfire/engine/procedures/energy_channels.h"
#include "gridfire/reaction/reaction.h"

#include "fourdst/atomic/species.h"

#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <vector>

using fourdst::atomic::Species;
using gridfire::engine::classify_energy_channels;
using gridfire::engine::EnergyChannels;
using gridfire::reaction::ReaclibReaction;
using gridfire::reaction::ReactionSet;

namespace {
    const Species& sp(const char* name) { return fourdst::atomic::species.at(name); }

    std::vector<Species> standard_fuels() {
        return {sp("H-1"), sp("He-4"), sp("C-12"), sp("Ne-20"), sp("O-16"), sp("Si-28")};
    }

    /// Builds a ReactionSet from (id, reactants, products) triples. Rates are irrelevant to classification.
    ReactionSet make_set(const std::vector<std::tuple<std::string, std::vector<Species>, std::vector<Species>>>& defs) {
        std::vector<std::unique_ptr<gridfire::reaction::Reaction>> reactions;
        for (const auto& [id, reactants, products] : defs) {
            reactions.push_back(std::make_unique<ReaclibReaction>(id, id, 4, reactants, products, 1.0, "test", gridfire::reaction::RateCoefficientSet{}));
        }
        return ReactionSet(std::move(reactions));
    }

    std::string channel_of(const ReactionSet& set, const EnergyChannels& ch, const std::string& id) {
        for (size_t i = 0; i < set.size(); ++i) {
            if (std::string(set[i].id()) == id) return ch.names[ch.reaction_channel[i]];
        }
        return "<missing>";
    }

    /// A small synthetic network exercising lineage, the lightest-reactant rule, and the pp/CNO split.
    ReactionSet synthetic_network() {
        return make_set({
            {"pp", {sp("H-1"), sp("H-1")}, {sp("H-2")}},
            {"dp", {sp("H-2"), sp("H-1")}, {sp("He-3")}},
            {"he3he3", {sp("He-3"), sp("He-3")}, {sp("He-4"), sp("H-1"), sp("H-1")}},
            {"he3he4", {sp("He-3"), sp("He-4")}, {sp("Be-7")}},
            {"c12p", {sp("C-12"), sp("H-1")}, {sp("N-13")}},
            {"n13decay", {sp("N-13")}, {sp("C-13")}},
            {"n14a", {sp("N-14"), sp("He-4")}, {sp("F-18")}},
            {"c12a", {sp("C-12"), sp("He-4")}, {sp("O-16")}},
            {"c12c12", {sp("C-12"), sp("C-12")}, {sp("Ne-20"), sp("He-4")}},
            {"o16o16", {sp("O-16"), sp("O-16")}, {sp("Si-28"), sp("He-4")}},
            {"li7li7", {sp("Li-7"), sp("Li-7")}, {sp("C-14")}},  // No fuel lineage for Li-7
        });
    }
}

TEST(EnergyChannelClassifier, ChannelNamesFollowFuelOrderWithOtherLast) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels ch = classify_energy_channels(set, standard_fuels(), true);
    const std::vector<std::string> expected = {"H_pp", "H_CNO", "He", "C", "Ne", "O", "Si", "other"};
    EXPECT_EQ(ch.names, expected);
}

TEST(EnergyChannelClassifier, EveryReactionGetsExactlyOneValidChannel) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels ch = classify_energy_channels(set, standard_fuels(), true);
    ASSERT_EQ(ch.reaction_channel.size(), set.size());
    for (const size_t c : ch.reaction_channel) {
        EXPECT_LT(c, ch.names.size());
    }
}

TEST(EnergyChannelClassifier, SyntheticNetworkLabels) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels ch = classify_energy_channels(set, standard_fuels(), true);
    EXPECT_EQ(channel_of(set, ch, "pp"), "H_pp");
    EXPECT_EQ(channel_of(set, ch, "dp"), "H_pp");           // d inherits H lineage
    EXPECT_EQ(channel_of(set, ch, "he3he3"), "H_pp");       // no H reactant, but He-3 is an H intermediate
    EXPECT_EQ(channel_of(set, ch, "he3he4"), "H_pp");       // lightest reactant He-3 wins over the He-4 fuel
    EXPECT_EQ(channel_of(set, ch, "c12p"), "H_CNO");
    EXPECT_EQ(channel_of(set, ch, "n13decay"), "H_CNO");    // N-13 inherits H lineage, heavy nucleus -> CNO
    EXPECT_EQ(channel_of(set, ch, "n14a"), "He");           // CNO ash + alpha is helium burning
    EXPECT_EQ(channel_of(set, ch, "c12a"), "He");
    EXPECT_EQ(channel_of(set, ch, "c12c12"), "C");
    EXPECT_EQ(channel_of(set, ch, "o16o16"), "O");
    EXPECT_EQ(channel_of(set, ch, "li7li7"), "other");
}

TEST(EnergyChannelClassifier, WithoutPPCNOSplitHydrogenIsOneChannel) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels ch = classify_energy_channels(set, standard_fuels(), false);
    EXPECT_EQ(ch.names.front(), "H");
    EXPECT_EQ(channel_of(set, ch, "pp"), "H");
    EXPECT_EQ(channel_of(set, ch, "c12p"), "H");
}

TEST(EnergyChannelClassifier, CustomFuelListChangesChannels) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels ch = classify_energy_channels(set, {sp("H-1"), sp("He-4")}, true);
    const std::vector<std::string> expected = {"H_pp", "H_CNO", "He", "other"};
    EXPECT_EQ(ch.names, expected);
    EXPECT_EQ(channel_of(set, ch, "c12c12"), "other");      // C-12 is no longer a fuel and nothing here produces it
    EXPECT_EQ(channel_of(set, ch, "o16o16"), "He");         // O-16 is no longer a fuel; it inherits He lineage from c12a
}

TEST(EnergyChannelClassifier, IsDeterministic) {
    const ReactionSet set = synthetic_network();
    const EnergyChannels a = classify_energy_channels(set, standard_fuels(), true);
    const EnergyChannels b = classify_energy_channels(set, standard_fuels(), true);
    EXPECT_EQ(a.names, b.names);
    EXPECT_EQ(a.reaction_channel, b.reaction_channel);
}

// --- The real main sequence network ---

class EnergyChannelMainSequence : public ::testing::Test {
protected:
    static void SetUpTestSuite() { gridfire::test::quiet_logging(); }
};

TEST_F(EnergyChannelMainSequence, CanonicalReactionsLandInTheirChannels) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto ctx = shared.fresh_context();
    const auto& reactions = shared.engine.getNetworkReactions(*ctx);
    const EnergyChannels ch = classify_energy_channels(reactions, standard_fuels(), true);

    const std::vector<std::pair<std::string, std::string>> expected = {
        {"p(p,e+)d", "H_pp"}, {"d(p,g)he3", "H_pp"}, {"he3(he3,2p)he4", "H_pp"}, {"he4(he3,g)be7", "H_pp"},
        {"be7(p,g)b8", "H_pp"}, {"b8(,e+ a)he4", "H_pp"},
        {"c12(p,g)n13", "H_CNO"}, {"n13(,e+)c13", "H_CNO"}, {"n15(p,a)c12", "H_CNO"}, {"o17(p,a)n14", "H_CNO"},
        {"he4(2a,)c12", "He"}, {"c12(a,g)o16", "He"}, {"n14(a,g)f18", "He"}, {"o16(a,g)ne20", "He"},
        {"c12(c12,a)ne20", "C"}, {"c12(c12,p)na23", "C"}, {"o16(o16,a)si28", "O"},
    };
    for (const auto& [id, want] : expected) {
        EXPECT_EQ(channel_of(reactions, ch, id), want) << "reaction " << id;
    }
}
