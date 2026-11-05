// Unit tests for electron screening (gridfire::screening).
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/screening/screening_bare.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/screening/screening_weak.h"

#include "fourdst/atomic/species.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

using fourdst::atomic::Species;

namespace {
    const Species& sp(const char* name) { return fourdst::atomic::species.at(name); }

    gridfire::reaction::ReactionSet two_body_set(const std::vector<std::pair<Species, Species>>& pairs) {
        std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
        int i = 0;
        for (const auto& [a, b] : pairs) {
            const std::string id = "pair" + std::to_string(i++);
            list.push_back(std::make_unique<gridfire::reaction::ReaclibReaction>(
                id, id, 4, std::vector{a, b}, std::vector<Species>{}, 0.0, "test", gridfire::reaction::RateCoefficientSet{}));
        }
        return gridfire::reaction::ReactionSet(std::move(list));
    }

    /// Debye-Hueckel (weak) screening from first principles, in CGS: H_12 = Z_1 Z_2 e^2 / (k T lambda_D) with
    /// lambda_D^-2 = 4 pi e^2 rho N_A zeta / (k T) and zeta = sum_i (Z_i^2 + Z_i) Y_i (ions plus unscreened electrons).
    double debye_huckel_H(const double Z1Z2, const double T, const double rho, const double zeta) {
        constexpr double e = 4.80320471e-10, k = 1.380649e-16, Na = 6.02214076e23, pi = 3.14159265358979;
        const double inv_lambda = std::sqrt(4.0 * pi * e * e * rho * Na * zeta / (k * T));
        return Z1Z2 * e * e * inv_lambda / (k * T);
    }
}

TEST(WeakScreening, MatchesDebyeHuckelFromFirstPrinciples) {
    const std::vector<Species> species = {sp("H-1"), sp("He-4"), sp("C-12"), sp("N-14")};
    const std::vector<double> Y = {0.70 / 1.007825, 0.28 / 4.002603, 0.012 / 12.0, 0.008 / 14.003074};
    double zeta = 0.0;
    for (size_t i = 0; i < species.size(); ++i) zeta += (species[i].z() * species[i].z() + species[i].z()) * Y[i];

    const auto reactions = two_body_set({{sp("H-1"), sp("H-1")}, {sp("He-4"), sp("He-4")}, {sp("C-12"), sp("H-1")}, {sp("N-14"), sp("H-1")}});
    const std::vector<double> Z1Z2 = {1, 4, 6, 7};
    const gridfire::screening::WeakScreeningModel model;
    // Solar core, a hotter/denser H-burning state, and helium burning: all in the weak regime (H << 1)
    for (const auto& [T, rho] : std::vector<std::pair<double, double>>{{1.5e7, 150.0}, {3e7, 50.0}, {1.2e8, 1e4}}) {
        const auto factors = model.calculateScreeningFactors(reactions, species, Y, T / 1e9, rho);
        for (size_t r = 0; r < Z1Z2.size(); ++r) {
            const double H = debye_huckel_H(Z1Z2[r], T, rho, zeta);
            ASSERT_LT(H, 0.5) << "test state should be weakly screened";
            // The Salpeter coefficient 0.188 is rounded, so allow 0.5% in H
            EXPECT_NEAR(std::log(factors[r]), H, 5e-3 * H) << "Z1Z2 = " << Z1Z2[r] << " at T = " << T << ", rho = " << rho;
        }
    }
}

TEST(WeakScreening, SolarCorePPEnhancementIsAFewPercent) {
    // Classic result (e.g. Bahcall): weak screening raises the solar-centre pp rate by ~5%
    const std::vector<Species> species = {sp("H-1"), sp("He-4")};
    const std::vector<double> Y = {0.35 / 1.007825, 0.65 / 4.002603};
    const auto factors = gridfire::screening::WeakScreeningModel().calculateScreeningFactors(
        two_body_set({{sp("H-1"), sp("H-1")}}), species, Y, 0.0157, 150.0);
    EXPECT_GT(factors[0], 1.03);
    EXPECT_LT(factors[0], 1.08);
}

TEST(Screening, BareFactorsAreOne) {
    const std::vector<Species> species = {sp("H-1"), sp("He-4")};
    const auto factors = gridfire::screening::BareScreeningModel().calculateScreeningFactors(
        two_body_set({{sp("H-1"), sp("H-1")}, {sp("He-4"), sp("He-4")}}), species, {0.7, 0.07}, 0.015, 150.0);
    for (const double f : factors) EXPECT_EQ(f, 1.0);
}

TEST(Screening, ConfigDefaultsToBlended) {
    EXPECT_EQ(gridfire::config::GridFireConfig{}.engine.screening, gridfire::screening::ScreeningType::BLENDED);
}

TEST(Screening, GraphEngineDefaultsToBlended) {
    const auto& shared = gridfire::test::main_sequence_engine();
    EXPECT_EQ(shared.graph.getScreeningModel(*shared.fresh_context()), gridfire::screening::ScreeningType::BLENDED);
}

TEST(WeakScreening, LogDerivativesMatchFiniteDifferences) {
    const std::vector<Species> species = {sp("H-1"), sp("He-4"), sp("C-12"), sp("N-14")};
    const std::vector<double> Y = {0.70 / 1.007825, 0.28 / 4.002603, 0.012 / 12.0, 0.008 / 14.003074};
    const auto reactions = two_body_set({{sp("H-1"), sp("H-1")}, {sp("C-12"), sp("H-1")}, {sp("N-14"), sp("H-1")}});
    const gridfire::screening::WeakScreeningModel model;
    const double T9 = 0.015, rho = 150.0;
    const auto dlnf = model.calculateScreeningLogDerivatives(reactions, species, Y, T9, rho);
    ASSERT_TRUE(dlnf.has_value());
    for (size_t k = 0; k < species.size(); ++k) {
        const double h = 1e-4 * Y[k];
        std::vector<double> up = Y, down = Y;
        up[k] += h;
        down[k] -= h;
        const auto f_up = model.calculateScreeningFactors(reactions, species, up, T9, rho);
        const auto f_down = model.calculateScreeningFactors(reactions, species, down, T9, rho);
        for (size_t r = 0; r < reactions.size(); ++r) {
            const double fd = (std::log(f_up[r]) - std::log(f_down[r])) / (2.0 * h);
            const double analytic = (*dlnf)[r * species.size() + k];
            EXPECT_NEAR(analytic, fd, 1e-6 * std::abs(fd)) << "reaction " << r << ", species " << species[k].name();
        }
    }
}

TEST(Screening, BareLogDerivativesAreZero) {
    const std::vector<Species> species = {sp("H-1"), sp("He-4")};
    const auto dlnf = gridfire::screening::BareScreeningModel().calculateScreeningLogDerivatives(
        two_body_set({{sp("H-1"), sp("H-1")}}), species, {0.7, 0.07}, 0.015, 150.0);
    ASSERT_TRUE(dlnf.has_value());
    for (const double d : *dlnf) EXPECT_EQ(d, 0.0);
}
