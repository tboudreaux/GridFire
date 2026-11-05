// Unit tests for CD09ScreeningModel (Chugunov & DeWitt 2009, PRC 80, 014611) against an independent implementation of
// the published equations, its limits, its analytic derivatives, and its use through the engine.
#include "gf_test_support.h"

#include "gridfire/config/config.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/screening/screening_CD09.h"
#include "gridfire/screening/screening_blended.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/screening/screening_weak.h"

#include "fourdst/atomic/species.h"
#include "fourdst/constants/const.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using fourdst::atomic::Species;
using gridfire::screening::CD09ScreeningModel;

namespace {
    const Species& sp(const char* name) { return fourdst::atomic::species.at(name); }

    gridfire::reaction::ReactionSet make_set(const std::vector<std::vector<Species>>& reactant_lists, const std::string& prefix = "r") {
        std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
        int i = 0;
        for (const auto& reactants : reactant_lists) {
            const std::string id = prefix + std::to_string(i++);
            list.push_back(std::make_unique<gridfire::reaction::ReaclibReaction>(
                id, id, 4, reactants, std::vector<Species>{}, 0.0, "test", gridfire::reaction::RateCoefficientSet{}));
        }
        return gridfire::reaction::ReactionSet(std::move(list));
    }

    // ---------------- Independent reference: CD09 as published, with an exact electron degeneracy factor ----------
    struct Constants {
        double e2, kB, hbar, h, me, NA, u;
    };
    const Constants& k() {
        static const Constants c = [] {
            const auto& t = fourdst::constant::Constants::getInstance();
            return Constants{t.get("alpha").value * t.get("hbar").value * t.get("c").value, t.get("kB").value,
                             t.get("hbar").value, t.get("h").value, t.get("mE").value, t.get("N_a").value, t.get("u").value};
        }();
        return c;
    }

    /// CD09 eq. (24): Potekhin & Chabrier (2000)
    double f0(const double g) {
        constexpr double A1 = -0.907, A2 = 0.62954, B1 = 0.00456, B2 = 211.6, B3 = -1e-4, B4 = 0.00462;
        const double A3 = -std::sqrt(3.0) / 2 - A1 / std::sqrt(A2);
        return A1 * (std::sqrt(g * (A2 + g)) - A2 * std::log(std::sqrt(g / A2) + std::sqrt(1 + g / A2)))
             + 2 * A3 * (std::sqrt(g) - std::atan(std::sqrt(g))) + B1 * (g - B2 * std::log(1 + g / B2))
             + B3 / 2 * std::log(1 + g * g / B4);
    }

    /// Normalised Fermi-Dirac integrals F_{1/2}, F_{-1/2} by quadrature (x = t^2)
    void fermi_dirac(const double eta, double& f12, double& fm12) {
        const int n = 4000;
        const double tmax = std::sqrt(std::max(eta, 0.0) + 60.0), dt = tmax / n;
        double a = 0, b = 0;
        for (int i = 0; i <= n; ++i) {
            const double t = i * dt, w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2), x = t * t - eta;
            if (x > 700) continue;
            const double occupation = 1.0 / (std::exp(x) + 1.0);
            a += w * t * t * occupation;
            b += w * occupation;
        }
        f12 = a * dt / 3 * 4 / std::sqrt(std::numbers::pi);
        fm12 = b * dt / 3 * 2 / std::sqrt(std::numbers::pi);
    }

    /// Exact theta_e = F_{-1/2}(eta) / F_{1/2}(eta) at F_{1/2}(eta) = psi = n_e lambda_e^3 / 2
    double theta_exact(const double rho, const double Ye, const double T) {
        const double psi = 0.5 * rho * k().NA * Ye * std::pow(k().h, 3) / std::pow(2 * std::numbers::pi * k().me * k().kB * T, 1.5);
        double lo = -60, hi = 3000, f12 = 0, fm12 = 0;
        for (int it = 0; it < 70; ++it) {
            const double mid = 0.5 * (lo + hi);
            fermi_dirac(mid, f12, fm12);
            (f12 < psi ? lo : hi) = mid;
        }
        fermi_dirac(0.5 * (lo + hi), f12, fm12);
        return fm12 / f12;
    }

    struct Comp { std::vector<Species> species; std::vector<double> Y; };

    /// h for the pair (Zi, Ai) + (Zj, Aj). theta_e < 0 selects CD09's own ion-only weak limit (eq. A3).
    double h_reference(const double Zi, const double Ai, const double Zj, const double Aj, const Comp& c, const double T,
                       const double rho, const double theta_e) {
        double SZ = 0, SZ2 = 0;
        for (size_t i = 0; i < c.species.size(); ++i) {
            SZ += c.species[i].z() * c.Y[i];
            SZ2 += c.species[i].z() * c.species[i].z() * c.Y[i];
        }
        const double ae = std::cbrt(3.0 / (4 * std::numbers::pi * rho * k().NA * SZ));
        const double Ge = k().e2 / (ae * k().kB * T);
        const double Gij = Ge * 2 * Zi * Zj / (std::cbrt(Zi) + std::cbrt(Zj));
        const double mu = Ai * Aj / (Ai + Aj) * k().u;
        const double tau = std::cbrt(27 * std::numbers::pi * std::numbers::pi * mu * std::pow(Zi * Zj * k().e2, 2) /
                                     (2 * k().kB * T * k().hbar * k().hbar));
        const double zeta = 3 * Gij / tau;
        const double y = 4 * Zi * Zj / ((Zi + Zj) * (Zi + Zj));
        const double c1 = 0.013 * y * y, c2 = 0.406 * std::pow(y, 0.14), c3 = 0.062 * std::pow(y, 0.19) + 1.8 / Gij;
        const double t = std::cbrt(1 + c1 * zeta + c2 * zeta * zeta + c3 * zeta * zeta * zeta);
        const double S = f0(Ge * std::pow(Zi, 5.0 / 3) / t) + f0(Ge * std::pow(Zj, 5.0 / 3) / t) - f0(Ge * std::pow(Zi + Zj, 5.0 / 3) / t);
        const double eta = std::pow(Zi + Zj, 2.5) - std::pow(Zi, 2.5) - std::pow(Zj, 2.5);
        double C;
        if (theta_e < 0) {
            C = 3 * Zi * Zj * std::sqrt(SZ2 / SZ) / eta; // CD09 eq. (A3)
        } else {
            const double hWeak = Zi * Zj * 0.188 * std::sqrt(rho * (SZ2 + theta_e * SZ) / std::pow(T / 1e6, 3));
            C = hWeak / (std::pow(Ge, 1.5) * eta / std::sqrt(3.0));
        }
        return (C + Gij * Gij) / (1 + Gij * Gij) * S;
    }

    double ln_reference(const std::vector<Species>& reactants, const Comp& c, const double T, const double rho, const double theta_e) {
        if (reactants.size() == 3) {
            const double mHe = sp("He-4").mass();
            return h_reference(2, mHe, 2, mHe, c, T, rho, theta_e) + h_reference(2, mHe, 4, sp("Be-8").mass(), c, T, rho, theta_e);
        }
        return h_reference(reactants[0].z(), reactants[0].mass(), reactants[1].z(), reactants[1].mass(), c, T, rho, theta_e);
    }

    Comp composition(const std::vector<std::pair<const char*, double>>& X) {
        Comp c;
        for (const auto& [name, x] : X) {
            c.species.push_back(sp(name));
            c.Y.push_back(x / sp(name).mass());
        }
        return c;
    }

    struct State {
        const char* label;
        double T, rho;
        std::vector<std::pair<const char*, double>> X;
        std::vector<std::vector<const char*>> reactions;
    };

    const std::vector<State>& states() {
        static const std::vector<State> s = {
            {"solar core", 1.57e7, 150, {{"H-1", 0.34}, {"He-4", 0.64}, {"C-12", 0.0005}, {"N-14", 0.006}, {"O-16", 0.0095}},
             {{"H-1", "H-1"}, {"H-1", "C-12"}, {"H-1", "N-14"}, {"He-3", "He-4"}}},
            {"0.3 Msun core", 8e6, 1e3, {{"H-1", 0.6}, {"He-4", 0.39}, {"N-14", 0.01}}, {{"H-1", "H-1"}, {"H-1", "N-14"}}},
            {"HB core", 1.2e8, 1e4, {{"He-4", 0.9}, {"C-12", 0.07}, {"O-16", 0.03}}, {{"He-4", "He-4", "He-4"}, {"He-4", "C-12"}}},
            {"RGB-tip flash", 1e8, 1e6, {{"He-4", 0.98}, {"C-12", 0.01}, {"N-14", 0.01}}, {{"He-4", "He-4", "He-4"}, {"He-4", "N-14"}}},
            {"C burning", 8e8, 2e5, {{"C-12", 0.25}, {"O-16", 0.73}, {"Ne-20", 0.02}}, {{"C-12", "C-12"}, {"O-16", "O-16"}}},
            {"WD C ignition", 3e8, 3e9, {{"C-12", 0.5}, {"O-16", 0.5}}, {{"C-12", "C-12"}, {"C-12", "O-16"}}},
        };
        return s;
    }

    std::vector<std::vector<Species>> reactant_lists(const State& s) {
        std::vector<std::vector<Species>> lists;
        for (const auto& r : s.reactions) {
            std::vector<Species> v;
            for (const char* n : r) v.push_back(sp(n));
            lists.push_back(v);
        }
        return lists;
    }
}

// The model is CD09 eqs. (23)-(25) and (A4) with a weak limit that includes the electrons, weighted by their
// degeneracy. Against the same equations with an exact (quadrature) theta_e, the only difference is the closed-form
// theta_e (< 0.7%), which moves ln f by well under 1%.
TEST(CD09Screening, MatchesPublishedEquationsWithExactElectronDegeneracy) {
    gridfire::test::quiet_logging();
    const CD09ScreeningModel model;
    for (const auto& s : states()) {
        const Comp c = composition(s.X);
        double Ye = 0;
        for (size_t i = 0; i < c.species.size(); ++i) Ye += c.species[i].z() * c.Y[i];
        const double theta = theta_exact(s.rho, Ye, s.T);
        const auto lists = reactant_lists(s);
        const auto factors = model.calculateScreeningFactors(make_set(lists), c.species, c.Y, s.T / 1e9, s.rho);
        for (size_t r = 0; r < lists.size(); ++r) {
            const double expected = ln_reference(lists[r], c, s.T, s.rho, theta);
            EXPECT_NEAR(std::log(factors[r]), expected, 5e-3 * std::abs(expected)) << s.label << ", reaction " << r;
        }
    }
}

// With strongly degenerate electrons (theta_e -> 0) the model is CD09 exactly as published, including eq. (A3)
TEST(CD09Screening, DegenerateElectronsReproduceThePaper) {
    gridfire::test::quiet_logging();
    const CD09ScreeningModel model;
    const State& wd = states().back();
    const Comp c = composition(wd.X);
    const auto lists = reactant_lists(wd);
    const auto factors = model.calculateScreeningFactors(make_set(lists), c.species, c.Y, wd.T / 1e9, wd.rho);
    for (size_t r = 0; r < lists.size(); ++r) {
        const double paper = ln_reference(lists[r], c, wd.T, wd.rho, -1.0);
        EXPECT_NEAR(std::log(factors[r]), paper, 1e-4 * paper) << "reaction " << r;
    }
}

// Weakly coupled, non-degenerate plasma: CD09 reduces to the Salpeter weak screening of WeakScreeningModel
TEST(CD09Screening, WeakLimitIsSalpeter) {
    gridfire::test::quiet_logging();
    const CD09ScreeningModel cd09;
    const gridfire::screening::WeakScreeningModel weak;
    const Comp c = composition({{"H-1", 0.7}, {"He-4", 0.28}, {"C-12", 0.01}, {"O-16", 0.01}});
    const auto set = make_set({{sp("H-1"), sp("H-1")}, {sp("H-1"), sp("C-12")}, {sp("He-4"), sp("C-12")}});
    const auto f_cd09 = cd09.calculateScreeningFactors(set, c.species, c.Y, 0.05, 1e-2);
    const auto f_weak = weak.calculateScreeningFactors(set, c.species, c.Y, 0.05, 1e-2);
    for (size_t r = 0; r < f_cd09.size(); ++r) {
        ASSERT_LT(std::log(f_weak[r]), 1e-3) << "test state should be weakly coupled";
        EXPECT_NEAR(std::log(f_cd09[r]), std::log(f_weak[r]), 0.02 * std::log(f_weak[r])) << "reaction " << r;
    }
}

TEST(CD09Screening, AnalyticDerivativesMatchFiniteDifferences) {
    gridfire::test::quiet_logging();
    const CD09ScreeningModel model;
    for (const auto& s : states()) {
        const Comp c = composition(s.X);
        std::vector<std::vector<Species>> lists;
        for (size_t i = 0; i < c.species.size(); ++i)
            for (size_t j = i; j < c.species.size(); ++j) lists.push_back({c.species[i], c.species[j]});
        lists.push_back({sp("He-4"), sp("He-4"), sp("He-4")});
        const auto set = make_set(lists);
        const double T9 = s.T / 1e9;
        const auto lnf = [&](const std::vector<double>& Y, const double t9, const double rho) {
            auto f = model.calculateScreeningFactors(set, c.species, Y, t9, rho);
            for (auto& v : f) v = std::log(v);
            return f;
        };
        const auto dY = *model.calculateScreeningLogDerivatives(set, c.species, c.Y, T9, s.rho);
        const auto dState = *model.calculateScreeningStateLogDerivatives(set, c.species, c.Y, T9, s.rho);
        for (size_t k = 0; k < c.species.size(); ++k) {
            const double h = 1e-5 * c.Y[k];
            auto up = c.Y, down = c.Y;
            up[k] += h;
            down[k] -= h;
            const auto a = lnf(up, T9, s.rho), b = lnf(down, T9, s.rho);
            for (size_t r = 0; r < lists.size(); ++r) {
                const double fd = (a[r] - b[r]) / (2 * h);
                EXPECT_NEAR(dY[r * c.species.size() + k], fd, 1e-5 * std::abs(fd) + 1e-10) << s.label << " r" << r << " Y" << k;
            }
        }
        const double eps = 1e-5;
        const auto tp = lnf(c.Y, T9 * (1 + eps), s.rho), tm = lnf(c.Y, T9 * (1 - eps), s.rho);
        const auto rp = lnf(c.Y, T9, s.rho * (1 + eps)), rm = lnf(c.Y, T9, s.rho * (1 - eps));
        for (size_t r = 0; r < lists.size(); ++r) {
            const double fdT = (tp[r] - tm[r]) / (2 * eps), fdRho = (rp[r] - rm[r]) / (2 * eps);
            EXPECT_NEAR(dState[r].dlnT, fdT, 1e-5 * std::abs(fdT) + 1e-10) << s.label << " r" << r;
            EXPECT_NEAR(dState[r].dlnRho, fdRho, 1e-5 * std::abs(fdRho) + 1e-10) << s.label << " r" << r;
        }
        // The CppAD path (used to record the engine's tape) computes the same factors
        const std::vector<CppAD::AD<double>> Yad(c.Y.begin(), c.Y.end());
        const auto f_ad = model.calculateScreeningFactors(set, c.species, Yad, CppAD::AD<double>(T9), CppAD::AD<double>(s.rho));
        const auto f = model.calculateScreeningFactors(set, c.species, c.Y, T9, s.rho);
        for (size_t r = 0; r < lists.size(); ++r) {
            EXPECT_NEAR(CppAD::Value(f_ad[r]), f[r], 1e-13 * f[r]) << s.label << " r" << r;
        }
    }
}

// The pair cache is shared across threads and keyed by reaction address, but a new reaction at a reused address must
// never pick up the old reaction's pairs
TEST(CD09Screening, CacheIsThreadSafeAndSurvivesReusedAddresses) {
    gridfire::test::quiet_logging();
    const CD09ScreeningModel model;
    const Comp c = composition({{"H-1", 0.7}, {"He-4", 0.28}, {"C-12", 0.01}, {"O-16", 0.01}});
    std::vector<double> expected;
    {
        const auto set = make_set({{sp("H-1"), sp("C-12")}, {sp("He-4"), sp("O-16")}});
        expected = CD09ScreeningModel().calculateScreeningFactors(set, c.species, c.Y, 0.02, 100.0);
        (void)model.calculateScreeningFactors(make_set({{sp("H-1"), sp("H-1")}, {sp("He-4"), sp("He-4")}}), c.species, c.Y, 0.02, 100.0);
    }
    for (int round = 0; round < 50; ++round) {
        const auto set = make_set({{sp("H-1"), sp("C-12")}, {sp("He-4"), sp("O-16")}});
        EXPECT_EQ(model.calculateScreeningFactors(set, c.species, c.Y, 0.02, 100.0), expected) << "round " << round;
    }

    const auto shared_set = make_set({{sp("H-1"), sp("H-1")}, {sp("H-1"), sp("N-14")}, {sp("He-4"), sp("C-12")}, {sp("C-12"), sp("C-12")}});
    const auto serial = model.calculateScreeningFactors(shared_set, c.species, c.Y, 0.02, 100.0);
    std::vector<std::thread> threads;
    std::vector<int> mismatches(8, 0);
    for (int t = 0; t < 8; ++t) {
        threads.emplace_back([&, t] {
            const CD09ScreeningModel* m = &model;
            for (int i = 0; i < 200; ++i) {
                const auto local = make_set({{sp("He-4"), sp("N-14")}, {sp("H-1"), sp("O-16")}}, "t" + std::to_string(t));
                (void)m->calculateScreeningFactors(local, c.species, c.Y, 0.02, 100.0);
                if (m->calculateScreeningFactors(shared_set, c.species, c.Y, 0.02, 100.0) != serial) ++mismatches[t];
            }
        });
    }
    for (auto& thread : threads) thread.join();
    for (int t = 0; t < 8; ++t) EXPECT_EQ(mismatches[t], 0) << "thread " << t;
}

TEST(CD09Screening, SelectableThroughTheConfiguration) {
    gridfire::test::quiet_logging();
    EXPECT_NE(dynamic_cast<CD09ScreeningModel*>(gridfire::screening::selectScreeningModel(gridfire::screening::ScreeningType::CD09).get()), nullptr);

    gridfire::config::GridFireConfig config;
    config.engine.screening = gridfire::screening::ScreeningType::CD09;
    gridfire::policy::MainSequencePolicy policy(gridfire::test::solar_composition(), config);
    const auto built = policy.construct();
    const auto& graph = dynamic_cast<const gridfire::engine::GraphEngine&>(built.engine);
    EXPECT_EQ(graph.getScreeningModel(*built.scratch_blob), gridfire::screening::ScreeningType::CD09);

    // The engine's analytic Jacobian and eps derivatives use CD09's analytic screening derivatives
    const auto& species = graph.getNetworkSpecies(*built.scratch_blob);
    const auto comp = gridfire::test::make_composition(species, gridfire::test::abundances_in_order(gridfire::test::solar_composition(), species));
    ASSERT_TRUE(graph.calculateAnalyticEpsDerivatives(comp, 0.0157, 150.0).has_value());
    const auto analytic = graph.generateAnalyticJacobianMatrix(*built.scratch_blob, comp, 0.0157, 150.0);
    ASSERT_TRUE(analytic.has_value());
    const auto autodiff = graph.generateAutoDiffJacobianMatrix(*built.scratch_blob, comp, 0.0157, 150.0);
    double worst = 0, largest = 0;
    for (const auto& row : species) {
        for (const auto& col : species) {
            const double a = (*analytic)(row, col), b = autodiff(row, col);
            worst = std::max(worst, std::abs(a - b));
            largest = std::max(largest, std::abs(b));
        }
    }
    EXPECT_LT(worst, 1e-10 * largest) << "analytic and CppAD Jacobians disagree under CD09 screening";
    const auto epsA = *graph.calculateAnalyticEpsDerivatives(comp, 0.0157, 150.0);
    const auto epsAD = graph.calculateAutoDiffEpsDerivatives(*built.scratch_blob, comp, 0.0157, 150.0);
    EXPECT_NEAR(epsA.dEps_dT, epsAD.dEps_dT, 1e-9 * std::abs(epsAD.dEps_dT));
    EXPECT_NEAR(epsA.dEps_dRho, epsAD.dEps_dRho, 1e-9 * std::abs(epsAD.dEps_dRho));
}

// ---------------------------------------------------------------------------------------------------------------------
// BlendedScreeningModel: Salpeter where weak screening holds, CD09 where it does not
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    std::unique_ptr<gridfire::screening::ScreeningModel> blended(const gridfire::screening::ScreeningBlendWindow& w = {}) {
        return gridfire::screening::selectScreeningModel(gridfire::screening::ScreeningType::BLENDED, w);
    }
}

// In the Sun every pp-chain and CNO pair has H < 0.5, so the blend is exactly Salpeter (as the literature supports)
TEST(BlendedScreening, SolarCoreIsSalpeter) {
    gridfire::test::quiet_logging();
    const State& sun = states().front();
    const Comp c = composition(sun.X);
    const auto set = make_set(reactant_lists(sun));
    const auto f = blended()->calculateScreeningFactors(set, c.species, c.Y, sun.T / 1e9, sun.rho);
    const auto weak = gridfire::screening::WeakScreeningModel().calculateScreeningFactors(set, c.species, c.Y, sun.T / 1e9, sun.rho);
    EXPECT_EQ(f, weak);
}

// Strongly coupled, degenerate matter: the blend is exactly CD09
TEST(BlendedScreening, WhiteDwarfIsCD09) {
    gridfire::test::quiet_logging();
    const State& wd = states().back();
    const Comp c = composition(wd.X);
    const auto set = make_set(reactant_lists(wd));
    const auto f = blended()->calculateScreeningFactors(set, c.species, c.Y, wd.T / 1e9, wd.rho);
    const auto cd09 = CD09ScreeningModel().calculateScreeningFactors(set, c.species, c.Y, wd.T / 1e9, wd.rho);
    for (size_t r = 0; r < f.size(); ++r) EXPECT_NEAR(f[r], cd09[r], 1e-12 * cd09[r]) << "reaction " << r;
}

// Inside the window the exponents are blended with the smoothstep weight in ln H
TEST(BlendedScreening, BlendsExponentsInsideTheWindow) {
    gridfire::test::quiet_logging();
    const gridfire::screening::BlendedScreeningModel model(std::make_unique<gridfire::screening::WeakScreeningModel>(),
                                                           std::make_unique<CD09ScreeningModel>());
    EXPECT_EQ(model.weight(0.3), 0.0);
    EXPECT_EQ(model.weight(0.5), 0.0);
    EXPECT_EQ(model.weight(1.5), 1.0);
    EXPECT_EQ(model.weight(1.9), 1.0);
    EXPECT_NEAR(model.weight(std::sqrt(0.5 * 1.5)), 0.5, 1e-12) << "halfway in ln H";
    int blended_reactions = 0;
    for (const auto& s : states()) {
        const Comp c = composition(s.X);
        std::vector<std::vector<Species>> pairs; // every pair of species present, so some land inside the window
        for (size_t i = 0; i < c.species.size(); ++i)
            for (size_t j = i; j < c.species.size(); ++j) pairs.push_back({c.species[i], c.species[j]});
        const auto set = make_set(pairs);
        const auto f = model.calculateScreeningFactors(set, c.species, c.Y, s.T / 1e9, s.rho);
        const auto weak = gridfire::screening::WeakScreeningModel().calculateScreeningFactors(set, c.species, c.Y, s.T / 1e9, s.rho);
        const auto cd09 = CD09ScreeningModel().calculateScreeningFactors(set, c.species, c.Y, s.T / 1e9, s.rho);
        for (size_t r = 0; r < f.size(); ++r) {
            const double H = std::log(weak[r]), w = model.weight(H);
            if (w > 0.0 && w < 1.0) ++blended_reactions;
            EXPECT_NEAR(std::log(f[r]), (1 - w) * H + w * std::log(cd09[r]), 1e-12) << s.label << " r" << r;
        }
    }
    EXPECT_GT(blended_reactions, 0) << "no test state lands inside the window";
}

// ln f and its derivative are continuous across the window edges: sweep density through the window for p + N-14
TEST(BlendedScreening, ContinuousAcrossTheWindow) {
    gridfire::test::quiet_logging();
    const auto model = blended();
    const Comp c = composition({{"H-1", 0.6}, {"He-4", 0.39}, {"N-14", 0.01}});
    const auto set = make_set({{sp("H-1"), sp("N-14")}});
    double previous = std::nan(""), previous_slope = std::nan("");
    const double step = 0.002; // in ln rho
    double worst_jump = 0, worst_slope_jump = 0;
    for (double lnRho = std::log(1.0); lnRho < std::log(3e4); lnRho += step) {
        const double value = std::log(model->calculateScreeningFactors(set, c.species, c.Y, 0.01, std::exp(lnRho))[0]);
        const double slope = (*model->calculateScreeningStateLogDerivatives(set, c.species, c.Y, 0.01, std::exp(lnRho)))[0].dlnRho;
        if (!std::isnan(previous)) {
            worst_jump = std::max(worst_jump, std::abs(value - previous - slope * step) / step);
            worst_slope_jump = std::max(worst_slope_jump, std::abs(slope - previous_slope));
        }
        previous = value;
        previous_slope = slope;
    }
    EXPECT_LT(worst_jump, 0.05) << "ln f jumps (or its derivative disagrees with it) across the window";
    EXPECT_LT(worst_slope_jump, 0.05) << "d ln f / d ln rho jumps across the window";
}

TEST(BlendedScreening, AnalyticDerivativesMatchFiniteDifferences) {
    gridfire::test::quiet_logging();
    const auto model = blended();
    for (const auto& s : states()) {
        const Comp c = composition(s.X);
        std::vector<std::vector<Species>> lists;
        for (size_t i = 0; i < c.species.size(); ++i)
            for (size_t j = i; j < c.species.size(); ++j) lists.push_back({c.species[i], c.species[j]});
        lists.push_back({sp("He-4"), sp("He-4"), sp("He-4")});
        const auto set = make_set(lists);
        const double T9 = s.T / 1e9;
        const auto lnf = [&](const std::vector<double>& Y, const double t9, const double rho) {
            auto f = model->calculateScreeningFactors(set, c.species, Y, t9, rho);
            for (auto& v : f) v = std::log(v);
            return f;
        };
        const auto dY = *model->calculateScreeningLogDerivatives(set, c.species, c.Y, T9, s.rho);
        const auto dState = *model->calculateScreeningStateLogDerivatives(set, c.species, c.Y, T9, s.rho);
        for (size_t k = 0; k < c.species.size(); ++k) {
            const double h = 1e-6 * c.Y[k];
            auto up = c.Y, down = c.Y;
            up[k] += h;
            down[k] -= h;
            const auto a = lnf(up, T9, s.rho), b = lnf(down, T9, s.rho);
            for (size_t r = 0; r < lists.size(); ++r) {
                const double fd = (a[r] - b[r]) / (2 * h);
                EXPECT_NEAR(dY[r * c.species.size() + k], fd, 1e-5 * std::abs(fd) + 1e-9) << s.label << " r" << r << " Y" << k;
            }
        }
        const double eps = 1e-6;
        const auto tp = lnf(c.Y, T9 * (1 + eps), s.rho), tm = lnf(c.Y, T9 * (1 - eps), s.rho);
        const auto rp = lnf(c.Y, T9, s.rho * (1 + eps)), rm = lnf(c.Y, T9, s.rho * (1 - eps));
        for (size_t r = 0; r < lists.size(); ++r) {
            const double fdT = (tp[r] - tm[r]) / (2 * eps), fdRho = (rp[r] - rm[r]) / (2 * eps);
            EXPECT_NEAR(dState[r].dlnT, fdT, 1e-5 * std::abs(fdT) + 1e-9) << s.label << " r" << r;
            EXPECT_NEAR(dState[r].dlnRho, fdRho, 1e-5 * std::abs(fdRho) + 1e-9) << s.label << " r" << r;
        }
        const std::vector<CppAD::AD<double>> Yad(c.Y.begin(), c.Y.end());
        const auto f_ad = model->calculateScreeningFactors(set, c.species, Yad, CppAD::AD<double>(T9), CppAD::AD<double>(s.rho));
        const auto f = model->calculateScreeningFactors(set, c.species, c.Y, T9, s.rho);
        for (size_t r = 0; r < lists.size(); ++r) EXPECT_NEAR(CppAD::Value(f_ad[r]), f[r], 1e-13 * f[r]) << s.label << " r" << r;
    }
}

TEST(BlendedScreening, SelectableAndValidated) {
    gridfire::test::quiet_logging();
    EXPECT_THROW(blended({.weakExponent = 0.5, .strongExponent = 2.5}), std::invalid_argument) << "window above the weak cap";
    EXPECT_THROW(blended({.weakExponent = 1.0, .strongExponent = 0.5}), std::invalid_argument) << "inverted window";
    EXPECT_THROW(blended({.weakExponent = 0.0, .strongExponent = 1.0}), std::invalid_argument) << "zero lower edge";

    gridfire::config::GridFireConfig config;
    config.engine.screening = gridfire::screening::ScreeningType::BLENDED;
    config.engine.screeningBlend = {.weakExponent = 0.05, .strongExponent = 0.3}; // so the solar network sits in the window
    gridfire::policy::MainSequencePolicy policy(gridfire::test::solar_composition(), config);
    const auto built = policy.construct();
    const auto& graph = dynamic_cast<const gridfire::engine::GraphEngine&>(built.engine);
    EXPECT_EQ(graph.getScreeningModel(*built.scratch_blob), gridfire::screening::ScreeningType::BLENDED);
    const auto& species = graph.getNetworkSpecies(*built.scratch_blob);
    const auto comp = gridfire::test::make_composition(species, gridfire::test::abundances_in_order(gridfire::test::solar_composition(), species));
    const auto analytic = graph.generateAnalyticJacobianMatrix(*built.scratch_blob, comp, 0.0157, 150.0);
    ASSERT_TRUE(analytic.has_value());
    const auto autodiff = graph.generateAutoDiffJacobianMatrix(*built.scratch_blob, comp, 0.0157, 150.0);
    double worst = 0, largest = 0;
    for (const auto& row : species) {
        for (const auto& col : species) {
            worst = std::max(worst, std::abs((*analytic)(row, col) - autodiff(row, col)));
            largest = std::max(largest, std::abs(autodiff(row, col)));
        }
    }
    EXPECT_LT(worst, 1e-10 * largest) << "analytic and CppAD Jacobians disagree under BLENDED screening";
    const auto epsA = *graph.calculateAnalyticEpsDerivatives(comp, 0.0157, 150.0);
    const auto epsAD = graph.calculateAutoDiffEpsDerivatives(*built.scratch_blob, comp, 0.0157, 150.0);
    EXPECT_NEAR(epsA.dEps_dT, epsAD.dEps_dT, 1e-9 * std::abs(epsAD.dEps_dT));
    EXPECT_NEAR(epsA.dEps_dRho, epsAD.dEps_dRho, 1e-9 * std::abs(epsAD.dEps_dRho));
}
