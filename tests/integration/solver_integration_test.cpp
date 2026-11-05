// Integration tests: end-to-end runs of PointSolver and GridSolver.
//
// Three kinds of checks:
//   * Exact answers. Tiny networks whose reactions have constant rates (REACLIB a0 only) have closed-form solutions
//     (exponential decay, Bateman chains, relaxation to equilibrium), so they measure the integrator directly.
//   * Bookkeeping that needs no reference data: NetOut.energy equals the mass defect of the composition change minus
//     the neutrino losses; baryon number is conserved; chained runs reproduce a single run; NetOut's instantaneous
//     fields describe its own final state.
//   * Loose physics sanity checks and GridSolver consistency on the main sequence network.
#include "gf_test_support.h"

#include "gridfire/engine/engine_graph.h"
#include "gridfire/reaction/reaclib.h"
#include "gridfire/reaction/reaction.h"
#include "gridfire/screening/screening_types.h"
#include "gridfire/solver/strategies/GridSolver.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/config/config.h"
#include "gridfire/utils/gf_omp.h"

#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils/utils.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <cmath>
#include <format>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using fourdst::atomic::Species;
using gridfire::NetIn;
using gridfire::NetOut;
using gridfire::engine::GraphEngine;
using gridfire::reaction::ReaclibReaction;
using gridfire::reaction::ReactionSet;
using gridfire::solver::PointSolver;
using gridfire::solver::PointSolverContext;

namespace {
    const Species& sp(const char* name) { return fourdst::atomic::species.at(name); }

    /// Most failures SolverRobustness.StressGridFailureBudget accepts (measured; lower it as robustness improves)
    constexpr size_t STRESS_FAILURE_BUDGET = 0;

    /// A reaction with a temperature- and density-independent rate (all REACLIB coefficients but a0 are zero).
    /// Q = 0 so the engine never attaches a reverse rate, and the type is plain REACLIB so there are no neutrino losses.
    std::unique_ptr<gridfire::reaction::Reaction> constant_rate(
        const std::string& id,
        const std::vector<Species>& reactants,
        const std::vector<Species>& products,
        const double rate
    ) {
        return std::make_unique<ReaclibReaction>(id, id, 1, reactants, products, 0.0, "test",
                                                 gridfire::reaction::RateCoefficientSet{std::log(rate), 0, 0, 0, 0, 0, 0});
    }

    /// A small GraphEngine together with a state blob for it.
    struct SmallNetwork {
        std::unique_ptr<GraphEngine> engine;
        std::unique_ptr<gridfire::engine::scratch::StateBlob> context;

        explicit SmallNetwork(const ReactionSet& reactions) {
            gridfire::test::quiet_logging();
            GF_PAR_INIT();
            engine = std::make_unique<GraphEngine>(reactions);
            context = engine->constructStateBlob(nullptr);
        }
    };

    ReactionSet make_set(std::vector<std::unique_ptr<gridfire::reaction::Reaction>> reactions) {
        return ReactionSet(std::move(reactions));
    }

    struct RunOptions {
        double rel_tol = 1e-10;
        double abs_tol = 1e-12;
    };

    NetOut run(const gridfire::engine::DynamicEngine& engine, const gridfire::engine::scratch::StateBlob& context,
               const fourdst::composition::Composition& composition, const double T, const double rho,
               const double tMax, const RunOptions& options = {}) {
        PointSolverContext sctx(context);
        sctx.set_stdout_logging(false);
        sctx.rel_tol = options.rel_tol;
        sctx.abs_tol = options.abs_tol;
        const PointSolver solver(engine);
        NetIn in;
        in.composition = composition;
        in.temperature = T;
        in.density = rho;
        in.tMax = tMax;
        in.dt0 = 1e-12;
        in.energy = 0.0;
        return solver.evaluate(sctx, in);
    }

    double Y(const NetOut& out, const char* species) {
        return out.composition.contains(sp(species)) ? out.composition.getMolarAbundance(sp(species)) : 0.0;
    }

    /// GridFire's definition of the energy released by a change of composition: -sum_i dY_i m_i u N_A c^2 [erg/g].
    double mass_defect_energy(const fourdst::composition::CompositionAbstract& before,
                              const fourdst::composition::CompositionAbstract& after) {
        const auto& k = gridfire::test::constants();
        double dm = 0.0;
        for (const auto& species : after.getRegisteredSpecies()) {
            const double y0 = before.contains(species) ? before.getMolarAbundance(species) : 0.0;
            dm += (after.getMolarAbundance(species) - y0) * species.mass();
        }
        return -dm * k.u * k.Na * k.c * k.c;
    }

    double baryon_number(const fourdst::composition::CompositionAbstract& comp) {
        double total = 0.0;
        for (const auto& species : comp.getRegisteredSpecies()) total += species.a() * comp.getMolarAbundance(species);
        return total;
    }

    constexpr double MeV_to_erg = 1.602176634e-6;
}

// ---------------------------------------------------------------------------------------------------------------------
// Exact answers
// ---------------------------------------------------------------------------------------------------------------------

TEST(ExactSolutions, SingleBetaDecayFollowsExponential) {
    // The real REACLIB 13N(,e+)13C rate is a constant (half-life 598 s).
    auto reaction = gridfire::reaclib::get_all_reaclib_reactions().get("n13(,e+)c13");
    ASSERT_TRUE(reaction.has_value());
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(std::move(*reaction));
    const SmallNetwork net(make_set(std::move(list)));

    const auto initial = gridfire::test::make_composition({sp("N-13"), sp("C-13")}, {1.0 / 13.0, 0.0});
    const auto rhs = net.engine->calculateRHSAndEnergy(*net.context, initial, 0.015, 100.0, false);
    ASSERT_TRUE(rhs.has_value());
    const double lambda = -rhs->dydt.at(sp("N-13")) / initial.getMolarAbundance(sp("N-13"));
    EXPECT_NEAR(std::log(2.0) / lambda, 598.0, 2.0) << "13N half-life";

    for (const double t : {60.0, 600.0, 3000.0}) {
        const NetOut out = run(*net.engine, *net.context, initial, 1.5e7, 100.0, t);
        const double exact = (1.0 / 13.0) * std::exp(-lambda * t);
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "N-13"), exact), 1e-6) << "t = " << t;
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "C-13"), 1.0 / 13.0 - exact), 1e-6) << "t = " << t;
    }
}

TEST(ExactSolutions, TwoStepChainFollowsBateman) {
    const double l1 = 1e-2, l2 = 3e-3;
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("o20(,e-)f20", {sp("O-20")}, {sp("F-20")}, l1));
    list.push_back(constant_rate("f20(,e-)ne20", {sp("F-20")}, {sp("Ne-20")}, l2));
    const SmallNetwork net(make_set(std::move(list)));

    const double y0 = 1.0 / 20.0;
    const auto initial = gridfire::test::make_composition({sp("O-20"), sp("F-20"), sp("Ne-20")}, {y0, 0.0, 0.0});
    for (const double t : {10.0, 100.0, 500.0, 2000.0}) {
        const NetOut out = run(*net.engine, *net.context, initial, 1e8, 1.0, t);
        const double a = y0 * std::exp(-l1 * t);
        const double b = y0 * l1 / (l2 - l1) * (std::exp(-l1 * t) - std::exp(-l2 * t));
        const double c = y0 - a - b;
        // rtol = 1e-10 per step; allow the global error to grow to 1e-6 of the abundance, or of 1e-8 y0 once decayed away
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "O-20"), a, 1e-8 * y0), 1e-6) << "t = " << t;
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "F-20"), b, 1e-8 * y0), 1e-6) << "t = " << t;
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "Ne-20"), c, 1e-8 * y0), 1e-6) << "t = " << t;
    }
}

TEST(ExactSolutions, ReversiblePairRelaxesToEquilibrium) {
    const double kf = 2e-2, kr = 5e-3;
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("n13(,e+)c13", {sp("N-13")}, {sp("C-13")}, kf));
    list.push_back(constant_rate("c13(,e-)n13", {sp("C-13")}, {sp("N-13")}, kr));
    const SmallNetwork net(make_set(std::move(list)));

    const double total = 1.0 / 13.0;
    const double a_eq = total * kr / (kf + kr);
    const auto initial = gridfire::test::make_composition({sp("N-13"), sp("C-13")}, {total, 0.0});
    for (const double t : {5.0, 50.0, 500.0, 5000.0}) {
        const NetOut out = run(*net.engine, *net.context, initial, 1e8, 1.0, t);
        const double a = a_eq + (total - a_eq) * std::exp(-(kf + kr) * t);
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "N-13"), a), 1e-6) << "t = " << t;
        EXPECT_LT(gridfire::test::rel_diff(Y(out, "C-13"), total - a), 1e-6) << "t = " << t;
    }
}

TEST(ExactSolutions, ErrorShrinksAsToleranceTightens) {
    const double l1 = 1e-2, l2 = 3e-3, t = 500.0;
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("o20(,e-)f20", {sp("O-20")}, {sp("F-20")}, l1));
    list.push_back(constant_rate("f20(,e-)ne20", {sp("F-20")}, {sp("Ne-20")}, l2));
    const SmallNetwork net(make_set(std::move(list)));

    const double y0 = 1.0 / 20.0;
    const auto initial = gridfire::test::make_composition({sp("O-20"), sp("F-20"), sp("Ne-20")}, {y0, 0.0, 0.0});
    const double exact = y0 * l1 / (l2 - l1) * (std::exp(-l1 * t) - std::exp(-l2 * t));

    std::vector<double> errors;
    for (const double rtol : {1e-4, 1e-6, 1e-8, 1e-10}) {
        const NetOut out = run(*net.engine, *net.context, initial, 1e8, 1.0, t, {rtol, 1e-12});
        errors.push_back(gridfire::test::rel_diff(Y(out, "F-20"), exact));
        // CVODE controls local error; the global error may exceed rtol by a modest factor, never by orders of magnitude
        EXPECT_LT(errors.back(), 100.0 * rtol) << "rtol = " << rtol;
    }
    EXPECT_LT(errors.back(), errors.front() * 1e-3) << "tightening rtol by 1e6 should reduce the error by well over 1e3";
}

TEST(ExactSolutions, PublicForwardFlowsReproduceScreenedRHS) {
    // Strong two-body reactions with Q = 0 have no reverse rate, so dY/dt = sum_r nu_r * (forward flow)_r exactly. At
    // solar-core density weak screening changes these flows by several percent, so this checks that the public
    // per-reaction flow (used for culling and QSE timescales) includes the same screening as the RHS.
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("c12_p", {sp("C-12"), sp("H-1")}, {sp("N-13")}, 1e-10));
    list.push_back(constant_rate("c12_a", {sp("C-12"), sp("He-4")}, {sp("O-16")}, 1e-12));
    const SmallNetwork net(make_set(std::move(list)));
    ASSERT_EQ(net.engine->getScreeningModel(*net.context), gridfire::screening::ScreeningType::WEAK);

    const auto comp = gridfire::test::make_composition({sp("H-1"), sp("He-4"), sp("C-12"), sp("N-13"), sp("O-16")},
                                                       {0.7, 0.07, 1e-3, 0.0, 0.0});
    const double T9 = 0.015, rho = 150.0;
    const auto rhs = net.engine->calculateRHSAndEnergy(*net.context, comp, T9, rho, false);
    ASSERT_TRUE(rhs.has_value());
    std::map<Species, double> from_flows;
    for (const auto& reaction : net.engine->getNetworkReactions(*net.context)) {
        const double flow = net.engine->calculateMolarReactionFlow(*net.context, *reaction, comp, T9, rho);
        for (const auto& s : net.engine->getNetworkSpecies(*net.context)) from_flows[s] += reaction->stoichiometry(s) * flow;
    }
    for (const auto& [s, dydt] : rhs->dydt) {
        EXPECT_LT(gridfire::test::rel_diff(from_flows[s], dydt, 1e-300), 1e-12) << s.name();
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Energy bookkeeping
// ---------------------------------------------------------------------------------------------------------------------

TEST(EnergyBookkeeping, BetaDecayEnergyIsMassDefectMinusNeutrinos) {
    auto reaction = gridfire::reaclib::get_all_reaclib_reactions().get("n13(,e+)c13");
    ASSERT_TRUE(reaction.has_value());
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(std::move(*reaction));
    const SmallNetwork net(make_set(std::move(list)));

    const auto initial = gridfire::test::make_composition({sp("N-13"), sp("C-13")}, {1.0 / 13.0, 0.0});
    const NetOut out = run(*net.engine, *net.context, initial, 1.5e7, 100.0, 1000.0);

    const double decayed = 1.0 / 13.0 - Y(out, "N-13");
    const double neutrinos = decayed * 0.706 * gridfire::test::constants().Na * MeV_to_erg; // Bahcall <E_nu> for 13N
    const double expected = mass_defect_energy(initial, out.composition) - neutrinos;
    EXPECT_LT(gridfire::test::rel_diff(out.energy, expected), 1e-6)
        << std::format("energy {:.10e} vs mass defect {:.10e} - neutrinos {:.10e}", out.energy,
                       mass_defect_energy(initial, out.composition), neutrinos);

    // The reported neutrino loss is accumulated as (rate at last RHS call) x (step), so it is only first-order accurate
    EXPECT_LT(gridfire::test::rel_diff(out.specific_neutrino_energy_loss, neutrinos), 0.05)
        << std::format("accumulated neutrino loss {:.6e} vs exact {:.6e}", out.specific_neutrino_energy_loss, neutrinos);
}

TEST(EnergyBookkeeping, StrongOnlyEnergyIsMassDefect) {
    // A constant-rate strong chain has no neutrino losses, so the released energy is exactly the mass defect
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("c12_alpha", {sp("C-12"), sp("He-4")}, {sp("O-16")}, 1e-3));
    const SmallNetwork net(make_set(std::move(list)));

    const auto initial = gridfire::test::make_composition({sp("He-4"), sp("C-12"), sp("O-16")}, {0.5 / 4.0, 0.5 / 12.0, 0.0});
    const NetOut out = run(*net.engine, *net.context, initial, 2e8, 1e3, 1.0);
    ASSERT_GT(1.0 - Y(out, "C-12") / (0.5 / 12.0), 1e-3) << "test should burn a measurable amount";
    EXPECT_LT(gridfire::test::rel_diff(out.energy, mass_defect_energy(initial, out.composition)), 1e-7);
    EXPECT_EQ(out.specific_neutrino_energy_loss, 0.0);
}

TEST(EnergyBookkeeping, SolarRunEnergyIsMassDefectMinusNeutrinos) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const auto initial = gridfire::test::solar_composition();
    const NetOut out = run(shared.engine, *context, initial, 1.5e7, 150.0, 3e15, {1e-8, 1e-14});

    const double defect = mass_defect_energy(initial, out.composition);
    ASSERT_GT(defect, 0.0);
    // Neutrinos carry ~2-3% of the energy; their accumulated total is first-order accurate
    const double nu_fraction = out.specific_neutrino_energy_loss / defect;
    EXPECT_GT(nu_fraction, 0.015);
    EXPECT_LT(nu_fraction, 0.05);
    EXPECT_LT(gridfire::test::rel_diff(out.energy + out.specific_neutrino_energy_loss, defect), 2e-3)
        << std::format("energy {:.8e} + neutrinos {:.8e} vs mass defect {:.8e}", out.energy,
                       out.specific_neutrino_energy_loss, defect);
    RecordProperty("neutrino_fraction", std::format("{:.4f}", nu_fraction));
}

// ---------------------------------------------------------------------------------------------------------------------
// Conservation, restarts and final-state consistency on the main sequence network
// ---------------------------------------------------------------------------------------------------------------------

TEST(SolverConsistency, BaryonNumberConservedAndStepsCounted) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const auto initial = gridfire::test::solar_composition();
    // sum_i A_i Y_i is not 1 (Y is per gram, and masses are not A amu), but reactions must leave it unchanged.
    const double B0 = baryon_number(initial);
    // The 4e7 K run burns hydrogen to exhaustion. Default tolerances: tighter ones currently fail there (see KnownBugs).
    for (const auto& [T, tMax] : std::vector<std::pair<double, double>>{{1.5e7, 1e17}, {2.5e7, 1e14}, {4e7, 1e13}}) {
        const NetOut out = run(shared.engine, *context, initial, T, 100.0, tMax, {1e-5, 1e-8});
        ASSERT_LT(Y(out, "H-1"), 0.9 * initial.getMolarAbundance(sp("H-1"))) << "T = " << T << ": should burn substantially";
        EXPECT_LT(std::abs(baryon_number(out.composition) / B0 - 1.0), 1e-10) << "T = " << T;
        for (const auto& species : out.composition.getRegisteredSpecies()) {
            EXPECT_GE(out.composition.getMolarAbundance(species), 0.0) << species.name() << " at T = " << T;
        }
        EXPECT_GT(out.num_steps, 10) << "T = " << T << ": num_steps should count every CVODE step, across reinitialisations";
    }
}

TEST(SolverConsistency, ChainedRunsReproduceSingleRun) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const auto initial = gridfire::test::solar_composition();
    const double T = 2e7, rho = 100.0, tMax = 1e14;
    const RunOptions options{1e-9, 1e-15};

    const NetOut whole = run(shared.engine, *context, initial, T, rho, tMax, options);
    const NetOut first = run(shared.engine, *context, initial, T, rho, tMax / 2, options);
    const NetOut second = run(shared.engine, *context, first.composition, T, rho, tMax / 2, options);

    for (const auto& species : whole.composition.getRegisteredSpecies()) {
        const double X_whole = whole.composition.getMassFraction(species);
        if (X_whole < 1e-10) continue;
        EXPECT_LT(gridfire::test::rel_diff(second.composition.getMassFraction(species), X_whole), 1e-5) << species.name();
    }
    EXPECT_LT(gridfire::test::rel_diff(first.energy + second.energy, whole.energy), 1e-5);
}

TEST(SolverConsistency, InstantaneousOutputsDescribeFinalState) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const double T = 1.8e7, rho = 120.0;
    const NetOut out = run(shared.engine, *context, gridfire::test::solar_composition(), T, rho, 1e13, {1e-8, 1e-14});

    const auto probe = shared.fresh_context();
    const auto rhs = shared.engine.calculateRHSAndEnergy(*probe, out.composition, T / 1e9, rho, false);
    ASSERT_TRUE(rhs.has_value());
    const auto [dEps_dT, dEps_dRho] = shared.engine.calculateEpsDerivatives(*probe, out.composition, T / 1e9, rho);
    const auto channels = shared.engine.calculateEnergyChannels(*probe, out.composition, T / 1e9, rho);

    EXPECT_LT(gridfire::test::rel_diff(out.eps_nuc, rhs->nuclearEnergyGenerationRate), 1e-12);
    EXPECT_LT(gridfire::test::rel_diff(out.neutrino_energy_loss_rate, rhs->neutrinoEnergyLossRate), 1e-12);
    EXPECT_LT(gridfire::test::rel_diff(out.dEps_dT, dEps_dT), 1e-12);
    EXPECT_LT(gridfire::test::rel_diff(out.dEps_dRho, dEps_dRho), 1e-12);
    EXPECT_LT(gridfire::test::rel_diff(out.dlnEps_dlnT, T * dEps_dT / rhs->nuclearEnergyGenerationRate), 1e-12);
    EXPECT_LT(gridfire::test::rel_diff(out.dlnEps_dlnRho, rho * dEps_dRho / rhs->nuclearEnergyGenerationRate), 1e-12);

    ASSERT_EQ(out.eps_by_channel.size(), channels.size());
    ASSERT_EQ(out.energy_channel_names.size(), channels.size());
    for (size_t i = 0; i < channels.size(); ++i) {
        EXPECT_LT(gridfire::test::rel_diff(out.eps_by_channel[i], channels[i], 1e-30), 1e-12) << out.energy_channel_names[i];
    }
    const double sum = std::accumulate(out.eps_by_channel.begin(), out.eps_by_channel.end(), 0.0);
    EXPECT_LT(gridfire::test::rel_diff(sum, out.eps_nuc), 1e-10);

    // Hydrogen burning: eps ~ T^4 (pp) to T^20 (CNO), and eps ~ rho for two-body reactions
    EXPECT_GT(out.dlnEps_dlnT, 3.0);
    EXPECT_LT(out.dlnEps_dlnT, 25.0);
    EXPECT_NEAR(out.dlnEps_dlnRho, 1.0, 0.05);
}

// ---------------------------------------------------------------------------------------------------------------------
// Physics sanity (loose bounds)
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    double channel(const NetOut& out, const std::string& name) {
        for (size_t i = 0; i < out.energy_channel_names.size(); ++i) {
            if (out.energy_channel_names[i] == name) return out.eps_by_channel[i];
        }
        ADD_FAILURE() << "no energy channel named " << name;
        return 0.0;
    }

    /// Solar composition burned at T long enough for the CN cycle to approach equilibrium, as in a main sequence core.
    NetOut burned_solar(const double T) {
        const auto& shared = gridfire::test::main_sequence_engine();
        const auto context = shared.fresh_context();
        // ~100 Myr: long enough to convert the initial 12C to 14N (the CN cycle then runs at its slow equilibrium rate),
        // but ~1% of the hydrogen burning time
        return run(shared.engine, *context, gridfire::test::solar_composition(), T, 150.0, 3e15, {1e-8, 1e-14});
    }
}

TEST(PhysicsSanity, SolarCoreIsPPDominated) {
    const NetOut out = burned_solar(1.5e7);
    const double pp = channel(out, "H_pp"), cno = channel(out, "H_CNO");
    RecordProperty("pp_fraction", std::format("{:.4f}", pp / out.eps_nuc));
    EXPECT_GT(pp / out.eps_nuc, 0.95) << std::format("pp {:.4e}, CNO {:.4e}, total {:.4e}", pp, cno, out.eps_nuc);
    EXPECT_GT(out.eps_nuc, 1.0);
    EXPECT_LT(out.eps_nuc, 1e2);
}

TEST(PhysicsSanity, PPAndCNOCrossNearSeventeenToTwentyMillionKelvin) {
    const NetOut cool = burned_solar(1.6e7);
    const NetOut hot = burned_solar(2.3e7);
    EXPECT_GT(channel(cool, "H_pp"), channel(cool, "H_CNO"));
    EXPECT_LT(channel(hot, "H_pp"), channel(hot, "H_CNO"));
}

TEST(PhysicsSanity, TripleAlphaTurnsOnNearOneHundredMillionKelvin) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto species = shared.engine.getNetworkSpecies(*shared.fresh_context());
    std::vector<double> Y(species.size(), 0.0);
    for (size_t i = 0; i < species.size(); ++i) {
        if (species[i] == sp("He-4")) Y[i] = 0.99 / 4.0;
        if (species[i] == sp("C-12")) Y[i] = 0.01 / 12.0;
    }
    const auto helium = gridfire::test::make_composition(species, Y);

    auto he_eps = [&](const double T) {
        const auto context = shared.fresh_context();
        const NetOut out = run(shared.engine, *context, helium, T, 1e4, 1.0, {1e-8, 1e-14});
        return channel(out, "He");
    };
    const double cold = he_eps(5e7), warm = he_eps(1e8), hot = he_eps(2e8);
    RecordProperty("eps_He_1e8", std::format("{:.4e}", warm));
    EXPECT_LT(cold, 1e-3 * warm);
    EXPECT_GT(warm, 1.0);
    EXPECT_GT(hot, 1e3 * warm);
}

// ---------------------------------------------------------------------------------------------------------------------
// GridSolver
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    std::vector<NetOut> run_grid(const std::vector<NetIn>& inputs, const int threads) {
        const auto& shared = gridfire::test::main_sequence_engine();
        const auto context = shared.fresh_context();
        const PointSolver point(shared.engine);
        const gridfire::solver::GridSolver grid(shared.engine, point);
        gridfire::solver::GridSolverContext gctx(*context);
        gctx.zone_completion_logging = false;
        gctx.set_stdout_logging(false);
#ifdef _OPENMP
        const int previous = omp_get_max_threads();
        omp_set_num_threads(threads);
#else
        (void)threads;
#endif
        auto results = grid.evaluate(gctx, inputs);
#ifdef _OPENMP
        omp_set_num_threads(previous);
#endif
        return results;
    }

    std::vector<NetIn> zone_inputs(const std::vector<double>& temperatures) {
        std::vector<NetIn> inputs;
        for (const double T : temperatures) {
            NetIn in;
            in.composition = gridfire::test::solar_composition();
            in.temperature = T;
            in.density = 100.0;
            in.tMax = 1e12;
            in.dt0 = 1e-12;
            in.energy = 0.0;
            inputs.push_back(in);
        }
        return inputs;
    }

    void expect_identical(const NetOut& a, const NetOut& b, const std::string& what) {
        EXPECT_EQ(a.energy, b.energy) << what;
        EXPECT_EQ(a.eps_nuc, b.eps_nuc) << what;
        EXPECT_EQ(a.dEps_dT, b.dEps_dT) << what;
        EXPECT_EQ(a.eps_by_channel, b.eps_by_channel) << what;
        for (const auto& species : a.composition.getRegisteredSpecies()) {
            EXPECT_EQ(a.composition.getMolarAbundance(species), b.composition.getMolarAbundance(species)) << what << " " << species.name();
        }
    }
}

TEST(GridSolver, IdenticalZonesGiveIdenticalResults) {
    const auto results = run_grid(zone_inputs({1.6e7, 1.6e7, 1.6e7, 1.6e7}), 4);
    ASSERT_EQ(results.size(), 4u);
    for (size_t z = 1; z < results.size(); ++z) expect_identical(results[0], results[z], std::format("zone {}", z));
}

TEST(GridSolver, ZonesMatchStandalonePointSolver) {
    const auto inputs = zone_inputs({1.4e7, 1.8e7, 2.2e7, 3.0e7});
    const auto results = run_grid(inputs, 4);
    const auto& shared = gridfire::test::main_sequence_engine();
    for (size_t z = 0; z < inputs.size(); ++z) {
        const auto context = shared.fresh_context();
        PointSolverContext sctx(*context);
        sctx.set_stdout_logging(false);
        const NetOut alone = PointSolver(shared.engine).evaluate(sctx, inputs[z]);
        expect_identical(results[z], alone, std::format("zone {} (T = {:.2e})", z, inputs[z].temperature));
    }
}

TEST(GridSolver, ThreadCountDoesNotChangeResults) {
    const auto inputs = zone_inputs({1.4e7, 1.8e7, 2.2e7, 3.0e7, 1.5e7, 2.5e7});
    const auto serial = run_grid(inputs, 1);
    const auto parallel = run_grid(inputs, 4);
    ASSERT_EQ(serial.size(), parallel.size());
    for (size_t z = 0; z < serial.size(); ++z) expect_identical(serial[z], parallel[z], std::format("zone {}", z));
}

// ---------------------------------------------------------------------------------------------------------------------
// Robustness
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    struct StressOutcome {
        int runs = 0;
        std::vector<std::string> failures;
    };

    /// Default-tolerance burns over a grid of hydrogen-burning conditions, like those a stellar evolution code asks
    /// for: T = 1e7-5e7 K, rho = 1-1e4 g/cm^3, solar H or 97% of it moved to He-4, and burns to tMax = 1e12-1e17 s
    /// (the longer ones exhaust hydrogen). Modelled on the stress grid in SOLVER_FINDINGS_v1.2.0.
    StressOutcome run_stress_grid(const gridfire::config::GridFireConfig& config) {
        const auto& shared = gridfire::test::main_sequence_engine();
        const auto solar = gridfire::test::solar_composition();
        StressOutcome outcome;
        for (const double T : {1e7, 1.5e7, 3e7, 5e7}) {
            for (const double rho : {1.0, 150.0, 1e4}) {
                for (const double hydrogenScale : {1.0, 0.03}) {
                    std::vector<std::string> symbols;
                    std::vector<double> X;
                    for (const auto& s : solar.getRegisteredSpecies()) {
                        symbols.emplace_back(s.name());
                        X.push_back(solar.getMassFraction(s));
                    }
                    const auto h = std::ranges::find(symbols, "H-1") - symbols.begin();
                    const auto he = std::ranges::find(symbols, "He-4") - symbols.begin();
                    X[he] += X[h] * (1.0 - hydrogenScale);
                    X[h] *= hydrogenScale;
                    const auto composition = fourdst::composition::buildCompositionFromMassFractions(symbols, X);

                    for (const double tMax : {1e12, 1.6e15, 1e17}) {
                        ++outcome.runs;
                        const auto context = shared.fresh_context();
                        PointSolverContext sctx(*context);
                        sctx.set_stdout_logging(false);
                        NetIn in;
                        in.composition = composition;
                        in.temperature = T;
                        in.density = rho;
                        in.tMax = tMax;
                        in.dt0 = 1e-12;
                        in.energy = 0.0;
                        try {
                            (void)PointSolver(shared.engine, config).evaluate(sctx, in);
                        } catch (const std::exception& e) {
                            outcome.failures.push_back(std::format("T={:.1e} rho={:.0e} X_H x{} tMax={:.1e}: {}", T, rho, hydrogenScale, tMax, e.what()));
                        }
                    }
                }
            }
        }
        return outcome;
    }
}

TEST(SolverRobustness, StressGridFailureBudget) {
    const StressOutcome outcome = run_stress_grid(gridfire::config::GridFireConfig{});
    std::string list;
    for (const auto& failure : outcome.failures) list += "\n  " + failure;
    RecordProperty("failures", static_cast<int>(outcome.failures.size()));
    std::cout << std::format("stress grid: {} of {} runs failed{}\n", outcome.failures.size(), outcome.runs, list);
    EXPECT_LE(outcome.failures.size(), STRESS_FAILURE_BUDGET) << "failures:" << list;
}

// ---------------------------------------------------------------------------------------------------------------------
// Failure handling
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    NetIn valid_input() {
        NetIn in;
        in.composition = gridfire::test::solar_composition();
        in.temperature = 1.5e7;
        in.density = 100.0;
        in.tMax = 1e10;
        in.dt0 = 1e-12;
        in.energy = 0.0;
        return in;
    }

    void evaluate_point(const NetIn& in) {
        const auto& shared = gridfire::test::main_sequence_engine();
        const auto context = shared.fresh_context();
        PointSolverContext sctx(*context);
        sctx.set_stdout_logging(false);
        (void)PointSolver(shared.engine).evaluate(sctx, in);
    }
}

TEST(FailureHandling, ValidInputRuns) {
    EXPECT_NO_THROW(evaluate_point(valid_input()));
}

TEST(FailureHandling, NonPositiveTemperatureIsRejected) {
    for (const double T : {0.0, -1e7}) {
        NetIn in = valid_input();
        in.temperature = T;
        EXPECT_ANY_THROW(evaluate_point(in)) << "T = " << T;
    }
}

TEST(FailureHandling, NonPositiveDensityIsRejected) {
    for (const double rho : {0.0, -1.0}) {
        NetIn in = valid_input();
        in.density = rho;
        EXPECT_ANY_THROW(evaluate_point(in)) << "rho = " << rho;
    }
}

TEST(FailureHandling, NonFiniteInputsAreRejected) {
    NetIn in = valid_input();
    in.temperature = std::nan("");
    EXPECT_ANY_THROW(evaluate_point(in)) << "T = NaN";
    in = valid_input();
    in.density = std::numeric_limits<double>::infinity();
    EXPECT_ANY_THROW(evaluate_point(in)) << "rho = inf";
}

TEST(FailureHandling, SpeciesOutsideNetworkIsRejected) {
    NetIn in = valid_input();
    in.composition = fourdst::composition::buildCompositionFromMassFractions(
        std::vector<std::string>{"H-1", "He-4", "U-238"}, std::vector<double>{0.7, 0.29, 0.01});
    EXPECT_ANY_THROW(evaluate_point(in));
}


TEST(FailureHandling, NonPositiveTMaxDoesNotIntegrate) {
    for (const double tMax : {0.0, -1.0}) {
        NetIn in = valid_input();
        in.tMax = tMax;
        // Either rejecting the input or returning the unchanged composition with zero energy is acceptable
        try {
            const auto& shared = gridfire::test::main_sequence_engine();
            const auto context = shared.fresh_context();
            PointSolverContext sctx(*context);
            sctx.set_stdout_logging(false);
            const NetOut out = PointSolver(shared.engine).evaluate(sctx, in);
            EXPECT_EQ(out.energy, 0.0) << "tMax = " << tMax;
            EXPECT_DOUBLE_EQ(out.composition.getMolarAbundance(sp("H-1")), in.composition.getMolarAbundance(sp("H-1")));
        } catch (const std::exception&) {
            SUCCEED();
        }
    }
}

// Known bug: the energy accumulator's error weight is rel_tol * |E| + abs_tol, with abs_tol (a mass-fraction tolerance)
// used directly in erg/g. At t = 0, E = 0, so a tight abs_tol demands e.g. 1e-16 erg/g of an energy growing at
// ~1e15 erg/g/s, and CVODE fails on its first step. Making the accumulator a CVODE quadrature with its own (relative)
// tolerance would fix it. This test skips while the bug is present and fails once it is fixed, as a reminder to
// turn it into a normal assertion.
TEST(KnownBugs, TightAbsoluteToleranceFailsOnEnergyAccumulator) {
    std::vector<std::unique_ptr<gridfire::reaction::Reaction>> list;
    list.push_back(constant_rate("o20(,e-)f20", {sp("O-20")}, {sp("F-20")}, 1e-2));
    const SmallNetwork net(make_set(std::move(list)));
    const auto initial = gridfire::test::make_composition({sp("O-20"), sp("F-20")}, {0.05, 0.0});
    bool failed = false;
    try {
        (void)run(*net.engine, *net.context, initial, 1e8, 1.0, 100.0, {1e-8, 1e-16});
    } catch (const std::exception&) {
        failed = true;
    }
    if (failed) {
        GTEST_SKIP() << "Known bug: abs_tol is applied to the energy accumulator in erg/g; abs_tol = 1e-16 fails at t = 0";
    }
    FAIL() << "abs_tol = 1e-16 now integrates: the energy-tolerance bug appears fixed; make this a normal test";
}


// Tight tolerances used to fail with CV_CONSTR_FAIL as hydrogen ran out (H-1 ~ 1e-8): clamping negative iterates
// flattened the RHS of destroyed species while the exact Jacobian still reported their slope. Evaluating raw iterates,
// plus restarts after a failed step, fixed it; the result must also agree with the default-tolerance run.
TEST(SolverRobustness, TightTolerancesIntegrateThroughHydrogenExhaustion) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const auto initial = gridfire::test::solar_composition();
    const NetOut tight = run(shared.engine, *context, initial, 4e7, 100.0, 1e13, {1e-8, 1e-12});
    const NetOut loose = run(shared.engine, *context, initial, 4e7, 100.0, 1e13, {1e-5, 1e-8});
    EXPECT_LT(tight.composition.getMassFraction(sp("H-1")), 1e-6) << "the burn should exhaust hydrogen";
    EXPECT_LT(std::abs(baryon_number(tight.composition) / baryon_number(initial) - 1.0), 1e-10);
    EXPECT_LT(gridfire::test::rel_diff(tight.energy, loose.energy), 1e-4);
    EXPECT_LT(gridfire::test::rel_diff(tight.composition.getMassFraction(sp("He-4")), loose.composition.getMassFraction(sp("He-4"))), 1e-5);
}


// Known bug: GridSolver catches a zone's GridFireError, prints it, and returns a default-constructed NetOut for that
// zone (with uninitialised energy and num_steps), so the caller cannot tell the zone failed. It should rethrow after
// the parallel loop, or report per-zone status. Skips while the bug is present; fails once fixed.
TEST(KnownBugs, GridSolverHidesZoneFailures) {
    auto inputs = zone_inputs({1.6e7, 1.6e7});
    inputs[1].temperature = -1.0; // Rejected by PointSolver::evaluate
    bool threw = false;
    try {
        (void)run_grid(inputs, 2);
    } catch (const std::exception&) {
        threw = true;
    }
    if (!threw) {
        GTEST_SKIP() << "Known bug: GridSolver swallows a failing zone's exception and returns a default NetOut";
    }
    FAIL() << "GridSolver now reports zone failures: the bug appears fixed; make this a normal test";
}

// ---------------------------------------------------------------------------------------------------------------------
// evaluateRates: instantaneous rates through the solver stack, without integrating
// ---------------------------------------------------------------------------------------------------------------------

namespace {
    NetIn state(const fourdst::composition::Composition& composition, const double T, const double rho) {
        NetIn in;
        in.composition = composition;
        in.temperature = T;
        in.density = rho;
        in.tMax = 0.0; // ignored by evaluateRates
        in.dt0 = 0.0;
        in.energy = 0.0;
        return in;
    }
}

// evaluateRates at evolve()'s final state reproduces the final-state values evolve() reports
TEST(SolverEvaluateRates, PointSolverMatchesEvolveFinalState) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    PointSolverContext sctx(*context);
    sctx.set_stdout_logging(false);
    const PointSolver solver(shared.engine);

    NetIn in = state(gridfire::test::solar_composition(), 1.5e7, 150.0);
    in.tMax = 3e15;
    in.dt0 = 1e-12;
    const NetOut out = solver.evaluate(sctx, in);

    const gridfire::InstantaneousRates rates = solver.evaluateRates(sctx, state(out.composition, 1.5e7, 150.0));
    EXPECT_NEAR(rates.eps_nuc, out.eps_nuc, 1e-12 * std::abs(out.eps_nuc));
    EXPECT_NEAR(rates.dEps_dT, out.dEps_dT, 1e-12 * std::abs(out.dEps_dT));
    EXPECT_NEAR(rates.dEps_dRho, out.dEps_dRho, 1e-12 * std::abs(out.dEps_dRho));
    EXPECT_NEAR(rates.neutrino_energy_loss_rate, out.neutrino_energy_loss_rate, 1e-12 * out.neutrino_energy_loss_rate);
    EXPECT_EQ(rates.energy_channel_names, out.energy_channel_names);
    ASSERT_EQ(rates.eps_by_channel.size(), out.eps_by_channel.size());
    for (size_t c = 0; c < rates.eps_by_channel.size(); ++c) {
        EXPECT_NEAR(rates.eps_by_channel[c], out.eps_by_channel[c], 1e-12 * std::abs(out.eps_nuc)) << out.energy_channel_names[c];
    }
    EXPECT_EQ(rates.dydt.size(), shared.engine.getNetworkSpecies(*context).size()) << "dY/dt covers every network species";
    EXPECT_LT(rates.dydt.at(sp("H-1")), 0.0);
    EXPECT_GT(rates.neutrino_flux, 0.0);
}

// The result is the engine's own RHS and derivatives at that state, and calling it leaves integration unaffected
TEST(SolverEvaluateRates, IsAPassThroughThatLeavesIntegrationUnchanged) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const PointSolver solver(shared.engine);
    const auto comp = gridfire::test::solar_composition();

    PointSolverContext sctx(*context);
    sctx.set_stdout_logging(false);
    const auto rates = solver.evaluateRates(sctx, state(comp, 2e7, 80.0));
    const auto& network = shared.engine.getNetworkSpecies(*context);
    const auto engine_comp = gridfire::test::make_composition(network, gridfire::test::abundances_in_order(comp, network));
    const auto rhs = shared.engine.calculateRHSAndEnergy(*context, engine_comp, 0.02, 80.0, false);
    ASSERT_TRUE(rhs.has_value());
    EXPECT_EQ(rates.eps_nuc, rhs->nuclearEnergyGenerationRate);
    EXPECT_EQ(rates.dydt, rhs->dydt);

    NetIn in = state(comp, 1.5e7, 150.0);
    in.tMax = 1e14;
    in.dt0 = 1e-12;
    PointSolverContext plain(*shared.fresh_context());
    plain.set_stdout_logging(false);
    const NetOut reference = solver.evaluate(plain, in);
    (void)solver.evaluateRates(sctx, state(comp, 3e7, 10.0));
    const NetOut after_rates = solver.evaluate(sctx, in);
    EXPECT_EQ(after_rates.energy, reference.energy);
    EXPECT_EQ(after_rates.eps_nuc, reference.eps_nuc);
    EXPECT_EQ(after_rates.num_steps, reference.num_steps);
}

TEST(SolverEvaluateRates, GridSolverMatchesPointSolverPerZone) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const PointSolver point(shared.engine);
    const gridfire::solver::GridSolver grid(shared.engine, point);
    gridfire::solver::GridSolverContext gctx(*context);
    gctx.zone_completion_logging = false;

    std::vector<NetIn> inputs;
    for (const double T : {1.0e7, 1.5e7, 2.5e7, 4e7}) {
        inputs.push_back(state(gridfire::test::solar_composition(), T, 100.0));
    }
    const auto grid_rates = grid.evaluateRates(gctx, inputs);
    ASSERT_EQ(grid_rates.size(), inputs.size());
    PointSolverContext sctx(*context);
    for (size_t z = 0; z < inputs.size(); ++z) {
        const auto single = point.evaluateRates(sctx, inputs[z]);
        EXPECT_EQ(grid_rates[z].eps_nuc, single.eps_nuc) << "zone " << z;
        EXPECT_EQ(grid_rates[z].dEps_dT, single.dEps_dT) << "zone " << z;
        EXPECT_EQ(grid_rates[z].dEps_dRho, single.dEps_dRho) << "zone " << z;
    }
}

TEST(SolverEvaluateRates, RejectsInvalidStates) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const PointSolver solver(shared.engine);
    PointSolverContext sctx(*context);
    const auto comp = gridfire::test::solar_composition();
    EXPECT_THROW((void)solver.evaluateRates(sctx, state(comp, -1.0, 100.0)), gridfire::exceptions::SolverError);
    EXPECT_THROW((void)solver.evaluateRates(sctx, state(comp, 1.5e7, std::nan(""))), gridfire::exceptions::SolverError);

    // A species outside the network
    auto outside = comp;
    outside.registerSpecies(sp("U-238"));
    outside.setMolarAbundance(sp("U-238"), 1e-8);
    EXPECT_THROW((void)solver.evaluateRates(sctx, state(outside, 1.5e7, 100.0)), gridfire::exceptions::SolverError);

    // The wrong context type
    gridfire::solver::GridSolverContext wrong(*context);
    EXPECT_THROW((void)solver.evaluateRates(wrong, state(comp, 1.5e7, 100.0)), gridfire::exceptions::SolverError);
}

// ---------------------------------------------------------------------------------------------------------------------
// Per-zone failure reporting
// ---------------------------------------------------------------------------------------------------------------------

// A zone that cannot be integrated (T = 1e10 K with no step-failure restarts) is reported in its NetOut and in
// zone_failures; the other zones complete as if it were not there
TEST(GridSolverFailures, FailedZoneIsReportedAndOthersComplete) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    gridfire::config::GridFireConfig config;
    config.solver.pointSolver.maxStepFailureRestarts = 0;
    const PointSolver point(shared.engine, config);
    const gridfire::solver::GridSolver grid(shared.engine, point);
    gridfire::solver::GridSolverContext gctx(*context);
    gctx.zone_completion_logging = false;

    std::vector<NetIn> inputs = zone_inputs({1.5e7, 1e10, 2e7});
    inputs[1].density = 1e4;
    const auto results = grid.evaluate(gctx, inputs);
    ASSERT_EQ(results.size(), 3u);
    ASSERT_EQ(gctx.zone_failures.size(), 3u);

    EXPECT_FALSE(results[1].converged);
    EXPECT_NE(gctx.zone_failures[1], nullptr);
    EXPECT_FALSE(results[1].failure_message.empty());
    EXPECT_LT(results[1].time_reached, inputs[1].tMax);
    EXPECT_EQ(results[1].energy, 0.0);
    EXPECT_EQ(results[1].composition.getMolarAbundance(sp("H-1")), inputs[1].composition.getMolarAbundance(sp("H-1")))
        << "a failed zone carries its input composition";
    EXPECT_THROW(std::rethrow_exception(gctx.zone_failures[1]), gridfire::exceptions::CVODESolverFailureError);

    for (const size_t z : {0u, 2u}) {
        EXPECT_TRUE(results[z].converged) << "zone " << z;
        EXPECT_EQ(gctx.zone_failures[z], nullptr) << "zone " << z;
        EXPECT_EQ(results[z].time_reached, inputs[z].tMax) << "zone " << z;
        const NetOut alone = run(shared.engine, *context, inputs[z].composition, inputs[z].temperature,
                                 inputs[z].density, inputs[z].tMax, {.rel_tol = 1e-5, .abs_tol = 1e-8});
        EXPECT_EQ(results[z].energy, alone.energy) << "zone " << z << " is affected by the failed zone";
    }
}

TEST(GridSolverFailures, PointSolverReportsTimeReached) {
    const auto& shared = gridfire::test::main_sequence_engine();
    const auto context = shared.fresh_context();
    const NetOut out = run(shared.engine, *context, gridfire::test::solar_composition(), 1.5e7, 150.0, 1e13);
    EXPECT_TRUE(out.converged);
    EXPECT_EQ(out.time_reached, 1e13);
    EXPECT_TRUE(out.failure_message.empty());
}
