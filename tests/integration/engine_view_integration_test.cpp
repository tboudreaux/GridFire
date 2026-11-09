// Integration tests: PointSolver driven through hand-built engine views.
//
// The shipped policies hand the solver a bare GraphEngine, so nothing else exercises solving through a view. Each test
// builds one view on top of the main sequence GraphEngine, integrates a solar-core point with it and compares the
// specific energy (eps_nuc averaged over the run) and final H-1 with the same run on the bare engine.
// MultiscalePartitioningEngineView is not supported by PointSolver; the test checks that it is rejected.
#include "gf_test_support.h"

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/views/engine_adaptive.h"
#include "gridfire/engine/views/engine_defined.h"
#include "gridfire/engine/views/engine_multiscale.h"
#include "gridfire/exceptions/error_solver.h"
#include "gridfire/policy/stellar_policy.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/utils/gf_omp.h"

#include "fourdst/atomic/species.h"

#include <gtest/gtest.h>

#include <format>
#include <memory>
#include <string>
#include <vector>

using gridfire::NetIn;
using gridfire::NetOut;
using gridfire::engine::DynamicEngine;
using gridfire::engine::GraphEngine;
using gridfire::engine::scratch::StateBlob;
using gridfire::solver::PointSolver;
using gridfire::solver::PointSolverContext;

namespace {
    constexpr double T = 1.5e7;
    constexpr double RHO = 150.0;
    constexpr double T_MAX = 3.15e13; // ~1 Myr

    /// A main sequence policy of its own: the views need a non-const base engine.
    struct Network {
        std::unique_ptr<gridfire::policy::MainSequencePolicy> policy;
        std::unique_ptr<gridfire::policy::ConstructionResults> results;
        GraphEngine* graph = nullptr;

        Network() {
            gridfire::test::quiet_logging();
            GF_PAR_INIT();
            policy = std::make_unique<gridfire::policy::MainSequencePolicy>(gridfire::test::solar_composition());
            results = std::make_unique<gridfire::policy::ConstructionResults>(policy->construct());
            graph = dynamic_cast<GraphEngine*>(policy->get_engine_stack().front().get());
        }
    };

    Network& network() {
        static Network n;
        return n;
    }

    NetOut run(const DynamicEngine& engine, const StateBlob& context) {
        PointSolverContext sctx(context);
        sctx.set_stdout_logging(false);
        sctx.rel_tol = 1e-8;
        sctx.abs_tol = 1e-14;
        const PointSolver solver(engine);
        NetIn in;
        in.composition = gridfire::test::solar_composition();
        in.temperature = T;
        in.density = RHO;
        in.tMax = T_MAX;
        in.dt0 = 1e-12;
        in.energy = 0.0;
        return solver.evaluate(sctx, in);
    }

    const NetOut& bare_result() {
        static const NetOut out = [] {
            Network& net = network();
            const auto ctx = net.results->scratch_blob->clone_structure();
            return run(*net.graph, *ctx);
        }();
        return out;
    }

    void expect_energy_matches_bare(const NetOut& view_out, const double tolerance, const std::string& view) {
        const NetOut& bare = bare_result();
        ASSERT_GT(bare.energy, 0.0);
        const double eps_view = view_out.energy / T_MAX;
        const double eps_bare = bare.energy / T_MAX;
        EXPECT_LT(gridfire::test::rel_diff(eps_view, eps_bare), tolerance)
            << std::format("{}: mean eps_nuc {:.8e} erg/g/s vs bare engine {:.8e}", view, eps_view, eps_bare);
        const auto& H1 = fourdst::atomic::species.at("H-1");
        EXPECT_LT(gridfire::test::rel_diff(view_out.composition.getMolarAbundance(H1), bare.composition.getMolarAbundance(H1)),
                  tolerance) << view << ": final H-1";
    }
}

TEST(EngineViewSolve, AdaptiveEngineView) {
    Network& net = network();
    gridfire::engine::AdaptiveEngineView view(*net.graph);
    const auto ctx = view.constructStateBlob(net.results->scratch_blob.get());
    const NetOut out = run(view, *ctx);
    // The view culls to the species that carry flow at the start, so the solver integrates a smaller system; what it
    // drops is negligible here (eps agrees with the bare engine to ~4e-9)
    // (the solver works on its own copy of ctx, so project a copy here to see the culled set)
    const auto probe = ctx->clone_structure();
    NetIn in;
    in.composition = gridfire::test::solar_composition();
    in.temperature = T;
    in.density = RHO;
    (void)view.project(*probe, in);
    EXPECT_LT(view.getNetworkSpecies(*probe).size(), net.graph->getNetworkSpecies(*probe).size());
    expect_energy_matches_bare(out, 1e-6, "AdaptiveEngineView");
}

TEST(EngineViewSolve, DefinedEngineViewOverWholeNetwork) {
    Network& net = network();
    std::vector<std::string> ids;
    const auto base_ctx = net.results->scratch_blob->clone_structure();
    for (const auto& reaction : net.graph->getNetworkReactions(*base_ctx)) {
        ids.emplace_back(reaction->id());
    }
    gridfire::engine::DefinedEngineView view(ids, *net.graph);
    const auto ctx = view.constructStateBlob(net.results->scratch_blob.get());
    const NetOut out = run(view, *ctx);
    // Same reactions as the bare engine, so the same answer up to integration error (measured ~1e-8)
    expect_energy_matches_bare(out, 1e-6, "DefinedEngineView");
}

TEST(EngineViewSolve, MultiscalePartitioningEngineViewIsRejected) {
    Network& net = network();
    // Solving through it does not converge with tolerance and does not conserve CNO nuclei, so PointSolver refuses
    // it up front, also when it sits below another view
    gridfire::engine::MultiscalePartitioningEngineView multiscale(*net.graph);
    EXPECT_THROW(PointSolver{multiscale}, gridfire::exceptions::UnsupportedEngineError);
    gridfire::engine::AdaptiveEngineView adaptive(multiscale);
    EXPECT_THROW(PointSolver{adaptive}, gridfire::exceptions::UnsupportedEngineError);
}

TEST(EngineViewSolve, ViewBlobFromNullptrBuildsBaseState) {
    Network& net = network();
    gridfire::engine::AdaptiveEngineView view(*net.graph);
    // With no blob given, a view starts from a fresh blob for its base engine
    const auto ctx = view.constructStateBlob(nullptr);
    const NetOut out = run(view, *ctx);
    expect_energy_matches_bare(out, 1e-6, "AdaptiveEngineView (nullptr blob)");
}
