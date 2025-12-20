#include "gridfire/solver/strategies/GridSolver.h"

#include "gridfire/exceptions/error_solver.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/utils/macros.h"
#include "gridfire/utils/gf_omp.h"

#include <cstdio>
#include <print>

namespace gridfire::solver {
    void GridSolverContext::init() {}
    void GridSolverContext::reset() {
        solver_workspaces.clear();
        timestep_callbacks.clear();
    }

    void GridSolverContext::set_callback(const std::function<void(const TimestepContextBase &)> &callback) {
        for (auto &cb : timestep_callbacks) {
            cb = callback;
        }
    }

    void GridSolverContext::set_callback(const std::function<void(const TimestepContextBase &)> &callback, const size_t zone_idx) {
        if (zone_idx >= timestep_callbacks.size()) {
            throw exceptions::SolverError("GridSolverContext::set_callback: zone_idx out of range.");
        }
        timestep_callbacks[zone_idx] = callback;
    }

    void GridSolverContext::clear_callback() {
        for (auto &cb : timestep_callbacks) {
            cb = nullptr;
        }
    }

    void GridSolverContext::clear_callback(const size_t zone_idx) {
        if (zone_idx >= timestep_callbacks.size()) {
            throw exceptions::SolverError("GridSolverContext::clear_callback: zone_idx out of range.");
        }
        timestep_callbacks[zone_idx] = nullptr;
    }

    void GridSolverContext::set_stdout_logging(const bool enable) {
        zone_stdout_logging = enable;
    }

    void GridSolverContext::set_detailed_logging(const bool enable) {
        zone_detailed_logging = enable;
    }

    GridSolverContext::GridSolverContext(
        const engine::scratch::StateBlob &ctx_template
    ) :
    ctx_template(ctx_template) {}


    GridSolver::GridSolver(
        const engine::DynamicEngine &engine,
        const SingleZoneDynamicNetworkSolver &solver
    ) :
    MultiZoneNetworkSolver(engine, solver) {
        GF_PAR_INIT();
    }

    std::vector<NetOut> GridSolver::evaluate(
        SolverContextBase& ctx,
        const std::vector<NetIn>& netIns
    ) const {
        auto* sctx_p = dynamic_cast<GridSolverContext*>(&ctx);
        if (!sctx_p) {
            throw exceptions::SolverError("GridSolver::evaluate: SolverContextBase is not of type GridSolverContext.");
        }

        const size_t n_zones = netIns.size();
        if (n_zones == 0) { return {}; }

        std::vector<NetOut> results(n_zones);

        sctx_p->solver_workspaces.resize(n_zones);

        GF_OMP(
            parallel for default(none) shared(sctx_p, n_zones),
        for (size_t zone_idx = 0; zone_idx < n_zones; ++zone_idx)) {
            sctx_p->solver_workspaces[zone_idx] = std::make_unique<PointSolverContext>(sctx_p->ctx_template);
            sctx_p->solver_workspaces[zone_idx]->set_stdout_logging(sctx_p->zone_stdout_logging);
            sctx_p->solver_workspaces[zone_idx]->set_detailed_logging(sctx_p->zone_detailed_logging);
        }

        GF_OMP(
            parallel for default(none) shared(results, sctx_p, netIns, n_zones),
        for (size_t zone_idx = 0; zone_idx < n_zones; ++zone_idx)) {
            try {
                results[zone_idx] = m_solver.evaluate(
                    *sctx_p->solver_workspaces[zone_idx],
                    netIns[zone_idx]
                );
            } catch (exceptions::GridFireError& e) {
                std::println("CVODE Solver Failure in zone {}: {}", zone_idx, e.what());
            }
            if (sctx_p->zone_completion_logging) {
                std::println("Thread {} completed zone {}", GF_OMP_THREAD_NUM, zone_idx);
            }
        }
        return results;
    }
}
