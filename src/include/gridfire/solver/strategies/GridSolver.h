#pragma once

#include "gridfire/solver/strategies/strategy_abstract.h"

#include <functional>

namespace gridfire::solver {
    struct GridSolverContext final : SolverContextBase {
        std::vector<std::unique_ptr<SolverContextBase>> solver_workspaces;
        std::vector<std::function<void(const TimestepContextBase&)>> timestep_callbacks;
        const engine::scratch::StateBlob& ctx_template;

        bool zone_completion_logging = true;
        bool zone_stdout_logging = false;
        bool zone_detailed_logging = false;

        void init() override;
        void reset();

        void set_callback(const std::function<void(const TimestepContextBase&)> &callback);
        void set_callback(const std::function<void(const TimestepContextBase&)> &callback, size_t zone_idx);

        void clear_callback();
        void clear_callback(size_t zone_idx);

        void set_stdout_logging(bool enable) override;
        void set_detailed_logging(bool enable) override;

        explicit GridSolverContext(const engine::scratch::StateBlob& ctx_template);
    };

    class GridSolver final : public MultiZoneDynamicNetworkSolver {
    public:
        GridSolver(
            const engine::DynamicEngine& engine,
            const SingleZoneDynamicNetworkSolver& solver
        );

        ~GridSolver() override = default;

        std::vector<NetOut> evaluate(
            SolverContextBase& ctx,
            const std::vector<NetIn>& netIns
        ) const override;
    };
}