#pragma once

#include "gridfire/solver/solver.h"

#include <vector>
#include <tuple>
#include <string>
#include <any>

class PySingleZoneDynamicNetworkSolver final : public gridfire::solver::SingleZoneDynamicNetworkSolver {
public:
    explicit PySingleZoneDynamicNetworkSolver(const gridfire::engine::DynamicEngine &engine) : gridfire::solver::SingleZoneDynamicNetworkSolver(engine) {}

    gridfire::NetOut evaluate(
        gridfire::solver::SolverContextBase &solver_ctx,
        const gridfire::NetIn &netIn
    ) const override;
};

class PyMultiZoneDynamicNetworkSolver final : public gridfire::solver::MultiZoneDynamicNetworkSolver {
public:
    explicit PyMultiZoneDynamicNetworkSolver(
        const gridfire::engine::DynamicEngine &engine,
        const gridfire::solver::SingleZoneDynamicNetworkSolver &local_solver
    ) : gridfire::solver::MultiZoneDynamicNetworkSolver(engine, local_solver) {}

    std::vector<gridfire::NetOut> evaluate(
        gridfire::solver::SolverContextBase &solver_ctx,
        const std::vector<gridfire::NetIn> &netIns
    ) const override;
};

class PyTimestepContextBase final : public gridfire::solver::TimestepContextBase {
public:
    [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;
};

class PySolverContextBase final : public gridfire::solver::SolverContextBase {
public:
    void init() override;
    void set_stdout_logging(bool enable) override;
    void set_detailed_logging(bool enable) override;
};
