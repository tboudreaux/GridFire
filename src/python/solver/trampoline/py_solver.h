#pragma once

#include "gridfire/solver/solver.h"

#include <vector>
#include <tuple>
#include <string>
#include <any>

class PyDynamicNetworkSolverStrategy final : public gridfire::solver::DynamicNetworkSolverStrategy {
    explicit PyDynamicNetworkSolverStrategy(gridfire::engine::DynamicEngine &engine) : gridfire::solver::DynamicNetworkSolverStrategy(engine) {}
    gridfire::NetOut evaluate(const gridfire::NetIn &netIn) override;
    void set_callback(const std::any &callback) override;
    [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;
};

class PySolverContextBase final : public gridfire::solver::SolverContextBase {
public:
    [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;
};