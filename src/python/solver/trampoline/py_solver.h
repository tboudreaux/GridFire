#pragma once

#include "gridfire/solver/solver.h"

#include <vector>
#include <tuple>
#include <string>
#include <any>

class PyDynamicNetworkSolverStrategy final : public gridfire::solver::DynamicNetworkSolverStrategy {
    explicit PyDynamicNetworkSolverStrategy(gridfire::DynamicEngine &engine) : gridfire::solver::DynamicNetworkSolverStrategy(engine) {}
    gridfire::NetOut evaluate(const gridfire::NetIn &netIn) override;
    void set_callback(const std::any &callback) override;
    std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;
};