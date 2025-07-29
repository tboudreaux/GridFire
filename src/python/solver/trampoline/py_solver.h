#pragma once

#include "gridfire/solver/solver.h"

#include <vector>

class PyDynamicNetworkSolverStrategy final : public gridfire::solver::DynamicNetworkSolverStrategy {
    explicit PyDynamicNetworkSolverStrategy(gridfire::DynamicEngine &engine) : gridfire::solver::DynamicNetworkSolverStrategy(engine) {}
    gridfire::NetOut evaluate(const gridfire::NetIn &netIn) override;
};