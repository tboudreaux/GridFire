#pragma once

#include "fourdst/config/config.h"

#include <string>
#include <vector>

namespace gridfire::config {
    struct BoundaryFluxConfig {
        double relativeThreshold = 3e-8;
        double absoluteThreshold = 1e-24;
    };

    struct TriggerConfig {
        double offDiagonalThreshold = 1e10;
        double timestepCollapseRatio = 0.5;
        double maxConvergenceFailures  = 2;
        BoundaryFluxConfig boundaryFlux;
    };
    struct PointSolverConfig {
        double absTol = 1.0e-8;
        double relTol = 1.0e-5;
        TriggerConfig trigger;
    };

    struct SolverConfig {
        PointSolverConfig pointSolver;
    };

    struct AdaptiveEngineViewConfig {
        double relativeCullingThreshold = 1.0e-75;
    };

    struct EngineViewConfig {
        AdaptiveEngineViewConfig adaptiveEngineView;
    };

    struct EnergyChannelConfig {
        std::vector<std::string> fuels = {"H-1", "He-4", "C-12", "Ne-20", "O-16", "Si-28"}; ///< Fuels ordered by burning stage
        bool splitHydrogenPPCNO = true; ///< Report hydrogen burning as separate pp and CNO channels
    };

    struct EngineConfig {
        EngineViewConfig views;
        EnergyChannelConfig energyChannels;

    };

    struct GridFireConfig {
        SolverConfig solver;
        EngineConfig engine;
    };


}