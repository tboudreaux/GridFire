#pragma once

#include "fourdst/config/config.h"

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

    struct EngineConfig {
        EngineViewConfig views;

    };

    struct GridFireConfig {
        SolverConfig solver;
        EngineConfig engine;
    };


}