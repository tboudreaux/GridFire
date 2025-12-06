#pragma once

#include "fourdst/config/config.h"

namespace gridfire::config {
    struct CVODESolverConfig {
        double absTol = 1.0e-8;
        double relTol = 1.0e-5;
    };

    struct SolverConfig {
        CVODESolverConfig cvode;
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