#pragma once

#include "fourdst/config/config.h"

namespace gridfire::config {
    struct CVODESolverConfig {
        double absTol = 1.0e-8;
        double relTol = 1.0e-5;
    };


    struct SpectralSolverConfig {
        struct Trigger {
            double timestepCollapseRatio = 0.5;
            size_t maxConvergenceFailures = 2;
            double relativeFailureRate = 0.5;
            size_t windowSize = 10;
        };
        struct MonitorFunctionConfig {
            double structure_weight = 1.0;
            double abundance_weight = 10.0;
            double alpha = 0.2;
            double beta = 0.8;
        };
        struct BasisConfig {
            size_t num_elements = 50;
        };
        double absTol = 1.0e-8;
        double relTol = 1.0e-5;
        size_t degree = 3;
        MonitorFunctionConfig monitorFunction;
        BasisConfig basis;
        Trigger trigger;
    };

    struct SolverConfig {
        CVODESolverConfig cvode;
        SpectralSolverConfig spectral;
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