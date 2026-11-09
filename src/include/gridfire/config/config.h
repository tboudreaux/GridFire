#pragma once

#include "fourdst/config/config.h"
#include "gridfire/engine/types/jacobian_method.h"
#include "gridfire/screening/screening_types.h"

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
        /// Clamp negative CVODE iterates to zero before evaluating the RHS and Jacobian (the v1.2.0 behaviour).
        /// Off by default: clamping flattens the RHS of a species below zero, removing the restoring term of its
        /// destruction reactions while the (exact) Jacobian still reports it, so Newton converges to negative
        /// abundances and steps fail the positivity constraint.
        bool clampNegativeIterates = false;
        /// When CVODE gives up on a step (repeated error-test, convergence or positivity-constraint failures), restart
        /// BDF from the last accepted state (order 1, fresh history) up to this many times per evaluate call before
        /// failing. The integration interval and tolerances are unchanged.
        int maxStepFailureRestarts = 3;
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
        screening::ScreeningType screening = screening::ScreeningType::BLENDED; ///< Electron screening of reaction rates (default BLENDED) (WEAK: Salpeter 1954; CD09: Chugunov & DeWitt 2009, weak to strong coupling; BLENDED: Salpeter handing over to CD09 as coupling grows; BARE: none)
        screening::ScreeningBlendWindow screeningBlend; ///< Where BLENDED screening hands over from Salpeter to CD09 (weak-screening exponent H)
        /// Species carried in the composition without any reactions (e.g. the metals a stellar code lumps into Z,
        /// represented by one isotope): they count towards Y_e and screening, and their abundances never change.
        /// A passive species may not take part in any reaction of the network.
        std::vector<std::string> passiveSpecies;
        engine::JacobianMethod jacobian = engine::JacobianMethod::ANALYTIC; ///< How GraphEngine computes the species Jacobian
        /// Include reverse (photodisintegration) rates, from detailed balance, for strong reactions with two reactants
        /// and two products. Negligible for hydrogen and helium burning, so off by default.
        bool useReverseReactions = false;

    };

    struct GridFireConfig {
        SolverConfig solver;
        EngineConfig engine;
    };


}