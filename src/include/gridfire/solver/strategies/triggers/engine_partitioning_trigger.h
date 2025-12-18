#pragma once

#include "gridfire/trigger/trigger_abstract.h"
#include "gridfire/trigger/trigger_result.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "fourdst/logging/logging.h"

#include <string>
#include <deque>
#include <memory>

/**
 * @file engine_partitioning_trigger.h
 * @brief CVODE-specific triggers that decide when to (re)partition the reaction network engine.
 *
 * @details
 * This header provides three concrete Trigger<CVODESolverStrategy::TimestepContext> implementations:
 *  - SimulationTimeTrigger: fires when the simulated time advances by a fixed interval.
 *  - OffDiagonalTrigger: fires when any off-diagonal Jacobian entry exceeds a threshold.
 *  - TimestepCollapseTrigger: fires when the timestep changes sharply relative to a moving average.
 *
 * It also provides a convenience factory (makeEnginePartitioningTrigger) composing these triggers
 * with logical combinators defined in trigger_logical.h.
 *
 * See the implementation for details: src/lib/solver/strategies/triggers/engine_partitioning_trigger.cpp
 *
 * @par Related headers:
 *  - trigger_abstract.h: base Trigger interface and lifecycle semantics
 *  - trigger_logical.h: AND/OR/NOT/EveryNth composition utilities
 */
namespace gridfire::trigger::solver::CVODE {
    /**
     * @class SimulationTimeTrigger
     * @brief Triggers when the current simulation time advances by at least a fixed interval.
     *
     * @details
     *  - check(ctx) returns true when (ctx.t - last_trigger_time) >= interval.
     *  - update(ctx) will, if check(ctx) is true, record ctx.t as the new last_trigger_time and
     *    store a small delta relative to the configured interval (for diagnostics/logging).
     *  - Counters (hits/misses/updates/resets) are maintained for diagnostics; they are
     *    mutable to allow updates from const check().
     *
     * @par Constraints/Errors:
     *  - Constructing with a non-positive interval throws std::invalid_argument.
     *
     * @note Thread-safety: not thread-safe; intended for single-threaded trigger evaluation.
     *
     * See also: engine_partitioning_trigger.cpp for the concrete logic and logging.
     */
    class SimulationTimeTrigger final : public Trigger<gridfire::solver::PointSolverTimestepContext> {
    public:
        /**
         * @brief Construct with a positive time interval between firings.
         * @param interval Strictly positive time interval (simulation units) between triggers.
         * @throws std::invalid_argument if interval <= 0.
         */
        explicit SimulationTimeTrigger(double interval);
        /**
         * @brief Evaluate whether enough simulated time has elapsed since the last trigger.
         * @param ctx CVODE timestep context providing the current simulation time (ctx.t).
         * @return true if (ctx.t - last_trigger_time) >= interval; false otherwise.
         *
         * @post increments hit/miss counters and may emit trace logs.
         */
        bool check(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /**
         * @brief Update internal state; if check(ctx) is true, advance last_trigger_time.
         * @param ctx CVODE timestep context.
         *
         * @note update() calls check(ctx) and, on success, records the overshoot delta
         * (ctx.t - last_trigger_time) - interval for diagnostics.
         */
        void update(const gridfire::solver::PointSolverTimestepContext &ctx) override;

        void step(const gridfire::solver::PointSolverTimestepContext &ctx) override;
        /**
         * @brief Reset counters and last trigger bookkeeping (time and delta) to zero.
         */
        void reset() override;

        /** @brief Stable human-readable name. */
        std::string name() const override;
        /**
         * @brief Structured explanation of the evaluation outcome.
         * @param ctx CVODE timestep context.
         * @return TriggerResult including name, value, and description.
         */
        TriggerResult why(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /** @brief Textual description including configured interval. */
        std::string describe() const override;
        /** @brief Number of true evaluations since last reset. */
        size_t numTriggers() const override;
        /** @brief Number of false evaluations since last reset. */
        size_t numMisses() const override;
    private:
        /** @brief Logger used for trace/error diagnostics. */
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /** @name Diagnostics counters */
        ///@{
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
        ///@}

        /** @brief Required time interval between successive triggers. */
        double m_interval;

        /** @brief Time at which the trigger last fired; initialized to 0. */
        mutable double m_last_trigger_time = 0.0;
        /** @brief Overshoot relative to interval at the last firing; for diagnostics. */
        mutable double m_last_trigger_time_delta = 0.0;
    };

    /**
     * @class OffDiagonalTrigger
     * @brief Triggers when any off-diagonal Jacobian entry magnitude exceeds a threshold.
     *
     * Semantics:
     *  - Iterates over all species pairs (row != col) and queries the engine's Jacobian
     *    via ctx.engine.getJacobianMatrixEntry(row, col). If any |entry| > threshold,
     *    check(ctx) returns true (short-circuits on first exceedance).
     *  - update(ctx) only records an update counter; it does not cache Jacobian values.
     *
     * @note Complexity: O(S^2) per check for S species (due to dense scan).
     *
     * @par Constraints/Errors:
     *  - Constructing with threshold < 0 throws std::invalid_argument.
     *
     * @par See also
     *  - engine_partitioning_trigger.cpp for concrete logic and trace logging.
     */
    class OffDiagonalTrigger final : public Trigger<gridfire::solver::PointSolverTimestepContext> {
    public:
        /**
         * @brief Construct with a non-negative magnitude threshold.
         * @param threshold Off-diagonal Jacobian magnitude threshold (>= 0).
         * @throws std::invalid_argument if threshold < 0.
         */
        explicit OffDiagonalTrigger(double threshold);
        /**
         * @brief Check if any off-diagonal Jacobian entry exceeds the threshold.
         * @param ctx CVODE timestep context providing access to engine species and Jacobian.
         * @return true if max_{i!=j} |J(i,j)| > threshold; false otherwise.
         *
         * @post increments hit/miss counters and may emit trace logs.
         */
        bool check(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /**
         * @brief Record an update; does not mutate any Jacobian-related state.
         * @param ctx CVODE timestep context (unused except for symmetry with interface).
         */
        void update(const gridfire::solver::PointSolverTimestepContext &ctx) override;
        void step(const gridfire::solver::PointSolverTimestepContext &ctx) override;
        /** @brief Reset counters to zero. */
        void reset() override;

        /** @brief Stable human-readable name. */
        std::string name() const override;
        /**
         * @brief Structured explanation of the evaluation outcome.
         * @param ctx CVODE timestep context.
         */
        TriggerResult why(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /** @brief Textual description including configured threshold. */
        std::string describe() const override;
        /** @brief Number of true evaluations since last reset. */
        size_t numTriggers() const override;
        /** @brief Number of false evaluations since last reset. */
        size_t numMisses() const override;
    private:
        /** @brief Logger used for trace/error diagnostics. */
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /** @name Diagnostics counters */
        ///@{
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
        ///@}

        /** @brief Off-diagonal magnitude threshold (>= 0). */
        double m_threshold;
    };

    /**
     * @class TimestepCollapseTrigger
     * @brief Triggers when the timestep deviates from its recent average beyond a threshold.
     *
     * @details
     *  - Maintains a sliding window of recent dt values (size = windowSize).
     *  - check(ctx):
     *      - If the window is empty, returns false.
     *      - Computes the arithmetic mean of values in the window and compares either:
     *          - relative: |dt - mean| / mean >= threshold
     *          - absolute: |dt - mean| >= threshold
     *  - update(ctx): pushes ctx.dt into the fixed-size window (dropping oldest when full).
     *
     * @par Constraints/Errors:
     *  - threshold must be >= 0.
     *  - If relative==true, threshold must be in [0, 1]. Violations throw std::invalid_argument.
     *
     * @note
     *  - With windowSize==1, the mean is the most recent prior dt.
     *  - Counter fields are mutable to allow updates during const check().
     *
     * See also: engine_partitioning_trigger.cpp for exact logic and logging.
     */
    class TimestepCollapseTrigger final : public Trigger<gridfire::solver::PointSolverTimestepContext> {
    public:
        /**
         * @brief Construct with threshold and relative/absolute mode; window size defaults to 1.
         * @param threshold Non-negative threshold; if relative, must be in [0, 1].
         * @param relative If true, use relative deviation; otherwise use absolute deviation.
         * @throws std::invalid_argument on invalid threshold constraints.
         */
        explicit TimestepCollapseTrigger(double threshold, bool relative);
        /**
         * @brief Construct with threshold, mode, and window size.
         * @param threshold Non-negative threshold; if relative, must be in [0, 1].
         * @param relative If true, use relative deviation; otherwise use absolute deviation.
         * @param windowSize Number of dt samples to average over (>= 1 recommended).
         * @throws std::invalid_argument on invalid threshold constraints.
         */
        explicit TimestepCollapseTrigger(double threshold, bool relative, size_t windowSize);
        /**
         * @brief Evaluate whether the current dt deviates sufficiently from recent average.
         * @param ctx CVODE timestep context providing current dt.
         * @return true if deviation exceeds the configured threshold; false otherwise.
         *
         * @post increments hit/miss counters and may emit trace logs.
         */
        bool check(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /**
         * @brief Update sliding window with the most recent dt and increment update counter.
         * @param ctx CVODE timestep context.
         */
        void update(const gridfire::solver::PointSolverTimestepContext &ctx) override;
        void step(const gridfire::solver::PointSolverTimestepContext &ctx) override;
        /** @brief Reset counters and clear the dt window. */
        void reset() override;

        /** @brief Stable human-readable name. */
        std::string name() const override;
        /** @brief Structured explanation of the evaluation outcome. */
        TriggerResult why(const gridfire::solver::PointSolverTimestepContext &ctx) const override;
        /** @brief Textual description including threshold, mode, and window size. */
        std::string describe() const override;
        /** @brief Number of true evaluations since last reset. */
        size_t numTriggers() const override;
        /** @brief Number of false evaluations since last reset. */
        size_t numMisses() const override;
    private:
        /** @brief Logger used for trace/error diagnostics. */
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /** @name Diagnostics counters */
        ///@{
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;
        ///@}

        /** @brief Threshold for absolute or relative deviation. */
        double m_threshold;
        /** @brief When true, use relative deviation; otherwise absolute deviation. */
        bool m_relative;
        /** @brief Number of dt samples to maintain in the moving window. */
        size_t m_windowSize;

        /** @brief Sliding window of recent timesteps (most recent at back). */
        std::deque<double> m_timestep_window;
    };

    class ConvergenceFailureTrigger final : public Trigger<gridfire::solver::PointSolverTimestepContext> {
    public:
        explicit ConvergenceFailureTrigger(size_t totalFailures, float relativeFailureRate, size_t windowSize);

        bool check(const gridfire::solver::PointSolverTimestepContext &ctx) const override;

        void update(const gridfire::solver::PointSolverTimestepContext &ctx) override;

        void step(const gridfire::solver::PointSolverTimestepContext &ctx) override;

        void reset() override;

        [[nodiscard]] std::string name() const override;

        [[nodiscard]] std::string describe() const override;

        [[nodiscard]] TriggerResult why(const gridfire::solver::PointSolverTimestepContext &ctx) const override;

        [[nodiscard]] size_t numTriggers() const override;

        [[nodiscard]] size_t numMisses() const override;

    private:
        /** @brief Logger used for trace/error diagnostics. */
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /** @name Diagnostics counters */
        ///@{
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;

        size_t m_totalFailures;
        float m_relativeFailureRate;
        size_t m_windowSize;

        std::deque<size_t> m_window;

    private:
        float current_mean() const;
        bool abs_failure(const gridfire::solver::PointSolverTimestepContext& ctx) const;
        bool rel_failure(const gridfire::solver::PointSolverTimestepContext& ctx) const;
    };

    /**
     * @brief Compose a trigger suitable for deciding engine re-partitioning during CVODE solves.
     *
     * Policy (as of implementation):
     *  - OR of the following three conditions:
     *    1) Every 1000th firing of SimulationTimeTrigger(simulationTimeInterval)
     *    2) OffDiagonalTrigger(offDiagonalThreshold)
     *    3) Every 10th firing of TimestepCollapseTrigger(timestepGrowthThreshold,
     *       timestepGrowthRelative, timestepGrowthWindowSize)
     *
     * See engine_partitioning_trigger.cpp for construction details using OrTrigger and
     * EveryNthTrigger from trigger_logical.h.
     *
     * @param simulationTimeInterval Interval used by SimulationTimeTrigger (> 0).
     * @param offDiagonalThreshold Off-diagonal Jacobian magnitude threshold (>= 0).
     * @param timestepCollapseRatio Threshold for timestep deviation (>= 0, and <= 1 when relative).
     * @param maxConvergenceFailures Window size for timestep averaging (>= 1 recommended).
     * @return A unique_ptr to a composed Trigger<TimestepContext> implementing the policy above.
     *
     * @note The exact policy is subject to change; this function centralizes that decision.
     */
    std::unique_ptr<Trigger<gridfire::solver::PointSolverTimestepContext>> makeEnginePartitioningTrigger(
        double simulationTimeInterval,
        double offDiagonalThreshold,
        double timestepCollapseRatio,
        size_t maxConvergenceFailures
    );
}
