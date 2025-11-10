#include "gridfire/solver/strategies/triggers/engine_partitioning_trigger.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"

#include "gridfire/trigger/trigger_logical.h"
#include "gridfire/trigger/trigger_abstract.h"

#include "quill/LogMacros.h"

#include <memory>
#include <deque>
#include <string>

namespace {
    template <typename T>
    void push_to_fixed_deque(std::deque<T>& dq, T value, size_t max_size) {
        dq.push_back(value);
        if (dq.size() > max_size) {
            dq.pop_front();
        }
    }
}

namespace gridfire::trigger::solver::CVODE {
    SimulationTimeTrigger::SimulationTimeTrigger(double interval) : m_interval(interval) {
        if (interval <= 0.0) {
            LOG_ERROR(m_logger, "Interval must be positive, currently it is {}", interval);
            throw std::invalid_argument("Interval must be positive, currently it is " + std::to_string(interval));
        }
    }

    bool SimulationTimeTrigger::check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        if (ctx.t - m_last_trigger_time >= m_interval) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "SimulationTimeTrigger triggered at t = {}, last trigger time was {}, delta = {}", ctx.t, m_last_trigger_time, m_last_trigger_time_delta);
            return true;
        }
        m_misses++;
        return false;
    }

    void SimulationTimeTrigger::update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) {
        if (check(ctx)) {
            m_last_trigger_time_delta = (ctx.t - m_last_trigger_time) - m_interval;
            m_last_trigger_time = ctx.t;
            m_updates++;
        }
    }

    void SimulationTimeTrigger::step(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) {
        // --- SimulationTimeTrigger::step does nothing and is intentionally left blank --- //
    }

    void SimulationTimeTrigger::reset() {
        m_misses = 0;
        m_hits = 0;
        m_updates = 0;
        m_last_trigger_time = 0.0;
        m_last_trigger_time_delta = 0.0;
        m_resets++;
    }

    std::string SimulationTimeTrigger::name() const {
        return "Simulation Time Trigger";
    }

    TriggerResult SimulationTimeTrigger::why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        TriggerResult result;
        result.name = name();
        if (check(ctx)) {
            result.value = true;
            result.description = "Triggered because current time " + std::to_string(ctx.t) + " - last trigger time " + std::to_string(m_last_trigger_time - m_last_trigger_time_delta) + " >= interval " + std::to_string(m_interval);
        } else {
            result.value = false;
            result.description = "Not triggered because current time " + std::to_string(ctx.t) + " - last trigger time " + std::to_string(m_last_trigger_time) + " < interval " + std::to_string(m_interval);
        }
        return result;
    }

    std::string SimulationTimeTrigger::describe() const {
        return "SimulationTimeTrigger(interval=" + std::to_string(m_interval) + ")";
    }

    size_t SimulationTimeTrigger::numTriggers() const {
        return m_hits;
    }

    size_t SimulationTimeTrigger::numMisses() const {
        return m_misses;
    }

    OffDiagonalTrigger::OffDiagonalTrigger(
        double threshold
    ) : m_threshold(threshold) {
        if (threshold < 0.0) {
            LOG_ERROR(m_logger, "Threshold must be non-negative, currently it is {}", threshold);
            throw std::invalid_argument("Threshold must be non-negative, currently it is " + std::to_string(threshold));
        }
    }

    bool OffDiagonalTrigger::check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        for (const auto& rowSpecies : ctx.engine.getNetworkSpecies()) {
            for (const auto& colSpecies : ctx.engine.getNetworkSpecies()) {
                double DRowDCol = std::abs(ctx.engine.getJacobianMatrixEntry(rowSpecies, colSpecies));
                if (rowSpecies != colSpecies && DRowDCol > m_threshold) {
                    m_hits++;
                    LOG_TRACE_L2(m_logger, "OffDiagonalTrigger triggered at t = {} due to entry ({}, {}) = {}", ctx.t, rowSpecies.name(), colSpecies.name(), DRowDCol);
                    return true;
                }
            }
        }
        m_misses++;
        return false;
    }


    void OffDiagonalTrigger::update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) {
        m_updates++;
    }

    void OffDiagonalTrigger::step(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) {
        // --- OffDiagonalTrigger::step does nothing and is intentionally left blank --- //
    }

    void OffDiagonalTrigger::reset() {
        m_misses = 0;
        m_hits = 0;
        m_updates = 0;
        m_resets++;
    }

    std::string OffDiagonalTrigger::name() const {
        return "Off-Diagonal Trigger";
    }

    TriggerResult OffDiagonalTrigger::why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        TriggerResult result;
        result.name = name();

        if (check(ctx)) {
            result.value = true;
            result.description = "Triggered because an off-diagonal Jacobian entry exceeded the threshold " + std::to_string(m_threshold);
        } else {
            result.value = false;
            result.description = "Not triggered because no off-diagonal Jacobian entry exceeded the threshold " + std::to_string(m_threshold);
        }
        return result;
    }

    std::string OffDiagonalTrigger::describe() const {
        return "OffDiagonalTrigger(threshold=" + std::to_string(m_threshold) + ")";
    }

    size_t OffDiagonalTrigger::numTriggers() const {
        return m_hits;
    }

    size_t OffDiagonalTrigger::numMisses() const {
        return m_misses;
    }

    TimestepCollapseTrigger::TimestepCollapseTrigger(
        const double threshold,
        const bool relative
    ) : TimestepCollapseTrigger(threshold, relative, 1){}


    TimestepCollapseTrigger::TimestepCollapseTrigger(
        double threshold,
        const bool relative,
        const size_t windowSize
    ) : m_threshold(threshold), m_relative(relative), m_windowSize(windowSize) {
        if (threshold < 0.0) {
            LOG_ERROR(m_logger, "Threshold must be non-negative, currently it is {}", threshold);
            throw std::invalid_argument("Threshold must be non-negative, currently it is " + std::to_string(threshold));
        }
        if (relative && threshold > 1.0) {
            LOG_ERROR(m_logger, "Relative threshold must be between 0 and 1, currently it is {}", threshold);
            throw std::invalid_argument("Relative threshold must be between 0 and 1, currently it is " + std::to_string(threshold));
        }
    }

    bool TimestepCollapseTrigger::check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        if (m_timestep_window.empty()) {
            m_misses++;
            return false;
        }
        double averageTimestep = 0.0;
        for (const auto& dt : m_timestep_window) {
            averageTimestep += dt;
        }
        averageTimestep /= static_cast<double>(m_timestep_window.size());
        if (m_relative && (std::abs(ctx.dt - averageTimestep) / averageTimestep) >= m_threshold) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "TimestepCollapseTrigger triggered at t = {} due to relative growth: dt = {}, average dt = {}, threshold = {}", ctx.t, ctx.dt, averageTimestep, m_threshold);
            return true;
        } else if (!m_relative && std::abs(ctx.dt - averageTimestep) >= m_threshold) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "TimestepCollapseTrigger triggered at t = {} due to absolute growth: dt = {}, average dt = {}, threshold = {}", ctx.t, ctx.dt, averageTimestep, m_threshold);
            return true;
        }
        m_misses++;
        return false;
    }

    void TimestepCollapseTrigger::update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) {
        push_to_fixed_deque(m_timestep_window, ctx.dt, m_windowSize);
        m_updates++;
    }

    void TimestepCollapseTrigger::step(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) {
        // --- TimestepCollapseTrigger::step does nothing and is intentionally left blank --- //
    }

    void TimestepCollapseTrigger::reset() {
        m_misses = 0;
        m_hits = 0;
        m_updates = 0;
        m_resets++;
        m_timestep_window.clear();
    }

    std::string TimestepCollapseTrigger::name() const {
        return "TimestepCollapseTrigger";
    }

    TriggerResult TimestepCollapseTrigger::why(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) const {
        TriggerResult result;
        result.name = name();

        if (check(ctx)) {
            result.value = true;
            result.description = "Triggered because timestep change exceeded the threshold " + std::to_string(m_threshold);
        } else {
            result.value = false;
            result.description = "Not triggered because timestep change did not exceed the threshold " + std::to_string(m_threshold);
        }
        return result;
    }

    std::string TimestepCollapseTrigger::describe() const {
        return "TimestepCollapseTrigger(threshold=" + std::to_string(m_threshold) + ", relative=" + (m_relative ? "true" : "false") + ", windowSize=" + std::to_string(m_windowSize) + ")";
    }

    size_t TimestepCollapseTrigger::numTriggers() const {
        return m_hits;
    }

    size_t TimestepCollapseTrigger::numMisses() const {
        return m_misses;
    }

    ConvergenceFailureTrigger::ConvergenceFailureTrigger(
        const size_t totalFailures,
        const float relativeFailureRate,
        const size_t windowSize
    ) :
    m_totalFailures(totalFailures),
    m_relativeFailureRate(relativeFailureRate),
    m_windowSize(windowSize) {}

    bool ConvergenceFailureTrigger::check(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) const {
        if (m_window.size() != m_windowSize) {
            m_misses++;
            return false; // Short circuit if not enough data has been seen yet.
        }
        if (abs_failure(ctx) || rel_failure(ctx)) {
            m_hits++;
            return true;
        }
        m_misses++;
        return false;
    }

    void ConvergenceFailureTrigger::update(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) {
        // --- ConvergenceFailureTrigger::update does nothing and is intentionally left blank --- //
    }

    void ConvergenceFailureTrigger::step(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) {
        push_to_fixed_deque(m_window, ctx.currentConvergenceFailures, m_windowSize);
        m_updates++;
    }

    void ConvergenceFailureTrigger::reset() {
        m_window.clear();
        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
        m_resets++;
    }

    std::string ConvergenceFailureTrigger::name() const {
        return "ConvergenceFailureTrigger";
    }

    std::string ConvergenceFailureTrigger::describe() const {
        return "ConvergenceFailureTrigger(abs_failure_threshold=" + std::to_string(m_totalFailures) + ", rel_failure_threshold=" + std::to_string(m_relativeFailureRate) + ", windowSize=" + std::to_string(m_windowSize) + ")";
    }

    TriggerResult ConvergenceFailureTrigger::why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const {
        TriggerResult result;
        result.name = name();

        if (m_window.size() != m_windowSize) {
            result.value = false;
            result.description = "Not triggered because trigger has not seen sufficient data to build up window of size " + std::to_string(m_windowSize);
            return result;
        }
        if (abs_failure(ctx)) {
            result.value = true;
            result.description = std::format("Triggered because number of convergence failures ({}) exceeded absolute tolerances", ctx.currentConvergenceFailures);
            return result;
        }
        if (rel_failure(ctx)) {
            result.value = true;
            result.description = std::format("Triggered because number of convergence failures - the mean ({} - {}) exceeded tolerances relative to mean ({} * {})", ctx.currentConvergenceFailures, current_mean(), current_mean(), m_relativeFailureRate);
            return result;
        }

        result.value = false;
        result.description = "Not triggered because total number of convergence failures and relative number of convergence triggers did not grow sufficiently";
        return result;
    }

    size_t ConvergenceFailureTrigger::numTriggers() const {
        return m_hits;
    }

    size_t ConvergenceFailureTrigger::numMisses() const {
        return m_misses;
    }

    float ConvergenceFailureTrigger::current_mean() const {
        size_t acc = 0;
        for (const auto nlcfails: m_window) {
            acc += nlcfails;
        }
        return static_cast<float>(acc) / static_cast<float>(m_windowSize);
    }

    bool ConvergenceFailureTrigger::abs_failure(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) const {
        if (ctx.currentConvergenceFailures > m_totalFailures) {
            return true;
        }
        return false;
    }

    bool ConvergenceFailureTrigger::rel_failure(
        const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx
    ) const {
        const float mean = current_mean();
        if (static_cast<float>(ctx.currentConvergenceFailures) - mean  > m_relativeFailureRate * mean) {
            return true;
        }
        return false;
    }

    std::unique_ptr<Trigger<gridfire::solver::CVODESolverStrategy::TimestepContext>> makeEnginePartitioningTrigger(
        const double simulationTimeInterval,
        const double offDiagonalThreshold,
        const double timestepGrowthThreshold,
        const bool timestepGrowthRelative,
        const size_t timestepGrowthWindowSize
    ) {
        using ctx_t = gridfire::solver::CVODESolverStrategy::TimestepContext;

        // Create the individual conditions that can trigger a repartitioning
        // The current trigger logic is as follows
        // 1. Trigger every 1000th time that the simulation time exceeds the simulationTimeInterval
        // 2. OR if any off-diagonal Jacobian entry exceeds the offDiagonalThreshold
        // 3. OR every 10th time that the timestep growth exceeds the timestepGrowthThreshold (relative or absolute)
        // 4. OR if the number of convergence failures begins to grow

        // TODO: This logic likely needs to be revisited; however, for now it is easy enough to change and test and it works reasonably well
        auto simulationTimeTrigger = std::make_unique<EveryNthTrigger<ctx_t>>(std::make_unique<SimulationTimeTrigger>(simulationTimeInterval), 1000);
        auto offDiagTrigger = std::make_unique<OffDiagonalTrigger>(offDiagonalThreshold);
        auto timestepGrowthTrigger = std::make_unique<EveryNthTrigger<ctx_t>>(std::make_unique<TimestepCollapseTrigger>(timestepGrowthThreshold, timestepGrowthRelative, timestepGrowthWindowSize), 10);
        auto convergenceFailureTrigger = std::make_unique<ConvergenceFailureTrigger>(5, 1.0f, 10);

        // Combine the triggers using logical OR
        auto orTriggerA = std::make_unique<OrTrigger<ctx_t>>(std::move(simulationTimeTrigger), std::move(offDiagTrigger));
        auto orTriggerB = std::make_unique<OrTrigger<ctx_t>>(std::move(orTriggerA), std::move(timestepGrowthTrigger));
        // auto orTriggerC = std::make_unique<OrTrigger<ctx_t>>(std::move(orTriggerB), std::move(convergenceFailureTrigger));

        return convergenceFailureTrigger;
    }
}