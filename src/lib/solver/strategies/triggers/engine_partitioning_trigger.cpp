#include "gridfire/solver/strategies/triggers/engine_partitioning_trigger.h"
#include "gridfire/solver/strategies/PointSolver.h"

#include "gridfire/trigger/trigger_logical.h"
#include "gridfire/trigger/trigger_abstract.h"

#include "sundials/sundials_nvector.h"

#include "quill/LogMacros.h"

#include <memory>
#include <deque>
#include <string>

#include "gridfire/utils/utils.h"

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

    bool SimulationTimeTrigger::check(const gridfire::solver::PointSolverTimestepContext &ctx) const {
        if (ctx.t - m_last_trigger_time >= m_interval) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "SimulationTimeTrigger triggered at t = {}, last trigger time was {}, delta = {}", ctx.t, m_last_trigger_time, m_last_trigger_time_delta);
            return true;
        }
        m_misses++;
        return false;
    }

    void SimulationTimeTrigger::update(const gridfire::solver::PointSolverTimestepContext &ctx) {
        if (check(ctx)) {
            m_last_trigger_time_delta = (ctx.t - m_last_trigger_time) - m_interval;
            m_last_trigger_time = ctx.t;
            m_updates++;
        }
    }

    void SimulationTimeTrigger::step(
        const gridfire::solver::PointSolverTimestepContext &ctx
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

    TriggerResult SimulationTimeTrigger::why(const gridfire::solver::PointSolverTimestepContext &ctx) const {
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

    bool OffDiagonalTrigger::check(const gridfire::solver::PointSolverTimestepContext &ctx) const {
        //TODO : This currently does nothing
        return false;
    }


    void OffDiagonalTrigger::update(const gridfire::solver::PointSolverTimestepContext &ctx) {
        m_updates++;
    }

    void OffDiagonalTrigger::step(
        const gridfire::solver::PointSolverTimestepContext &ctx
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

    TriggerResult OffDiagonalTrigger::why(const gridfire::solver::PointSolverTimestepContext &ctx) const {
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

    bool TimestepCollapseTrigger::check(const gridfire::solver::PointSolverTimestepContext &ctx) const {
        if (m_timestep_window.size() < m_windowSize) {
            m_misses++;
            return false;
        }
        double averageTimestep = 0.0;
        for (const auto& dt : m_timestep_window) {
            averageTimestep += dt;
        }
        averageTimestep /= static_cast<double>(m_timestep_window.size());
        if (ctx.dt > averageTimestep) {
            m_misses++;
            return false; // Only trigger on timestep collapse (i.e., decrease in dt)
        }
        const double diff = std::abs(ctx.dt - averageTimestep);
        if (m_relative &&  1-(diff / averageTimestep) <= m_threshold) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "TimestepCollapseTrigger triggered at t = {} due to relative growth: dt = {}, average dt = {}, threshold = {}", ctx.t, ctx.dt, averageTimestep, m_threshold);
            return true;
        } else if (!m_relative && diff >= m_threshold) {
            m_hits++;
            LOG_TRACE_L2(m_logger, "TimestepCollapseTrigger triggered at t = {} due to absolute growth: dt = {}, average dt = {}, threshold = {}", ctx.t, ctx.dt, averageTimestep, m_threshold);
            return true;
        }
        m_misses++;
        return false;
    }

    void TimestepCollapseTrigger::update(const gridfire::solver::PointSolverTimestepContext &ctx) {
        m_updates++;
        m_timestep_window.clear();
    }

    void TimestepCollapseTrigger::step(
        const gridfire::solver::PointSolverTimestepContext &ctx
    ) {
        push_to_fixed_deque(m_timestep_window, ctx.dt, m_windowSize);
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
        const gridfire::solver::PointSolverTimestepContext &ctx
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
        const gridfire::solver::PointSolverTimestepContext &ctx
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
        const gridfire::solver::PointSolverTimestepContext &ctx
    ) {
        m_window.clear();
    }

    void ConvergenceFailureTrigger::step(
        const gridfire::solver::PointSolverTimestepContext &ctx
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

    TriggerResult ConvergenceFailureTrigger::why(const gridfire::solver::PointSolverTimestepContext &ctx) const {
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
        const gridfire::solver::PointSolverTimestepContext &ctx
    ) const {
        if (ctx.currentConvergenceFailures > m_totalFailures) {
            return true;
        }
        return false;
    }

    bool ConvergenceFailureTrigger::rel_failure(
        const gridfire::solver::PointSolverTimestepContext &ctx
    ) const {
        const float mean = current_mean();
        if (mean < 10) {
            return false; // Avoid triggering on small means
        }
        if (static_cast<float>(ctx.currentConvergenceFailures) - mean  > m_relativeFailureRate * mean) {
            return true;
        }
        return false;
    }

    BoundaryFluxTrigger::BoundaryFluxTrigger(
        const double relativeThreshold,
        const double absoluteThreshold
    ) :
    m_relativeThreshold(relativeThreshold),
    m_absoluteThreshold(absoluteThreshold) {
        if (m_relativeThreshold <= 0.0) {
            throw exceptions::GridFireError(std::format("Relative threshold must be positive and non zero, currently it is {}", m_relativeThreshold));
        }
    }

    void BoundaryFluxTrigger::step(const gridfire::solver::PointSolverTimestepContext &ctx) {
        // Does nothing; not a stateful trigger
    }


    bool BoundaryFluxTrigger::check(const gridfire::solver::PointSolverTimestepContext &ctx) const {
        // First get the current total flow through all active reactions
        sunrealtype* y_data = N_VGetArrayPointer(ctx.state);
        std::vector<double> Y(y_data, y_data + ctx.networkSpecies.size());
        // Adjust any tiny negative abundances to zero using std::ranges
        std::ranges::transform(
            Y,
            Y.begin(),
            [](const double y) {
                if (y < 0 && y > -1e-16) {
                    return 0.0;
                }
                return y;
            }
        );
        const fourdst::composition::Composition comp(ctx.networkSpecies, Y);

        const double net_active_flow = get_reaction_set_flow(
            ctx.engine.getNetworkReactions(ctx.state_ctx),
            ctx,
            comp,
            ctx.T9,
            ctx.rho,
            ReactionSetType::ACTIVE
        );

        const reaction::ReactionSet inactiveReactions = ctx.engine.getInactiveNetworkReactions(ctx.state_ctx);
        if (inactiveReactions.empty()) {
            m_misses++;
            return false; // No inactive reactions to consider
        }

        const double net_boundary_flow = get_reaction_set_flow(
            inactiveReactions,
            ctx,
            comp,
            ctx.T9,
            ctx.rho,
            ReactionSetType::INACTIVE
        );


        if (net_boundary_flow > m_absoluteThreshold) {
            m_hits++;
            return true;
        }

        const double relative_boundary_flow = net_boundary_flow / (net_active_flow + 1e-300); // Avoid division by zero
        if (relative_boundary_flow >= m_relativeThreshold) {
            m_hits++;
            return true;
        }

        m_misses++;
        return false;

    }

    void BoundaryFluxTrigger::update(const gridfire::solver::PointSolverTimestepContext &ctx) {
        // No-op since this is a stateless trigger
        m_updates++;
    }

    void BoundaryFluxTrigger::reset() {
        m_hits = 0;
        m_misses = 0;
        m_updates = 0;
        m_resets++;
    }

    std::string BoundaryFluxTrigger::name() const {
        return "BoundaryFluxTrigger";
    }

    std::string BoundaryFluxTrigger::describe() const {
        return std::format("BoundaryFluxTrigger(rel={}, abs={})", m_relativeThreshold, m_absoluteThreshold);
    }

    TriggerResult BoundaryFluxTrigger::why(const gridfire::solver::PointSolverTimestepContext &ctx) const {
        sunrealtype* y_data = N_VGetArrayPointer(ctx.state);
        const std::vector<double> Y(y_data, y_data + ctx.networkSpecies.size());
        const fourdst::composition::Composition comp(ctx.networkSpecies, Y);

        const double net_active_flow = get_reaction_set_flow(
            ctx.engine.getNetworkReactions(ctx.state_ctx),
            ctx,
            comp,
            ctx.T9,
            ctx.rho,
            ReactionSetType::ACTIVE
        );
        const reaction::ReactionSet inactiveReactions = ctx.engine.getInactiveNetworkReactions(ctx.state_ctx);
        const double net_boundary_flow = get_reaction_set_flow(
            inactiveReactions,
            ctx,
            comp,
            ctx.T9,
            ctx.rho,
            ReactionSetType::INACTIVE
        );

        TriggerResult result;
        result.name = name();
        if (check(ctx)) {
            result.value = true;
            result.description = std::format(
                "Triggered because boundary flux ({} mol/s) exceeded thresholds: absolute threshold = {} mol/s, relative threshold = {} (boundary flow = {} mol/s, active flow = {} mol/s)",
                net_boundary_flow,
                m_absoluteThreshold,
                m_relativeThreshold,
                net_boundary_flow,
                net_active_flow
            );
        } else {
            result.value = false;
            result.description = std::format(
                "Not triggered because boundary flux ({} mol/g/s) did not exceed thresholds: absolute threshold = {} mol/g/s, relative threshold = {} (boundary flow = {} mol/g/s, active flow = {} mol/g/s)",
                net_boundary_flow,
                m_absoluteThreshold,
                m_relativeThreshold,
                net_boundary_flow,
                net_active_flow
            );
        }

        return result;
    }

    size_t BoundaryFluxTrigger::numMisses() const {
        return m_misses;
    }

    double BoundaryFluxTrigger::get_reaction_set_flow(
        const reaction::ReactionSet &reactions,
        const gridfire::solver::PointSolverTimestepContext &ctx,
        const fourdst::composition::Composition &comp,
        const double T9,
        const double rho,
        const ReactionSetType type
    ) {
        double flow = 0.0;
        for (const auto& reaction: reactions) {
            double rFlow = 0.0;
            if (type == ReactionSetType::ACTIVE) {
                rFlow = ctx.engine.calculateMolarReactionFlow(ctx.state_ctx, *reaction, comp, T9, rho);
            } else {
                rFlow = ctx.engine.getInactiveReactionMolarReactionFlow(ctx.state_ctx, *reaction, comp, T9, rho);
            }
            flow += std::abs(rFlow);
        }

        return flow;

    }

    size_t BoundaryFluxTrigger::numTriggers() const {
        return m_hits;
    }

    std::unique_ptr<Trigger<gridfire::solver::PointSolverTimestepContext>> makeEnginePartitioningTrigger(
        const config::TriggerConfig& cfg
    ) {
        using ctx_t = gridfire::solver::PointSolverTimestepContext;

        // 1. INSTABILITY TRIGGERS
        auto convergenceFailureTrigger = std::make_unique<ConvergenceFailureTrigger>(
            cfg.maxConvergenceFailures,
            1.0f,
            10
        );

        auto timestepCollapseTrigger = std::make_unique<TimestepCollapseTrigger>(
            cfg.timestepCollapseRatio,
            true, // relative
            5
        );

        auto instabilityGroup = std::make_unique<OrTrigger<ctx_t>>(
            std::move(convergenceFailureTrigger),
            std::move(timestepCollapseTrigger)
        );

        // 2. MAINTENANCE TRIGGERS
        auto offDiagTrigger = std::make_unique<OffDiagonalTrigger>(cfg.offDiagonalThreshold);

        // 3. PREDICTIVE TRIGGERS
        auto boundaryFluxTrigger = std::make_unique<BoundaryFluxTrigger>(
            cfg.boundaryFlux.relativeThreshold,
            cfg.boundaryFlux.absoluteThreshold
        );

        // Combine boundary flux into off-diagonal trigger
        auto nonInstabilityGroup = std::make_unique<OrTrigger<ctx_t>>(
            std::move(offDiagTrigger),
            std::move(boundaryFluxTrigger)
        );

        // Combine: (Instability) OR (Structure Change)
        return std::make_unique<OrTrigger<ctx_t>>(
            std::move(instabilityGroup),
            std::move(nonInstabilityGroup)
        );
    }


}