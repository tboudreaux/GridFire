#pragma once

#include "gridfire/trigger/trigger_abstract.h"
#include "gridfire/trigger/trigger_result.h"
#include "gridfire/solver/strategies/CVODE_solver_strategy.h"
#include "fourdst/logging/logging.h"

#include <string>
#include <deque>
#include <memory>

namespace gridfire::trigger::solver::CVODE {
    class SimulationTimeTrigger final : public Trigger<gridfire::solver::CVODESolverStrategy::TimestepContext> {
    public:
        explicit SimulationTimeTrigger(double interval);
        bool check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        void update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) override;
        void reset() override;

        std::string name() const override;
        TriggerResult why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;

        double m_interval;

        mutable double m_last_trigger_time = 0.0;
        mutable double m_last_trigger_time_delta = 0.0;
    };

    class OffDiagonalTrigger final : public Trigger<gridfire::solver::CVODESolverStrategy::TimestepContext> {
    public:
        explicit OffDiagonalTrigger(double threshold);
        bool check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        void update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) override;
        void reset() override;

        std::string name() const override;
        TriggerResult why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;

        double m_threshold;
    };

    class TimestepCollapseTrigger final : public Trigger<gridfire::solver::CVODESolverStrategy::TimestepContext> {
    public:
        explicit TimestepCollapseTrigger(double threshold, bool relative);
        explicit TimestepCollapseTrigger(double threshold, bool relative, size_t windowSize);
        bool check(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        void update(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) override;
        void reset() override;

        std::string name() const override;
        TriggerResult why(const gridfire::solver::CVODESolverStrategy::TimestepContext &ctx) const override;
        std::string describe() const override;
        size_t numTriggers() const override;
        size_t numMisses() const override;
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        mutable size_t m_hits = 0;
        mutable size_t m_misses = 0;
        mutable size_t m_updates = 0;
        mutable size_t m_resets = 0;

        double m_threshold;
        bool m_relative;
        size_t m_windowSize;

        std::deque<double> m_timestep_window;
    };

    std::unique_ptr<Trigger<gridfire::solver::CVODESolverStrategy::TimestepContext>> makeEnginePartitioningTrigger(
        const double simulationTimeInterval,
        const double offDiagonalThreshold,
        const double timestepGrowthThreshold,
        const bool timestepGrowthRelative,
        const size_t timestepGrowthWindowSize
    );
}