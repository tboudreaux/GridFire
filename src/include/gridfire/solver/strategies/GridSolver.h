#pragma once

#include "gridfire/solver/strategies/strategy_abstract.h"

#include <exception>
#include <functional>

namespace gridfire::solver {
    struct GridSolverContext final : SolverContextBase {
        std::vector<std::unique_ptr<SolverContextBase>> solver_workspaces;
        /// Exception of each zone of the latest evaluate() call (null for zones that converged)
        std::vector<std::exception_ptr> zone_failures;

        /// Shared-composition group of each zone (-1: burnt on its own); empty when no groups are set
        std::vector<int> shared_composition_groups;
        /// Mass of each zone, used to weight the shells of a group (any units; only ratios within a group matter)
        std::vector<double> zone_masses;

        /**
         * @brief Burn the zones of each group as one composition from the next evaluate() on.
         *
         * Zones with the same non-negative group id share one composition, the mass-weighted average of their inputs,
         * which burns at the mass-weighted mean of their rates (SingleZoneNetworkSolver::evaluateSharedComposition):
         * a convective zone under instantaneous mixing. Each zone still gets its own energy, neutrino losses and
         * instantaneous rates. Zones with group -1 burn on their own, exactly as without groups. The members of a group
         * use the first member's tMax and dt0.
         *
         * @param group_of_zone One entry per zone of the next evaluate() call: -1, or a group id >= 0
         * @param zone_mass One entry per zone: the zone's mass (positive for every zone in a group)
         * @throws exceptions::SolverError if the sizes differ, an id is below -1, or a grouped zone's mass is not
         *         positive and finite
         */
        void set_shared_composition_groups(std::vector<int> group_of_zone, std::vector<double> zone_mass);

        /// Burn every zone on its own again.
        void clear_shared_composition_groups();
        std::vector<std::function<void(const TimestepContextBase&)>> timestep_callbacks;
        const engine::scratch::StateBlob& ctx_template;

        bool zone_completion_logging = true;
        bool zone_stdout_logging = false;
        bool zone_detailed_logging = false;

        void init() override;
        void reset();

        void set_callback(const std::function<void(const TimestepContextBase&)> &callback);
        void set_callback(const std::function<void(const TimestepContextBase&)> &callback, size_t zone_idx);

        void clear_callback();
        void clear_callback(size_t zone_idx);

        void set_stdout_logging(bool enable) override;
        void set_detailed_logging(bool enable) override;

        explicit GridSolverContext(const engine::scratch::StateBlob& ctx_template);
    };

    /**
     * @brief Group ids for GridSolverContext::set_shared_composition_groups from a per-zone flag, for zones ordered
     *        in space: each contiguous run of flagged zones (e.g. convective) is one group, the rest are -1.
     */
    [[nodiscard]] std::vector<int> shared_composition_groups_from_mask(const std::vector<bool>& flagged);

    class GridSolver final : public MultiZoneDynamicNetworkSolver {
    public:
        GridSolver(
            const engine::DynamicEngine& engine,
            const SingleZoneDynamicNetworkSolver& solver
        );

        ~GridSolver() override = default;

        std::vector<NetOut> evaluate(
            SolverContextBase& ctx,
            const std::vector<NetIn>& netIns
        ) const override;

        /**
         * @brief The single-zone solver's evaluateRates for every zone (in parallel with OpenMP), without integrating.
         */
        [[nodiscard]] std::vector<InstantaneousRates> evaluateRates(
            SolverContextBase& ctx,
            const std::vector<NetIn>& netIns
        ) const override;
    private:
        /// Burn the zones of one shared-composition group and fill their results
        void evaluate_group(
            GridSolverContext& ctx,
            const std::vector<NetIn>& netIns,
            const std::vector<size_t>& zones,
            std::vector<NetOut>& results
        ) const;

        /// One single-zone context per zone, made from ctx_template (reused when the zone count is unchanged).
        static void prepare_workspaces(GridSolverContext& ctx, size_t n_zones);
    };
}