#include "gridfire/solver/strategies/GridSolver.h"

#include "gridfire/exceptions/error_solver.h"
#include "gridfire/solver/strategies/PointSolver.h"
#include "gridfire/utils/macros.h"
#include "gridfire/utils/gf_omp.h"

#include <algorithm>
#include <exception>
#include <format>
#include <map>
#include <cmath>

#include <cstdio>
#include <print>

namespace gridfire::solver {
    void GridSolverContext::init() {}
    void GridSolverContext::reset() {
        solver_workspaces.clear();
        timestep_callbacks.clear();
    }

    void GridSolverContext::set_callback(const std::function<void(const TimestepContextBase &)> &callback) {
        for (auto &cb : timestep_callbacks) {
            cb = callback;
        }
    }

    void GridSolverContext::set_callback(const std::function<void(const TimestepContextBase &)> &callback, const size_t zone_idx) {
        if (zone_idx >= timestep_callbacks.size()) {
            throw exceptions::SolverError("GridSolverContext::set_callback: zone_idx out of range.");
        }
        timestep_callbacks[zone_idx] = callback;
    }

    void GridSolverContext::clear_callback() {
        for (auto &cb : timestep_callbacks) {
            cb = nullptr;
        }
    }

    void GridSolverContext::clear_callback(const size_t zone_idx) {
        if (zone_idx >= timestep_callbacks.size()) {
            throw exceptions::SolverError("GridSolverContext::clear_callback: zone_idx out of range.");
        }
        timestep_callbacks[zone_idx] = nullptr;
    }

    void GridSolverContext::set_stdout_logging(const bool enable) {
        zone_stdout_logging = enable;
    }

    void GridSolverContext::set_detailed_logging(const bool enable) {
        zone_detailed_logging = enable;
    }

    void GridSolverContext::set_shared_composition_groups(std::vector<int> group_of_zone, std::vector<double> zone_mass) {
        if (group_of_zone.size() != zone_mass.size()) {
            throw exceptions::SolverError(std::format(
                "set_shared_composition_groups: {} group ids but {} zone masses.", group_of_zone.size(), zone_mass.size()));
        }
        for (size_t z = 0; z < group_of_zone.size(); ++z) {
            if (group_of_zone[z] < -1) {
                throw exceptions::SolverError(std::format(
                    "set_shared_composition_groups: zone {} has group id {}; use -1 for a zone burnt on its own.", z, group_of_zone[z]));
            }
            if (group_of_zone[z] >= 0 && (!std::isfinite(zone_mass[z]) || zone_mass[z] <= 0.0)) {
                throw exceptions::SolverError(std::format(
                    "set_shared_composition_groups: zone {} is in group {} but its mass is {}.", z, group_of_zone[z], zone_mass[z]));
            }
        }
        shared_composition_groups = std::move(group_of_zone);
        zone_masses = std::move(zone_mass);
    }

    void GridSolverContext::clear_shared_composition_groups() {
        shared_composition_groups.clear();
        zone_masses.clear();
    }

    std::vector<int> shared_composition_groups_from_mask(const std::vector<bool>& flagged) {
        std::vector<int> groups(flagged.size(), -1);
        int next = 0;
        for (size_t z = 0; z < flagged.size(); ++z) {
            if (!flagged[z]) { continue; }
            if (z > 0 && flagged[z - 1]) {
                groups[z] = groups[z - 1];
            } else {
                groups[z] = next++;
            }
        }
        return groups;
    }

    GridSolverContext::GridSolverContext(
        const engine::scratch::StateBlob &ctx_template
    ) :
    ctx_template(ctx_template) {}


    GridSolver::GridSolver(
        const engine::DynamicEngine &engine,
        const SingleZoneDynamicNetworkSolver &solver
    ) :
    MultiZoneNetworkSolver(engine, solver) {
        GF_PAR_INIT();
    }

    std::vector<NetOut> GridSolver::evaluate(
        SolverContextBase& ctx,
        const std::vector<NetIn>& netIns
    ) const {
        auto* sctx_p = dynamic_cast<GridSolverContext*>(&ctx);
        if (!sctx_p) {
            throw exceptions::SolverError("GridSolver::evaluate: SolverContextBase is not of type GridSolverContext.");
        }

        const size_t n_zones = netIns.size();
        if (n_zones == 0) { return {}; }

        std::vector<NetOut> results(n_zones);

        sctx_p->solver_workspaces.resize(n_zones);

        GF_OMP(
            parallel for default(none) shared(sctx_p, n_zones),
        for (size_t zone_idx = 0; zone_idx < n_zones; ++zone_idx)) {
            sctx_p->solver_workspaces[zone_idx] = std::make_unique<PointSolverContext>(sctx_p->ctx_template);
            sctx_p->solver_workspaces[zone_idx]->set_stdout_logging(sctx_p->zone_stdout_logging);
            sctx_p->solver_workspaces[zone_idx]->set_detailed_logging(sctx_p->zone_detailed_logging);
        }

        sctx_p->zone_failures.assign(n_zones, nullptr);

        // Tasks: one per zone, or one per shared-composition group
        const bool grouped = !sctx_p->shared_composition_groups.empty();
        if (grouped && sctx_p->shared_composition_groups.size() != n_zones) {
            throw exceptions::SolverError(std::format(
                "GridSolver::evaluate: shared-composition groups are set for {} zones but {} zones were given.",
                sctx_p->shared_composition_groups.size(), n_zones));
        }
        std::vector<std::vector<size_t>> tasks;
        std::map<int, size_t> task_of_group;
        for (size_t z = 0; z < n_zones; ++z) {
            const int group = grouped ? sctx_p->shared_composition_groups[z] : -1;
            if (group < 0) {
                tasks.push_back({z});
                continue;
            }
            const auto [it, inserted] = task_of_group.try_emplace(group, tasks.size());
            if (inserted) { tasks.emplace_back(); }
            tasks[it->second].push_back(z);
        }
        const size_t n_tasks = tasks.size();

        // A zone (or group) that fails is reported in its NetOut (converged = false, its input composition, the time
        // reached) and in zone_failures, so one bad zone neither stops the others nor needs an exception to leave the
        // parallel region
        const auto record_failure = [&](const size_t zone, const std::exception_ptr& error, const size_t workspace) {
            sctx_p->zone_failures[zone] = error;
            NetOut failed{};
            failed.composition = netIns[zone].composition;
            failed.converged = false;
            if (const auto* point = dynamic_cast<const PointSolverContext*>(sctx_p->solver_workspaces[workspace].get())) {
                failed.time_reached = point->last_accepted_time;
            }
            try {
                std::rethrow_exception(error);
            } catch (const std::exception& e) {
                failed.failure_message = e.what();
            } catch (...) {
                failed.failure_message = "unknown error";
            }
            results[zone] = std::move(failed);
        };

        GF_OMP(
            parallel for schedule(dynamic, 1) default(none) shared(results, sctx_p, netIns, tasks, n_tasks, record_failure),
        for (size_t task = 0; task < n_tasks; ++task)) {
            const std::vector<size_t>& zones = tasks[task];
            const size_t first = zones.front();
            try {
                if (zones.size() == 1) {
                    results[first] = m_solver.evaluate(*sctx_p->solver_workspaces[first], netIns[first]);
                } else {
                    evaluate_group(*sctx_p, netIns, zones, results);
                }
            } catch (...) {
                const std::exception_ptr error = std::current_exception();
                for (const size_t zone : zones) {
                    record_failure(zone, error, first);
                }
            }
            if (sctx_p->zone_completion_logging) {
                std::println("Thread {} completed {} zone(s) starting at {}", GF_OMP_THREAD_NUM, zones.size(), first);
            }
        }
        return results;
    }

    void GridSolver::evaluate_group(
        GridSolverContext& ctx,
        const std::vector<NetIn>& netIns,
        const std::vector<size_t>& zones,
        std::vector<NetOut>& results
    ) const {
        const size_t first = zones.front();
        double total_mass = 0.0;
        for (const size_t zone : zones) {
            total_mass += ctx.zone_masses[zone];
        }

        // The group's composition is the mass-weighted average of its members' (identical inputs give that input)
        std::map<fourdst::atomic::Species, double> mixed;
        std::vector<ShellState> shells;
        shells.reserve(zones.size());
        for (const size_t zone : zones) {
            const double weight = ctx.zone_masses[zone] / total_mass;
            for (const auto& species : netIns[zone].composition.getRegisteredSpecies()) {
                mixed[species] += weight * netIns[zone].composition.getMolarAbundance(species);
            }
            shells.push_back({.temperature = netIns[zone].temperature, .density = netIns[zone].density, .weight = weight});
        }
        std::vector<fourdst::atomic::Species> species;
        std::vector<double> abundances;
        for (const auto& [sp, Y] : mixed) {
            species.push_back(sp);
            abundances.push_back(Y);
        }
        NetIn groupIn = netIns[first];
        groupIn.composition = fourdst::composition::Composition(species, abundances);

        SolverContextBase& workspace = *ctx.solver_workspaces[first];
        const SharedCompositionResult burnt = m_solver.evaluateSharedComposition(workspace, groupIn, shells);

        // Each member: the shared composition, its own energy and neutrino losses, and its own instantaneous rates
        for (size_t i = 0; i < zones.size(); ++i) {
            const size_t zone = zones[i];
            NetOut out = burnt.group;
            out.energy = burnt.shell_energy[i];
            out.specific_neutrino_energy_loss = burnt.shell_neutrino_energy_loss[i];
            out.specific_neutrino_flux = burnt.shell_neutrino_flux[i];

            NetIn state = netIns[zone];
            state.composition = burnt.group.composition;
            const InstantaneousRates rates = m_solver.evaluateRates(workspace, state);
            out.eps_nuc = rates.eps_nuc;
            out.dEps_dT = rates.dEps_dT;
            out.dEps_dRho = rates.dEps_dRho;
            out.neutrino_energy_loss_rate = rates.neutrino_energy_loss_rate;
            out.energy_channel_names = rates.energy_channel_names;
            out.eps_by_channel = rates.eps_by_channel;
            out.dlnEps_dlnT = 0.0;
            out.dlnEps_dlnRho = 0.0;
            if (rates.eps_nuc != 0.0 && std::isfinite(rates.eps_nuc)) {
                out.dlnEps_dlnT = state.temperature * rates.dEps_dT / rates.eps_nuc;
                out.dlnEps_dlnRho = state.density * rates.dEps_dRho / rates.eps_nuc;
            }
            results[zone] = std::move(out);
        }
    }

    std::vector<InstantaneousRates> GridSolver::evaluateRates(
        SolverContextBase& ctx,
        const std::vector<NetIn>& netIns
    ) const {
        auto* sctx_p = dynamic_cast<GridSolverContext*>(&ctx);
        if (!sctx_p) {
            throw exceptions::SolverError("GridSolver::evaluateRates: SolverContextBase is not of type GridSolverContext.");
        }
        const size_t n_zones = netIns.size();
        std::vector<InstantaneousRates> results(n_zones);
        if (n_zones == 0) { return results; }
        prepare_workspaces(*sctx_p, n_zones);

        // Exceptions cannot leave an OpenMP region, so the first one is kept and rethrown afterwards
        std::exception_ptr failure = nullptr;
        GF_OMP(
            parallel for default(none) shared(results, sctx_p, netIns, n_zones, failure),
        for (size_t zone_idx = 0; zone_idx < n_zones; ++zone_idx)) {
            try {
                results[zone_idx] = m_solver.evaluateRates(*sctx_p->solver_workspaces[zone_idx], netIns[zone_idx]);
            } catch (...) {
#if defined(GF_USE_OPENMP)
#pragma omp critical(gf_grid_evaluate_rates_failure)
#endif
                {
                    if (!failure) { failure = std::current_exception(); }
                }
            }
        }
        if (failure) {
            std::rethrow_exception(failure);
        }
        return results;
    }

    void GridSolver::prepare_workspaces(GridSolverContext& ctx, const size_t n_zones) {
        if (ctx.solver_workspaces.size() == n_zones &&
            std::ranges::all_of(ctx.solver_workspaces, [](const auto& w) { return w != nullptr; })) {
            return;
        }
        ctx.solver_workspaces.resize(n_zones);
        for (auto& workspace : ctx.solver_workspaces) {
            if (!workspace) {
                workspace = std::make_unique<PointSolverContext>(ctx.ctx_template);
                workspace->set_stdout_logging(ctx.zone_stdout_logging);
                workspace->set_detailed_logging(ctx.zone_detailed_logging);
            }
        }
    }
}
