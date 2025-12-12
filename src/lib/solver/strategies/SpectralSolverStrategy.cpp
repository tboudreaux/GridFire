#include "gridfire/solver/strategies/SpectralSolverStrategy.h"
#include "gridfire/utils/macros.h"
#include "gridfire/utils/table_format.h"

#include <sunlinsol/sunlinsol_dense.h>

#include "gridfire/utils/sundials.h"

#include "quill/LogMacros.h"
#include "sunmatrix/sunmatrix_dense.h"
#include <mutex>

#include "gridfire/utils/gf_omp.h"
#include "gridfire/utils/formatters/jacobian_format.h"
#include "gridfire/utils/logging.h"

namespace {
    std::pair<size_t, std::vector<double>> evaluate_bspline(
        double x,
        const gridfire::solver::SpectralSolverStrategy::SplineBasis& basis
    ) {
        const int p = basis.degree;
        const std::vector<double>& t = basis.knots;

        auto it = std::ranges::upper_bound(t, x);
        size_t i = std::distance(t.begin(), it) - 1;

        if (i < static_cast<size_t>(p)) i = p;
        if (i >= t.size() - 1 - p) i = t.size() - 2 - p;

        if (x >= t.back()) {
            i = t.size() - p - 2;
        }

        // Cox-de Boor algorithm
        std::vector<double> N(p + 1);
        std::vector<double> left(p + 1);
        std::vector<double> right(p + 1);

        N[0] = 1.0;

        for (int j = 1; j <= p; ++j) {
            left[j] = x - t[i + 1 - j];
            right[j] = t[i + j] - x;
            double saved = 0.0;

            for (int r = 0; r < j; ++r) {
                double temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;

                saved = left[j - r] * temp;
            }
            N[j] = saved;
        }
        return {i - p, N};
    }
}



namespace gridfire::solver {

    SpectralSolverStrategy::SpectralSolverStrategy(
        const engine::DynamicEngine& engine
    ) : MultiZoneNetworkSolver<engine::DynamicEngine> (engine) {
        LOG_INFO(m_logger, "Initializing SpectralSolverStrategy");

        utils::check_sundials_flag(SUNContext_Create(SUN_COMM_NULL, &m_sun_ctx), "SUNContext_Create", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

        m_absTol = m_config->solver.spectral.absTol;
        m_relTol = m_config->solver.spectral.relTol;

        LOG_INFO(m_logger, "SpectralSolverStrategy initialized successfully");

        GF_PAR_INIT()
    }

    SpectralSolverStrategy::~SpectralSolverStrategy() {
        LOG_INFO(m_logger, "Destroying SpectralSolverStrategy");


        if (m_cvode_mem) {
            CVodeFree(&m_cvode_mem);
            m_cvode_mem = nullptr;
        }

        if (m_LS) SUNLinSolFree(m_LS);
        if (m_J) SUNMatDestroy(m_J);
        if (m_Y) N_VDestroy(m_Y);
        if (m_constraints) N_VDestroy(m_constraints);

        if (m_sun_ctx) {
            SUNContext_Free(&m_sun_ctx);
            m_sun_ctx = nullptr;
        }

        if (m_T_coeffs) N_VDestroy(m_T_coeffs);
        if (m_rho_coeffs) N_VDestroy(m_rho_coeffs);

        LOG_INFO(m_logger, "SpectralSolverStrategy destroyed successfully");
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Main Evaluation Loop
    /////////////////////////////////////////////////////////////////////////////////

    std::vector<NetOut> SpectralSolverStrategy::evaluate(
        const std::vector<NetIn>& netIns,
        const std::vector<double>& mass_coords,
        const engine::scratch::StateBlob &ctx_template
    ) {
        LOG_INFO(m_logger, "Starting spectral solver evaluation for {} zones", netIns.size());

        assert(std::ranges::all_of(netIns, [&netIns](const NetIn& in) { return in.tMax == netIns[0].tMax; }) && "All NetIn entries must have the same tMax for spectral solver evaluation.");

        std::vector<NetIn> updatedNetIns = netIns;

        std::vector<std::unique_ptr<engine::scratch::StateBlob>> workspaces;
        workspaces.resize(updatedNetIns.size());

        LOG_INFO(m_logger, "Building workspaces...");
        GF_OMP(parallel for,)
        for (size_t shellID = 0; shellID < updatedNetIns.size(); ++shellID) {
            workspaces[shellID] = ctx_template.clone_structure();
        }
        LOG_INFO(m_logger, "Workspaces built successfully.");

        LOG_INFO(m_logger, "Projecting initial conditions onto engine network...");
        GF_OMP(parallel for,)
        for (size_t shellID = 0; shellID < updatedNetIns.size(); ++shellID) {
            updatedNetIns[shellID].composition = m_engine.project(*workspaces[shellID], updatedNetIns[shellID]);
        }
        LOG_INFO(m_logger, "Initial conditions projected successfully.");

        /////////////////////////////////////
        /// Build the species union set   ///
        /////////////////////////////////////

        LOG_INFO(m_logger, "Collecting global species set from all zones...");
        std::set<fourdst::atomic::Species> global_active_species;
        for (const auto &updatedNetIn : updatedNetIns) {
            const auto& local_species = updatedNetIn.composition.getRegisteredSpecies();
            global_active_species.insert(local_species.begin(), local_species.end());
        }
        m_global_species_list.clear();
        m_global_species_list.reserve(global_active_species.size());
        m_global_species_list.insert(
            m_global_species_list.end(),
            global_active_species.begin(),
            global_active_species.end()
        );
        LOG_INFO(m_logger, "Done collecting global species set. Total unique species: {}", m_global_species_list.size());


        /////////////////////////////////////
        /// Evaluate the monitor function ///
        /////////////////////////////////////
        LOG_INFO(m_logger, "Evaluating monitor function...");
        const std::vector<double> monitor_function = evaluate_monitor_function(updatedNetIns);
        LOG_INFO(m_logger, "Monitor function evaluated successfully...");

        /////////////////////////////////////////////
        /// Determine number of quadratude nodes  ///
        /////////////////////////////////////////////
        const size_t num_nodes = nyquist_elements(
            m_config->solver.spectral.basis.num_elements,
            updatedNetIns.size()
        );
        LOG_INFO(m_logger, "Configuration requested {} quadrature nodes, actually using {} based on Nyquist criterion and number of zones.", m_config->solver.spectral.basis.num_elements, num_nodes);

        /////////////////////////////////////
        ///  Generate Quadrature Basis    ///
        /////////////////////////////////////
        LOG_INFO(m_logger, "Generating quadrature basis...");
        m_current_basis = generate_basis_from_monitor(monitor_function, mass_coords, num_nodes);
        LOG_INFO(m_logger, "Quadrature basis generated successfully with {} basis functions and {} quadrature nodes.", m_current_basis.knots.size() - m_current_basis.degree - 1, m_current_basis.quadrature_nodes.size());

        /////////////////////////////////////////
        /// Rebuild workspaces for only basis ///
        /////////////////////////////////////////
        LOG_INFO(m_logger, "Rebuilding workspaces for {} quadrature nodes...", m_current_basis.quadrature_nodes.size());
        workspaces.clear();
        workspaces.resize(m_current_basis.quadrature_nodes.size());
        GF_OMP(parallel for,)
        for (size_t shellID = 0; shellID < workspaces.size(); ++shellID) { // NOLINT(*-loop-convert)
            workspaces[shellID] = ctx_template.clone_structure();
        }

        LOG_INFO(m_logger, "Workspaces rebuilt successfully.");

        ////////////////////////////////////////
        ///  Project Initial Conditions      ///
        ////////////////////////////////////////
        LOG_INFO(m_logger, "Projecting initial conditions onto quadrature nodes...");
        const size_t num_basis_funcs = m_current_basis.knots.size() - m_current_basis.degree - 1;

        std::vector<BasisEval> shell_cache(updatedNetIns.size());

        GF_OMP(parallel for,)
        for (size_t shellID = 0; shellID < shell_cache.size(); ++shellID) {
            auto [start, phi] = evaluate_bspline(mass_coords[shellID], m_current_basis);
            shell_cache[shellID] = {.start_idx=start, .phi=phi};
        }

        const DenseLinearSolver proj_solver(num_basis_funcs, m_sun_ctx);
        proj_solver.init_from_cache(num_basis_funcs, shell_cache);

        if (m_T_coeffs) N_VDestroy(m_T_coeffs);
        m_T_coeffs = N_VNew_Serial(static_cast<sunindextype>(num_basis_funcs), m_sun_ctx);
        project_specific_variable(updatedNetIns, mass_coords, shell_cache, proj_solver, m_T_coeffs, 0, [](const NetIn& s) { return s.temperature; }, true);

        if (m_rho_coeffs) N_VDestroy(m_rho_coeffs);
        m_rho_coeffs = N_VNew_Serial(static_cast<sunindextype>(num_basis_funcs), m_sun_ctx);
        project_specific_variable(updatedNetIns, mass_coords, shell_cache, proj_solver, m_rho_coeffs, 0, [](const NetIn& s) { return s.density; }, true);

        const size_t num_species = m_global_species_list.size();
        size_t current_offset = 0;

        const size_t total_coefficients = num_basis_funcs * (num_species + 1);

        if (m_Y) N_VDestroy(m_Y);
        if (m_constraints) N_VDestroy(m_constraints);

        m_Y = N_VNew_Serial(static_cast<sunindextype>(total_coefficients), m_sun_ctx);
        N_VConst(0.0, m_Y);
        m_constraints = N_VClone(m_Y);
        N_VConst(0.0, m_constraints); // For now no constraints on coefficients

        for (const auto& sp : m_global_species_list) {
            project_specific_variable(
                updatedNetIns,
                mass_coords,
                shell_cache,
                proj_solver,
                m_Y,
                current_offset,
                [&sp](const NetIn& s) {
                    if (!s.composition.contains(sp)) return 0.0;
                    return s.composition.getMolarAbundance(sp);
                },
                false
            );
            current_offset += num_basis_funcs;
        }

        sunrealtype* y_data = N_VGetArrayPointer(m_Y);
        const size_t energy_offset = num_species * num_basis_funcs;

        assert(energy_offset == current_offset && "Energy offset calculation mismatch in spectral solver initialization.");

        for (size_t i = 0; i < num_basis_funcs; ++i) {
            y_data[energy_offset + i] = 0.0;
        }
        LOG_INFO(m_logger, "Done projecting initial conditions onto quadrature nodes.");

        LOG_INFO(m_logger, "Projecting quadrature conditions onto workspaces");
        for (size_t q = 0; q < m_current_basis.quadrature_nodes.size(); ++q) {
            // ReSharper disable once CppUseStructuredBinding
            GridPoint gp = reconstruct_at_quadrature(m_Y, q, m_current_basis);
            NetIn quad_netin;
            quad_netin.temperature = gp.T9 * 1e9;
            quad_netin.density = gp.rho;
            quad_netin.tMax = 1; // Not needed for projection
            quad_netin.composition = gp.composition;

            (void)m_engine.project(*workspaces[q], quad_netin); // We do not need to capture this output just do the projection
        }
        LOG_INFO(m_logger, "Done projecting initial conditions onto workspaces.");


        DenseLinearSolver mass_solver(num_basis_funcs, m_sun_ctx);
        mass_solver.init_from_basis(num_basis_funcs, m_current_basis);

        /////////////////////////////////////
        /// CVODE Initialization          ///
        /////////////////////////////////////

        LOG_INFO(m_logger, "Initializing CVODE resources...");
        CVODEUserData data;
        data.solver_instance = this;
        data.engine = &m_engine;
        data.mass_matrix_solver_instance = &mass_solver;
        data.basis = &m_current_basis;
        data.workspaces = &workspaces;

        const double absTol = m_absTol.value_or(1e-10);
        const double relTol = m_relTol.value_or(1e-6);

        const bool size_changed = m_last_size != total_coefficients;
        m_last_size = total_coefficients;

        if (m_cvode_mem == nullptr || size_changed) {
            if (m_cvode_mem) {
                CVodeFree(&m_cvode_mem);
                m_cvode_mem = nullptr;
            }
            if (m_LS) {
                SUNLinSolFree(m_LS);
                m_LS = nullptr;
            }
            if (m_J) {
                SUNMatDestroy(m_J);
                m_J = nullptr;
            }

            m_cvode_mem = CVodeCreate(CV_BDF, m_sun_ctx);
            utils::check_sundials_flag(m_cvode_mem == nullptr ? -1 : 0, "CVodeCreate", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

            utils::check_sundials_flag(CVodeInit(m_cvode_mem, cvode_rhs_wrapper, 0.0, m_Y), "CVodeInit", utils::SUNDIALS_RET_CODE_TYPES::CVODE);
            m_J = SUNDenseMatrix(static_cast<sunindextype>(total_coefficients), static_cast<sunindextype>(total_coefficients), m_sun_ctx);
            m_LS = SUNLinSol_Dense(m_Y, m_J, m_sun_ctx);
            utils::check_sundials_flag(CVodeSetLinearSolver(m_cvode_mem, m_LS, m_J), "CVodeSetLinearSolver", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

            utils::check_sundials_flag(CVodeSetJacFn(m_cvode_mem, cvode_jac_wrapper), "CVodeSetJacFn", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

        } else {
            utils::check_sundials_flag(CVodeReInit(m_cvode_mem, 0.0, m_Y), "CVodeReInit", utils::SUNDIALS_RET_CODE_TYPES::CVODE);
        }

        utils::check_sundials_flag(CVodeSStolerances(m_cvode_mem, relTol, absTol), "CVodeSStolerances", utils::SUNDIALS_RET_CODE_TYPES::CVODE);
        utils::check_sundials_flag(CVodeSetUserData(m_cvode_mem, &data), "CVodeSetUserData", utils::SUNDIALS_RET_CODE_TYPES::CVODE);
        LOG_INFO(m_logger, "CVODE resources initialized successfully.");

        /////////////////////////////////////
        /// Time Integration Loop         ///
        /////////////////////////////////////

        const double target_time = updatedNetIns[0].tMax;
        double current_time = 0.0;

        LOG_INFO(m_logger, "Beginning time integration to target time {:10.4e}...", target_time);
        while (current_time < target_time) {
            int flag = CVode(m_cvode_mem, target_time, m_Y, &current_time, CV_ONE_STEP);
            utils::check_sundials_flag(flag, "CVode", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

            std::println("Advanced to time: {:10.4e} / {:10.4e}", current_time, target_time);
            LOG_INFO(m_logger, "Advanced to time: {:10.4e} / {:10.4e}", current_time, target_time);
        }
        LOG_INFO(m_logger, "Time integration complete. Reconstructing solution...");

        std::vector<NetOut> results = reconstruct_solution(updatedNetIns, mass_coords, m_Y, m_current_basis, target_time);
        LOG_INFO(m_logger, "Spectral solver evaluation complete for all zones.");
        return results;
    }

    void SpectralSolverStrategy::set_callback(const std::any &callback) {
        m_callback = std::any_cast<TimestepCallback>(callback);
    }

    std::vector<std::tuple<std::string, std::string>> SpectralSolverStrategy::describe_callback_context() const {
        throw std::runtime_error("SpectralSolverStrategy does not yet implement describe_callback_context.");
    }

    bool SpectralSolverStrategy::get_stdout_logging_enabled() const {
        return m_stdout_logging_enabled;
    }

    void SpectralSolverStrategy::set_stdout_logging_enabled(bool logging_enabled) {
        m_stdout_logging_enabled = logging_enabled;
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Static Wrappers for SUNDIALS Callbacks
    ////////////////////////////////////////////////////////////////////////////////

    int SpectralSolverStrategy::cvode_rhs_wrapper(
        const sunrealtype t,
        const N_Vector y_coeffs,
        const N_Vector ydot_coeffs,
        void *user_data
    ) {
        auto *data = static_cast<CVODEUserData*>(user_data);
        const auto *instance = data->solver_instance;

        try {
            return instance -> calculate_rhs(t, y_coeffs, ydot_coeffs, data);
        } catch (const std::exception& e) {
            LOG_CRITICAL(instance->m_logger, "Uncaught exception in Spectral Solver RHS wrapper at time {}: {}", t, e.what());
            return -1;
        } catch (...) {
            LOG_CRITICAL(instance->m_logger, "Unknown uncaught exception in Spectral Solver RHS wrapper at time {}", t);
            return -1;
        }
    }

    int SpectralSolverStrategy::cvode_jac_wrapper(
        const sunrealtype t,
        const N_Vector y,
        const N_Vector ydot,
        const SUNMatrix J,
        void *user_data,
        const N_Vector tmp1,
        const N_Vector tmp2,
        const N_Vector tmp3
    ) {
        const auto *data = static_cast<CVODEUserData*>(user_data);
        const auto *instance = data->solver_instance;

        try {
            return instance->calculate_jacobian(t, y, ydot, J, data, tmp1, tmp2, tmp3);
        } catch (const std::exception& e) {
            LOG_CRITICAL(instance->m_logger, "Uncaught exception in Spectral Solver Jacobian wrapper at time {}: {}", t, e.what());
            return -1;
        } catch (...) {
            LOG_CRITICAL(instance->m_logger, "Unknown uncaught exception in Spectral Solver Jacobian wrapper at time {}", t);
            return -1;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// RHS implementation
    ////////////////////////////////////////////////////////////////////////////////

    // ReSharper disable once CppDFAUnreachableFunctionCall
    int SpectralSolverStrategy::calculate_rhs(
        sunrealtype t,
        N_Vector y_coeffs,
        N_Vector ydot_coeffs,
        CVODEUserData* data
    ) const {
        const auto& basis = m_current_basis;
        DenseLinearSolver* mass_solver = data->mass_matrix_solver_instance;

        const size_t num_basis_funcs = basis.knots.size() - basis.degree - 1;
        const size_t num_species = m_global_species_list.size();

        const size_t total_dofs = num_basis_funcs * (num_species + 1);

        sunrealtype* global_rhs_data = N_VGetArrayPointer(ydot_coeffs);
        N_VConst(0.0, ydot_coeffs);

        std::atomic<bool> failure_flag{false};

        GF_OMP(parallel,) {
            std::vector<double> local_rhs(total_dofs, 0.0);

            GF_OMP(for nowait,)
            for (size_t q = 0; q < basis.quadrature_nodes.size(); ++q) {
                if (failure_flag.load(std::memory_order_relaxed)) continue;

                double wq = basis.quadrature_weights[q];
                const auto& [start_idx, phi] = basis.quad_evals[q];

                GridPoint gp = reconstruct_at_quadrature(y_coeffs, q, basis);
                LOG_TRACE_L2(m_logger, "RHS Evaluation at time {}: Quad Node {}, T9 = {:10.4e}, rho = {:10.4e}", t, q, gp.T9, gp.rho);
                auto results = m_engine.calculateRHSAndEnergy(
                    *data->workspaces->at(q),
                    gp.composition,
                    gp.T9,
                    gp.rho,
                    false
                );
                if (!results) {
                    failure_flag.store(true);
                    LOG_CRITICAL(m_logger, "Engine failed to calculate RHS at time {}: {}", t, EngineStatus_to_string(results.error()));
                    GF_OMP(critical, {
                        throw std::runtime_error("Engine failure during RHS calculation.");
                    }) {
                        std::fill_n(global_rhs_data, total_dofs, -1.0);
                    }
                    continue;
                }

                const auto& [dydt, eps_nuc, contributions, nu_loss, nu_flux] = results.value();
                for (size_t s = 0; s < num_species; ++s) {
                    const auto& sp = m_global_species_list[s];
                    double rate = 0.0;
                    if (dydt.contains(sp)) {
                        rate = dydt.at(sp);
                    }

                    size_t species_offset = s * num_basis_funcs;

                    for (size_t k  = 0; k < phi.size(); ++k) {
                        size_t global_idx = species_offset + start_idx + k;
                        local_rhs[global_idx] += wq * phi[k] * rate;
                    }
                }

                size_t energy_offset = num_species * num_basis_funcs;

                for (size_t k = 0; k < phi.size(); ++k) {
                    size_t global_idx = energy_offset + start_idx + k;
                    local_rhs[global_idx] += eps_nuc * wq * phi[k];
                }
            }

            GF_OMP(critical,) {
                if (!failure_flag.load(std::memory_order_relaxed)) {
                    for (size_t i = 0; i < total_dofs; ++i) {
                        global_rhs_data[i] += local_rhs[i];
                    }
                }
            }
        }

        if (failure_flag.load()) {
            return -1;
        }

        size_t total_vars = num_species + 1;
        mass_solver -> solve_inplace(ydot_coeffs, total_vars, num_basis_funcs);
        return 0;
    }

    int SpectralSolverStrategy::calculate_jacobian(
        sunrealtype t,
        N_Vector y_coeffs,
        N_Vector ydot_coeffs,
        SUNMatrix J,
        const CVODEUserData *data,
        N_Vector tmp1,
        N_Vector tmp2,
        N_Vector tmp3
    ) const {
        const auto& basis = m_current_basis;
        DenseLinearSolver* mass_solver = data->mass_matrix_solver_instance;

        const size_t num_basis_funcs = basis.knots.size() - basis.degree - 1;
        const size_t num_species = m_global_species_list.size();
        const size_t total_dofs = num_basis_funcs * (num_species + 1);

        SUNMatZero(J);
        sunrealtype* J_global_data = SUNDenseMatrix_Data(J);

        std::atomic<bool> failure_flag(false);

        GF_OMP(parallel,) {
#if defined(GF_USE_OPENMP)
            int thread_id = omp_get_thread_num();
#else
            int thread_id = 0;
#endif
            GF_OMP(for schedule(static) nowait,)
            for (size_t q = 0; q < basis.quadrature_nodes.size(); ++q) {
                if (failure_flag.load(std::memory_order_relaxed)) continue;

                double wq = basis.quadrature_weights[q];
                const auto& [start_idx, phi] = basis.quad_evals[q];

                GridPoint gp = reconstruct_at_quadrature(y_coeffs, q, basis);
                try {
                    engine::NetworkJacobian jac = m_engine.generateJacobianMatrix(
                        *data->workspaces->at(q),
                        gp.composition,
                        gp.T9,
                        gp.rho
                    );

                    auto accumulate_term = [&](size_t row_offset, size_t col_offset, double value) {
                        if (std::abs(value) < 1e-100) return; // regularization

                        const double w_val = wq * value;

                        for (size_t k_col = 0; k_col < phi.size(); ++k_col) {
                            const size_t global_col = col_offset + start_idx + k_col;
                            const size_t col_data_idx = global_col * total_dofs;

                            for (size_t k_row = 0; k_row < phi.size(); ++k_row) {
                                const size_t global_row = row_offset + start_idx + k_row;
                                const double contribution = w_val * phi[k_row] * phi[k_col];

                                GF_OMP(atomic update,)
                                J_global_data[global_row + col_data_idx] += contribution;
                            }
                        }
                    };

                    for (size_t s_row = 0; s_row < num_species; ++s_row) {
                        const auto& sp_row = m_global_species_list[s_row];

                        if (!gp.composition.contains(sp_row)) continue;

                        for (size_t s_col = 0; s_col < num_species; ++s_col) {
                            const auto& sp_col = m_global_species_list[s_col];

                            if (!gp.composition.contains(sp_col)) continue;

                            double val = jac(sp_row, sp_col);

                            accumulate_term(s_row * num_basis_funcs, s_col * num_basis_funcs, val);
                        }
                    }
                } catch (const exceptions::GridFireError &e) {
                    failure_flag.store(true);
                    std::string error_msg = std::format("Engine failed to calculate Jacobian, due to known internal GridFire error, at time {}: {}", t, e.what());
                    LOG_CRITICAL(m_logger, "{}", error_msg);
                    GF_OMP(critical, {
                        throw std::runtime_error(error_msg);
                    }) {
                        SUNMatZero(J);
                    }
                    continue;
                } catch (const std::exception &e) {
                    failure_flag.store(true);
                    std::string error_msg = std::format("Engine failed to calculate Jacobian, due to some known yet non GridFire error, at time {}: {}", t, e.what());
                    LOG_CRITICAL(m_logger, "{}", error_msg);
                    GF_OMP(critical, {
                        throw std::runtime_error(error_msg);
                    }) {
                        SUNMatZero(J);
                    }
                    continue;
                } catch (...) {
                    failure_flag.store(true);
                    std::string error_msg = std::format("Engine failed to calculate Jacobian, due to unknown error, at time {}.", t);
                    LOG_CRITICAL(m_logger, "{}", error_msg);
                    GF_OMP(critical, {
                        throw std::runtime_error(error_msg);
                    }) {
                        SUNMatZero(J);
                    }
                    continue;
                }
            }
        }

        if (failure_flag.load()) { return -1; }
        inspect_jacobian(J, "Physics Assembly (Pre-Mass-Matrix)");

        GF_OMP(parallel for,)
        for (size_t col = 0; col < total_dofs; ++col) {
            sunrealtype* col_ptr = J_global_data + (col * total_dofs);
            mass_solver->solve_inplace_ptr(col_ptr, num_species + 1, num_basis_funcs);
        }
        inspect_jacobian(J, "Final Jacobian (Post-Mass-Matrix)");
        return 0;
    }

    size_t SpectralSolverStrategy::nyquist_elements(
        const size_t requested_elements,
        const size_t num_shells
    ) {
        size_t max_allowed_elements = std::max(size_t(1), num_shells/2);
        size_t actual_elements = std::min(requested_elements, max_allowed_elements);
        if (num_shells <= 5) {
            actual_elements = 1;
        }
        return actual_elements;
    }


    ////////////////////////////////////////////////////////////////////////////////
    /// Spectral Utilities
    ///     These include basis generation, monitor function evaluation
    ///     projection and reconstruction routines.
    ////////////////////////////////////////////////////////////////////////////////

    std::vector<double> SpectralSolverStrategy::evaluate_monitor_function(
        const std::vector<NetIn>& current_shells
    ) const {
        assert(!m_global_species_list.empty() && "Global species list must be initialized before evaluating monitor function.");
        assert(!current_shells.empty() && "Current shells list must not be empty when evaluating monitor function.");

        const size_t n_shells = current_shells.size();
        if (n_shells < 3) {
            return std::vector<double>(n_shells, 1.0); // NOLINT(*-return-braced-init-list)
        }

        std::vector<double> M(n_shells, 1.0);
        std::vector<double> data(n_shells);

        auto accumulate_variable = [&](auto getter, double weight, bool use_log) {
            double min_val = std::numeric_limits<double>::max();
            double max_val = std::numeric_limits<double>::lowest();

            for (size_t i = 0 ; i < n_shells; ++i) {
                double val = getter(current_shells[i]);
                if (use_log) {
                    val = std::log10(std::max(val, 1e-100));
                }

                data[i] = val;

                if (val < min_val) min_val = val;
                if (val > max_val) max_val = val;
            }

            const double scale = max_val - min_val;
            if (scale < 1e-10) return;

            for (size_t i = 1; i < n_shells - 1; ++i) {
                const double v_prev = data[i-1];
                const double v_curr = data[i];
                const double v_next = data[i+1];

                // Finite difference estimates for first and second derivatives
                double d1 = std::abs(v_next - v_prev) / 2.0;
                double d2 = std::abs(v_next - 2.0 * v_curr + v_prev);

                d1 /= scale;
                d2 /= scale;

                const double alpha = m_config->solver.spectral.monitorFunction.alpha;
                const double beta = m_config->solver.spectral.monitorFunction.beta;

                M[i] += weight * (alpha * d1 + beta * d2);
            }
        };

        auto safe_get_abundance = [](const NetIn& netIn, const fourdst::atomic::Species& sp) -> double {
            if (netIn.composition.contains(sp)) {
                return netIn.composition.getMolarAbundance(sp);
            }
            return 0.0;
        };

        const double structure_weight = m_config->solver.spectral.monitorFunction.structure_weight;
        double abundance_weight = m_config->solver.spectral.monitorFunction.abundance_weight;
        accumulate_variable([](const NetIn& s) { return s.temperature; }, structure_weight, true);
        accumulate_variable([](const NetIn& s) { return s.density; }, structure_weight, true);

        for (const auto& sp : m_global_species_list) {
            accumulate_variable([&sp, &safe_get_abundance](const NetIn& s) { return safe_get_abundance(s, sp); }, abundance_weight, false);
        }

        //////////////////////////////
        /// Smoothing the Monitor  ///
        //////////////////////////////

        std::vector<double> M_smooth = M;
        for (size_t i = 1; i < n_shells - 1; ++i) {
            M_smooth[i] = (M[i-1] + 2.0 * M[i] + M[i+1]) / 4.0;
        }

        M_smooth[0] = M_smooth[1];
        M_smooth[n_shells-1] = M_smooth[n_shells-2];

        return M_smooth;

    }

    SpectralSolverStrategy::SplineBasis SpectralSolverStrategy::generate_basis_from_monitor(
        const std::vector<double>& monitor_values,
        const std::vector<double>& mass_coordinates,
        const size_t actual_elements
    ) const {
        SplineBasis basis;
        basis.degree = 3; // Cubic Spline

        const size_t n_shells = monitor_values.size();

        std::vector<double> I(n_shells, 0.0);
        double current_integral = 0.0;

        for (size_t i = 1; i < n_shells; ++i) {
            const double dx = mass_coordinates[i] - mass_coordinates[i-1];
            double dI = 0.5 * (monitor_values[i] + monitor_values[i-1]) * dx;

            dI = std::max(dI, 1e-30);
            current_integral += dI;
            I[i] = current_integral;
        }

        const double total_integral = I.back();
        for (size_t i = 0; i < n_shells; ++i) {
            I[i] /= total_integral;
        }

        basis.knots.reserve(actual_elements + 1 + 2 * basis.degree);

        // Note that these imply that mass_coordinates must be sorted in increasing order
        double min_mass = mass_coordinates.front();
        double max_mass = mass_coordinates.back();

        for (int i = 0; i < basis.degree; ++i) {
            basis.knots.push_back(min_mass);
        }

        for (size_t k = 1; k < actual_elements; ++k) {
            double target_I = static_cast<double>(k) / static_cast<double>(actual_elements);

            auto it = std::ranges::lower_bound(I, target_I);
            size_t idx = std::distance(I.begin(), it);

            if (idx == 0) idx = 1;
            if (idx >= n_shells) idx = n_shells - 1;

            double I0 = I[idx-1];
            double I1 = I[idx];
            double m0 = mass_coordinates[idx-1];
            double m1 = mass_coordinates[idx];

            double fraction = (target_I - I0) / (I1 - I0);
            double knot_location = m0 + fraction * (m1 - m0);

            basis.knots.push_back(knot_location);
        }

        for (int i = 0; i < basis.degree; ++i) {
            basis.knots.push_back(max_mass);
        }

        constexpr double sqrt_3_over_5 = 0.77459666924;
        constexpr double five_over_nine = 5.0 / 9.0;
        constexpr double eight_over_nine = 8.0 / 9.0;
        static constexpr std::array<double, 3> gl_nodes = {-sqrt_3_over_5, 0.0, sqrt_3_over_5};
        static constexpr std::array<double, 3> gl_weights = {five_over_nine, eight_over_nine, five_over_nine};

        basis.quadrature_nodes.clear();
        basis.quadrature_weights.clear();

        for (size_t i = basis.degree; i < basis.knots.size() - basis.degree - 1; ++i) {
            double a = basis.knots[i];
            double b = basis.knots[i+1];

            if ( b - a < 1e-14) continue;

            double mid = 0.5 * (a + b);
            double half_width = 0.5 * (b - a);

            for (size_t j = 0; j < gl_nodes.size(); ++j) {
                double phys_node = mid + gl_nodes[j] * half_width;
                double phys_weight = gl_weights[j] * half_width;

                basis.quadrature_nodes.push_back(phys_node);
                basis.quadrature_weights.push_back(phys_weight);

                auto [start, phi] = evaluate_bspline(phys_node, basis);
                basis.quad_evals.push_back({start, phi});
            }
        }

        return basis;
    }


    SpectralSolverStrategy::GridPoint SpectralSolverStrategy::reconstruct_at_quadrature(
        const N_Vector y_coeffs,
        const size_t quad_index,
        const SplineBasis &basis
    ) const {
        auto [start_idx, vals] = basis.quad_evals[quad_index];

        const sunrealtype* T_ptr = N_VGetArrayPointer(m_T_coeffs);
        const sunrealtype* rho_ptr = N_VGetArrayPointer(m_rho_coeffs);
        const sunrealtype* y_data = N_VGetArrayPointer(y_coeffs);

        const size_t num_basis_funcs = basis.knots.size() - basis.degree - 1;
        const size_t num_species = m_global_species_list.size();

        double logT = 0.0;
        double logRho = 0.0;

        for (size_t k = 0; k < vals.size(); ++k) {
            size_t idx = start_idx + k;
            logT += T_ptr[idx] * vals[k];
            logRho += rho_ptr[idx] * vals[k];
        }

        GridPoint result;
        result.T9 = std::pow(10.0, logT) / 1e9;
        result.rho = std::pow(10.0, logRho);

        for (size_t s = 0; s < num_species; ++s) {
            const auto& species = m_global_species_list[s];

            double abundance = 0.0;
            const size_t offset = s * num_basis_funcs;

            for (size_t k = 0; k < vals.size(); ++k) {
                abundance += y_data[offset + start_idx + k] * vals[k];
            }

            if (abundance < 0.0) abundance = 0.0;

            result.composition.registerSpecies(species);
            result.composition.setMolarAbundance(species, abundance);
        }

        return result;
    }

    std::vector<NetOut> SpectralSolverStrategy::reconstruct_solution(
        const std::vector<NetIn>& original_inputs,
        const std::vector<double>& mass_coordinates,
        const N_Vector final_coeffs,
        const SplineBasis& basis,
        const double dt
    ) const {
        const size_t n_shells = original_inputs.size();
        const size_t num_basis_funcs = basis.knots.size() - basis.degree - 1;
        const size_t num_species = m_global_species_list.size();

        // Pre-allocate output vector so threads can write directly
        std::vector<NetOut> outputs(n_shells);

        const sunrealtype* c_data = N_VGetArrayPointer(final_coeffs);

        GF_OMP(parallel for,)
        for (size_t shellID = 0; shellID < n_shells; ++shellID) {
            const double x = mass_coordinates[shellID];

            auto [start_idx, vals] = evaluate_bspline(x, basis);

            auto reconstruct_var = [&](const size_t coeff_offset) -> double {
                double result = 0.0;
                for (size_t i = 0; i < vals.size(); ++i) {
                    result += c_data[coeff_offset + start_idx + i] * vals[i];
                }
                return result;
            };

            fourdst::composition::Composition comp_new;

            for (size_t s_idx = 0; s_idx < num_species; ++s_idx) {
                const auto& sp = m_global_species_list[s_idx];

                const size_t current_offset = s_idx * num_basis_funcs;
                double Y_val = reconstruct_var(current_offset);

                if (Y_val < 0.0 && Y_val >= -1e-16) {
                     Y_val = 0.0;
                }

                comp_new.registerSpecies(sp);
                comp_new.setMolarAbundance(sp, Y_val);
            }

            const double energy = reconstruct_var(num_species * num_basis_funcs);

            NetOut netOut;
            netOut.composition = comp_new;
            netOut.energy = energy;
            netOut.num_steps = -1; // Not tracked in spectral solver

            outputs[shellID] = std::move(netOut);
        }

        return outputs;
    }


    void SpectralSolverStrategy::project_specific_variable(
        const std::vector<NetIn> &current_shells,
        const std::vector<double> &mass_coordinates,
        const std::vector<BasisEval> &shell_cache,
        const DenseLinearSolver &linear_solver,
        N_Vector output_vec,
        size_t output_offset,
        const std::function<double(const NetIn &)> &getter,
        bool use_log
    ) {
        const size_t n_shells = current_shells.size();

        sunrealtype* out_ptr = N_VGetArrayPointer(output_vec);
        size_t basis_size = N_VGetLength(linear_solver.temp_vector);

        for (size_t i = 0; i < basis_size; ++i ) {
            out_ptr[output_offset + i] = 0.0;
        }

        for (size_t shellID = 0; shellID < n_shells; ++shellID) {
            double val = getter(current_shells[shellID]);
            if (use_log) val = std::log10(std::max(val, 1e-100));

            const auto& eval = shell_cache[shellID];

            for (size_t i = 0; i < eval.phi.size(); ++i) {
                out_ptr[output_offset + eval.start_idx + i] += val * eval.phi[i];
            }
        }

        sunrealtype* tmp_data = N_VGetArrayPointer(linear_solver.temp_vector);
        for (size_t i = 0; i < basis_size; ++i) tmp_data[i] = out_ptr[output_offset + i];

        utils::check_sundials_flag(SUNLinSolSolve(linear_solver.LS, linear_solver.A, linear_solver.temp_vector, linear_solver.temp_vector, 0.0), "SUNLinSolSolve - Projection Solver", utils::SUNDIALS_RET_CODE_TYPES::CVODE);

        for (size_t i = 0; i < basis_size; ++i) out_ptr[output_offset + i] = tmp_data[i];
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// Debugging Utilities
    ///////////////////////////////////////////////////////////////////////////////////
    void SpectralSolverStrategy::inspect_jacobian(SUNMatrix J, const std::string &context) const {
        sunrealtype* data = SUNDenseMatrix_Data(J);
        sunindextype rows = SUNDenseMatrix_Rows(J);
        sunindextype cols = SUNDenseMatrix_Columns(J);

        // --- 1. Gather Statistics ---
        size_t nan_count = 0;
        size_t inf_count = 0;
        size_t zero_diag_count = 0;

        double max_val = 0.0;
        double min_val = 0.0;
        double min_diag_val = std::numeric_limits<double>::max();
        double max_diag_val = std::numeric_limits<double>::lowest();

        bool matrix_is_empty = (rows == 0 || cols == 0);
        bool non_zero_elements_found = false;

        // Iterate Column-Major (standard for SUNDIALS)
        for (sunindextype j = 0; j < cols; ++j) {
            // Diagonal Check
            if (j < rows) {
                double diag = data[j * rows + j];
                double abs_diag = std::abs(diag);

                if (abs_diag < 1e-100) zero_diag_count++;

                if (abs_diag < min_diag_val) min_diag_val = abs_diag;
                if (abs_diag > max_diag_val) max_diag_val = abs_diag;
            }

            for (sunindextype i = 0; i < rows; ++i) {
                double val = data[j * rows + i];

                if (std::isnan(val)) {
                    nan_count++;
                } else if (std::isinf(val)) {
                    inf_count++;
                } else {
                    double abs_val = std::abs(val);
                    if (abs_val > 0.0) {
                        if (!non_zero_elements_found) {
                            // First non-zero element initializes the range
                            max_val = abs_val;
                            min_val = abs_val;
                            non_zero_elements_found = true;
                        } else {
                            if (abs_val > max_val) max_val = abs_val;
                            if (abs_val < min_val) min_val = abs_val;
                        }
                    }
                }
            }
        }

        if (!non_zero_elements_found) {
            min_diag_val = 0.0;
            max_diag_val = 0.0;
        }

        // --- 2. Build Data Vectors ---
        std::vector<std::string> metrics = {
            "Dimensions",
            "NaN Count",
            "Inf Count",
            "Zero Diagonals",
            "Global Range (abs)",
            "Diagonal Range (abs)",
            "Status"
        };

        std::vector<std::string> values;
        values.reserve(metrics.size());

        // Dimensions
        values.push_back(std::format("{} x {}", rows, cols));

        // Errors
        values.push_back(std::to_string(nan_count));
        values.push_back(std::to_string(inf_count));
        values.push_back(std::to_string(zero_diag_count));

        // Ranges
        values.push_back(std::format("[{:.2e}, {:.2e}]", min_val, max_val));
        values.push_back(std::format("[{:.2e}, {:.2e}]", min_diag_val, max_diag_val));

        // Status
        bool failed = (nan_count > 0 || inf_count > 0 || zero_diag_count > 0 || matrix_is_empty);
        values.push_back(failed ? "FAIL" : "OK");

        // --- 3. Construct Columns Manually ---
        std::vector<std::unique_ptr<gridfire::utils::ColumnBase>> columns;

        columns.push_back(std::make_unique<gridfire::utils::Column<std::string>>("Metric", metrics));
        columns.push_back(std::make_unique<gridfire::utils::Column<std::string>>("Value", values));

        // --- 4. Format and Log ---
        std::string table_name = std::format("Jacobian Inspection: {}", context);
        std::string report_str = gridfire::utils::format_table(table_name, columns);

        if (failed) {
            std::println("{}", report_str);
            LOG_CRITICAL(m_logger, "\n{}", report_str);
        } else {
            std::println("{}", report_str);
            LOG_INFO(m_logger, "\n{}", report_str);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// SpectralSolverStrategy::MassMatrixSolver Implementation
    ///////////////////////////////////////////////////////////////////////////////

    SpectralSolverStrategy::DenseLinearSolver::DenseLinearSolver(
        size_t size,
        SUNContext sun_ctx
    ) : ctx(sun_ctx) {
        A = SUNDenseMatrix(size, size, sun_ctx);
        temp_vector = N_VNew_Serial(size, sun_ctx);

        LS = SUNLinSol_Dense(temp_vector, A, sun_ctx);

        if (!A || !temp_vector || !LS) {
            throw std::runtime_error("Failed to create MassMatrixSolver components.");
        }

        zero();
    }

    SpectralSolverStrategy::DenseLinearSolver::~DenseLinearSolver() {
        if (LS) SUNLinSolFree(LS);
        if (A) SUNMatDestroy(A);
        if (temp_vector) N_VDestroy(temp_vector);
    }

    void SpectralSolverStrategy::DenseLinearSolver::zero() const {
        SUNMatZero(A);
    }

    void SpectralSolverStrategy::DenseLinearSolver::init_from_cache(
        const size_t num_basis_funcs,
        const std::vector<BasisEval> &shell_cache
    ) const {
        sunrealtype* a_data = SUNDenseMatrix_Data(A);
        for (const auto&[start_idx, phi] : shell_cache) {
            for (size_t i = 0; i < phi.size(); ++i) {
                const size_t row = start_idx + i;

                for (size_t j = 0; j < phi.size(); ++j) {
                    const size_t col = start_idx + j;

                    a_data[col * num_basis_funcs + row] += phi[i] * phi[j];
                }
            }
        }

        // Apply a small diagonal perturbation for numerical stability
        for (size_t i = 0; i < num_basis_funcs; ++i) {
            constexpr double epsilon = 1e-14;
            a_data[i * num_basis_funcs + i] += epsilon;
        }

        setup();
    }

    void SpectralSolverStrategy::DenseLinearSolver::init_from_basis(
        const size_t num_basis_funcs,
        const SplineBasis &basis
    ) const {
        sunrealtype* m_data = SUNDenseMatrix_Data(A);
        for (size_t q = 0; q < basis.quadrature_nodes.size(); ++q) {
            double w_q = basis.quadrature_weights[q];
            const auto& eval = basis.quad_evals[q];

            for (size_t i = 0; i < eval.phi.size(); ++i) {
                size_t row = eval.start_idx + i;
                for (size_t j = 0; j < eval.phi.size(); ++j) {
                    size_t col = eval.start_idx + j;
                    m_data[col * num_basis_funcs + row] += w_q * eval.phi[j] * eval.phi[i];
                }
            }
        }
        setup();
    }

    void SpectralSolverStrategy::DenseLinearSolver::setup() const {
        utils::check_sundials_flag(SUNLinSolSetup(LS, A), "SUNLinSolSetup - Mass Matrix Solver", utils::SUNDIALS_RET_CODE_TYPES::CVODE);
    }

    // ReSharper disable once CppMemberFunctionMayBeConst
    void SpectralSolverStrategy::DenseLinearSolver::solve_inplace(const N_Vector x, const size_t num_vars, const size_t basis_size) const {
        sunrealtype* x_data = N_VGetArrayPointer(x);
        sunrealtype* tmp_data = N_VGetArrayPointer(temp_vector);

        for (size_t v = 0; v < num_vars; ++v) {
            const size_t offset = v * basis_size;

            for (size_t i = 0; i < basis_size; ++i) {
                tmp_data[i] = x_data[offset + i];
            }
            SUNLinSolSolve(LS, A, temp_vector, temp_vector, 0.0);

            for (size_t i = 0; i < basis_size; ++i) {
                x_data[offset + i] = tmp_data[i];
            }
        }
    }

    void SpectralSolverStrategy::DenseLinearSolver::solve_inplace_ptr(
        sunrealtype *data_ptr,
        const size_t num_vars,
        const size_t basis_size
    ) const {
        sunrealtype* tmp_data = N_VGetArrayPointer(temp_vector);
        for (size_t v = 0; v < num_vars; ++v) {
            const size_t offset = v * basis_size;
            sunrealtype* var_segment = data_ptr + offset;

            for (size_t i = 0; i < basis_size; ++i) {
                tmp_data[i] = var_segment[i];
            }

            const int flag = SUNLinSolSolve(LS, A, temp_vector, temp_vector, 0.0);

            if (flag != 0) {
                GF_OMP(critical,) {
                    std::cerr << "Mass matrix inversion failed in SpectralSolverStrategy::DenseLinearSolver::solve_inplace_ptr. This is a critical error which cannot, at present, be recovered from. If you see this message discard any results reported after it and please report this to the GridFire developers." << flag << std::endl;
                    std::abort();
                }
                // TODO: We cannot trivially throw here if this function fails since it may be run in parallel. For now we trust init_from_cache to not
                //       generate a singular matrix; however, in the future we may want to implement a more robust error handling strategy.
            }

            for (size_t i = 0; i < basis_size; ++i) {
                var_segment[i] = tmp_data[i];
            }
        }
    }
}
