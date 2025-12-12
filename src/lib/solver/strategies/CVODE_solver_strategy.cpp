#include "gridfire/solver/strategies/CVODE_solver_strategy.h"

#include "gridfire/types/types.h"
#include "gridfire/utils/table_format.h"
#include "gridfire/engine/diagnostics/dynamic_engine_diagnostics.h"

#include "quill/LogMacros.h"

#include "fourdst/composition/composition.h"

// ReSharper disable once CppUnusedIncludeDirective
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

#include "fourdst/atomic/species.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/solver/strategies/triggers/engine_partitioning_trigger.h"
#include "gridfire/trigger/procedures/trigger_pprint.h"
#include "gridfire/exceptions/error_solver.h"
#include "gridfire/utils/sundials.h"


namespace gridfire::solver {
     using namespace gridfire::engine;

    CVODESolverStrategy::TimestepContext::TimestepContext(
        const double t,
        const N_Vector &state,
        const double dt,
        const double last_step_time,
        const double t9,
        const double rho,
        const size_t num_steps,
        const DynamicEngine &engine,
        const std::vector<fourdst::atomic::Species> &networkSpecies,
        const size_t currentConvergenceFailure,
        const size_t currentNonlinearIterations,
        const std::map<fourdst::atomic::Species, std::unordered_map<std::string, double>> &reactionContributionMap,
        scratch::StateBlob& ctx
    ) :
    t(t),
    state(state),
    dt(dt),
    last_step_time(last_step_time),
    T9(t9),
    rho(rho),
    num_steps(num_steps),
    engine(engine),
    networkSpecies(networkSpecies),
    currentConvergenceFailures(currentConvergenceFailure),
    currentNonlinearIterations(currentNonlinearIterations),
    reactionContributionMap(reactionContributionMap),
    state_ctx(ctx)
    {}

    std::vector<std::tuple<std::string, std::string>> CVODESolverStrategy::TimestepContext::describe() const {
        std::vector<std::tuple<std::string, std::string>> description;
        description.emplace_back("t", "Current Time");
        description.emplace_back("state", "Current State Vector (N_Vector)");
        description.emplace_back("dt", "Last Timestep Size");
        description.emplace_back("last_step_time", "Time at Last Step");
        description.emplace_back("T9", "Temperature in GK");
        description.emplace_back("rho", "Density in g/cm^3");
        description.emplace_back("num_steps", "Number of Steps Taken So Far");
        description.emplace_back("engine", "Reference to the DynamicEngine");
        description.emplace_back("networkSpecies", "Reference to the list of network species");
        description.emplace_back("currentConvergenceFailures", "Number of convergence failures for the current step");
        description.emplace_back("currentNonLinearIterations", "Number of nonlinear iterations for the current step");
        return description;
    }


    CVODESolverStrategy::CVODESolverStrategy(
        const DynamicEngine &engine,
        const scratch::StateBlob& ctx
    ): SingleZoneNetworkSolver<DynamicEngine>(engine, ctx) {
        // PERF: In order to support MPI this function must be changed
        const int flag = SUNContext_Create(SUN_COMM_NULL, &m_sun_ctx);
        if (flag < 0) {
            throw std::runtime_error("Failed to create SUNDIALS context (SUNDIALS Errno: " + std::to_string(flag) + ")");
        }
    }

    CVODESolverStrategy::~CVODESolverStrategy() {
        LOG_TRACE_L1(m_logger, "Cleaning up CVODE resources...");
        cleanup_cvode_resources(true);

        if (m_sun_ctx) {
            SUNContext_Free(&m_sun_ctx);
        }
    }

    NetOut CVODESolverStrategy::evaluate(const NetIn& netIn) {
        return evaluate(netIn, false);
    }

    NetOut CVODESolverStrategy::evaluate(
        const NetIn &netIn,
        bool displayTrigger,
        bool forceReinitialize
    ) {
        LOG_TRACE_L1(m_logger, "Starting solver evaluation with T9: {} and rho: {}", netIn.temperature/1e9, netIn.density);
        LOG_TRACE_L1(m_logger, "Building engine update trigger....");
        auto trigger = trigger::solver::CVODE::makeEnginePartitioningTrigger(1e12, 1e10, 0.5, 2);
        LOG_TRACE_L1(m_logger, "Engine update trigger built!");


        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)

        // The tolerance selection algorithm is that
        //  1. Default tolerances are taken from the config file
        //  2. If the user has set tolerances in code, those override the config
        //  3. If the user has not set tolerances in code and the config does not have them, use hardcoded defaults

        auto absTol = m_config->solver.cvode.absTol;
        auto relTol = m_config->solver.cvode.relTol;

        if (m_absTol) {
            absTol = *m_absTol;
        }
        if (m_relTol) {
            relTol = *m_relTol;
        }

        bool resourcesExist = (m_cvode_mem != nullptr) && (m_Y != nullptr);

        bool inconsistentComposition = netIn.composition.hash() != m_last_composition_hash;
        fourdst::composition::Composition equilibratedComposition;

        if (forceReinitialize || !resourcesExist || inconsistentComposition) {
            cleanup_cvode_resources(true);
            LOG_INFO(
                m_logger,
                "Preforming full CVODE initialization (Reason: {})",
                forceReinitialize ? "Forced reinitialization" :
                        (!resourcesExist ? "CVODE resources do not exist" :
                            "Input composition inconsistent with previous state"));
            LOG_TRACE_L1(m_logger, "Starting engine update chain...");
            equilibratedComposition = m_engine.project(*m_scratch_blob, netIn);
            LOG_TRACE_L1(m_logger, "Engine updated and equilibrated composition found!");

            size_t numSpecies = m_engine.getNetworkSpecies(*m_scratch_blob).size();
            uint64_t N = numSpecies + 1;

            LOG_TRACE_L1(m_logger, "Number of species: {} ({} independent variables)", numSpecies, N);
            LOG_TRACE_L1(m_logger, "Initializing CVODE resources");
            m_cvode_mem = CVodeCreate(CV_BDF, m_sun_ctx);
            utils::check_cvode_flag(m_cvode_mem == nullptr ? -1 : 0, "CVodeCreate");

            initialize_cvode_integration_resources(N, numSpecies, 0.0, equilibratedComposition, absTol, relTol, 0.0);
            m_last_size = N;
        } else {
            LOG_INFO(m_logger, "Reusing existing CVODE resources (size: {})", m_last_size);

            const size_t numSpecies = m_engine.getNetworkSpecies(*m_scratch_blob).size();
            sunrealtype *y_data = N_VGetArrayPointer(m_Y);
            for (size_t i = 0; i < numSpecies; i++) {
                const auto& species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
                if (netIn.composition.contains(species)) {
                    y_data[i] = netIn.composition.getMolarAbundance(species);
                } else {
                    y_data[i] = std::numeric_limits<double>::min();
                }
            }
            y_data[numSpecies] = 0.0; // Reset energy accumulator
            utils::check_cvode_flag(CVodeSStolerances(m_cvode_mem, relTol, absTol), "CVodeSStolerances");
            utils::check_cvode_flag(CVodeReInit(m_cvode_mem, 0.0, m_Y), "CVodeReInit");

            equilibratedComposition = netIn.composition; // Use the provided composition as-is if we already have validated CVODE resources and that the composition is consistent with the previous state
        }

        size_t numSpecies = m_engine.getNetworkSpecies(*m_scratch_blob).size();
        CVODEUserData user_data {
            .solver_instance = this,
            .ctx = *m_scratch_blob,
            .engine = &m_engine,
        };
        LOG_TRACE_L1(m_logger, "CVODE resources successfully initialized!");

        double current_time = 0;
        // ReSharper disable once CppTooWideScope
        [[maybe_unused]] double last_callback_time = 0;
        m_num_steps = 0;
        double accumulated_energy = 0.0;

        double accumulated_neutrino_energy_loss = 0.0;
        double accumulated_total_neutrino_flux = 0.0;

        size_t total_convergence_failures = 0;
        size_t total_nonlinear_iterations = 0;
        size_t total_update_stages_triggered = 0;

        size_t prev_nonlinear_iterations = 0;
        size_t prev_convergence_failures = 0;

        size_t total_steps = 0;

        LOG_TRACE_L1(m_logger, "Starting CVODE iteration");
        fourdst::composition::Composition postStep = equilibratedComposition;
        while (current_time < netIn.tMax) {
            user_data.T9 = T9;
            user_data.rho = netIn.density;
            user_data.networkSpecies = &m_engine.getNetworkSpecies(*m_scratch_blob);
            user_data.captured_exception.reset();

            utils::check_cvode_flag(CVodeSetUserData(m_cvode_mem, &user_data), "CVodeSetUserData");

            LOG_TRACE_L2(m_logger, "Taking one CVODE step...");
            int flag = CVode(m_cvode_mem, netIn.tMax, m_Y, &current_time, CV_ONE_STEP);
            LOG_TRACE_L2(m_logger, "CVODE step complete. Current time: {}, step status: {}", current_time, utils::cvode_ret_code_map.at(flag));

            if (user_data.captured_exception){
                LOG_CRITICAL(m_logger, "An exception was captured during RHS evaluation ({}). Rethrowing...", user_data.captured_exception->what());
                std::rethrow_exception(std::make_exception_ptr(*user_data.captured_exception));
            }

            utils::check_cvode_flag(flag, "CVode");

            long int n_steps;
            double last_step_size;
            CVodeGetNumSteps(m_cvode_mem, &n_steps);
            CVodeGetLastStep(m_cvode_mem, &last_step_size);
            long int nliters, nlcfails;
            CVodeGetNumNonlinSolvIters(m_cvode_mem, &nliters);
            CVodeGetNumNonlinSolvConvFails(m_cvode_mem, &nlcfails);

            sunrealtype* y_data = N_VGetArrayPointer(m_Y);
            const double current_energy = y_data[numSpecies]; // Specific energy rate

            // TODO: Accumulate neutrino loss through the state vector directly which will allow CVODE to properly integrate it
            accumulated_neutrino_energy_loss += user_data.neutrino_energy_loss_rate * last_step_size;
            accumulated_total_neutrino_flux += user_data.total_neutrino_flux * last_step_size;

            size_t iter_diff = (total_nonlinear_iterations + nliters) - prev_nonlinear_iterations;
            size_t convFail_diff = (total_convergence_failures + nlcfails) - prev_convergence_failures;
            if (m_stdout_logging_enabled) {
                std::println(
                    "Step: {:6} | Updates: {:3} | Epoch Steps: {:4} | t: {:.3e} [s] | dt: {:15.6E} [s] | Iterations: {:6} (+{:2}) | Total Convergence Failures: {:2} (+{:2})",
                    total_steps + n_steps,
                    total_update_stages_triggered,
                    n_steps,
                    current_time,
                    last_step_size,
                    total_nonlinear_iterations + nliters,
                    iter_diff,
                    total_convergence_failures + nlcfails,
                    convFail_diff
                );
            }
            for (size_t i = 0; i < numSpecies; ++i) {
                const auto& species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
                if (y_data[i] > 0.0) {
                    postStep.setMolarAbundance(species, y_data[i]);
                }
            }
            // fourdst::composition::Composition collectedComposition = m_engine.collectComposition(postStep, netIn.temperature/1e9, netIn.density);
            // for (size_t i = 0; i < numSpecies; ++i) {
            //     y_data[i] = collectedComposition.getMolarAbundance(m_engine.getNetworkSpecies()[i]);
            // }
            LOG_INFO(m_logger, "Completed {:5} steps to time {:10.4E} [s] (dt = {:15.6E} [s]). Current specific energy: {:15.6E} [erg/g]", total_steps + n_steps, current_time, last_step_size, current_energy);
            LOG_DEBUG(m_logger, "Current composition (molar abundance): {}", [&]() -> std::string {
                std::stringstream ss;
                for (size_t i = 0; i < numSpecies; ++i) {
                    const auto& species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
                    ss << species.name() << ": (y_data = " << y_data[i] << ", collected = " << postStep.getMolarAbundance(species) << ")";
                    if (i < numSpecies - 1) {
                        ss << ", ";
                    }
                }
                return ss.str();
            }());

            static const std::map<fourdst::atomic::Species,
                      std::unordered_map<std::string,double>> kEmptyMap{};

            const auto& rcMap = user_data.reaction_contribution_map.has_value()
                                  ? user_data.reaction_contribution_map.value()
                                  : kEmptyMap;

            auto ctx = TimestepContext(
                current_time,
                m_Y,
                last_step_size,
                last_callback_time,
                T9,
                netIn.density,
                n_steps,
                m_engine,
                m_engine.getNetworkSpecies(*m_scratch_blob),
                convFail_diff,
                iter_diff,
                rcMap,
                *m_scratch_blob
            );

            prev_nonlinear_iterations = nliters + total_nonlinear_iterations;
            prev_convergence_failures = nlcfails + total_convergence_failures;

            if (m_callback.has_value()) {
                m_callback.value()(ctx);
            }
            trigger->step(ctx);

            if (m_detailed_step_logging) {
                log_step_diagnostics(*m_scratch_blob, user_data, true, true, true, "step_" + std::to_string(total_steps + n_steps) + ".json");
            }

            if (trigger->check(ctx)) {
                if (m_stdout_logging_enabled && displayTrigger) {
                    trigger::printWhy(trigger->why(ctx));
                }
                trigger->update(ctx);
                accumulated_energy += current_energy; // Add the specific energy rate to the accumulated energy
                total_nonlinear_iterations += nliters;
                total_convergence_failures += nlcfails;
                total_steps += n_steps;

                sunrealtype* end_of_step_abundances = N_VGetArrayPointer(ctx.state);

                LOG_INFO(
                    m_logger,
                    "Engine Update Triggered at time {} ({} update{} triggered). Current total specific energy {} [erg/g]",
                    current_time,
                    total_update_stages_triggered,
                    total_update_stages_triggered == 1 ? "" : "s",
                    accumulated_energy);
                total_update_stages_triggered++;

                fourdst::composition::Composition temp_comp;
                std::vector<double> mass_fractions;
                auto num_species_at_stop = static_cast<long int>(m_engine.getNetworkSpecies(*m_scratch_blob).size());

                if (num_species_at_stop > m_Y->ops->nvgetlength(m_Y) - 1) {
                    LOG_ERROR(
                        m_logger,
                        "Number of species at engine update ({}) exceeds the number of species in the CVODE solver ({}). This should never happen.",
                        num_species_at_stop,
                        m_Y->ops->nvgetlength(m_Y) - 1 // -1 due to energy in the last index
                    );
                    throw std::runtime_error("Number of species at engine update exceeds the number of species in the CVODE solver. This should never happen.");
                }

                for (const auto& species: m_engine.getNetworkSpecies(*m_scratch_blob)) {
                    const size_t sid = m_engine.getSpeciesIndex(*m_scratch_blob, species);
                    temp_comp.registerSpecies(species);
                    double y = end_of_step_abundances[sid];
                    if (y > 0.0) {
                        temp_comp.setMolarAbundance(species, y);
                    }
                }

#ifndef NDEBUG
                for (long int i = 0; i < num_species_at_stop; ++i) {
                    const auto& species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
                    if (std::abs(temp_comp.getMolarAbundance(species) - y_data[i]) > 1e-12) {
                        throw exceptions::UtilityError("Conversion from solver state to composition molar abundance failed verification.");
                    }
                }
#endif

                NetIn netInTemp;
                netInTemp.temperature = T9 * 1e9; // Convert back to Kelvin
                netInTemp.density = netIn.density;
                netInTemp.composition = temp_comp;

                LOG_DEBUG(
                    m_logger,
                    "Prior to Engine update composition is (molar abundance) {}",
                    [&]() -> std::string {
                        std::stringstream ss;
                        const size_t ns = temp_comp.size();
                        size_t count = 0;
                        for (const auto &symbol : temp_comp | std::views::keys) {
                            ss << symbol << ": " << temp_comp.getMolarAbundance(symbol);
                            if (count < ns - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }()
                );

                LOG_DEBUG(
                    m_logger,
                    "Prior to Engine Update active reactions are: {}",
                    [&]() -> std::string {
                        std::stringstream ss;
                        const gridfire::reaction::ReactionSet& reactions = m_engine.getNetworkReactions(*m_scratch_blob);
                        size_t count = 0;
                        for (const auto& reaction : reactions) {
                            ss << reaction -> id();
                            if (count < reactions.size() - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }()
                );
                fourdst::composition::Composition currentComposition = m_engine.project(*m_scratch_blob, netInTemp);
                LOG_DEBUG(
                    m_logger,
                    "After to Engine update composition is (molar abundance) {}",
                    [&]() -> std::string {
                        std::stringstream ss;
                        const size_t ns = currentComposition.size();
                        size_t count = 0;
                        for (const auto &symbol : currentComposition | std::views::keys) {
                            ss << symbol << ": " << currentComposition.getMolarAbundance(symbol);
                            if (count < ns - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }()
                );
                LOG_DEBUG(
                    m_logger,
                    "Fractional Abundance Changes: {}",
                    [&]() -> std::string {
                        std::stringstream ss;
                        const size_t ns = currentComposition.size();
                        size_t count = 0;
                        for (const auto &symbol : currentComposition | std::views::keys) {
                            if (!temp_comp.contains(symbol)) {
                                ss << symbol << ": New Species";
                            } else {
                                const double old_X = temp_comp.getMolarAbundance(symbol);
                                const double new_X = currentComposition.getMolarAbundance(symbol);
                                const double frac_change = (new_X - old_X) / (old_X + std::numeric_limits<double>::epsilon());
                                ss << symbol << ": " << frac_change;
                            }
                            if (count < ns - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }()
                );
                LOG_DEBUG(
                    m_logger,
                    "After Engine Update active reactions are: {}",
                    [&]() -> std::string {
                        std::stringstream ss;
                        const gridfire::reaction::ReactionSet& reactions = m_engine.getNetworkReactions(*m_scratch_blob);
                        size_t count = 0;
                        for (const auto& reaction : reactions) {
                            ss << reaction -> id();
                            if (count < reactions.size() - 1) {
                                ss << ", ";
                            }
                            count++;
                        }
                        return ss.str();
                    }()
                );
                LOG_INFO(
                    m_logger,
                    "Due to a triggered engine update the composition was updated from size {} to {} species.",
                    num_species_at_stop,
                    m_engine.getNetworkSpecies(*m_scratch_blob).size()
                );

                numSpecies = m_engine.getNetworkSpecies(*m_scratch_blob).size();
                size_t N = numSpecies + 1;

                LOG_INFO(m_logger, "Starting CVODE reinitialization after engine update...");
                cleanup_cvode_resources(true);

                m_cvode_mem = CVodeCreate(CV_BDF, m_sun_ctx);
                utils::check_cvode_flag(m_cvode_mem == nullptr ? -1 : 0, "CVodeCreate");

                initialize_cvode_integration_resources(N, numSpecies, current_time, currentComposition, absTol, relTol, accumulated_energy);

                utils::check_cvode_flag(CVodeReInit(m_cvode_mem, current_time, m_Y), "CVodeReInit");
                // throw exceptions::DebugException("Debug");
                LOG_INFO(m_logger, "Done reinitializing CVODE after engine update. The next log messages will be from the first step after reinitialization...");
            }

        }

        if (m_stdout_logging_enabled) { // Flush the buffer if standard out logging is enabled
            std::cout << std::flush;
        }

        LOG_INFO(m_logger, "CVODE iteration complete");

        sunrealtype* y_data = N_VGetArrayPointer(m_Y);
        accumulated_energy += y_data[numSpecies];
        std::vector<double> y_vec(y_data, y_data + numSpecies);

        for (double & i : y_vec) {
            if (i < 0 && std::abs(i) < 1e-16) {
                i = 0.0; // Regularize tiny negative abundances to zero
            }
        }

        LOG_INFO(m_logger, "Constructing final composition= with {} species", numSpecies);

        fourdst::composition::Composition topLevelComposition(m_engine.getNetworkSpecies(*m_scratch_blob), y_vec);
        LOG_INFO(m_logger, "Final composition constructed from solver state successfully! ({})", [&topLevelComposition]() -> std::string {
            std::ostringstream ss;
            size_t i = 0;
            for (const auto& [species, abundance] : topLevelComposition) {
                ss << species.name() << ": " << abundance;
                if (i < topLevelComposition.size() - 1) {
                    ss << ", ";
                }
                ++i;
            }
            return ss.str();
        }());

        LOG_INFO(m_logger, "Collecting final composition...");
        fourdst::composition::Composition outputComposition = m_engine.collectComposition(*m_scratch_blob, topLevelComposition, netIn.temperature/1e9, netIn.density);

        assert(outputComposition.getRegisteredSymbols().size() == equilibratedComposition.getRegisteredSymbols().size());

        LOG_INFO(m_logger, "Final composition constructed successfully! ({})", [&outputComposition]() -> std::string {
            std::ostringstream ss;
            size_t i = 0;
            for (const auto& [species, abundance] : outputComposition) {
                ss << species.name() << ": " << abundance;
                if (i < outputComposition.size() - 1) {
                    ss << ", ";
                }
                ++i;
            }
            return ss.str();
        }());

        LOG_INFO(m_logger, "Constructing output data...");
        NetOut netOut;
        netOut.composition = outputComposition;
        netOut.energy = accumulated_energy;
        utils::check_cvode_flag(CVodeGetNumSteps(m_cvode_mem, reinterpret_cast<long int *>(&netOut.num_steps)), "CVodeGetNumSteps");

        LOG_TRACE_L2(m_logger, "generating final nuclear energy generation rate derivatives...");
        auto [dEps_dT, dEps_dRho] = m_engine.calculateEpsDerivatives(
            *m_scratch_blob,
            outputComposition,
            T9,
            netIn.density
        );
        LOG_TRACE_L2(m_logger, "Found dEps/dT: {:0.3E} and dEps/dRho: {:0.3E}", dEps_dT, dEps_dRho);

        netOut.dEps_dT = dEps_dT;
        netOut.dEps_dRho = dEps_dRho;
        netOut.specific_neutrino_energy_loss = accumulated_neutrino_energy_loss;
        netOut.specific_neutrino_flux = accumulated_total_neutrino_flux;
        LOG_TRACE_L2(m_logger, "Output data built!");
        LOG_TRACE_L2(m_logger, "Solver evaluation complete!.");

        m_last_composition_hash = netOut.composition.hash();
        m_last_size = netOut.composition.size() + 1;
        CVodeGetLastStep(m_cvode_mem, &m_last_good_time_step);
        return netOut;
    }

    void CVODESolverStrategy::set_callback(const std::any &callback) {
        m_callback = std::any_cast<TimestepCallback>(callback);
    }

    bool CVODESolverStrategy::get_stdout_logging_enabled() const {
        return m_stdout_logging_enabled;
    }

    void CVODESolverStrategy::set_stdout_logging_enabled(const bool logging_enabled) {
        m_stdout_logging_enabled = logging_enabled;
    }

    void CVODESolverStrategy::set_absTol(double absTol) {
        m_absTol = absTol;
    }

    void CVODESolverStrategy::set_relTol(double relTol) {
        m_relTol = relTol;
    }

    double CVODESolverStrategy::get_absTol() const {
        if (m_absTol.has_value()) {
            return m_absTol.value();
        } else {
            return -1.0;
        }
    }

    double CVODESolverStrategy::get_relTol() const {
        if (m_relTol.has_value()) {
            return m_relTol.value();
        } else {
            return -1.0;
        }
    }

    std::vector<std::tuple<std::string, std::string>> CVODESolverStrategy::describe_callback_context() const {
        return {};
    }

    int CVODESolverStrategy::cvode_rhs_wrapper(
        const sunrealtype t,
        const N_Vector y,
        const N_Vector ydot,
        void *user_data
    ) {
        auto* data = static_cast<CVODEUserData*>(user_data);
        const auto* instance = data->solver_instance;

        try {
            LOG_TRACE_L2(instance->m_logger, "CVODE RHS wrapper called at time {}", t);
            // ReSharper disable once CppUseStructuredBinding
            const auto result = instance->calculate_rhs(t, y, ydot, data);
            data->reaction_contribution_map = result.reaction_contribution_map;
            data->neutrino_energy_loss_rate = result.neutrino_energy_loss_rate;
            data->total_neutrino_flux = result.total_neutrino_flux;
            LOG_TRACE_L2(instance->m_logger, "CVODE RHS wrapper completed successfully at time {}", t);
            return 0;
        } catch (const exceptions::EngineError& e) {
            LOG_ERROR(instance->m_logger, "EngineError caught in CVODE RHS wrapper at time {}: {}. Will attempt to recover...", t, e.what());
            return 1; // 1 Indicates a recoverable error, CVODE will retry the step
        } catch (const std::exception& e) {
            LOG_CRITICAL(instance->m_logger, "Unrecoverable and Unknown exception caught in CVODE RHS wrapper at time {} ({})", t, e.what());
            // data->captured_exception = std::make_unique<exceptions::GridFireError>(e.what());
            return -1; // Some unrecoverable error
        }
    }

    int CVODESolverStrategy::cvode_jac_wrapper(
        sunrealtype t,
        N_Vector y,
        N_Vector ydot,
        SUNMatrix J,
        void *user_data,
        N_Vector tmp1,
        N_Vector tmp2,
        N_Vector tmp3
    ) {
        const auto* data = static_cast<CVODEUserData*>(user_data);
        const auto* engine = data->engine;
        const auto* solver_instance = data->solver_instance;

        LOG_TRACE_L2(solver_instance->m_logger, "CVODE Jacobian wrapper starting");
        const size_t numSpecies = engine->getNetworkSpecies(data->ctx).size();
        sunrealtype* y_data = N_VGetArrayPointer(y);


        // Solver constraints should keep these values very close to 0 but floating point noise can still result in very
        // small negative numbers which can result in NaN's and more immediate crashes in the composition
        // finalization stage
        for (size_t i = 0; i < numSpecies; ++i) {
            if (y_data[i] < 0.0) {
                y_data[i] = 0.0;
            }
        }
        std::vector<double> y_vec(y_data, y_data + numSpecies);
        fourdst::composition::Composition composition(engine->getNetworkSpecies(data->ctx), y_vec);
        LOG_TRACE_L2(solver_instance->m_logger, "Generating Jacobian matrix at time {} with {} species in composition (mean molecular mass: {})", t, composition.size(), composition.getMeanParticleMass());
        LOG_TRACE_L2(solver_instance->m_logger, "Composition is {}", [&composition]() -> std::string {
            std::stringstream ss;
            size_t i = 0;
            for (const auto& [species, abundance] : composition) {
                ss << species.name() << ": " << abundance;
                if (i < composition.size() - 1) {
                    ss << ", ";
                }
                i++;
            }
            return ss.str();
        }());

        LOG_TRACE_L2(solver_instance->m_logger, "Generating Jacobian matrix at time {}", t);
        NetworkJacobian jac = engine->generateJacobianMatrix(data->ctx, composition, data->T9, data->rho);
        LOG_TRACE_L2(solver_instance->m_logger, "Regularizing Jacobian matrix at time {}", t);
        jac = regularize_jacobian(jac, composition, solver_instance->m_logger);
        LOG_TRACE_L2(solver_instance->m_logger, "Done regularizing Jacobian matrix at time {}", t);
        if (!jac.infs().empty() || !jac.nans().empty()) {
            auto infString = [&jac]() -> std::string {
                std::stringstream ss;
                size_t i = 0;
                std::vector<JacobianEntry> entries = jac.infs();
                for (const auto &[fst, snd]: entries | std::views::keys) {
                    ss << "J(" << fst.name() << ", " << snd.name() << ")";
                    if (i < entries.size() - 1) {
                        ss << ", ";
                    }
                    i++;
                }
                if (entries.empty()) {
                    ss << "None";
                }
                return ss.str();
            };
            auto nanString = [&jac]() -> std::string {
                std::stringstream ss;
                size_t i = 0;
                std::vector<JacobianEntry> entries = jac.nans();
                for (const auto &[fst, snd]: entries | std::views::keys) {
                    ss << "J(" << fst.name() << ", " << snd.name() << ")";
                    if (i < entries.size() - 1) {
                        ss << ", ";
                    }
                    i++;
                }
                if (entries.empty()) {
                    ss << "None";
                }
                return ss.str();
            };
            LOG_ERROR(
                solver_instance->m_logger,
                "Jacobian matrix generated at time {} contains {} infinite entries ({}) and {} NaN entries ({}). This will lead to a solver failure. GridFire will now halt.",
                t,
                jac.infs().size(),
                infString(),
                jac.nans().size(),
                nanString()
            );
            throw exceptions::IllConditionedJacobianError(std::format("Jacobian matrix generated at time {} contains {} infinite entries ({}) and {} NaN entries ({}). This will lead to a solver failure. In order to ensure tractability GridFire will not proceed. Focus on improving conditioning of the Jacobian matrix. If you believe this is an error please contact the GridFire developers.", t, jac.infs().size(), infString(), jac.nans().size(), nanString()));
        }
        LOG_TRACE_L2(solver_instance->m_logger, "Jacobian matrix created at time {} of shape ({} x {})", t, std::get<0>(jac.shape()), std::get<1>(jac.shape()));

        sunrealtype* J_data = SUNDenseMatrix_Data(J);
        const long int N = SUNDenseMatrix_Columns(J);

        LOG_TRACE_L2(solver_instance->m_logger, "Transferring Jacobian matrix data to SUNDenseMatrix format at time {}", t);
        for (size_t j = 0; j < numSpecies; ++j) {
            const fourdst::atomic::Species& species_j = engine->getNetworkSpecies(data->ctx)[j];
            for (size_t i = 0; i < numSpecies; ++i) {
                const fourdst::atomic::Species& species_i = engine->getNetworkSpecies(data->ctx)[i];
                // J(i,j) = d(f_i)/d(y_j)
                // Column-major order format for SUNDenseMatrix: J_data[j*N + i] indexes J(i,j)
                const double dYi_dt = jac(species_i, species_j);
                J_data[j * N + i] = dYi_dt;
            }
        }
        LOG_TRACE_L2(solver_instance->m_logger, "Done transferring Jacobian matrix data to SUNDenseMatrix format at time {}", t);

        // For now assume that the energy derivatives wrt. abundances are zero
        // TODO: Need a better way to build this part of the output jacobian so it properly pushes the solver
        //       in the right direction. Currently we effectively are doing a fixed point iteration in energy space.
        for (long int i = 1; i < N; ++i) {
            J_data[(N - 1) * N + i] = 0.0; // df(energy_dot)/df(y_i)
            J_data[i * N + (N - 1)] = 0.0; // df(f_i)/df(energy_dot)
        }

        return 0;
    }

    CVODESolverStrategy::CVODERHSOutputData CVODESolverStrategy::calculate_rhs(
        const sunrealtype t,
        N_Vector y,
        N_Vector ydot,
        const CVODEUserData *data
    ) const {
        const size_t numSpecies = m_engine.getNetworkSpecies(data->ctx).size();
        sunrealtype* y_data = N_VGetArrayPointer(y);

        // Solver constraints should keep these values very close to 0 but floating point noise can still result in very
        // small negative numbers which can result in NaN's and more immediate crashes in the composition
        // finalization stage
        for (size_t i = 0; i < numSpecies; ++i) {
            if (y_data[i] < 0.0) {
                y_data[i] = 0.0;
            }
        }
        std::vector<double> y_vec(y_data, y_data + numSpecies);
        fourdst::composition::Composition composition(m_engine.getNetworkSpecies(*m_scratch_blob), y_vec);

        LOG_TRACE_L2(m_logger, "Calculating RHS at time {} with {} species in composition", t, composition.size());
        const auto result = m_engine.calculateRHSAndEnergy(*m_scratch_blob, composition, data->T9, data->rho, false);
        if (!result) {
            LOG_CRITICAL(m_logger, "Failed to calculate RHS at time {}: {}", t, EngineStatus_to_string(result.error()));
            throw exceptions::BadRHSEngineError(std::format("Failed to calculate RHS at time {}: {}", t, EngineStatus_to_string(result.error())));
        }

        sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
        const auto& [
            dydt,
            nuclearEnergyGenerationRate,
            reactionContributions,
            neutrinoEnergyLossRate,
            totalNeutrinoFlux
        ] = result.value();

        LOG_TRACE_L2(m_logger, "Done calculating RHS at time {}, specific nuclear energy generation rate: {}", t, nuclearEnergyGenerationRate);
        LOG_TRACE_L2(m_logger, "RHS at time {} is {}", t, [&dydt]() -> std::string {
            std::stringstream ss;
            size_t i = 0;
            for (const auto& [species, rate] : dydt) {
                ss << "dY(" << species.name() << ")/dt" << ": " << rate;
                if (i < dydt.size() - 1) {
                    ss << ", ";
                }
                i++;
            }
            return ss.str();
        }());

        for (size_t i = 0; i < numSpecies; ++i) {
            fourdst::atomic::Species species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
            ydot_data[i] = dydt.at(species);
        }
        ydot_data[numSpecies] = nuclearEnergyGenerationRate; // Set the last element to the specific energy rate

        return {reactionContributions, result.value().neutrinoEnergyLossRate, result.value().totalNeutrinoFlux};
    }

    void CVODESolverStrategy::initialize_cvode_integration_resources(
        const uint64_t N,
        const size_t numSpecies,
        const double current_time,
        const fourdst::composition::Composition &composition,
        const double absTol,
        const double relTol,
        const double accumulatedEnergy
    ) {
        LOG_TRACE_L2(m_logger, "Initializing CVODE integration resources with N: {}, current_time: {}, absTol: {}, relTol: {}", N, current_time, absTol, relTol);
        cleanup_cvode_resources(false); // Cleanup any existing resources before initializing new ones

        m_Y = utils::init_sun_vector(N, m_sun_ctx);
        m_YErr = N_VClone(m_Y);

        sunrealtype *y_data = N_VGetArrayPointer(m_Y);
        for (size_t i = 0; i < numSpecies; i++) {
            const auto& species = m_engine.getNetworkSpecies(*m_scratch_blob)[i];
            if (composition.contains(species)) {
                y_data[i] = composition.getMolarAbundance(species);
            } else {
                y_data[i] = std::numeric_limits<double>::min(); // Species not in the composition, set to a small value
            }
        }
        y_data[numSpecies] = accumulatedEnergy; // Specific energy rate, initialized to zero


        utils::check_cvode_flag(CVodeInit(m_cvode_mem, cvode_rhs_wrapper, current_time, m_Y), "CVodeInit");
        utils::check_cvode_flag(CVodeSStolerances(m_cvode_mem, relTol, absTol), "CVodeSStolerances");

        // Constraints
        // We constrain the solution vector using CVODE's built in constraint flags as outlines on page 53 of the CVODE manual
        // https://computing.llnl.gov/sites/default/files/cv_guide-5.7.0.pdf
        // For completeness and redundancy against future dead links the flags have been copied here
        // 0.0: No constraint on the corresponding component of y.
        // 1.0: The corresponding component of y is constrained to be >= 0.
        // -1.0: The corresponding component of y is constrained to be <= 0.
        // 2.0: The corresponding component of y is constrained to be > 0.
        // -2.0: The corresponding component of y is constrained to be < 0
        // Here we use 1.0 for all species to ensure they remain non-negative.

        m_constraints = N_VClone(m_Y);
        if (m_constraints == nullptr) {
            LOG_ERROR(m_logger, "Failed to create constraints vector for CVODE");
            throw std::runtime_error("Failed to create constraints vector for CVODE");
        }
        N_VConst(1.0, m_constraints); // Set all constraints to >= 0 (note this is where the flag values are set)

        utils::check_cvode_flag(CVodeSetConstraints(m_cvode_mem, m_constraints), "CVodeSetConstraints");

        utils::check_cvode_flag(CVodeSetMaxStep(m_cvode_mem, 1.0e20), "CVodeSetMaxStep");

        m_J = SUNDenseMatrix(static_cast<sunindextype>(N), static_cast<sunindextype>(N), m_sun_ctx);
        utils::check_cvode_flag(m_J == nullptr ? -1 : 0, "SUNDenseMatrix");
        m_LS = SUNLinSol_Dense(m_Y, m_J, m_sun_ctx);
        utils::check_cvode_flag(m_LS == nullptr ? -1 : 0, "SUNLinSol_Dense");

        utils::check_cvode_flag(CVodeSetLinearSolver(m_cvode_mem, m_LS, m_J), "CVodeSetLinearSolver");
        utils::check_cvode_flag(CVodeSetJacFn(m_cvode_mem, cvode_jac_wrapper), "CVodeSetJacFn");
        LOG_TRACE_L2(m_logger, "CVODE solver initialized");
    }

    void CVODESolverStrategy::cleanup_cvode_resources(const bool memFree) {
        LOG_TRACE_L2(m_logger, "Cleaning up cvode resources");
        if (m_LS) SUNLinSolFree(m_LS);
        if (m_J) SUNMatDestroy(m_J);
        if (m_Y) N_VDestroy(m_Y);
        if (m_YErr) N_VDestroy(m_YErr);
        if (m_constraints) N_VDestroy(m_constraints);

        m_LS = nullptr;
        m_J = nullptr;
        m_Y = nullptr;
        m_YErr = nullptr;
        m_constraints = nullptr;

        if (memFree) {
            if (m_cvode_mem) CVodeFree(&m_cvode_mem);
            m_cvode_mem = nullptr;
        }
        LOG_TRACE_L2(m_logger, "Done Cleaning up cvode resources");
    }

    void CVODESolverStrategy::set_detailed_step_logging(const bool enabled) {
        m_detailed_step_logging = enabled;
    }

    void CVODESolverStrategy::log_step_diagnostics(
        scratch::StateBlob &ctx,
        const CVODEUserData &user_data,
        bool displayJacobianStiffness,
        bool displaySpeciesBalance,
        bool to_file, std::optional<std::string> filename
    ) const {
        if (to_file && !filename.has_value()) {
            LOG_ERROR(m_logger, "Filename must be provided when logging diagnostics to file.");
            throw exceptions::UtilityError("Filename must be provided when logging diagnostics to file.");
        }

        // --- 1. Get CVODE Step Statistics ---
        sunrealtype hlast, hcur, tcur;
        int qlast;

        utils::check_cvode_flag(CVodeGetLastStep(m_cvode_mem, &hlast), "CVodeGetLastStep");
        utils::check_cvode_flag(CVodeGetCurrentStep(m_cvode_mem, &hcur), "CVodeGetCurrentStep");
        utils::check_cvode_flag(CVodeGetLastOrder(m_cvode_mem, &qlast), "CVodeGetLastOrder");
        utils::check_cvode_flag(CVodeGetCurrentTime(m_cvode_mem, &tcur), "CVodeGetCurrentTime");

        nlohmann::json j;
        {
            std::vector<std::string> labels = {"Current Time (tcur)", "Last Step (hlast)", "Current Step (hcur)", "Last Order (qlast)"};
            std::vector<double> values = {static_cast<double>(tcur), static_cast<double>(hlast), static_cast<double>(hcur), static_cast<double>(qlast)};

            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Statistic", labels));
            columns.push_back(std::make_unique<utils::Column<double>>("Value", values));

            if (to_file) {
                j["CVODE_Step_Stats"] = utils::to_json(columns);
            } else {
                utils::print_table("CVODE Step Stats", columns);
            }
        }

        // --- 2. Get CVODE Cumulative Solver Statistics ---
        // These are the CRITICAL counters for diagnosing your problem
        long int nsteps, nfevals, nlinsetups, netfails, nniters, nconvfails, nsetfails;

        utils::check_cvode_flag(CVodeGetNumSteps(m_cvode_mem, &nsteps), "CVodeGetNumSteps");
        utils::check_cvode_flag(CVodeGetNumRhsEvals(m_cvode_mem, &nfevals), "CVodeGetNumRhsEvals");
        utils::check_cvode_flag(CVodeGetNumLinSolvSetups(m_cvode_mem, &nlinsetups), "CVodeGetNumLinSolvSetups");
        utils::check_cvode_flag(CVodeGetNumErrTestFails(m_cvode_mem, &netfails), "CVodeGetNumErrTestFails");
        utils::check_cvode_flag(CVodeGetNumNonlinSolvIters(m_cvode_mem, &nniters), "CVodeGetNumNonlinSolvIters");
        utils::check_cvode_flag(CVodeGetNumNonlinSolvConvFails(m_cvode_mem, &nconvfails), "CVodeGetNumNonlinSolvConvFails");
        utils::check_cvode_flag(CVodeGetNumLinConvFails(m_cvode_mem, &nsetfails), "CVodeGetNumLinConvFails");


        {
            std::vector<std::string> labels = {
                "Total Steps",
                "RHS Evals",
                "Linear Solver Setups (Jacobians)",
                "Total Newton Iters",
                "Error Test Fails",
                "Convergence Fails",
                "Linear Convergence Failures"
            };
            // --- ADDED nsetfails TO THIS LIST ---
            std::vector<long int> values = {nsteps, nfevals, nlinsetups, nniters, netfails, nconvfails, nsetfails};

            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<std::string>>("Counter", labels));
            columns.push_back(std::make_unique<utils::Column<long int>>("Count", values));

            if (to_file) {
                j["CVODE_Cumulative_Stats"] = utils::to_json(columns);
            } else {
                utils::print_table("CVODE Cumulative Stats", columns);
            }
        }

        // --- 3. Get Estimated Local Errors (Your Original Logic) ---
        utils::check_cvode_flag(CVodeGetEstLocalErrors(m_cvode_mem, m_YErr), "CVodeGetEstLocalErrors");

        sunrealtype *y_data = N_VGetArrayPointer(m_Y);
        sunrealtype *y_err_data = N_VGetArrayPointer(m_YErr);

        const auto absTol = m_config->solver.cvode.absTol;
        const auto relTol = m_config->solver.cvode.relTol;

        std::vector<double> err_ratios;
        const size_t num_components = N_VGetLength(m_Y);
        err_ratios.resize(num_components - 1); // Assuming -1 is for Energy or similar

        std::vector<double> Y_full(y_data, y_data + num_components - 1);
        std::vector<double> E_full(y_err_data, y_err_data + num_components - 1);

        auto result = diagnostics::report_limiting_species(ctx, *user_data.engine, Y_full, E_full, relTol, absTol, 10, to_file);
        if (to_file && result.has_value()) {
            j["Limiting_Species"] = result.value();
        }

        std::ranges::replace_if(
            Y_full,
            [](const double val) {
                return val < 0.0 && val > -1.0e-16;
            },
            0.0
        );

        for (size_t i = 0; i < num_components - 1; i++) {
            const double weight = relTol * std::abs(y_data[i]) + absTol;
            if (weight == 0.0) {
                err_ratios[i] = 0.0; // Avoid division by zero
                continue;
            }
            const double err_ratio = std::abs(y_err_data[i]) / weight;
            err_ratios[i] = err_ratio;
        }

        fourdst::composition::Composition composition(user_data.engine->getNetworkSpecies(*m_scratch_blob), Y_full);
        fourdst::composition::Composition collectedComposition = user_data.engine->collectComposition(*m_scratch_blob, composition, user_data.T9, user_data.rho);

        auto destructionTimescales = user_data.engine->getSpeciesDestructionTimescales(*m_scratch_blob, collectedComposition, user_data.T9, user_data.rho);
        auto netTimescales = user_data.engine->getSpeciesTimescales(*m_scratch_blob, collectedComposition, user_data.T9, user_data.rho);

        bool timescaleOkay = false;
        if (destructionTimescales && netTimescales) timescaleOkay = true;

        if (timescaleOkay) {
            std::vector<std::unique_ptr<utils::ColumnBase>> timescaleColumns;
            std::vector<fourdst::atomic::Species> species_list;
            std::vector<double> netTimescales_list;
            std::vector<double> destructionTimescales_list;
            std::vector<std::string> speciesStatus_list;

            for (const auto &sp: collectedComposition | std::views::keys) {
                species_list.emplace_back(sp);
                if (netTimescales.value().contains(sp)) netTimescales_list.emplace_back(netTimescales.value().at(sp));
                else netTimescales_list.emplace_back(std::numeric_limits<double>::infinity());

                if (destructionTimescales.value().contains(sp)) destructionTimescales_list.emplace_back(destructionTimescales.value().at(sp));
                else destructionTimescales_list.emplace_back(std::numeric_limits<double>::infinity());

                speciesStatus_list.push_back(SpeciesStatus_to_string(user_data.engine->getSpeciesStatus(*m_scratch_blob, sp)));
            }

            utils::Column<fourdst::atomic::Species> speciesColumn("Species", species_list);
            utils::Column<double> netTimescaleColumn("Net Timescale (s)", netTimescales_list);
            utils::Column<double> destructionTimescaleColumn("Destruction Timescale (s)", destructionTimescales_list);
            utils::Column<std::string> speciesStatusColumn("Species Status", speciesStatus_list);

            std::vector<std::unique_ptr<utils::ColumnBase>> table;
            table.push_back(std::make_unique<utils::Column<fourdst::atomic::Species>>(speciesColumn));
            table.push_back(std::make_unique<utils::Column<double>>(netTimescaleColumn));
            table.push_back(std::make_unique<utils::Column<double>>(destructionTimescaleColumn));
            table.push_back(std::make_unique<utils::Column<std::string>>(speciesStatusColumn));

            if (to_file) {
                j["Species_Timescales"] = utils::to_json(table);
            } else {
                utils::print_table("Species Timescales", table);
            }
        }


        if (err_ratios.empty()) {
            std::cout << "Error ratios vector is empty." << std::endl;
            return;
        }

        std::vector<size_t> indices(composition.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }

        std::ranges::sort(
            indices,
            [&err_ratios](const size_t i1, const size_t i2) {
                return err_ratios[i1] > err_ratios[i2];
            }
        );

        std::vector<fourdst::atomic::Species> sorted_species;
        std::vector<double> sorted_err_ratios;

        sorted_species.reserve(indices.size());
        sorted_err_ratios.reserve(indices.size());

        for (const auto idx: indices) {
            sorted_species.push_back(composition.getSpeciesAtIndex(idx));
            sorted_err_ratios.push_back(err_ratios[idx]);
        }

        {
            std::vector<std::unique_ptr<utils::ColumnBase>> columns;
            columns.push_back(std::make_unique<utils::Column<fourdst::atomic::Species>>("Species", sorted_species));
            columns.push_back(std::make_unique<utils::Column<double>>("Error Ratio", sorted_err_ratios));

            if (to_file) {
                j["Species_Error_Ratios_Linear"] = utils::to_json(columns);
            } else {
                utils::print_table("Species Error Ratios (Linear)", columns);
            }
        }

        // --- 4. Call Your Jacobian and Balance Diagnostics ---
        if (displayJacobianStiffness) {
            auto jStiff = diagnostics::inspect_jacobian_stiffness(ctx, *user_data.engine, composition, user_data.T9, user_data.rho, to_file);
            if (to_file && jStiff.has_value()) {
                j["Jacobian_Stiffness_Diagnostics"] = jStiff.value();
            }
        }
        if (displaySpeciesBalance) {
            // Limit this to the top N species to avoid spamming
            const size_t num_species_to_inspect = std::min(sorted_species.size(), static_cast<size_t>(5));
            for (size_t i = 0; i < num_species_to_inspect; ++i) {
                const auto& species = sorted_species[i];
                auto sbr = diagnostics::inspect_species_balance(ctx, *user_data.engine, std::string(species.name()), composition, user_data.T9, user_data.rho, to_file);
                if (to_file && sbr.has_value()) {
                    j[std::string("Species_Balance_Diagnostics_") + species.name().data()] = sbr.value();
                }
            }
        }

        if (to_file) {
            std::ofstream ofs(filename.value());
            if (!ofs.is_open()) {
                LOG_ERROR(m_logger, "Failed to open file {} for writing diagnostics.", filename.value());
                throw exceptions::UtilityError(std::format("Failed to open file {} for writing diagnostics.", filename.value()));
            }
            ofs << j.dump(4);
            ofs.close();
            LOG_TRACE_L2(m_logger, "Diagnostics written to file {}", filename.value());
        }
    }


}
