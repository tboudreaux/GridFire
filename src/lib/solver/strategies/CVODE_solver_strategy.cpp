#include "gridfire/solver/strategies/CVODE_solver_strategy.h"

#include "gridfire/network.h"
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

#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/solver/strategies/triggers/engine_partitioning_trigger.h"
#include "gridfire/trigger/procedures/trigger_pprint.h"


namespace {
    std::unordered_map<int, std::string> cvode_ret_code_map {
        {0, "CV_SUCCESS: The solver succeeded."},
        {1, "CV_TSTOP_RETURN: The solver reached the specified stopping time."},
        {2, "CV_ROOT_RETURN: A root was found."},
        {-99, "CV_WARNING: CVODE succeeded but in an unusual manner"},
        {-1, "CV_TOO_MUCH_WORK: The solver took too many internal steps."},
        {-2, "CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy requested."},
        {-3, "CV_ERR_FAILURE: The solver encountered a non-recoverable error."},
        {-4, "CV_CONV_FAILURE: The solver failed to converge."},
        {-5, "CV_LINIT_FAIL: The linear solver's initialization function failed."},
        {-6, "CV_LSETUP_FAIL: The linear solver's setup function failed."},
        {-7, "CV_LSOLVE_FAIL: The linear solver's solve function failed."},
        {-8, "CV_RHSFUNC_FAIL: The right-hand side function failed in an unrecoverable manner."},
        {-9, "CV_FIRST_RHSFUNC_ERR: The right-hand side function failed at the first call."},
        {-10, "CV_REPTD_RHSFUNC_ERR: The right-hand side function repeatedly failed recoverable."},
        {-11, "CV_UNREC_RHSFUNC_ERR: The right-hand side function failed unrecoverably."},
        {-12, "CV_RTFUNC_FAIL: The rootfinding function failed in an unrecoverable manner."},
        {-13, "CV_NLS_INIT_FAIL: The nonlinear solver's initialization function failed."},
        {-14, "CV_NLS_SETUP_FAIL: The nonlinear solver's setup function failed."},
        {-15, "CV_CONSTR_FAIL : The inequality constraint was violated and the solver was unable to recover."},
        {-16, "CV_NLS_FAIL: The nonlinear solver's solve function failed."},
        {-20, "CV_MEM_FAIL: Memory allocation failed."},
        {-21, "CV_MEM_NULL: The CVODE memory structure is NULL."},
        {-22, "CV_ILL_INPUT: An illegal input was detected."},
        {-23, "CV_NO_MALLOC: The CVODE memory structure has not been allocated."},
        {-24, "CV_BAD_K: The value of k is invalid."},
        {-25, "CV_BAD_T: The value of t is invalid."},
        {-26, "CV_BAD_DKY: The value of dky is invalid."},
        {-27, "CV_TOO_CLOSE: The time points are too close together."},
        {-28, "CV_VECTOROP_ERR: A vector operation failed."},
        {-29, "CV_PROJ_MEM_NULL: The projection memory structure is NULL."},
        {-30, "CV_PROJFUNC_FAIL: The projection function failed in an unrecoverable manner."},
        {-31, "CV_REPTD_PROJFUNC_ERR: THe projection function has repeated recoverable errors."}
    };
    void check_cvode_flag(const int flag, const std::string& func_name) {
        if (flag < 0) {
            if (!cvode_ret_code_map.contains(flag)) {
                throw std::runtime_error("CVODE error in " + func_name + ": Unknown error code: " + std::to_string(flag));
            }
            throw std::runtime_error("CVODE error in " + func_name + ": " + cvode_ret_code_map.at(flag));
        }
    }

    N_Vector init_sun_vector(uint64_t size, SUNContext sun_ctx) {
#ifdef SUNDIALS_HAVE_OPENMP
        N_Vector vec = N_VNew_OpenMP(size, 0, sun_ctx);
#elif SUNDIALS_HAVE_PTHREADS
        N_Vector vec = N_VNew_Pthreads(size, sun_ctx);
#else
        N_Vector vec = N_VNew_Serial(static_cast<long long>(size), sun_ctx);
#endif
        check_cvode_flag(vec == nullptr ? -1 : 0, "N_VNew");
        return vec;
    }
}

namespace gridfire::solver {

    CVODESolverStrategy::TimestepContext::TimestepContext(
        const double t,
        const N_Vector &state,
        const double dt,
        const double last_step_time,
        const double t9,
        const double rho,
        const size_t num_steps,
        const DynamicEngine &engine,
        const std::vector<fourdst::atomic::Species> &networkSpecies
    ) :
    t(t),
    state(state),
    dt(dt),
    last_step_time(last_step_time),
    T9(t9),
    rho(rho),
    num_steps(num_steps),
    engine(engine),
    networkSpecies(networkSpecies)
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
        return description;
    }


    CVODESolverStrategy::CVODESolverStrategy(DynamicEngine &engine): NetworkSolverStrategy<DynamicEngine>(engine) {
        // TODO: In order to support MPI this function must be changed
        const int flag = SUNContext_Create(SUN_COMM_NULL, &m_sun_ctx);
        if (flag < 0) {
            throw std::runtime_error("Failed to create SUNDIALS context (SUNDIALS Errno: " + std::to_string(flag) + ")");
        }
    }

    CVODESolverStrategy::~CVODESolverStrategy() {
        std::cout << "Cleaning up CVODE resources..." << std::endl;
        cleanup_cvode_resources(true);

        if (m_sun_ctx) {
            SUNContext_Free(&m_sun_ctx);
        }
    }

    NetOut CVODESolverStrategy::evaluate(const NetIn& netIn) {
        auto trigger = trigger::solver::CVODE::makeEnginePartitioningTrigger(1e12, 1e10, 1, true, 10);

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)

        const auto absTol = m_config.get<double>("gridfire:solver:CVODESolverStrategy:absTol", 1.0e-8);
        const auto relTol = m_config.get<double>("gridfire:solver:CVODESolverStrategy:relTol", 1.0e-8);

        fourdst::composition::Composition equilibratedComposition = m_engine.update(netIn);
        std::cout << "EXITED AT EXPECTED TESTING POINT" << std::endl;
        exit(0);

        size_t numSpecies = m_engine.getNetworkSpecies().size();
        uint64_t N = numSpecies + 1;

        m_cvode_mem = CVodeCreate(CV_BDF, m_sun_ctx);
        check_cvode_flag(m_cvode_mem == nullptr ? -1 : 0, "CVodeCreate");

        initialize_cvode_integration_resources(N, numSpecies, 0.0, equilibratedComposition, absTol, relTol, 0.0);

        CVODEUserData user_data;
        user_data.solver_instance = this;
        user_data.engine = &m_engine;

        double current_time = 0;
        // ReSharper disable once CppTooWideScope
        [[maybe_unused]] double last_callback_time = 0;
        m_num_steps = 0;
        double accumulated_energy = 0.0;
        int total_update_stages_triggered = 0;
        while (current_time < netIn.tMax) {
            user_data.T9 = T9;
            user_data.rho = netIn.density;
            user_data.networkSpecies = &m_engine.getNetworkSpecies();
            user_data.captured_exception.reset();

            check_cvode_flag(CVodeSetUserData(m_cvode_mem, &user_data), "CVodeSetUserData");

            int flag{};
            if (m_stdout_logging_enabled) {
                flag = CVode(m_cvode_mem, netIn.tMax, m_Y, &current_time, CV_ONE_STEP);
            } else {
                flag = CVode(m_cvode_mem, netIn.tMax, m_Y, &current_time, CV_NORMAL);
            }

            if (user_data.captured_exception){
                std::rethrow_exception(std::make_exception_ptr(*user_data.captured_exception));
            }

            check_cvode_flag(flag, "CVode");

            long int n_steps;
            double last_step_size;
            CVodeGetNumSteps(m_cvode_mem, &n_steps);
            CVodeGetLastStep(m_cvode_mem, &last_step_size);
            long int nliters, nlcfails;
            CVodeGetNumNonlinSolvIters(m_cvode_mem, &nliters);
            CVodeGetNumNonlinSolvConvFails(m_cvode_mem, &nlcfails);

            sunrealtype* y_data = N_VGetArrayPointer(m_Y);
            const double current_energy = y_data[numSpecies]; // Specific energy rate
            std::cout << std::scientific << std::setprecision(3)
                << "Step: " << std::setw(6) << n_steps
                << " | Time: " << current_time << " [s]"
                << " | Last Step Size: " << last_step_size
                << " | Accumulated Energy: " << current_energy << " [erg/g]"
                << " | NonLinIters: " << std::setw(2) << nliters
                << " | ConvFails: " << std::setw(2) << nlcfails
                << std::endl;

            auto ctx = TimestepContext(
                current_time,
                reinterpret_cast<N_Vector>(y_data),
                last_step_size,
                last_callback_time,
                T9,
                netIn.density,
                n_steps,
                m_engine,
                m_engine.getNetworkSpecies());

            if (trigger->check(ctx)) {
                trigger::printWhy(trigger->why(ctx));
                trigger->update(ctx);
                accumulated_energy += current_energy; // Add the specific energy rate to the accumulated energy
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
                size_t num_species_at_stop = m_engine.getNetworkSpecies().size();

                if (num_species_at_stop > m_Y->ops->nvgetlength(m_Y) - 1) {
                    LOG_ERROR(
                        m_logger,
                        "Number of species at engine update ({}) exceeds the number of species in the CVODE solver ({}). This should never happen.",
                        num_species_at_stop,
                        m_Y->ops->nvgetlength(m_Y) - 1 // -1 due to energy in the last index
                    );
                    throw std::runtime_error("Number of species at engine update exceeds the number of species in the CVODE solver. This should never happen.");
                }

                mass_fractions.reserve(num_species_at_stop);

                for (size_t i = 0; i < num_species_at_stop; ++i) {
                    const auto& species = m_engine.getNetworkSpecies()[i];
                    temp_comp.registerSpecies(species);
                    mass_fractions.push_back(y_data[i] * species.mass()); // Convert from molar abundance to mass fraction
                }
                temp_comp.setMassFraction(m_engine.getNetworkSpecies(), mass_fractions);
                temp_comp.finalize(true);

                NetIn netInTemp = netIn;
                netInTemp.temperature = T9 * 1e9; // Convert back to Kelvin
                netInTemp.density = netIn.density;
                netInTemp.composition = temp_comp;

                fourdst::composition::Composition currentComposition = m_engine.update(netInTemp);
                LOG_INFO(
                    m_logger,
                    "Due to a triggered stale engine the composition was updated from size {} to {} species.",
                    num_species_at_stop,
                    m_engine.getNetworkSpecies().size()
                );

                numSpecies = m_engine.getNetworkSpecies().size();
                N = numSpecies + 1;

                cleanup_cvode_resources(true);

                m_cvode_mem = CVodeCreate(CV_BDF, m_sun_ctx);
                check_cvode_flag(m_cvode_mem == nullptr ? -1 : 0, "CVodeCreate");

                initialize_cvode_integration_resources(N, numSpecies, current_time, currentComposition, absTol, relTol, accumulated_energy);

                check_cvode_flag(CVodeReInit(m_cvode_mem, current_time, m_Y), "CVodeReInit");
            }

        }

        sunrealtype* y_data = N_VGetArrayPointer(m_Y);
        accumulated_energy += y_data[numSpecies];

        std::vector<double> finalMassFractions(numSpecies);
        for (size_t i = 0; i < numSpecies; ++i) {
            const double molarMass = m_engine.getNetworkSpecies()[i].mass();
            finalMassFractions[i] = y_data[i] * molarMass; // Convert from molar abundance to mass fraction
            if (finalMassFractions[i] < MIN_ABUNDANCE_THRESHOLD) {
                finalMassFractions[i] = 0.0;
            }
        }

        std::vector<std::string> speciesNames;
        speciesNames.reserve(numSpecies);
        for (const auto& species : m_engine.getNetworkSpecies()) {
            speciesNames.emplace_back(species.name());
        }

        fourdst::composition::Composition outputComposition(speciesNames);
        outputComposition.setMassFraction(speciesNames, finalMassFractions);
        outputComposition.finalize(true);

        NetOut netOut;
        netOut.composition = outputComposition;
        netOut.energy = accumulated_energy;
        check_cvode_flag(CVodeGetNumSteps(m_cvode_mem, reinterpret_cast<long int *>(&netOut.num_steps)), "CVodeGetNumSteps");

        auto [dEps_dT, dEps_dRho] = m_engine.calculateEpsDerivatives(
            outputComposition,
            T9,
            netIn.density
        );

        netOut.dEps_dT = dEps_dT;
        netOut.dEps_dRho = dEps_dRho;

        return netOut;
    }

    void CVODESolverStrategy::set_callback(const std::any &callback) {
        m_callback = std::any_cast<TimestepCallback>(callback);
    }

    bool CVODESolverStrategy::get_stdout_logging_enabled() const {
        return m_stdout_logging_enabled;
    }

    void CVODESolverStrategy::set_stdout_logging_enabled(const bool value) {
        m_stdout_logging_enabled = value;
    }

    std::vector<std::tuple<std::string, std::string>> CVODESolverStrategy::describe_callback_context() const {
        return {};
    }

    int CVODESolverStrategy::cvode_rhs_wrapper(
        sunrealtype t,
        N_Vector y,
        N_Vector ydot,
        void *user_data
    ) {
        auto* data = static_cast<CVODEUserData*>(user_data);
        const auto* instance = data->solver_instance;

        try {
            instance->calculate_rhs(t, y, ydot, data);
            return 0;
        } catch (const exceptions::StaleEngineTrigger& e) {
            data->captured_exception = std::make_unique<exceptions::StaleEngineTrigger>(e);
            return 1; // 1 Indicates a recoverable error, CVODE will retry the step
        } catch (...) {
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

        const size_t numSpecies = engine->getNetworkSpecies().size();

        sunrealtype* J_data = SUNDenseMatrix_Data(J);
        const long int N = SUNDenseMatrix_Columns(J);

        for (size_t j = 0; j < numSpecies; ++j) {
            for (size_t i = 0; i < numSpecies; ++i) {
                const fourdst::atomic::Species& species_j = engine->getNetworkSpecies()[j];
                const fourdst::atomic::Species& species_i = engine->getNetworkSpecies()[i];
                // J(i,j) = d(f_i)/d(y_j)
                // Column-major order format for SUNDenseMatrix: J_data[j*N + i]
                J_data[j * N + i] = engine->getJacobianMatrixEntry(species_i, species_j);
            }
        }

        // For now assume that the energy derivatives wrt. abundances are zero
        for (size_t i = 0; i < N; ++i) {
            J_data[(N - 1) * N + i] = 0.0; // df(energy_dot)/df(y_i)
            J_data[i * N + (N - 1)] = 0.0; // df(f_i)/df(energy_dot)
        }

        return 0;
    }

    void CVODESolverStrategy::calculate_rhs(
        const sunrealtype t,
        N_Vector y,
        N_Vector ydot,
        const CVODEUserData *data
    ) const {
        const size_t numSpecies = m_engine.getNetworkSpecies().size();
        sunrealtype* y_data = N_VGetArrayPointer(y);

        // PERF: The trade off of ensured index consistency is some performance here. If this becomes a bottleneck we can revisit.
        //       The specific trade off is that we have decided to enforce that all interfaces accept composition objects rather
        //       than raw vectors of molar abundances. This then lets any method lookup the species by name rather than relying on
        //       the index in the vector being consistent. The trade off is that sometimes we need to construct a composition object
        //       which, at the moment, requires a somewhat expensive set of operations. Perhaps in the future we could enforce
        //       some consistent memory layout for the composition object to make this cheeper. That optimization would need to be
        //       done in the libcomposition library though...
        std::vector<double> y_vec(y_data, y_data + numSpecies);
        std::vector<std::string> symbols;
        symbols.reserve(numSpecies);
        for (const auto& species : m_engine.getNetworkSpecies()) {
            symbols.emplace_back(species.name());
        }
        std::vector<double> X;
        X.reserve(numSpecies);
        for (size_t i = 0; i < numSpecies; ++i) {
            const double molarMass = m_engine.getNetworkSpecies()[i].mass();
            X.push_back(y_vec[i] * molarMass); // Convert from molar abundance to mass fraction
        }
        fourdst::composition::Composition composition(symbols, X);

        std::ranges::replace_if(y_vec, [](const double val) { return val < 0.0; }, 0.0);

        const auto result = m_engine.calculateRHSAndEnergy(composition, data->T9, data->rho);
        if (!result) {
            throw exceptions::StaleEngineTrigger({data->T9, data->rho, y_vec, t, m_num_steps, y_data[numSpecies]});
        }

        sunrealtype* ydot_data = N_VGetArrayPointer(ydot);
        const auto& [dydt, nuclearEnergyGenerationRate] = result.value();

        for (size_t i = 0; i < numSpecies; ++i) {
            ydot_data[i] = dydt.at(m_engine.getNetworkSpecies()[i]);
        }
        ydot_data[numSpecies] = nuclearEnergyGenerationRate; // Set the last element to the specific energy rate
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
        cleanup_cvode_resources(false); // Cleanup any existing resources before initializing new ones

        m_Y = init_sun_vector(N, m_sun_ctx);
        m_YErr = N_VClone(m_Y);

        sunrealtype *y_data = N_VGetArrayPointer(m_Y);
        for (size_t i = 0; i < numSpecies; i++) {
            const auto& species = m_engine.getNetworkSpecies()[i];
            if (composition.contains(species)) {
                y_data[i] = composition.getMolarAbundance(species);
            } else {
                y_data[i] = std::numeric_limits<double>::min(); // Species not in the composition, set to a small value
            }
        }
        y_data[numSpecies] = accumulatedEnergy; // Specific energy rate, initialized to zero


        check_cvode_flag(CVodeInit(m_cvode_mem, cvode_rhs_wrapper, current_time, m_Y), "CVodeInit");
        check_cvode_flag(CVodeSStolerances(m_cvode_mem, relTol, absTol), "CVodeSStolerances");

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

        check_cvode_flag(CVodeSetConstraints(m_cvode_mem, m_constraints), "CVodeSetConstraints");

        check_cvode_flag(CVodeSetMaxStep(m_cvode_mem, 1.0e20), "CVodeSetMaxStep");

        m_J = SUNDenseMatrix(static_cast<sunindextype>(N), static_cast<sunindextype>(N), m_sun_ctx);
        check_cvode_flag(m_J == nullptr ? -1 : 0, "SUNDenseMatrix");
        m_LS = SUNLinSol_Dense(m_Y, m_J, m_sun_ctx);
        check_cvode_flag(m_LS == nullptr ? -1 : 0, "SUNLinSol_Dense");

        check_cvode_flag(CVodeSetLinearSolver(m_cvode_mem, m_LS, m_J), "CVodeSetLinearSolver");
        check_cvode_flag(CVodeSetJacFn(m_cvode_mem, cvode_jac_wrapper), "CVodeSetJacFn");
    }

    void CVODESolverStrategy::cleanup_cvode_resources(const bool memFree) {
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
    }

    void CVODESolverStrategy::log_step_diagnostics(const CVODEUserData &user_data) const {
        check_cvode_flag(CVodeGetEstLocalErrors(m_cvode_mem, m_YErr), "CVodeGetEstLocalErrors");

        sunrealtype *y_data = N_VGetArrayPointer(m_Y);
        sunrealtype *y_err_data = N_VGetArrayPointer(m_YErr);


        std::vector<double> err_ratios;
        std::vector<std::string> speciesNames;

        const auto absTol = m_config.get<double>("gridfire:solver:CVODESolverStrategy:absTol", 1.0e-8);
        const auto relTol = m_config.get<double>("gridfire:solver:CVODESolverStrategy:relTol", 1.0e-8);

        const size_t num_components = N_VGetLength(m_Y);
        err_ratios.resize(num_components - 1);

        std::vector<double> Y_full(y_data, y_data + num_components - 1);


        std::ranges::replace_if(
            Y_full,
            [](const double val) {
                return val < 0.0 && val > -1.0e-16;
            },
            0.0
        );

        for (size_t i = 0; i < num_components - 1; i++) {
            const double weight = relTol * std::abs(y_data[i]) + absTol;
            if (weight == 0.0) continue; // Skip components with zero weight

            const double err_ratio = std::abs(y_err_data[i]) / weight;

            err_ratios[i] = err_ratio;
            speciesNames.emplace_back(user_data.networkSpecies->at(i).name());
        }
        fourdst::composition::Composition composition(speciesNames, Y_full);

        if (err_ratios.empty()) {
            return;
        }

        std::vector<size_t> indices(speciesNames.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }

        std::ranges::sort(
            indices,
            [&err_ratios](const size_t i1, const size_t i2) {
                return err_ratios[i1] > err_ratios[i2];
            }
        );

        std::vector<std::string> sorted_speciesNames;
        std::vector<double> sorted_err_ratios;

        sorted_speciesNames.reserve(indices.size());
        sorted_err_ratios.reserve(indices.size());

        for (const auto idx: indices) {
            sorted_speciesNames.push_back(speciesNames[idx]);
            sorted_err_ratios.push_back(err_ratios[idx]);
        }



        std::vector<std::unique_ptr<utils::ColumnBase>> columns;
        columns.push_back(std::make_unique<utils::Column<std::string>>("Species", sorted_speciesNames));
        columns.push_back(std::make_unique<utils::Column<double>>("Error Ratio", sorted_err_ratios));

        std::cout << utils::format_table("Species Error Ratios", columns) << std::endl;
        diagnostics::inspect_jacobian_stiffness(*user_data.engine, composition, user_data.T9, user_data.rho);
        diagnostics::inspect_species_balance(*user_data.engine, "N-14", composition, user_data.T9, user_data.rho);
        diagnostics::inspect_species_balance(*user_data.engine, "n-1", composition, user_data.T9, user_data.rho);

    }
}
