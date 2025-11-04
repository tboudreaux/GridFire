#pragma once

#include "gridfire/solver/solver.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"
#include "gridfire/exceptions/exceptions.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/config/config.h"


#include <functional>
#include <any>
#include <string>
#include <vector>
#include <tuple>

// SUNDIALS/CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_types.h>

// Include headers for linear solvers and N_Vectors
// We will use preprocessor directives to select the correct ones
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

// These are the possible N_Vector implementations. We use the compiler defines to select the most appropriate one for the build.
// If none are defined, we default to serial.

// For precompiled binaries we will need to ensure that we have versions built for all three types (ideally with some runtime
// checks that will fail gracefully if the user tries to use an unsupported type).
#ifdef SUNDIALS_HAVE_OPENMP
    #include <nvector/nvector_openmp.h>
#endif
#ifdef SUNDIALS_HAVE_PTHREADS
    #include <nvector/nvector_pthreads.hh>
#endif
// Default to serial if no parallelism is enabled
#ifndef SUNDIALS_HAVE_OPENMP
    #ifndef SUNDIALS_HAVE_PTHREADS
        #include <nvector/nvector_serial.h>
    #endif
#endif

namespace gridfire::solver {
    /**
     * @class CVODESolverStrategy
     * @brief Stiff ODE integrator backed by SUNDIALS CVODE (BDF) for network + energy.
     *
     * Integrates the nuclear network abundances along with a final accumulator entry for specific
     * energy using CVODE's BDF method and a dense linear solver. The state vector layout is:
     *   [y_0, y_1, ..., y_{N-1}, eps], where eps is the accumulated specific energy (erg/g).
     *
     * Implementation summary:
     *  - Creates a SUNContext and CVODE memory; initializes the state from a Composition.
     *  - Enforces non-negativity on species via CVodeSetConstraints (>= 0 for all species slots).
     *  - Uses a user-provided DynamicEngine to compute RHS and to fill the dense Jacobian.
     *  - The Jacobian is assembled column-major into a SUNDenseMatrix; the energy row/column is
     *    currently set to zero (decoupled from abundances in the linearization).
     *  - An internal trigger can rebuild the engine/network; when triggered, CVODE resources are
     *    torn down and recreated with the new network size, preserving the energy accumulator.
     *  - The CVODE RHS wrapper captures exceptions::StaleEngineTrigger from the engine evaluation
     *    path as recoverable (return code 1) and stores a copy in user-data for the driver loop.
     *
     * @par Example
     * @code
     * using gridfire::solver::CVODESolverStrategy;
     * using gridfire::solver::NetIn;
     *
     * CVODESolverStrategy solver(engine);
     * NetIn in;
     * in.temperature = 1.0e9; // K
     * in.density = 1.0e6;     // g/cm^3
     * in.tMax = 1.0;          // s
     * in.composition = initialComposition;
     * auto out = solver.evaluate(in);
     * std::cout << "Final energy: " << out.energy << " erg/g\n";
     * @endcode
     */
    class CVODESolverStrategy final : public DynamicNetworkSolverStrategy {
    public:
        /**
         * @brief Construct the CVODE strategy and create a SUNDIALS context.
         * @param engine DynamicEngine used for RHS/Jacobian evaluation and network access.
         * @throws std::runtime_error If SUNContext_Create fails.
         */
        explicit CVODESolverStrategy(DynamicEngine& engine);
        /**
         * @brief Destructor: cleans CVODE/SUNDIALS resources and frees SUNContext.
         */
        ~CVODESolverStrategy() override;

        // Make the class non-copyable and non-movable to prevent shallow copies of CVODE pointers
        CVODESolverStrategy(const CVODESolverStrategy&) = delete;
        CVODESolverStrategy& operator=(const CVODESolverStrategy&) = delete;
        CVODESolverStrategy(CVODESolverStrategy&&) = delete;
        CVODESolverStrategy& operator=(CVODESolverStrategy&&) = delete;

        /**
         * @brief Integrate from t=0 to netIn.tMax and return final composition and energy.
         *
         * Implementation summary:
         *  - Converts temperature to T9, initializes CVODE memory and state (size = numSpecies + 1).
         *  - Repeatedly calls CVode in single-step or normal mode depending on stdout logging.
         *  - Wraps RHS to capture exceptions::StaleEngineTrigger as a recoverable step failure;
         *    if present after a step, it is rethrown for upstream handling.
         *  - Prints/collects diagnostics per step (step size, energy, solver iterations).
         *  - On trigger activation, rebuilds CVODE resources to reflect a changed network and
         *    reinitialized the state using the latest engine composition, preserving energy.
         *  - At the end, converts molar abundances to mass fractions and assembles NetOut,
         *    including derivatives of energy w.r.t. T and rho from the engine.
         *
         * @param netIn Inputs: temperature [K], density [g cm^-3], tMax [s], composition.
         * @return NetOut containing final Composition, accumulated energy [erg/g], step count,
         *         and dEps/dT, dEps/dRho.
         * @throws std::runtime_error If any CVODE or SUNDIALS call fails (negative return codes),
         *         or if internal consistency checks fail during engine updates.
         * @throws exceptions::StaleEngineTrigger Propagated if the engine signals a stale state
         *         during RHS evaluation (captured in the wrapper then rethrown here).
         */
        NetOut evaluate(const NetIn& netIn) override;

        /**
         * @brief Call to evaluate which will let the user control if the trigger reasoning is displayed
         * @param netIn Inputs: temperature [K], density [g cm^-3], tMax [s], composition.
         * @param displayTrigger Boolean flag to control if trigger reasoning is displayed
         * @return NetOut containing final Composition, accumulated energy [erg/g], step count,
         *         and dEps/dT, dEps/dRho.
         * @throws std::runtime_error If any CVODE or SUNDIALS call fails (negative return codes),
         *         or if internal consistency checks fail during engine updates.
         * @throws exceptions::StaleEngineTrigger Propagated if the engine signals a stale state
         *         during RHS evaluation (captured in the wrapper then rethrown here).
         */
        NetOut evaluate(const NetIn& netIn, bool displayTrigger);

        /**
         * @brief Install a timestep callback.
         * @param callback std::any containing TimestepCallback (std::function<void(const TimestepContext&)>).
         * @throws std::bad_any_cast If callback is not of the expected type.
         */
        void set_callback(const std::any &callback) override;

        /**
         * @brief Whether per-step logs are printed to stdout and CVode is stepped with CV_ONE_STEP.
         */
        [[nodiscard]] bool get_stdout_logging_enabled() const;

        /**
         * @brief Enable/disable per-step stdout logging.
         * @param logging_enabled Flag to control if a timestep summary is written to standard output or not
         */
        void set_stdout_logging_enabled(bool logging_enabled);

        /**
         * @brief Schema of fields exposed to the timestep callback context.
         */
        [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;

        /**
         * @struct TimestepContext
         * @brief Immutable view of the current integration state passed to callbacks.
         *
         * Fields capture CVODE time/state, step size, thermodynamic state, the engine reference,
         * and the list of network species used to interpret the state vector layout.
         */
        struct TimestepContext final : public SolverContextBase {
            // This struct can be identical to the one in DirectNetworkSolver
            const double t;                 ///< Current integration time [s].
            const N_Vector& state;          ///< Current CVODE state vector (N_Vector).
            const double dt;                ///< Last step size [s].
            const double last_step_time;    ///< Time at last callback [s].
            const double T9;                ///< Temperature in GK.
            const double rho;               ///< Density [g cm^-3].
            const size_t num_steps;         ///< Number of CVODE steps taken so far.
            const DynamicEngine& engine;    ///< Reference to the engine.
            const std::vector<fourdst::atomic::Species>& networkSpecies; ///< Species layout.
            const size_t currentConvergenceFailures; ///< Total number of convergence failures
            const size_t currentNonlinearIterations; ///< Total number of non linear iterations

            /**
             * @brief Construct a context snapshot.
             */
            TimestepContext(
                double t,
                const N_Vector& state,
                double dt,
                double last_step_time,
                double t9,
                double rho,
                size_t num_steps,
                const DynamicEngine& engine,
                const std::vector<fourdst::atomic::Species>& networkSpecies,
                size_t currentConvergenceFailure,
                size_t currentNonlinearIterations
            );

            /**
             * @brief Human-readable description of the context fields.
             */
            [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;
        };

        /**
         * @brief Type alias for a timestep callback.
         */
        using TimestepCallback = std::function<void(const TimestepContext& context)>; ///< Type alias for a timestep callback function.
    private:
        /**
         * @struct CVODEUserData
         * @brief A helper struct to pass C++ context to C-style CVODE callbacks.
         *
         * Carries pointers back to the solver instance and engine, the current thermodynamic
         * state, energy accumulator, and a slot to capture a copy of exceptions::StaleEngineTrigger
         * from RHS evaluation. The RHS wrapper treats this as a recoverable failure and returns 1
         * to CVODE, then the driver loop inspects and rethrows.
         */
        struct CVODEUserData {
            CVODESolverStrategy* solver_instance; // Pointer back to the class instance
            DynamicEngine* engine;
            double T9;
            double rho;
            double energy;
            const std::vector<fourdst::atomic::Species>* networkSpecies;
            std::unique_ptr<exceptions::StaleEngineTrigger> captured_exception = nullptr;
        };

    private:
        fourdst::config::Config& m_config = fourdst::config::Config::getInstance();
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /**
         * @brief CVODE RHS C-wrapper that delegates to calculate_rhs and captures exceptions.
         * @return 0 on success; 1 on recoverable StaleEngineTrigger; -1 on other failures.
         */
        static int cvode_rhs_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
        /**
         * @brief CVODE dense Jacobian C-wrapper that fills SUNDenseMatrix using the engine.
         *
         * Assembles J(i,j) = d(f_i)/d(y_j) for all species using engine->getJacobianMatrixEntry,
         * then zeros the last row and column corresponding to the energy variable.
         */
        static int cvode_jac_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        /**
         * @brief Compute RHS into ydot at time t from the engine and current state y.
         *
         * Converts the CVODE state to a Composition (mass fractions) and calls
         * engine.calculateRHSAndEnergy(T9, rho). Negative small abundances are clamped to zero
         * before constructing Composition. On stale engine, throws exceptions::StaleEngineTrigger.
         */
        void calculate_rhs(sunrealtype t, N_Vector y, N_Vector ydot, const CVODEUserData* data) const;

        /**
         * @brief Allocate and initialize CVODE vectors, linear algebra, tolerances, and constraints.
         *
         * State vector m_Y is sized to N (numSpecies + 1). Species slots are initialized from Composition
         * molar abundances when present, otherwise a tiny positive value; the last slot is set to
         * accumulatedEnergy. Sets scalar tolerances, non-negativity constraints for species, maximum
         * step size, creates a dense matrix and dense linear solver, and registers the Jacobian.
         */
        void initialize_cvode_integration_resources(
            uint64_t N,
            size_t numSpecies,
            double current_time,
            const fourdst::composition::Composition& composition,
            double absTol,
            double relTol,
            double accumulatedEnergy
        );

        /**
         * @brief Destroy CVODE vectors/linear algebra and optionally the CVODE memory block.
         * @param memFree If true, also calls CVodeFree on m_cvode_mem.
         */
        void cleanup_cvode_resources(bool memFree);

        /**
         * @brief Compute and print per-component error ratios; run diagnostic helpers.
         *
         * Gathers CVODE's estimated local errors, converts the state to a Composition, and prints a
         * sorted table of species with the highest error ratios; then invokes diagnostic routines to
         * inspect Jacobian stiffness and species balance.
         */
        void log_step_diagnostics(const CVODEUserData& user_data, bool displayJacobianStiffness) const;
    private:
        SUNContext m_sun_ctx = nullptr;   ///< SUNDIALS context (lifetime of the solver).
        void* m_cvode_mem = nullptr;      ///< CVODE memory block.
        N_Vector m_Y = nullptr;           ///< CVODE state vector (species + energy accumulator).
        N_Vector m_YErr = nullptr;        ///< Estimated local errors.
        SUNMatrix m_J = nullptr;          ///< Dense Jacobian matrix.
        SUNLinearSolver m_LS = nullptr;   ///< Dense linear solver.


        std::optional<TimestepCallback> m_callback;      ///< Optional per-step callback.
        int m_num_steps = 0;              ///< CVODE step counter (used for diagnostics and triggers).

        bool m_stdout_logging_enabled = true; ///< If true, print per-step logs and use CV_ONE_STEP.

        N_Vector m_constraints = nullptr; ///< CVODE constraints vector (>= 0 for species entries).
    };
}