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
#include <cvode/cvode.h>      // For CVDls (serial dense linear solver)
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

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
    class CVODESolverStrategy final : public DynamicNetworkSolverStrategy {
    public:
        explicit CVODESolverStrategy(DynamicEngine& engine);
        ~CVODESolverStrategy() override;

        // Make the class non-copyable and non-movable to prevent shallow copies of CVODE pointers
        CVODESolverStrategy(const CVODESolverStrategy&) = delete;
        CVODESolverStrategy& operator=(const CVODESolverStrategy&) = delete;
        CVODESolverStrategy(CVODESolverStrategy&&) = delete;
        CVODESolverStrategy& operator=(CVODESolverStrategy&&) = delete;

        NetOut evaluate(const NetIn& netIn) override;

        void set_callback(const std::any &callback) override;

        bool get_stdout_logging_enabled() const;

        void set_stdout_logging_enabled(const bool value);

        [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;

        struct TimestepContext final : public SolverContextBase {
            // This struct can be identical to the one in DirectNetworkSolver
            const double t;
            const N_Vector& state; // Note: state is now an N_Vector
            const double dt;
            const double last_step_time;
            const double T9;
            const double rho;
            const int num_steps;
            const DynamicEngine& engine;
            const std::vector<fourdst::atomic::Species>& networkSpecies;

            // Constructor
            TimestepContext(
                double t, const N_Vector& state, double dt, double last_step_time,
                double t9, double rho, int num_steps, const DynamicEngine& engine,
                const std::vector<fourdst::atomic::Species>& networkSpecies
            );

            [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;
        };

        using TimestepCallback = std::function<void(const TimestepContext& context)>; ///< Type alias for a timestep callback function.
    private:
        /**
         * @struct CVODEUserData
         * @brief A helper struct to pass C++ context to C-style CVODE callbacks.
         *
         * CVODE callbacks are C functions and use a void* pointer to pass user data.
         * This struct bundles all the necessary C++ objects (like 'this', engine references, etc.)
         * to be accessed safely within those static C wrappers.
         */
        struct CVODEUserData {
            CVODESolverStrategy* solver_instance; // Pointer back to the class instance
            DynamicEngine* engine;
            double T9;
            double rho;
            const std::vector<fourdst::atomic::Species>* networkSpecies;
            std::unique_ptr<exceptions::StaleEngineTrigger> captured_exception = nullptr;
        };

    private:
        Config& m_config = Config::getInstance();
        static int cvode_rhs_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);
        static int cvode_jac_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        int calculate_rhs(sunrealtype t, N_Vector y, N_Vector ydot, const CVODEUserData* data) const;

        void initialize_cvode_integration_resources(
            uint64_t N,
            size_t numSpecies,
            double current_time,
            const fourdst::composition::Composition& composition,
            double absTol,
            double relTol,
            double accumulatedEnergy
        );

        void cleanup_cvode_resources(bool memFree);
    private:
        SUNContext m_sun_ctx = nullptr;
        void* m_cvode_mem = nullptr;
        N_Vector m_Y = nullptr;
        SUNMatrix m_J = nullptr;
        SUNLinearSolver m_LS = nullptr;


        TimestepCallback m_callback;
        int m_num_steps = 0;

        bool m_stdout_logging_enabled = true;
    };
}