#pragma once

#include "gridfire/solver/strategies/strategy_abstract.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/types/types.h"
#include "gridfire/config/config.h"

#include "fourdst/logging/logging.h"
#include "fourdst/constants/const.h"

#include <vector>
#include <functional>
#include <cvode/cvode.h>
#include <sundials/sundials_types.h>

#include "gridfire/exceptions/error_engine.h"

#ifdef SUNDIALS_HAVE_OPENMP
    #include <nvector/nvector_openmp.h>
#endif
#ifdef SUNDIALS_HAVE_PTHREADS
    #include <nvector/nvector_pthreads.hh>
#endif
#ifndef SUNDIALS_HAVE_OPENMP
    #ifndef SUNDIALS_HAVE_PTHREADS
        #include <nvector/nvector_serial.h>
    #endif
#endif

namespace gridfire::solver {
    class SpectralSolverStrategy final : public MultiZoneDynamicNetworkSolver {
    public:
        explicit SpectralSolverStrategy(const engine::DynamicEngine& engine);
        ~SpectralSolverStrategy() override;

        std::vector<NetOut> evaluate(
            const std::vector<NetIn> &netIns,
            const std::vector<double>& mass_coords, const engine::scratch::StateBlob &ctx_template
        ) override;

        void set_callback(const std::any &callback) override;
        [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;

        [[nodiscard]] bool get_stdout_logging_enabled() const;
        void set_stdout_logging_enabled(bool logging_enabled);

    public:
        struct TimestepContext final : public SolverContextBase {
            TimestepContext(
                const double t,
                const N_Vector &state,
                const double dt,
                const double last_time_step,
                const engine::DynamicEngine &engine
            ) :
            t(t),
            state(state),
            dt(dt),
            last_time_step(last_time_step),
            engine(engine) {}

            [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;

            const double t;
            const N_Vector& state;
            const double dt;
            const double last_time_step;
            const engine::DynamicEngine& engine;
        };

        struct BasisEval {
            size_t start_idx;
            std::vector<double> phi;
        };

        struct SplineBasis {
            std::vector<double> knots;
            std::vector<double> quadrature_nodes;
            std::vector<double> quadrature_weights;
            int degree = 3;


            std::vector<BasisEval> quad_evals;
        };
    public:
        using TimestepCallback = std::function<void(const TimestepContext&)>;
    private:

        enum class ParallelInitializationResult : uint8_t {
            SUCCESS,
            FAILURE
        };

        struct SpectralCoefficients {
            size_t num_sets;
            size_t num_coefficients;
            std::vector<double> coefficients;

            double operator()(size_t i, size_t j) const;
        };

        struct GridPoint {
            double T9;
            double rho;
            fourdst::composition::Composition composition;
        };

        struct Constants {
            const double c = fourdst::constant::Constants::getInstance().get("c").value;
            const double N_a = fourdst::constant::Constants::getInstance().get("N_a").value;
            const double c2 = c * c;
        };

        struct DenseLinearSolver {
            SUNMatrix A;
            SUNLinearSolver LS;
            N_Vector temp_vector;
            SUNContext ctx;

            DenseLinearSolver(size_t size, SUNContext sun_ctx);
            ~DenseLinearSolver();

            DenseLinearSolver(const DenseLinearSolver&) = delete;
            DenseLinearSolver& operator=(const DenseLinearSolver&) = delete;

            void setup() const;
            void zero() const;

            void init_from_cache(size_t num_basis_funcs, const std::vector<BasisEval>& shell_cache) const;
            void init_from_basis(size_t num_basis_funcs, const SplineBasis& basis) const;

            void solve_inplace(N_Vector x, size_t num_vars, size_t basis_size) const;
            void solve_inplace_ptr(sunrealtype* data_ptr, size_t num_vars, size_t basis_size) const;
        };

        struct CVODEUserData {
            SpectralSolverStrategy* solver_instance{};
            std::vector<std::reference_wrapper<engine::scratch::StateBlob>> workspaces;
            const engine::DynamicEngine* engine{};
            std::unique_ptr<exceptions::EngineError> captured_exception{};

            std::vector<double> T9{};
            std::vector<double> rho{};
            double energy{};

            double neutrino_energy_loss_rate = 0.0;
            double total_neutrino_flux = 0.0;

            DenseLinearSolver* mass_matrix_solver_instance{};
            const SplineBasis* basis{};
        };

    private:
        fourdst::config::Config<config::GridFireConfig> m_config;
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");

        SUNContext m_sun_ctx = nullptr;   ///< SUNDIALS context (lifetime of the solver).
        void* m_cvode_mem = nullptr;      ///< CVODE memory block.
        N_Vector m_Y = nullptr;           ///< CVODE state vector (species + energy accumulator).
        SUNMatrix m_J = nullptr;          ///< Dense Jacobian matrix.
        SUNLinearSolver m_LS = nullptr;   ///< Dense linear solver.


        std::optional<TimestepCallback> m_callback;      ///< Optional per-step callback.
        int m_num_steps = 0;              ///< CVODE step counter (used for diagnostics and triggers).

        bool m_stdout_logging_enabled = true; ///< If true, print per-step logs and use CV_ONE_STEP.

        N_Vector m_constraints = nullptr; ///< CVODE constraints vector (>= 0 for species entries).

        std::optional<double> m_absTol;        ///< User-specified absolute tolerance.
        std::optional<double> m_relTol;        ///< User-specified relative tolerance.

        bool m_detailed_step_logging = false;    ///< If true, log detailed step diagnostics (error ratios, Jacobian, species balance).

        mutable size_t m_last_size = 0;
        mutable size_t m_last_composition_hash = 0ULL;
        mutable sunrealtype m_last_good_time_step = 0ULL;

        SplineBasis m_current_basis;

        Constants m_constants;

        N_Vector m_T_coeffs = nullptr;
        N_Vector m_rho_coeffs = nullptr;

        std::vector<fourdst::atomic::Species> m_global_species_list;


    private:
        std::vector<double> evaluate_monitor_function(const std::vector<NetIn>& current_shells) const;

        static SplineBasis generate_basis_from_monitor(const std::vector<double>& monitor_values, const std::vector<double>& mass_coordinates, size_t actual_elements);

        GridPoint reconstruct_at_quadrature(const N_Vector y_coeffs, size_t quad_index, const SplineBasis &basis) const;

        std::vector<NetOut> reconstruct_solution(const std::vector<NetIn>& original_inputs, const std::vector<double>& mass_coordinates, const N_Vector final_coeffs, const SplineBasis& basis, double dt) const;

        static int cvode_rhs_wrapper(sunrealtype t, N_Vector y, N_Vector, void* user_data);
        static int cvode_jac_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        int calculate_rhs(sunrealtype t, N_Vector y_coeffs, N_Vector ydot_coeffs, CVODEUserData* data) const;
        int calculate_jacobian(sunrealtype t, N_Vector y_coeffs, N_Vector ydot_coeffs, SUNMatrix J, const CVODEUserData *data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) const;

        static size_t nyquist_elements(size_t requested_elements, size_t num_shells) ;

        static void project_specific_variable(
            const std::vector<NetIn>& current_shells,
            const std::vector<double>& mass_coordinates,
            const std::vector<BasisEval>& shell_cache,
            const DenseLinearSolver& linear_solver,
            N_Vector output_vec,
            size_t output_offset,
            const std::function<double(const NetIn&)> &getter,
            bool use_log
        );

        void inspect_jacobian(SUNMatrix J, const std::string& context) const;
    };

}
