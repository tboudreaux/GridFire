#pragma once

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "quill/Logger.h"

#include <functional>
#include <any>
#include <vector>
#include <tuple>
#include <string>

namespace gridfire::solver {
    /**
     * @struct SolverContextBase
     * @brief Base class for solver callback contexts.
     *
     * This struct serves as a base class for contexts that can be passed to solver callbacks, it enforces
     * that derived classes implement a `describe` method that returns a vector of tuples describing
     * the context that a callback will receive when called.
     */
    struct SolverContextBase {
        virtual ~SolverContextBase() = default;

        /**
         * @brief Describe the context for callback functions.
         * @return A vector of tuples, each containing a string for the parameters name and a string for its type.
         *
         * This method should be overridden by derived classes to provide a description of the context
         * that will be passed to the callback function. The intent of this method is that an end user can investigate
         * the context that will be passed to the callback function, and use this information to craft their own
         * callback function.
         */
        [[nodiscard]] virtual std::vector<std::tuple<std::string, std::string>> describe() const = 0;
    };
    /**
     * @class NetworkSolverStrategy
     * @brief Abstract base class for network solver strategies.
     *
     * This class defines the interface for network solver strategies, which are responsible
     * for integrating the reaction network ODEs over a given timestep.  It is templated on the
     * engine type to allow for different engine implementations to be used with the same solver.
     *
     * @tparam EngineT The type of engine to use with this solver strategy.  Must inherit from Engine.
     */
    template <typename EngineT>
    class NetworkSolverStrategy {
    public:
        /**
         * @brief Constructor for the NetworkSolverStrategy.
         * @param engine The engine to use for evaluating the network.
         */
        explicit NetworkSolverStrategy(EngineT& engine) : m_engine(engine) {};

        /**
         * @brief Virtual destructor.
         */
        virtual ~NetworkSolverStrategy() = default;

        /**
         * @brief Evaluates the network for a given timestep.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         */
        virtual NetOut evaluate(const NetIn& netIn) = 0;

        /**
         * @brief set the callback function to be called at the end of each timestep.
         *
         * This function allows the user to set a callback function that will be called at the end of each timestep.
         * The callback function will receive a gridfire::solver::<SOMESOLVER>::TimestepContext object. Note that
         * depending on the solver, this context may contain different information. Further, the exact
         * signature of the callback function is left up to each solver. Every solver should provide a type or type alias
         * TimestepCallback that defines the signature of the callback function so that the user can easily
         * get that type information.
         *
         * @param callback The callback function to be called at the end of each timestep.
         */
        virtual void set_callback(const std::any& callback) = 0;

        /**
         * @brief Describe the context that will be passed to the callback function.
         * @return A vector of tuples, each containing a string for the parameter's name and a string for its type.
         *
         * This method should be overridden by derived classes to provide a description of the context
         * that will be passed to the callback function. The intent of this method is that an end user can investigate
         * the context that will be passed to the callback function, and use this information to craft their own
         * callback function.
         */
        [[nodiscard]] virtual std::vector<std::tuple<std::string, std::string>> describe_callback_context() const = 0;
    protected:
        EngineT& m_engine; ///< The engine used by this solver strategy.
    };

    /**
     * @brief Type alias for a network solver strategy that uses a DynamicEngine.
     */
    using DynamicNetworkSolverStrategy = NetworkSolverStrategy<DynamicEngine>;

    /**
     * @class DirectNetworkSolver
     * @brief A network solver that directly integrates the reaction network ODEs.
     *
     * This solver uses a Runge-Kutta method to directly integrate the reaction network
     * ODEs. It is simpler than the QSENetworkSolver, but it can be less efficient for
     * stiff networks with disparate timescales.
     *
     * @implements DynamicNetworkSolverStrategy
     */
    class DirectNetworkSolver final : public DynamicNetworkSolverStrategy {
    public:
        /**
         * @brief Constructor for the DirectNetworkSolver.
         * @param engine The dynamic engine to use for evaluating the network.
         */
        using DynamicNetworkSolverStrategy::DynamicNetworkSolverStrategy;

        /**
         * @struct TimestepContext
         * @brief Context for the timestep callback function for the DirectNetworkSolver.
         *
         * This struct contains the context that will be passed to the callback function at the end of each timestep.
         * It includes the current time, state, timestep size, cached results, and other relevant information.
         *
         * This type should be used when defining a callback function
         *
         * **Example:**
         * @code
         *  #include "gridfire/solver/solver.h"
         *
         *  #include <ofstream>
         *  #include <ranges>
         *
         *  static std::ofstream consumptionFile("consumption.txt");
         *  void callback(const gridfire::solver::DirectNetworkSolver::TimestepContext& context) {
         *      int H1Index = context.engine.getSpeciesIndex(fourdst::atomic::H_1);
         *      int He4Index = context.engine.getSpeciesIndex(fourdst::atomic::He_4);
         *
     *          consumptionFile << context.t << "," << context.state(H1Index) << "," << context.state(He4Index) << "\n";
         *  }
         *
         *  int main() {
         *      ... // Code to set up engine and solvers...
         *      solver.set_callback(callback);
         *      solver.evaluate(netIn);
         *      consumptionFile.close();
         *  }
         * @endcode
         */
        struct TimestepContext final : public SolverContextBase {
            const double t; ///< Current time.
            const boost::numeric::ublas::vector<double>& state; ///< Current state of the system.
            const double dt; ///< Time step size.
            const double cached_time; ///< Cached time for the last observed state.
            const double last_observed_time; ///< Last time the state was observed.
            const double last_step_time; ///< Last step time.
            const double T9; ///< Temperature in units of 10^9 K.
            const double rho; ///< Density in g/cm^3.
            const std::optional<StepDerivatives<double>>& cached_result; ///< Cached result of the step derivatives.
            const int num_steps; ///< Total number of steps taken.
            const DynamicEngine& engine; ///< Reference to the dynamic engine.
            const std::vector<fourdst::atomic::Species>& networkSpecies;

            TimestepContext(
                const double t,
                const boost::numeric::ublas::vector<double> &state,
                const double dt,
                const double cached_time,
                const double last_observed_time,
                const double last_step_time,
                const double t9,
                const double rho,
                const std::optional<StepDerivatives<double>> &cached_result,
                const int num_steps,
                const DynamicEngine &engine,
                const std::vector<fourdst::atomic::Species>& networkSpecies
            );

            /**
             * @brief Describe the context for callback functions.
             * @return A vector of tuples, each containing a string for the parameter's name and a string for its type.
             *
             * This method provides a description of the context that will be passed to the callback function.
             * The intent is that an end user can investigate the context and use this information to craft their own
             * callback function.
             *
             * @implements SolverContextBase::describe
             */
            [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;
        };

        /**
         * @brief Type alias for a timestep callback function.
         *
         * @brief The type alias for the callback function that will be called at the end of each timestep.
         *
         */
        using TimestepCallback = std::function<void(const TimestepContext& context)>; ///< Type alias for a timestep callback function.

        /**
         * @brief Evaluates the network for a given timestep using direct integration.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         */
        NetOut evaluate(const NetIn& netIn) override;

        /**
         * @brief Sets the callback function to be called at the end of each timestep.
         * @param callback The callback function to be called at the end of each timestep.
         *
         * This function allows the user to set a callback function that will be called at the end of each timestep.
         * The callback function will receive a gridfire::solver::DirectNetworkSolver::TimestepContext object.
         */
        void set_callback(const std::any &callback) override;

        /**
         * @brief Describe the context that will be passed to the callback function.
         * @return A vector of tuples, each containing a string for the parameter's name and a string for its type.
         *
         * This method provides a description of the context that will be passed to the callback function.
         * The intent is that an end user can investigate the context and use this information to craft their own
         * callback function.
         *
         * @implements SolverContextBase::describe
         */
        [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe_callback_context() const override;


    private:
        /**
         * @struct RHSManager
         * @brief Functor for calculating the right-hand side of the ODEs.
         *
         * This functor is used by the ODE solver to calculate the time derivatives of the
         * species abundances.  It takes the current abundances as input and returns the
         * time derivatives.
         */
        struct RHSManager {
            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.

            mutable double m_cached_time;
            mutable std::optional<StepDerivatives<double>> m_cached_result;

            mutable double m_last_observed_time = 0.0; ///< Last time the state was observed.


            quill::Logger* m_logger = LogManager::getInstance().newFileLogger("integration.log", "GridFire"); ///< Logger instance.
            mutable int m_num_steps = 0;
            mutable double m_last_step_time = 1e-20;

            TimestepCallback& m_callback;
            const std::vector<fourdst::atomic::Species>& m_networkSpecies;

            /**
             * @brief Constructor for the RHSFunctor.
             * @param engine The engine used to evaluate the network.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             * @param callback callback function to be called at the end of each timestep.
             * @param networkSpecies vector of species in the network in the correct order.
             */
            RHSManager(
                DynamicEngine& engine,
                const double T9,
                const double rho,
                TimestepCallback& callback,
                const std::vector<fourdst::atomic::Species>& networkSpecies
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho),
            m_cached_time(0),
            m_callback(callback),
            m_networkSpecies(networkSpecies){}

            /**
             * @brief Calculates the time derivatives of the species abundances.
             * @param Y Vector of current abundances.
             * @param dYdt Vector to store the time derivatives.
             * @param t Current time.
             */
            void operator()(
                const boost::numeric::ublas::vector<double>& Y,
                boost::numeric::ublas::vector<double>& dYdt,
                double t
            ) const;

            void observe(const boost::numeric::ublas::vector<double>& state, double t) const;
            void compute_and_cache(const boost::numeric::ublas::vector<double>& state, double t) const;

        };

        /**
         * @struct JacobianFunctor
         * @brief Functor for calculating the Jacobian matrix.
         *
         * This functor is used by the ODE solver to calculate the Jacobian matrix of the
         * ODEs. It takes the current abundances as input and returns the Jacobian matrix.
         */
        struct JacobianFunctor {
            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.

            /**
             * @brief Constructor for the JacobianFunctor.
             * @param engine The engine used to evaluate the network.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            JacobianFunctor(
                DynamicEngine& engine,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho) {}

            /**
             * @brief Calculates the Jacobian matrix.
             * @param Y Vector of current abundances.
             * @param J Matrix to store the Jacobian matrix.
             * @param t Current time.
             * @param dfdt Vector to store the time derivatives (not used).
             */
            void operator()(
                const boost::numeric::ublas::vector<double>& Y,
                boost::numeric::ublas::matrix<double>& J,
                double t,
                boost::numeric::ublas::vector<double>& dfdt
            ) const;

        };

    private:
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log"); ///< Logger instance.
        Config& m_config = Config::getInstance(); ///< Configuration instance.

        TimestepCallback m_callback;
    };
}