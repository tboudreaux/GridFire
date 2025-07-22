#pragma once

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_abstract.h"
#include "../engine/views/engine_adaptive.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "quill/Logger.h"

#include <vector>

namespace gridfire::solver {
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
         * @brief Evaluates the network for a given timestep using direct integration.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         */
        NetOut evaluate(const NetIn& netIn) override;
    private:
        /**
         * @struct RHSFunctor
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

            /**
             * @brief Constructor for the RHSFunctor.
             * @param engine The engine used to evaluate the network.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            RHSManager(
                DynamicEngine& engine,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho),
            m_cached_time(0) {}

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
    };
}