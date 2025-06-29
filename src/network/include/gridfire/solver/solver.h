#pragma once

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/engine_adaptive.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "quill/Logger.h"

#include "Eigen/Dense"

#include <vector>

namespace gridfire::solver {

    /**
     * @struct dynamicQSESpeciesIndices
     * @brief Structure to hold indices of dynamic and QSE species.
     *
     * This structure is used by the QSENetworkSolver to store the indices of species
     * that are treated dynamically and those that are assumed to be in Quasi-Steady-State
     * Equilibrium (QSE).
     */
    struct dynamicQSESpeciesIndices {
        std::vector<size_t> dynamicSpeciesIndices; ///< Indices of slow species that are not in QSE.
        std::vector<size_t> QSESpeciesIndices;  ///< Indices of fast species that are in QSE.
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
    protected:
        EngineT& m_engine; ///< The engine used by this solver strategy.
    };

    /**
     * @brief Type alias for a network solver strategy that uses a DynamicEngine.
     */
    using DynamicNetworkSolverStrategy = NetworkSolverStrategy<DynamicEngine>;

    /**
     * @brief Type alias for a network solver strategy that uses an AdaptiveEngineView.
     */
    using AdaptiveNetworkSolverStrategy = NetworkSolverStrategy<AdaptiveEngineView>;

    /**
     * @brief Type alias for a network solver strategy that uses a static Engine.
     */
    using StaticNetworkSolverStrategy = NetworkSolverStrategy<Engine>;

    /**
     * @class QSENetworkSolver
     * @brief A network solver that uses a Quasi-Steady-State Equilibrium (QSE) approach.
     *
     * This solver partitions the network into "fast" species in QSE and "slow" (dynamic) species.
     * The abundances of the fast species are determined by solving a system of algebraic
     * equations, while the abundances of the slow species are integrated using an ODE solver.
     * This hybrid approach is highly effective for stiff networks with disparate timescales.
     *
     * The QSE solver uses an AdaptiveEngineView to dynamically cull unimportant species and
     * reactions, which significantly improves performance for large networks.
     *
     * @implements AdaptiveNetworkSolverStrategy
     *
     * @see AdaptiveEngineView
     * @see DynamicEngine::getSpeciesTimescales()
     */
    class QSENetworkSolver final : public AdaptiveNetworkSolverStrategy {
    public:
        /**
         * @brief Constructor for the QSENetworkSolver.
         * @param engine The adaptive engine view to use for evaluating the network.
         */
        using AdaptiveNetworkSolverStrategy::AdaptiveNetworkSolverStrategy;

        /**
         * @brief Evaluates the network for a given timestep using the QSE approach.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the timestep.
         *
         * This method performs the following steps:
         *   1. Updates the adaptive engine view (if necessary).
         *   2. Partitions the species into dynamic and QSE species based on their timescales.
         *   3. Calculates the steady-state abundances of the QSE species.
         *   4. Integrates the ODEs for the dynamic species using a Runge-Kutta solver.
         *   5. Marshals the output variables into a NetOut struct.
         *
         * @throws std::runtime_error If the steady-state abundances cannot be calculated.
         *
         * @see AdaptiveEngineView::update()
         * @see packSpeciesTypeIndexVectors()
         * @see calculateSteadyStateAbundances()
         */
        NetOut evaluate(const NetIn& netIn) override;
    private: // methods
        /**
         * @brief Packs the species indices into vectors based on their type (dynamic or QSE).
         * @param Y Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @return A dynamicQSESpeciesIndices struct containing the indices of the dynamic and QSE species.
         *
         * This method determines whether each species should be treated dynamically or as
         * being in QSE based on its timescale and abundance.  Species with short timescales
         * or low abundances are assumed to be in QSE.
         *
         * @see DynamicEngine::getSpeciesTimescales()
         */
        dynamicQSESpeciesIndices packSpeciesTypeIndexVectors(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const;

        /**
         * @brief Calculates the steady-state abundances of the QSE species.
         * @param Y Vector of current abundances for all species.
         * @param T9 Temperature in units of 10^9 K.
         * @param rho Density in g/cm^3.
         * @param indices A dynamicQSESpeciesIndices struct containing the indices of the dynamic and QSE species.
         * @return An Eigen::VectorXd containing the steady-state abundances of the QSE species.
         *
         * This method solves a system of algebraic equations to determine the steady-state
         * abundances of the QSE species.
         *
         * @throws std::runtime_error If the steady-state abundances cannot be calculated.
         */
        Eigen::VectorXd calculateSteadyStateAbundances(
            const std::vector<double>& Y,
            const double T9,
            const double rho,
            const dynamicQSESpeciesIndices& indices
        ) const;

        /**
         * @brief Initializes the network with a short ignition phase.
         * @param netIn The input conditions for the network.
         * @return The output conditions after the ignition phase.
         *
         * This method performs a short integration of the network at a high temperature and
         * density to ignite the network and bring it closer to equilibrium.  This can improve
         * the convergence of the QSE solver.
         *
         * @see DirectNetworkSolver::evaluate()
         */
        NetOut initializeNetworkWithShortIgnition(
            const NetIn& netIn
        ) const;

        /**
         * @brief Determines whether the adaptive engine view should be updated.
         * @param conditions The current input conditions.
         * @return True if the view should be updated, false otherwise.
         *
         * This method implements a policy for determining when the adaptive engine view
         * should be updated.  The view is updated if the temperature or density has changed
         * significantly, or if a primary fuel source has been depleted.
         *
         * @see AdaptiveEngineView::update()
         */
        bool shouldUpdateView(const NetIn& conditions) const;
    private: // Nested functors for ODE integration
        /**
         * @struct RHSFunctor
         * @brief Functor for calculating the right-hand side of the ODEs for the dynamic species.
         *
         * This functor is used by the ODE solver to calculate the time derivatives of the
         * dynamic species.  It takes the current abundances of the dynamic species as input
         * and returns the time derivatives of those abundances.
         */
        struct RHSFunctor {
            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const std::vector<size_t>& m_dynamicSpeciesIndices; ///< Indices of the dynamic species.
            const std::vector<size_t>& m_QSESpeciesIndices; ///< Indices of the QSE species.
            const Eigen::VectorXd& m_Y_QSE; ///< Steady-state abundances of the QSE species.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.

            bool m_isViewInitialized = false;

            /**
             * @brief Constructor for the RHSFunctor.
             * @param engine The engine used to evaluate the network.
             * @param dynamicSpeciesIndices Indices of the dynamic species.
             * @param QSESpeciesIndices Indices of the QSE species.
             * @param Y_QSE Steady-state abundances of the QSE species.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            RHSFunctor(
                DynamicEngine& engine,
                const std::vector<size_t>& dynamicSpeciesIndices,
                const std::vector<size_t>& QSESpeciesIndices,
                const Eigen::VectorXd& Y_QSE,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_dynamicSpeciesIndices(dynamicSpeciesIndices),
            m_QSESpeciesIndices(QSESpeciesIndices),
            m_Y_QSE(Y_QSE),
            m_T9(T9),
            m_rho(rho) {}

            /**
             * @brief Calculates the time derivatives of the dynamic species.
             * @param YDynamic Vector of current abundances for the dynamic species.
             * @param dYdtDynamic Vector to store the time derivatives of the dynamic species.
             * @param t Current time.
             */
            void operator()(
                const boost::numeric::ublas::vector<double>& YDynamic,
                boost::numeric::ublas::vector<double>& dYdtDynamic,
                double t
            ) const;

        };

        /**
         * @struct JacobianFunctor
         * @brief Functor for calculating the Jacobian matrix of the ODEs for the dynamic species.
         *
         * This functor is used by the ODE solver to calculate the Jacobian matrix of the
         * ODEs for the dynamic species.  It takes the current abundances of the dynamic
         * species as input and returns the Jacobian matrix.
         *
         * @todo Implement the Jacobian functor.
         */
        struct JacobianFunctor {
            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const std::vector<size_t>& m_dynamicSpeciesIndices; ///< Indices of the dynamic species.
            const std::vector<size_t>& m_QSESpeciesIndices; ///< Indices of the QSE species.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.

            /**
             * @brief Constructor for the JacobianFunctor.
             * @param engine The engine used to evaluate the network.
             * @param dynamicSpeciesIndices Indices of the dynamic species.
             * @param QSESpeciesIndices Indices of the QSE species.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            JacobianFunctor(
                DynamicEngine& engine,
                const std::vector<size_t>& dynamicSpeciesIndices,
                const std::vector<size_t>& QSESpeciesIndices,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_dynamicSpeciesIndices(dynamicSpeciesIndices),
            m_QSESpeciesIndices(QSESpeciesIndices),
            m_T9(T9),
            m_rho(rho) {}

            /**
             * @brief Calculates the Jacobian matrix of the ODEs for the dynamic species.
             * @param YDynamic Vector of current abundances for the dynamic species.
             * @param JDynamic Matrix to store the Jacobian matrix.
             * @param t Current time.
             * @param dfdt Vector to store the time derivatives of the dynamic species (not used).
             */
            void operator()(
                const boost::numeric::ublas::vector<double>& YDynamic,
                boost::numeric::ublas::matrix<double>& JDynamic,
                double t,
                boost::numeric::ublas::vector<double>& dfdt
            ) const;
        };

        /**
         * @struct EigenFunctor
         * @brief Functor for calculating the residual and Jacobian for the QSE species using Eigen.
         *
         * @tparam T The numeric type to use for the calculation (double or ADDouble).
         */
        template<typename T>
        struct EigenFunctor {
            using InputType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            using OutputType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            using JacobianType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const std::vector<double>& m_YDynamic; ///< Abundances of the dynamic species.
            const std::vector<size_t>& m_dynamicSpeciesIndices; ///< Indices of the dynamic species.
            const std::vector<size_t>& m_QSESpeciesIndices; ///< Indices of the QSE species.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.

            /**
             * @brief Constructor for the EigenFunctor.
             * @param engine The engine used to evaluate the network.
             * @param YDynamic Abundances of the dynamic species.
             * @param dynamicSpeciesIndices Indices of the dynamic species.
             * @param QSESpeciesIndices Indices of the QSE species.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            EigenFunctor(
                DynamicEngine& engine,
                const std::vector<double>& YDynamic,
                const std::vector<size_t>& dynamicSpeciesIndices,
                const std::vector<size_t>& QSESpeciesIndices,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_YDynamic(YDynamic),
            m_dynamicSpeciesIndices(dynamicSpeciesIndices),
            m_QSESpeciesIndices(QSESpeciesIndices),
            m_T9(T9),
            m_rho(rho) {}

            /**
             * @brief Calculates the residual vector for the QSE species.
             * @param v_QSE Input vector of QSE species abundances (logarithmic).
             * @param f_QSE Output vector of residuals.
             * @return 0 for success.
             */
            int operator()(const InputType& v_QSE, OutputType& f_QSE) const;

            /**
             * @brief Calculates the Jacobian matrix for the QSE species.
             * @param v_QSE Input vector of QSE species abundances (logarithmic).
             * @param J_QSE Output Jacobian matrix.
             * @return 0 for success.
             */
            int df(const InputType& v_QSE, JacobianType& J_QSE) const;
        };
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log"); ///< Logger instance.
        fourdst::config::Config& m_config = fourdst::config::Config::getInstance(); ///< Configuration instance.

        bool m_isViewInitialized = false; ///< Flag indicating whether the adaptive engine view has been initialized.
        NetIn m_lastSeenConditions; ///< The last seen input conditions.
    };

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
        struct RHSFunctor {
            DynamicEngine& m_engine; ///< The engine used to evaluate the network.
            const double m_T9; ///< Temperature in units of 10^9 K.
            const double m_rho; ///< Density in g/cm^3.
            const size_t m_numSpecies; ///< The number of species in the network.

            /**
             * @brief Constructor for the RHSFunctor.
             * @param engine The engine used to evaluate the network.
             * @param T9 Temperature in units of 10^9 K.
             * @param rho Density in g/cm^3.
             */
            RHSFunctor(
                DynamicEngine& engine,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho),
            m_numSpecies(engine.getNetworkSpecies().size()) {}

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
            const size_t m_numSpecies; ///< The number of species in the network.

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
            m_rho(rho),
            m_numSpecies(engine.getNetworkSpecies().size()) {}

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
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log"); ///< Logger instance.
        fourdst::config::Config& m_config = fourdst::config::Config::getInstance(); ///< Configuration instance.
    };

    template<typename T>
    int QSENetworkSolver::EigenFunctor<T>::operator()(const InputType &v_QSE, OutputType &f_QSE) const {
        std::vector<double> YFull(m_engine.getNetworkSpecies().size(), 0.0);
        Eigen::VectorXd Y_QSE(v_QSE.array().exp());
        for (size_t i = 0; i < m_dynamicSpeciesIndices.size(); ++i) {
            YFull[m_dynamicSpeciesIndices[i]] = m_YDynamic[i];
        }
        for (size_t i = 0; i < m_QSESpeciesIndices.size(); ++i) {
            YFull[m_QSESpeciesIndices[i]] = Y_QSE(i);
        }
        const auto [full_dYdt, specificEnergyGenerationRate] = m_engine.calculateRHSAndEnergy(YFull, m_T9, m_rho);
        f_QSE.resize(m_QSESpeciesIndices.size());
        for (size_t i = 0; i < m_QSESpeciesIndices.size(); ++i) {
            f_QSE(i) = full_dYdt[m_QSESpeciesIndices[i]];
        }

        return 0; // Success
    }

    template <typename T>
    int QSENetworkSolver::EigenFunctor<T>::df(const InputType& v_QSE, JacobianType& J_QSE) const {
        std::vector<double> YFull(m_engine.getNetworkSpecies().size(), 0.0);
        Eigen::VectorXd Y_QSE(v_QSE.array().exp());
        for (size_t i = 0; i < m_dynamicSpeciesIndices.size(); ++i) {
            YFull[m_dynamicSpeciesIndices[i]] = m_YDynamic[i];
        }
        for (size_t i = 0; i < m_QSESpeciesIndices.size(); ++i) {
            YFull[m_QSESpeciesIndices[i]] = Y_QSE(i);
        }

        m_engine.generateJacobianMatrix(YFull, m_T9, m_rho);

        Eigen::MatrixXd J_orig(m_QSESpeciesIndices.size(), m_QSESpeciesIndices.size());
        for (size_t i = 0; i < m_QSESpeciesIndices.size(); ++i) {
            for (size_t j = 0; j < m_QSESpeciesIndices.size(); ++j) {
                J_orig(i, j) = m_engine.getJacobianMatrixEntry(m_QSESpeciesIndices[i], m_QSESpeciesIndices[j]);
            }
        }

        J_QSE = J_orig;
        for (long j = 0; j < J_QSE.cols(); ++j) {
            J_QSE.col(j) *= Y_QSE(j); // Chain rule for log space
        }
        return 0; // Success
    }


}