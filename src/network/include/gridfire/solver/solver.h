#pragma once

#include "gridfire/engine/engine_graph.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/network.h"

#include "fourdst/logging/logging.h"
#include "fourdst/config/config.h"

#include "quill/Logger.h"

#include "Eigen/Dense"

#include <vector>

namespace gridfire::solver {

    struct dynamicQSESpeciesIndices {
        std::vector<size_t> dynamicSpeciesIndices; // Slow species that are not in QSE
        std::vector<size_t> QSESpeciesIndices;  // Fast species that are in QSE
    };

    template <typename EngineT>
    class NetworkSolverStrategy {
    public:
        explicit NetworkSolverStrategy(EngineT& engine) : m_engine(engine) {};
        virtual ~NetworkSolverStrategy() = default;
        virtual NetOut evaluate(const NetIn& netIn) = 0;
    protected:
        EngineT& m_engine;
    };

    using DynamicNetworkSolverStrategy = NetworkSolverStrategy<DynamicEngine>;
    using StaticNetworkSolverStrategy = NetworkSolverStrategy<Engine>;


    class QSENetworkSolver final : public DynamicNetworkSolverStrategy {
    public:
        using DynamicNetworkSolverStrategy::DynamicNetworkSolverStrategy;
        NetOut evaluate(const NetIn& netIn) override;
    private: // methods
        dynamicQSESpeciesIndices packSpeciesTypeIndexVectors(
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const;
        Eigen::VectorXd calculateSteadyStateAbundances(
            const std::vector<double>& Y,
            const double T9,
            const double rho,
            const dynamicQSESpeciesIndices& indices
        ) const;
        NetOut initializeNetworkWithShortIgnition(
            const NetIn& netIn
        ) const;
    private: // Nested functors for ODE integration
        struct RHSFunctor {
            DynamicEngine& m_engine;
            const std::vector<size_t>& m_dynamicSpeciesIndices;
            const std::vector<size_t>& m_QSESpeciesIndices;
            const Eigen::VectorXd& m_Y_QSE;
            const double m_T9;
            const double m_rho;

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

            void operator()(
                const boost::numeric::ublas::vector<double>& YDynamic,
                boost::numeric::ublas::vector<double>& dYdtDynamic,
                double t
            ) const;

        };

        struct JacobianFunctor {
            DynamicEngine& m_engine;
            const std::vector<size_t>& m_dynamicSpeciesIndices;
            const std::vector<size_t>& m_QSESpeciesIndices;
            const double m_T9;
            const double m_rho;

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

            void operator()(
                const boost::numeric::ublas::vector<double>& YDynamic,
                boost::numeric::ublas::matrix<double>& JDynamic,
                double t,
                boost::numeric::ublas::vector<double>& dfdt
            ) const;
        };

        template<typename T>
        struct EigenFunctor {
            using InputType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            using OutputType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
            using JacobianType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

            DynamicEngine& m_engine;
            const std::vector<double>& m_YDynamic;
            const std::vector<size_t>& m_dynamicSpeciesIndices;
            const std::vector<size_t>& m_QSESpeciesIndices;
            const double m_T9;
            const double m_rho;

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

            int operator()(const InputType& v_QSE, OutputType& f_QSE) const;
            int df(const InputType& v_QSE, JacobianType& J_QSE) const;
        };
    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        fourdst::config::Config& m_config = fourdst::config::Config::getInstance();
    };

    class DirectNetworkSolver final : public DynamicNetworkSolverStrategy {
    public:
        using DynamicNetworkSolverStrategy::DynamicNetworkSolverStrategy;
        NetOut evaluate(const NetIn& netIn) override;
    private:
        struct RHSFunctor {
            DynamicEngine& m_engine;
            const double m_T9;
            const double m_rho;
            const size_t m_numSpecies;

            RHSFunctor(
                DynamicEngine& engine,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho),
            m_numSpecies(engine.getNetworkSpecies().size()) {}

            void operator()(
                const boost::numeric::ublas::vector<double>& Y,
                boost::numeric::ublas::vector<double>& dYdt,
                double t
            ) const;
        };
        struct JacobianFunctor {
            DynamicEngine& m_engine;
            const double m_T9;
            const double m_rho;
            const size_t m_numSpecies;

            JacobianFunctor(
                DynamicEngine& engine,
                const double T9,
                const double rho
            ) :
            m_engine(engine),
            m_T9(T9),
            m_rho(rho),
            m_numSpecies(engine.getNetworkSpecies().size()) {}

            void operator()(
                const boost::numeric::ublas::vector<double>& Y,
                boost::numeric::ublas::matrix<double>& J,
                double t,
                boost::numeric::ublas::vector<double>& dfdt
            ) const;

        };

    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        fourdst::config::Config& m_config = fourdst::config::Config::getInstance();
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