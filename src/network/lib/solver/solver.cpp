#include "gridfire/solver/solver.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/network.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"
#include "fourdst/config/config.h"

#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"

#include <boost/numeric/odeint.hpp>

#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>

#include "quill/LogMacros.h"

namespace gridfire::solver {

    NetOut QSENetworkSolver::evaluate(const NetIn &netIn) {
        using state_type = boost::numeric::ublas::vector<double>;
        using namespace boost::numeric::odeint;

        NetOut postIgnition = initializeNetworkWithShortIgnition(netIn);

        constexpr double abundance_floor = 1.0e-30;
        std::vector<double>Y_sanitized_initial;
        Y_sanitized_initial.reserve(m_engine.getNetworkSpecies().size());

        LOG_DEBUG(m_logger, "Sanitizing initial abundances with a floor of {:0.3E}...", abundance_floor);
        for (const auto& species : m_engine.getNetworkSpecies()) {
            double molar_abundance = 0.0;
            if (postIgnition.composition.contains(species)) {
                molar_abundance = postIgnition.composition.getMolarAbundance(std::string(species.name()));
            }

            if (molar_abundance < abundance_floor) {
                molar_abundance = abundance_floor;
            }
            Y_sanitized_initial.push_back(molar_abundance);
        }
        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const double rho = netIn.density; // Density in g/cm^3

        const auto indices = packSpeciesTypeIndexVectors(Y_sanitized_initial, T9, rho);
        Eigen::VectorXd Y_QSE;
        try {
            Y_QSE = calculateSteadyStateAbundances(Y_sanitized_initial, T9, rho, indices);
        } catch (const std::runtime_error& e) {
            LOG_ERROR(m_logger, "Failed to calculate steady state abundances. Aborting QSE evaluation.");
            m_logger->flush_log();
            throw std::runtime_error("Failed to calculate steady state abundances: " + std::string(e.what()));
        }

        state_type YDynamic_ublas(indices.dynamicSpeciesIndices.size() + 1);
        for (size_t i = 0; i < indices.dynamicSpeciesIndices.size(); ++i) {
            YDynamic_ublas(i) = Y_sanitized_initial[indices.dynamicSpeciesIndices[i]];
        }
        YDynamic_ublas(indices.dynamicSpeciesIndices.size()) = 0.0; // Placeholder for specific energy rate

        const RHSFunctor rhs_functor(m_engine, indices.dynamicSpeciesIndices, indices.QSESpeciesIndices, Y_QSE, T9, rho);
        const auto stepper = make_controlled<runge_kutta_dopri5<state_type>>(1.0e-8, 1.0e-8);

        size_t stepCount = integrate_adaptive(
            stepper,
            rhs_functor,
            YDynamic_ublas,
            0.0, // Start time
            netIn.tMax,
            netIn.dt0
        );

        std::vector<double> YFinal = Y_sanitized_initial;
        for (size_t i = 0; i <indices.dynamicSpeciesIndices.size(); ++i) {
            YFinal[indices.dynamicSpeciesIndices[i]] = YDynamic_ublas(i);
        }
        for (size_t i = 0; i < indices.QSESpeciesIndices.size(); ++i) {
            YFinal[indices.QSESpeciesIndices[i]] = Y_QSE(i);
        }

        const double finalSpecificEnergyRate = YDynamic_ublas(indices.dynamicSpeciesIndices.size());

        // --- Marshal output variables ---
        std::vector<std::string> speciesNames(m_engine.getNetworkSpecies().size());
        std::vector<double> finalMassFractions(m_engine.getNetworkSpecies().size());
        double massFractionSum = 0.0;

        for (size_t i = 0; i < speciesNames.size(); ++i) {
            const auto& species = m_engine.getNetworkSpecies()[i];
            speciesNames[i] = species.name();
            finalMassFractions[i] = YFinal[i] * species.mass(); // Convert from molar abundance to mass fraction
            massFractionSum += finalMassFractions[i];
        }
        for (auto& mf : finalMassFractions) {
            mf /= massFractionSum; // Normalize to get mass fractions
        }

        fourdst::composition::Composition outputComposition(speciesNames, finalMassFractions);
        NetOut netOut;
        netOut.composition = outputComposition;
        netOut.energy = finalSpecificEnergyRate; // Specific energy rate
        netOut.num_steps = stepCount;
        return netOut;
    }

    dynamicQSESpeciesIndices QSENetworkSolver::packSpeciesTypeIndexVectors(
        const std::vector<double>& Y,
        const double T9,
        const double rho
    ) const {
        constexpr double timescaleCutoff = 1.0e-5;
        constexpr double abundanceCutoff = 1.0e-15;

        LOG_INFO(m_logger, "Partitioning species using T9={:0.2f} and ρ={:0.2e}", T9, rho);
        LOG_INFO(m_logger, "Timescale Cutoff: {:.1e} s, Abundance Cutoff: {:.1e}", timescaleCutoff, abundanceCutoff);

        std::vector<size_t>dynamicSpeciesIndices; // Slow species that are not in QSE
        std::vector<size_t>QSESpeciesIndices;  // Fast species that are in QSE

        std::unordered_map<fourdst::atomic::Species, double> speciesTimescale = m_engine.getSpeciesTimescales(Y, T9, rho);

        for (size_t i = 0; i < m_engine.getNetworkSpecies().size(); ++i) {
            const auto& species = m_engine.getNetworkSpecies()[i];
            const double timescale = speciesTimescale[species];
            const double abundance = Y[i];

            if (std::isinf(timescale) || abundance < abundanceCutoff || timescale <= timescaleCutoff) {
                QSESpeciesIndices.push_back(i);
            } else {
                dynamicSpeciesIndices.push_back(i);
            }
        }
        LOG_INFO(m_logger, "Partitioning complete. Dynamical species: {}, QSE species: {}.", dynamicSpeciesIndices.size(), QSESpeciesIndices.size());
        LOG_INFO(m_logger, "Dynamic species: {}", [dynamicSpeciesIndices](const DynamicEngine& engine_wrapper) -> std::string {
            std::string result;
            int count = 0;
            for (const auto& i : dynamicSpeciesIndices) {
                result += std::string(engine_wrapper.getNetworkSpecies()[i].name());
                if (count < dynamicSpeciesIndices.size() - 2) {
                    result += ", ";
                } else if (count == dynamicSpeciesIndices.size() - 2) {
                    result += " and ";
                }
                count++;
            }
            return result;
        }(m_engine));
        LOG_INFO(m_logger, "QSE species: {}", [QSESpeciesIndices](const DynamicEngine& engine_wrapper) -> std::string {
            std::string result;
            int count = 0;
            for (const auto& i : QSESpeciesIndices) {
                result += std::string(engine_wrapper.getNetworkSpecies()[i].name());
                if (count < QSESpeciesIndices.size() - 2) {
                    result += ", ";
                } else if (count == QSESpeciesIndices.size() - 2) {
                    result += " and ";
                }
                count++;
            }
            return result;
        }(m_engine));
        return {dynamicSpeciesIndices, QSESpeciesIndices};
    }

    Eigen::VectorXd QSENetworkSolver::calculateSteadyStateAbundances(
        const std::vector<double> &Y,
        const double T9,
        const double rho,
        const dynamicQSESpeciesIndices &indices
    ) const {
        std::vector<double> Y_dyn;
        Eigen::VectorXd Y_qse_initial(indices.QSESpeciesIndices.size());
        for (const auto& i : indices.dynamicSpeciesIndices) { Y_dyn.push_back(Y[i]); }
        for (size_t i = 0; i < indices.QSESpeciesIndices.size(); ++i) {
            Y_qse_initial(i) = Y[indices.QSESpeciesIndices[i]];
            if (Y_qse_initial(i) < 1.0e-99) { Y_qse_initial(i) = 1.0e-99; }
        }

        Eigen::VectorXd v_qse = Y_qse_initial.array().log();

        EigenFunctor<double> qse_problem(m_engine, Y_dyn, indices.dynamicSpeciesIndices, indices.QSESpeciesIndices, T9, rho);
        LOG_INFO(m_logger, "--- QSE Pre-Solve Diagnostics ---");

        Eigen::VectorXd f_initial(indices.QSESpeciesIndices.size());
        qse_problem(v_qse, f_initial);
        LOG_INFO(m_logger, "Initial Guess ||f||: {:0.4e}", f_initial.norm());

        Eigen::MatrixXd J_initial(indices.QSESpeciesIndices.size(), indices.QSESpeciesIndices.size());
        qse_problem.df(v_qse, J_initial);
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(J_initial);
        double cond = svd.singularValues().maxCoeff() / svd.singularValues().minCoeff();
        LOG_INFO(m_logger, "Initial Jacobian Condition Number: {:0.4e}", cond);
        LOG_INFO(m_logger, "Starting QSE solve...");

        Eigen::HybridNonLinearSolver<EigenFunctor<double>> solver(qse_problem);
        solver.parameters.xtol = 1.0e-8; // Set tolerance

        // 5. Run the solver. It will modify v_qse in place.
        const int eigenStatus = solver.solve(v_qse);

        // 6. Check for convergence and return the result
        if(eigenStatus != Eigen::Success) {
            LOG_WARNING(m_logger, "--- QSE SOLVER FAILED ---");
            LOG_WARNING(m_logger, "Eigen status code: {}", eigenStatus);
            LOG_WARNING(m_logger, "Iterations performed: {}", solver.iter);

            // Log the final state that caused the failure
            Eigen::VectorXd Y_qse_final_fail = v_qse.array().exp();
            for(long i=0; i<v_qse.size(); ++i) {
                LOG_WARNING(m_logger, "Final v_qse[{}]: {:0.4e} -> Y_qse[{}]: {:0.4e}", i, v_qse(i), i, Y_qse_final_fail(i));
            }

            // Log the residual at the final state
            Eigen::VectorXd f_final(indices.QSESpeciesIndices.size());
            qse_problem(v_qse, f_final);
            LOG_WARNING(m_logger, "Final ||f||: {:0.4e}", f_final.norm());

            throw std::runtime_error("Eigen QSE solver did not converge.");
        }
        LOG_INFO(m_logger, "Eigen QSE solver converged in {} iterations.", solver.iter);

        return v_qse.array().exp();

    }

    NetOut QSENetworkSolver::initializeNetworkWithShortIgnition(const NetIn &netIn) const {
        const auto ignitionTemperature = m_config.get<double>(
            "gridfire:solver:QSE:ignition:temperature",
            2e8
        ); // 0.2 GK
        const auto ignitionDensity = m_config.get<double>(
            "gridfire:solver:QSE:ignition:density",
            1e6
        ); // 1e6 g/cm^3
        const auto ignitionTime = m_config.get<double>(
            "gridfire:solver:QSE:ignition:tMax",
            1e-7
        ); // 0.1 μs
        const auto ignitionStepSize = m_config.get<double>(
            "gridfire:solver:QSE:ignition:dt0",
            1e-15
        ); // 1e-15 seconds

        LOG_INFO(
            m_logger,
            "Igniting network with T={:<5.3E}, ρ={:<5.3E}, tMax={:<5.3E}, dt0={:<5.3E}...",
            ignitionTemperature,
            ignitionDensity,
            ignitionTime,
            ignitionStepSize
        );

        NetIn preIgnition = netIn;
        preIgnition.temperature = ignitionTemperature;
        preIgnition.density = ignitionDensity;
        preIgnition.tMax = ignitionTime;
        preIgnition.dt0 = ignitionStepSize;

        DirectNetworkSolver ignitionSolver(m_engine);
        NetOut postIgnition = ignitionSolver.evaluate(preIgnition);
        LOG_INFO(m_logger, "Network ignition completed in {} steps.", postIgnition.num_steps);
        return postIgnition;
    }

    void QSENetworkSolver::RHSFunctor::operator()(
        const boost::numeric::ublas::vector<double> &YDynamic,
        boost::numeric::ublas::vector<double> &dYdtDynamic,
        double t
    ) const {
        // --- Populate the slow / dynamic species vector ---
        std::vector<double> YFull(m_engine.getNetworkSpecies().size());
        for (size_t i = 0; i < m_dynamicSpeciesIndices.size(); ++i) {
            YFull[m_dynamicSpeciesIndices[i]] = YDynamic(i);
        }

        // --- Populate the QSE species vector ---
        for (size_t i = 0; i < m_QSESpeciesIndices.size(); ++i) {
            YFull[m_QSESpeciesIndices[i]] = m_Y_QSE(i);
        }

        auto [full_dYdt, specificEnergyRate] = m_engine.calculateRHSAndEnergy(YFull, m_T9, m_rho);

        dYdtDynamic.resize(m_dynamicSpeciesIndices.size() + 1);
        for (size_t i = 0; i < m_dynamicSpeciesIndices.size(); ++i) {
            dYdtDynamic(i) = full_dYdt[m_dynamicSpeciesIndices[i]];
        }
        dYdtDynamic[m_dynamicSpeciesIndices.size()] = specificEnergyRate;
    }

    NetOut DirectNetworkSolver::evaluate(const NetIn &netIn) {
        namespace ublas = boost::numeric::ublas;
        namespace odeint = boost::numeric::odeint;
        using fourdst::composition::Composition;

        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)
        const unsigned long numSpecies = m_engine.getNetworkSpecies().size();

        const auto absTol = m_config.get<double>("gridfire:solver:DirectNetworkSolver:absTol", 1.0e-8);
        const auto relTol = m_config.get<double>("gridfire:solver:DirectNetworkSolver:relTol", 1.0e-8);

        size_t stepCount = 0;

        RHSFunctor rhsFunctor(m_engine, T9, netIn.density);
        JacobianFunctor jacobianFunctor(m_engine, T9, netIn.density);

        ublas::vector<double> Y(numSpecies + 1);

        for (size_t i = 0; i < numSpecies; ++i) {
            const auto& species = m_engine.getNetworkSpecies()[i];
            try {
                Y(i) = netIn.composition.getMolarAbundance(std::string(species.name()));
            } catch (const std::runtime_error) {
                LOG_DEBUG(m_logger, "Species '{}' not found in composition. Setting abundance to 0.0.", species.name());
                Y(i) = 0.0;
            }
        }
        Y(numSpecies) = 0.0;

        const auto stepper = odeint::make_controlled<odeint::rosenbrock4<double>>(absTol, relTol);
        stepCount = odeint::integrate_adaptive(
            stepper,
            std::make_pair(rhsFunctor, jacobianFunctor),
            Y,
            0.0,
            netIn.tMax,
            netIn.dt0
        );

        std::vector<double> finalMassFractions(numSpecies);
        for (size_t i = 0; i < numSpecies; ++i) {
            const double molarMass = m_engine.getNetworkSpecies()[i].mass();
            finalMassFractions[i] = Y(i) * molarMass; // Convert from molar abundance to mass fraction
            if (finalMassFractions[i] < MIN_ABUNDANCE_THRESHOLD) {
                finalMassFractions[i] = 0.0;
            }
        }

        std::vector<std::string> speciesNames;
        speciesNames.reserve(numSpecies);
        for (const auto& species : m_engine.getNetworkSpecies()) {
            speciesNames.push_back(std::string(species.name()));
        }

        Composition outputComposition(speciesNames, finalMassFractions);
        outputComposition.finalize(true);

        NetOut netOut;
        netOut.composition = std::move(outputComposition);
        netOut.energy = Y(numSpecies); // Specific energy rate
        netOut.num_steps = stepCount;

        return netOut;
    }

    void DirectNetworkSolver::RHSFunctor::operator()(
        const boost::numeric::ublas::vector<double> &Y,
        boost::numeric::ublas::vector<double> &dYdt,
        double t
    ) const {
        const std::vector<double> y(Y.begin(), m_numSpecies + Y.begin());
        auto [dydt, eps] = m_engine.calculateRHSAndEnergy(y, m_T9, m_rho);
        dYdt.resize(m_numSpecies + 1);
        std::ranges::copy(dydt, dydt.begin());
        dYdt(m_numSpecies) = eps;
    }

    void DirectNetworkSolver::JacobianFunctor::operator()(
        const boost::numeric::ublas::vector<double> &Y,
        boost::numeric::ublas::matrix<double> &J,
        double t,
        boost::numeric::ublas::vector<double> &dfdt
    ) const {
        J.resize(m_numSpecies + 1, m_numSpecies + 1);
        J.clear();
        for (int i = 0; i < m_numSpecies + 1; ++i) {
            for (int j = 0; j < m_numSpecies + 1; ++j) {
                J(i, j) = m_engine.getJacobianMatrixEntry(i, j);
            }
        }
    }
}