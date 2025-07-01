#include "gridfire/solver/solver.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/network.h"

#include "gridfire/utils/logging.h"

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
#include <iomanip>

#include "quill/LogMacros.h"

namespace gridfire::solver {

    NetOut QSENetworkSolver::evaluate(const NetIn &netIn) {
        // --- Use the policy to decide whether to update the view ---
        if (shouldUpdateView(netIn)) {
            LOG_DEBUG(m_logger, "Solver update policy triggered, network view updating...");
            m_engine.update(netIn);
            LOG_DEBUG(m_logger, "Network view updated!");

            m_lastSeenConditions = netIn;
            m_isViewInitialized = true;
        }
        m_engine.generateJacobianMatrix(netIn.MolarAbundance(), netIn.temperature / 1e9, netIn.density);
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
            LOG_DEBUG(m_logger, "QSE Abundances: {}", [*this](const dynamicQSESpeciesIndices& indices, const Eigen::VectorXd& Y_QSE) -> std::string {
                std::stringstream ss;
                ss << std::scientific << std::setprecision(5);
                for (size_t i = 0; i < indices.QSESpeciesIndices.size(); ++i) {
                    ss << std::string(m_engine.getNetworkSpecies()[indices.QSESpeciesIndices[i]].name()) + ": ";
                    ss << Y_QSE(i);
                    if (i < indices.QSESpeciesIndices.size() - 2) {
                        ss << ", ";
                    } else if (i == indices.QSESpeciesIndices.size() - 2) {
                        ss << ", and ";
                    }

                }
                return ss.str();
            }(indices, Y_QSE));
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
            const double network_timescale = speciesTimescale.at(species);
            const double abundance = Y[i];

            double decay_timescale = std::numeric_limits<double>::infinity();
            const double half_life = species.halfLife();
            if (half_life > 0 && !std::isinf(half_life)) {
                constexpr double LN2 = 0.69314718056;
                decay_timescale = half_life / LN2;
            }

            const double final_timescale = std::min(network_timescale, decay_timescale);

            if (std::isinf(final_timescale) || abundance < abundanceCutoff || final_timescale <= timescaleCutoff) {
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
        LOG_TRACE_L1(m_logger, "Calculating steady state abundances for QSE species...");

        if (indices.QSESpeciesIndices.empty()) {
            LOG_DEBUG(m_logger, "No QSE species to solve for.");
            return Eigen::VectorXd(0);
        }
        // Use the EigenFunctor with Eigen's nonlinear solver
        EigenFunctor<double> functor(
            m_engine,
            Y,
            indices.dynamicSpeciesIndices,
            indices.QSESpeciesIndices,
            T9,
            rho
        );

        Eigen::VectorXd v_qse_log_initial(indices.QSESpeciesIndices.size());
        for (size_t i = 0; i < indices.QSESpeciesIndices.size(); ++i) {
            v_qse_log_initial(i) = std::log(std::max(Y[indices.QSESpeciesIndices[i]], 1e-99));
        }

        const static std::unordered_map<Eigen::LevenbergMarquardtSpace::Status, const char*> statusMessages = {
            {Eigen::LevenbergMarquardtSpace::NotStarted, "Not started"},
            {Eigen::LevenbergMarquardtSpace::Running, "Running"},
            {Eigen::LevenbergMarquardtSpace::ImproperInputParameters, "Improper input parameters"},
            {Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall, "Relative reduction too small"},
            {Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall, "Relative error too small"},
            {Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall, "Relative error and reduction too small"},
            {Eigen::LevenbergMarquardtSpace::CosinusTooSmall, "Cosine too small"},
            {Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation, "Too many function evaluations"},
            {Eigen::LevenbergMarquardtSpace::FtolTooSmall, "Function tolerance too small"},
            {Eigen::LevenbergMarquardtSpace::XtolTooSmall, "X tolerance too small"},
            {Eigen::LevenbergMarquardtSpace::GtolTooSmall, "Gradient tolerance too small"},
            {Eigen::LevenbergMarquardtSpace::UserAsked, "User asked to stop"}
        };

        Eigen::LevenbergMarquardt lm(functor);
        const Eigen::LevenbergMarquardtSpace::Status info = lm.minimize(v_qse_log_initial);

        if (info <= 0 || info >= 4) {
            LOG_ERROR(m_logger, "QSE species minimization failed with status: {} ({})",
                      static_cast<int>(info), statusMessages.at(info));
            throw std::runtime_error(
                "QSE species minimization failed with status: " + std::to_string(static_cast<int>(info)) +
                " (" + std::string(statusMessages.at(info)) + ")"
            );
        }
        LOG_DEBUG(m_logger, "QSE species minimization completed successfully with status: {} ({})",
                  static_cast<int>(info), statusMessages.at(info));
        return v_qse_log_initial.array().exp();

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

        const auto prevScreeningModel = m_engine.getScreeningModel();
        LOG_DEBUG(m_logger, "Setting screening model to BARE for high temperature and density ignition.");
        m_engine.setScreeningModel(screening::ScreeningType::BARE);
        DirectNetworkSolver ignitionSolver(m_engine);
        NetOut postIgnition = ignitionSolver.evaluate(preIgnition);
        LOG_INFO(m_logger, "Network ignition completed in {} steps.", postIgnition.num_steps);
        m_engine.setScreeningModel(prevScreeningModel);
        LOG_DEBUG(m_logger, "Restoring previous screening model: {}", static_cast<int>(prevScreeningModel));
        return postIgnition;
    }

    bool QSENetworkSolver::shouldUpdateView(const NetIn &conditions) const {
        // Policy 1: If the view has never been initialized, we must update.
        if (!m_isViewInitialized) {
            return true;
        }

        // Policy 2: Check for significant relative change in temperature.
        // Reaction rates are exponentially sensitive to temperature, so we use a tight threshold.
        const double temp_threshold = m_config.get<double>("gridfire:solver:policy:temp_threshold", 0.05); // 5%
        const double temp_relative_change = std::abs(conditions.temperature - m_lastSeenConditions.temperature) / m_lastSeenConditions.temperature;
        if (temp_relative_change > temp_threshold) {
            LOG_DEBUG(m_logger, "Temperature changed by {:.1f}%, triggering view update.", temp_relative_change * 100);
            return true;
        }

        // Policy 3: Check for significant relative change in density.
        const double rho_threshold = m_config.get<double>("gridfire:solver:policy:rho_threshold", 0.10); // 10%
        const double rho_relative_change = std::abs(conditions.density - m_lastSeenConditions.density) / m_lastSeenConditions.density;
        if (rho_relative_change > rho_threshold) {
            LOG_DEBUG(m_logger, "Density changed by {:.1f}%, triggering view update.", rho_relative_change * 100);
            return true;
        }

        // Policy 4: Check for fuel depletion.
        // If a primary fuel source changes significantly, the network structure may change.
        const double fuel_threshold = m_config.get<double>("gridfire:solver:policy:fuel_threshold", 0.15); // 15%

        // Example: Check hydrogen abundance
        const double h1_old = m_lastSeenConditions.composition.getMassFraction("H-1");
        const double h1_new = conditions.composition.getMassFraction("H-1");
        if (h1_old > 1e-12) { // Avoid division by zero
            const double h1_relative_change = std::abs(h1_new - h1_old) / h1_old;
            if (h1_relative_change > fuel_threshold) {
                LOG_DEBUG(m_logger, "H-1 mass fraction changed by {:.1f}%, triggering view update.", h1_relative_change * 100);
                return true;
            }
        }

        // If none of the above conditions are met, the current view is still good enough.
        return false;
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

        Composition outputComposition(speciesNames);
        outputComposition.setMassFraction(speciesNames, finalMassFractions);
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

        // std::string timescales = utils::formatNuclearTimescaleLogString(
        //     m_engine,
        //     y,
        //     m_T9,
        //     m_rho
        // );
        // LOG_TRACE_L2(m_logger, "{}", timescales);

        auto [dydt, eps] = m_engine.calculateRHSAndEnergy(y, m_T9, m_rho);
        dYdt.resize(m_numSpecies + 1);
        std::ranges::copy(dydt, dYdt.begin());
        dYdt(m_numSpecies) = eps;
    }

    void DirectNetworkSolver::JacobianFunctor::operator()(
        const boost::numeric::ublas::vector<double> &Y,
        boost::numeric::ublas::matrix<double> &J,
        double t,
        boost::numeric::ublas::vector<double> &dfdt
    ) const {
        J.resize(m_numSpecies+1, m_numSpecies+1);
        J.clear();
        for (int i = 0; i < m_numSpecies; ++i) {
            for (int j = 0; j < m_numSpecies; ++j) {
                J(i, j) = m_engine.getJacobianMatrixEntry(i, j);
            }
        }
    }
}