#include "gridfire/solver/solver.h"
#include "gridfire/engine/engine_graph.h"
#include "gridfire/network.h"
#include "gridfire/exceptions/error_engine.h"

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/composition.h"
#include "fourdst/config/config.h"

#include "unsupported/Eigen/NonLinearOptimization"

#include <boost/numeric/odeint.hpp>

#include <vector>
#include <string>
#include <stdexcept>
#include <iomanip>

#include "quill/LogMacros.h"

namespace gridfire::solver {
    NetOut DirectNetworkSolver::evaluate(const NetIn &netIn) {
        namespace ublas = boost::numeric::ublas;
        namespace odeint = boost::numeric::odeint;
        using fourdst::composition::Composition;


        const double T9 = netIn.temperature / 1e9; // Convert temperature from Kelvin to T9 (T9 = T / 1e9)

        const auto absTol = m_config.get<double>("gridfire:solver:DirectNetworkSolver:absTol", 1.0e-8);
        const auto relTol = m_config.get<double>("gridfire:solver:DirectNetworkSolver:relTol", 1.0e-8);

        Composition equilibratedComposition = m_engine.update(netIn);
        size_t numSpecies = m_engine.getNetworkSpecies().size();
        ublas::vector<double> Y(numSpecies + 1);

        RHSManager manager(m_engine, T9, netIn.density);
        JacobianFunctor jacobianFunctor(m_engine, T9, netIn.density);

        auto populateY = [&](const Composition& comp) {
            const size_t numSpeciesInternal = m_engine.getNetworkSpecies().size();
            Y.resize(numSpeciesInternal + 1);
            for (size_t i = 0; i < numSpeciesInternal; i++) {
                const auto& species = m_engine.getNetworkSpecies()[i];
                if (!comp.contains(species)) {
                    double lim = std::numeric_limits<double>::min();
                    LOG_DEBUG(m_logger, "Species '{}' not found in composition. Setting abundance to {:0.3E}.", species.name(), lim);
                    Y(i) = lim; // Species not in the composition, set to zero
                } else {
                    Y(i) = comp.getMolarAbundance(species);
                }
            }

            // TODO: a good starting point to make the temperature, density, and energy self consistent would be to turn this into an accumulator
            Y(numSpeciesInternal) = 0.0; // Specific energy rate, initialized to zero
        };

        // This is a quick debug that can be turned on. For solar code input parameters (T~1.5e7K, ρ~1.5e3 g/cm^3) this should be near 8e-17
        // std::cout << "D/H: " << equilibratedComposition.getMolarAbundance("H-2") / equilibratedComposition.getMolarAbundance("H-1") << std::endl;

        populateY(equilibratedComposition);
        const auto stepper = odeint::make_controlled<odeint::rosenbrock4<double>>(absTol, relTol);

        double current_time = 0.0;
        double current_initial_timestep = netIn.dt0;
        double accumulated_energy = 0.0;
        // size_t total_update_stages_triggered = 0;

        while (current_time < netIn.tMax) {
            try {
                odeint::integrate_adaptive(
                    stepper,
                    std::make_pair(manager, jacobianFunctor),
                    Y,
                    current_time,
                    netIn.tMax,
                    current_initial_timestep,
                    [&](const auto& state, double t) {
                        current_time = t;
                        manager.observe(state, t);
                    }
                );
                current_time = netIn.tMax;
            } catch (const exceptions::StaleEngineTrigger &e) {
                LOG_INFO(m_logger, "Catching StaleEngineTrigger at t = {:0.3E} with T9 = {:0.3E}, rho = {:0.3E}. Triggering update stage (last stage took {} steps).", current_time, T9, netIn.density, e.totalSteps());
                exceptions::StaleEngineTrigger::state staleState = e.getState();
                accumulated_energy += e.energy(); // Add the specific energy rate to the accumulated energy
                // total_update_stages_triggered++;

                Composition temp_comp;
                std::vector<double> mass_fractions;
                size_t num_species_at_stop = e.numSpecies();

                if (num_species_at_stop != m_engine.getNetworkSpecies().size()) {
                    throw std::runtime_error(
                        "StaleEngineError state has a different number of species than the engine. This should not happen."
                    );
                }
                mass_fractions.reserve(num_species_at_stop);

                for (size_t i = 0; i < num_species_at_stop; ++i) {
                    const auto& species = m_engine.getNetworkSpecies()[i];
                    temp_comp.registerSpecies(species);
                    mass_fractions.push_back(e.getMolarAbundance(i) * species.mass()); // Convert from molar abundance to mass fraction
                }
                temp_comp.setMassFraction(m_engine.getNetworkSpecies(), mass_fractions);
                temp_comp.finalize(true);

                NetIn netInTemp = netIn;
                netInTemp.temperature = e.temperature();
                netInTemp.density = e.density();
                netInTemp.composition = std::move(temp_comp);

                Composition currentComposition = m_engine.update(netInTemp);
                populateY(currentComposition);
                Y(Y.size() - 1) = e.energy(); // Set the specific energy rate from the stale state
                numSpecies = m_engine.getNetworkSpecies().size();

                // current_initial_timestep = 0.001 * manager.m_last_step_time; // set the new timestep to the last successful timestep before repartitioning
            }
        }

        accumulated_energy += Y(Y.size() - 1); // Add the specific energy rate to the accumulated energy

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
        netOut.energy = accumulated_energy; // Specific energy rate
        netOut.num_steps = manager.m_num_steps;

        return netOut;
    }

    void DirectNetworkSolver::RHSManager::operator()(
        const boost::numeric::ublas::vector<double> &Y,
        boost::numeric::ublas::vector<double> &dYdt,
        const double t
    ) const {
        const size_t numSpecies = m_engine.getNetworkSpecies().size();
        if (t != m_cached_time || !m_cached_result.has_value() || m_cached_result.value().dydt.size() != numSpecies + 1) {
            compute_and_cache(Y, t);
        }
        const auto&[dydt, nuclearEnergyGenerationRate] = m_cached_result.value();
        dYdt.resize(numSpecies + 1);
        std::ranges::copy(dydt, dYdt.begin());
        dYdt(numSpecies) = nuclearEnergyGenerationRate; // Set the last element to the specific energy rate
    }

    void DirectNetworkSolver::RHSManager::observe(
        const boost::numeric::ublas::vector<double> &state,
        const double t
    ) const {
        double dt = t - m_last_observed_time;
        compute_and_cache(state, t);
        LOG_INFO(
            m_logger,
            "(Step {}) Observed state at t = {:0.3E} (dt = {:0.3E})",
            m_num_steps,
            t,
            dt
        );
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(3);
        oss << "(Step: " << std::setw(10) << m_num_steps << ") t = " << t << " (dt = " << dt << ", eps_nuc: " << state(state.size() - 1) << " [erg])\n";
        std::cout << oss.str();
        m_last_observed_time = t;
        m_last_step_time = dt;

    }

    void DirectNetworkSolver::RHSManager::compute_and_cache(
        const boost::numeric::ublas::vector<double> &state,
        double t
    ) const {
        std::vector<double> y_vec(state.begin(), state.end() - 1);
        std::ranges::replace_if(
            y_vec,
            [](const double yi){
                return yi < 0.0;
            },
            0.0 // Avoid negative abundances
        );

        const auto result = m_engine.calculateRHSAndEnergy(y_vec, m_T9, m_rho);
        if (!result) {
            LOG_INFO(m_logger,
                "Triggering update stage due to stale engine detected at t = {:0.3E} with T9 = {:0.3E}, rho = {:0.3E}, y_vec (size: {})",
                t, m_T9, m_rho, y_vec.size());
            throw exceptions::StaleEngineTrigger({m_T9, m_rho, y_vec, t, m_num_steps, state(state.size() - 1)});
        }
        m_cached_result = result.value();
        m_cached_time = t;

        m_num_steps++;
    }

    void DirectNetworkSolver::JacobianFunctor::operator()(
        const boost::numeric::ublas::vector<double> &Y,
        boost::numeric::ublas::matrix<double> &J,
        double t,
        boost::numeric::ublas::vector<double> &dfdt
    ) const {
        size_t numSpecies = m_engine.getNetworkSpecies().size();
        J.resize(numSpecies+1, numSpecies+1);
        J.clear();
        for (size_t i = 0; i < numSpecies; ++i) {
            for (size_t j = 0; j < numSpecies; ++j) {
                J(i, j) = m_engine.getJacobianMatrixEntry(i, j);
            }
        }
    }

}