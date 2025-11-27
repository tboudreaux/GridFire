#ifndef GF_GRIDFIRE_CONTEXT_H
#define GF_GRIDFIRE_CONTEXT_H

#include "gridfire/gridfire.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <memory>
#include <vector>

struct GridFireContext {
    std::unique_ptr<gridfire::policy::NetworkPolicy> policy;
    gridfire::engine::DynamicEngine* engine;
    std::unique_ptr<gridfire::solver::DynamicNetworkSolverStrategy> solver;

    std::vector<fourdst::atomic::Species> speciesList;
    fourdst::composition::Composition working_comp;

    void init_species_map(const std::vector<std::string>& species_names);
    void init_engine_from_policy(const std::string& policy_name, const double *abundances, size_t num_species);
    void init_solver_from_engine(const std::string& solver_name);

    void init_composition_from_abundance_vector(const double* abundances, size_t num_species);

    int evolve(
        const double* Y_in,
        size_t num_species,
        double T,
        double rho,
        double tMax,
        double dt0,
        double* Y_out,
        double& energy_out,
        double& dEps_dT,
        double& dEps_dRho, double& mass_lost
    );

    std::string last_error;
};

#endif