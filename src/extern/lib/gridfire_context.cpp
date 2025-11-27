#include "gridfire/extern/gridfire_context.h"

#include "fourdst/atomic/species.h"

#include "fourdst/composition/exceptions/exceptions_composition.h"

void GridFireContext::init_species_map(const std::vector<std::string> &species_names) {
    for (const auto& name: species_names) {
        working_comp.registerSymbol(name);
    }
    this->speciesList.clear();
    this->speciesList.reserve(species_names.size());

    auto resolve_species_name = [](const std::string& name) -> fourdst::atomic::Species {
        if (fourdst::atomic::species.contains(name)) {
            return fourdst::atomic::species.at(name);
        }
        throw fourdst::composition::exceptions::UnknownSymbolError("Species " + name + " is not recognized in the atomic species database.");
    };

    for (const auto& name: species_names) {
        this->speciesList.push_back(resolve_species_name(name));
    }

}

void GridFireContext::init_engine_from_policy(const std::string &policy_name, const double *abundances, const size_t num_species) {
    init_composition_from_abundance_vector(abundances, num_species);

    enum class EnginePolicy {
        MAIN_SEQUENCE_POLICY
    };

    static const std::unordered_map<std::string, EnginePolicy> engine_map = {
        {"MAIN_SEQUENCE_POLICY", EnginePolicy::MAIN_SEQUENCE_POLICY}
    };

    if (!engine_map.contains(policy_name)) {
        throw gridfire::exceptions::PolicyError(
            std::format(
                "Engine Policy {} is not recognized. Valid policies are: {}",
                policy_name,
                gridfire::utils::iterable_to_delimited_string(engine_map, ", ", [](const auto& pair){ return pair.first; })
            )
        );
    }

    switch (engine_map.at(policy_name)) {
        case EnginePolicy::MAIN_SEQUENCE_POLICY: {
            this->policy = std::make_unique<gridfire::policy::MainSequencePolicy>(this->working_comp);
            this->engine = &policy->construct();
            break;
        }
        default:
            throw gridfire::exceptions::PolicyError(
                "Unhandled engine policy in GridFireContext::init_engine_from_policy"
            );
    }


}

void GridFireContext::init_solver_from_engine(const std::string &solver_name) {
    enum class SolverType {
        CVODE
    };

    static const std::unordered_map<std::string, SolverType> solver_map = {
        {"CVODE", SolverType::CVODE}
    };

    if (!solver_map.contains(solver_name)) {
        throw gridfire::exceptions::SolverError(
            std::format(
                "Solver {} is not recognized. Valid solvers are: {}",
                solver_name,
                gridfire::utils::iterable_to_delimited_string(solver_map, ", ", [](const auto& pair){ return pair.first; })
            )
        );
    }

    switch (solver_map.at(solver_name)) {
        case SolverType::CVODE: {
            this->solver = std::make_unique<gridfire::solver::CVODESolverStrategy>(*this->engine);
            break;
        }
        default:
            throw gridfire::exceptions::SolverError(
                "Unhandled solver type in GridFireContext::init_solver_from_engine"
            );
    }

}

void GridFireContext::init_composition_from_abundance_vector(const double *abundances, size_t num_species) {
    if (num_species == 0) {
        throw fourdst::composition::exceptions::InvalidCompositionError("Cannot initialize composition with zero species.");
    }
    if (num_species != working_comp.size()) {
        throw fourdst::composition::exceptions::InvalidCompositionError(
            std::format(
                "Number of species provided ({}) does not match the registered species count ({}).",
                num_species,
                working_comp.size()
            )
        );
    }

    for (size_t i = 0; i < num_species; i++) {
        this->working_comp.setMolarAbundance(this->speciesList[i], abundances[i]);
    }
}

int GridFireContext::evolve(
    const double* Y_in,
    const size_t num_species,
    const double T,
    const double rho,
    const double tMax,
    const double dt0,
    double* Y_out,
    double& energy_out,
    double& dEps_dT,
    double& dEps_dRho,
    double& specific_neutrino_energy_loss,
    double& specific_neutrino_flux,
    double& mass_lost
) {
    init_composition_from_abundance_vector(Y_in, num_species);

    gridfire::NetIn netIn;
    netIn.temperature = T;
    netIn.density = rho;
    netIn.dt0 = dt0;
    netIn.tMax = tMax;
    netIn.composition = this->working_comp;

    const gridfire::NetOut result = this->solver->evaluate(netIn);

    energy_out = result.energy;
    dEps_dT  = result.dEps_dT;
    dEps_dRho = result.dEps_dRho;
    specific_neutrino_energy_loss = result.specific_neutrino_energy_loss;
    specific_neutrino_flux = result.specific_neutrino_flux;

    std::set<fourdst::atomic::Species> seen_species;
    for (size_t i = 0; i < num_species; i++) {
        fourdst::atomic::Species species = this->speciesList[i];
        Y_out[i] = result.composition.getMolarAbundance(species);
        seen_species.insert(species);
    }

    mass_lost = 0.0;
    for (const auto& species : result.composition.getRegisteredSpecies()) {
        if (!seen_species.contains(species)) {
            mass_lost += species.mass() * result.composition.getMolarAbundance(species);
        }
    }

    return 0;
}