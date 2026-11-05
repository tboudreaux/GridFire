#include "gridfire/extern/gridfire_context.h"

#include "fourdst/atomic/species.h"

#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/exceptions/error_policy.h"
#include "gridfire/utils/logging.h"

#include <algorithm>
#include <cmath>

#include "fourdst/config/config.h"

void GFContext::init_species_map(const std::vector<std::string> &species_names) {
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

void GFContext::set_config(const gridfire::config::GridFireConfig &new_config) {
    if (engine != nullptr) {
        throw GFConfigError("The configuration must be set before the engine is constructed (gf_construct_engine_from_policy); the engine already built would not use it.");
    }
    config = new_config;
}

void GFContext::load_config(const std::string &path) {
    fourdst::config::Config<gridfire::config::GridFireConfig> loaded;
    loaded.load(path);
    set_config(*loaded);
}

void GFContext::save_config(const std::string &path) const {
    fourdst::config::Config<gridfire::config::GridFireConfig> current;
    current.mutate([this](auto& cfg) {
        cfg = config;
    });
    current.save(path);
}

void GFContext::init_engine_from_policy(const std::string &policy_name, const double *abundances, const size_t num_species) {
    const std::vector<double> Y_scratch(abundances, abundances + num_species);
    fourdst::composition::Composition comp = init_composition_from_abundance_vector(Y_scratch, num_species);

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
            this->policy = std::make_unique<gridfire::policy::MainSequencePolicy>(comp, config);
            const auto& [e, ctx] = policy->construct();
            this->engine = &e;
            this->engine_ctx = ctx->clone_structure();

            this->network_species_names.clear();
            for (const auto& species : e.getNetworkSpecies(*this->engine_ctx)) {
                this->network_species_names.emplace_back(species.name());
            }
            this->reaction_ids.clear();
            for (const auto& reaction : e.getNetworkReactions(*this->engine_ctx)) {
                this->reaction_ids.emplace_back(reaction->id());
            }

            break;
        }
        default:
            throw gridfire::exceptions::PolicyError(
                "Unhandled engine policy in GFPointContext::init_engine_from_policy"
            );
    }


}

fourdst::composition::Composition GFContext::init_composition_from_abundance_vector(const std::vector<double> &abundances, size_t num_species) const {
    if (num_species == 0) {
        throw fourdst::composition::exceptions::InvalidCompositionError("Cannot initialize composition with zero species.");
    }
    if (num_species != speciesList.size()) {
        throw fourdst::composition::exceptions::InvalidCompositionError(
            std::format(
                "Number of species provided ({}) does not match the registered species count ({}).",
                num_species,
                speciesList.size()
            )
        );
    }

    fourdst::composition::Composition comp;
    for (size_t i = 0; i < num_species; i++) {
        comp.registerSpecies(this->speciesList[i]);
        comp.setMolarAbundance(this->speciesList[i], abundances[i]);
    }

    return comp;
}


const std::vector<fourdst::atomic::Species>& GFContext::network_species() const {
    if (!engine || !engine_ctx) {
        throw gridfire::exceptions::GridFireError("The network is unavailable until an engine has been constructed (call gf_construct_engine_from_policy first).");
    }
    return engine->getNetworkSpecies(*engine_ctx);
}

const gridfire::reaction::ReactionSet& GFContext::network_reactions() const {
    if (!engine || !engine_ctx) {
        throw gridfire::exceptions::GridFireError("The network is unavailable until an engine has been constructed (call gf_construct_engine_from_policy first).");
    }
    return engine->getNetworkReactions(*engine_ctx);
}

void GFContext::register_network_species() {
    speciesList = network_species();
}

void GFContext::evaluate(
    const double* Y,
    const size_t num_species,
    const double* T,
    const double* rho,
    const size_t num_zones,
    const GFEvaluationOutputs& out
) const {
    const auto& network = network_species();
    if (num_species != speciesList.size()) {
        throw GFInvalidStateError(std::format(
            "Number of species provided ({}) does not match the registered species count ({}).", num_species, speciesList.size()));
    }

    // Where each registered species sits in the network (npos: not in the network)
    constexpr size_t npos = static_cast<size_t>(-1);
    std::vector<size_t> network_index(num_species, npos);
    for (size_t i = 0; i < num_species; ++i) {
        if (const auto it = std::ranges::find(network, speciesList[i]); it != network.end()) {
            network_index[i] = static_cast<size_t>(it - network.begin());
        }
    }

    const size_t num_channels = engine->getEnergyChannelNames(*engine_ctx).size();
    std::vector<double> Y_network(network.size());
    for (size_t zone = 0; zone < num_zones; ++zone) {
        const double* Y_zone = Y + zone * num_species;
        if (!std::isfinite(T[zone]) || T[zone] <= 0.0 || !std::isfinite(rho[zone]) || rho[zone] <= 0.0) {
            throw GFInvalidStateError(std::format(
                "Zone {}: temperature ({}) and density ({}) must be positive and finite.", zone, T[zone], rho[zone]));
        }

        // The engine needs every network species; those not registered are absent
        std::ranges::fill(Y_network, 0.0);
        for (size_t i = 0; i < num_species; ++i) {
            if (!std::isfinite(Y_zone[i]) || Y_zone[i] < 0.0) {
                throw GFInvalidStateError(std::format(
                    "Zone {}: molar abundance of {} ({}) must be finite and non-negative.", zone, speciesList[i].name(), Y_zone[i]));
            }
            if (network_index[i] == npos) {
                if (Y_zone[i] > 0.0) {
                    throw GFInvalidStateError(std::format(
                        "Zone {}: {} (molar abundance {:.3e}) is not in the network.", zone, speciesList[i].name(), Y_zone[i]));
                }
                continue;
            }
            Y_network[network_index[i]] = Y_zone[i];
        }
        const fourdst::composition::Composition comp(network, Y_network);

        const double T9 = T[zone] / 1e9;
        const auto rhs = engine->calculateRHSAndEnergy(*engine_ctx, comp, T9, rho[zone], false);
        if (!rhs) {
            throw gridfire::exceptions::BadRHSEngineError(std::format(
                "Zone {}: the engine could not evaluate the RHS ({}).", zone, gridfire::engine::EngineStatus_to_string(rhs.error())));
        }
        const auto [dEps_dT, dEps_dRho] = engine->calculateEpsDerivatives(*engine_ctx, comp, T9, rho[zone]);

        out.eps_nuc[zone] = rhs->nuclearEnergyGenerationRate;
        out.dEps_dT[zone] = dEps_dT;
        out.dEps_dRho[zone] = dEps_dRho;
        if (out.neutrino_energy_loss_rate) {
            out.neutrino_energy_loss_rate[zone] = rhs->neutrinoEnergyLossRate;
        }
        if (out.eps_by_channel) {
            const std::vector<double> channels = engine->calculateEnergyChannels(*engine_ctx, comp, T9, rho[zone]);
            std::ranges::copy(channels, out.eps_by_channel + zone * num_channels);
        }
        if (out.dYdt) {
            for (size_t i = 0; i < num_species; ++i) {
                const auto it = rhs->dydt.find(speciesList[i]);
                out.dYdt[zone * num_species + i] = it == rhs->dydt.end() ? 0.0 : it->second;
            }
        }
    }
}

void GFPointContext::init_solver_from_engine() {
    this->solver = std::make_unique<gridfire::solver::PointSolver>(*this->engine, config);
    this->solver_ctx = std::make_unique<gridfire::solver::PointSolverContext>(*engine_ctx);
}


int GFPointContext::evolve(
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
    double& mass_lost,
    const GFInstantaneousOutputs& instantaneous
) const {

    const std::vector<double> Y_scratch(Y_in, Y_in + num_species);
    const fourdst::composition::Composition comp = init_composition_from_abundance_vector(Y_scratch, num_species);

    gridfire::NetIn netIn;
    netIn.temperature = T;
    netIn.density = rho;
    netIn.dt0 = dt0;
    netIn.tMax = tMax;
    netIn.composition = comp;

    const gridfire::NetOut result = this->solver->evaluate(*solver_ctx, netIn);

    energy_out = result.energy;
    dEps_dT  = result.dEps_dT;
    dEps_dRho = result.dEps_dRho;
    specific_neutrino_energy_loss = result.specific_neutrino_energy_loss;
    specific_neutrino_flux = result.specific_neutrino_flux;
    instantaneous.write(result, 0);

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

void GFInstantaneousOutputs::write(const gridfire::NetOut& netOut, const size_t index) const {
    if (eps_nuc) eps_nuc[index] = netOut.eps_nuc;
    if (dlnEps_dlnT) dlnEps_dlnT[index] = netOut.dlnEps_dlnT;
    if (dlnEps_dlnRho) dlnEps_dlnRho[index] = netOut.dlnEps_dlnRho;
    if (neutrino_energy_loss_rate) neutrino_energy_loss_rate[index] = netOut.neutrino_energy_loss_rate;
    if (eps_by_channel) {
        const size_t num_channels = netOut.eps_by_channel.size();
        for (size_t channel = 0; channel < num_channels; ++channel) {
            eps_by_channel[index * num_channels + channel] = netOut.eps_by_channel[channel];
        }
    }
}

const std::vector<std::string>& GFContext::get_energy_channel_names() const {
    if (!engine || !engine_ctx) {
        throw gridfire::exceptions::GridFireError("Energy channels are unavailable until an engine has been constructed (call gf_construct_engine_from_policy first).");
    }
    return engine->getEnergyChannelNames(*engine_ctx);
}

void GFGridContext::init_solver_from_engine() {
    this->local_solver = std::make_unique<gridfire::solver::PointSolver>(*this->engine, config);
    this->solver = std::make_unique<gridfire::solver::GridSolver>(*this->engine, *this->local_solver);
    this->solver_ctx = std::make_unique<gridfire::solver::GridSolverContext>(*engine_ctx);
}

int GFGridContext::evolve(
    const double* Y_in,
    size_t num_species,
    const double *T,
    const double *rho,
    double tMax,
    double dt0,
    double *Y_out,
    double *energy_out,
    double *dEps_dT,
    double *dEps_dRho,
    double *specific_neutrino_energy_loss,
    double *specific_neutrino_flux,
    double *mass_lost,
    const GFInstantaneousOutputs& instantaneous
) const {
    if (this->get_zones() == 0) {
        throw gridfire::exceptions::GridFireError("GFGridContext has zero zones configured for evolution.");
    }

    if (Y_in == nullptr || T == nullptr || rho == nullptr) {
        throw gridfire::exceptions::GridFireError("Input abundance, temperature, or density pointers are null.");
    }

    if (Y_out == nullptr) {
        throw gridfire::exceptions::GridFireError("Output abundance pointer is null.");
    }

    std::vector<fourdst::composition::Composition> zone_compositions;

    zone_compositions.reserve(this->get_zones());

    std::vector<double> Y_scratch;
    Y_scratch.resize(num_species);

    for (size_t i = 0; i < this->get_zones(); ++i) {
        for (size_t j = 0; j < num_species; ++j) {
            Y_scratch[j] = Y_in[i * num_species + j];
        }
        zone_compositions.push_back(init_composition_from_abundance_vector(Y_scratch, num_species));
    }

    std::vector<gridfire::NetIn> netIns;
    netIns.reserve(this->get_zones());
    for (size_t i = 0; i < this->get_zones(); ++i) {
        gridfire::NetIn netIn;
        netIn.temperature = T[i];
        netIn.density = rho[i];
        netIn.dt0 = dt0;
        netIn.tMax = tMax;
        netIn.composition = zone_compositions[i];
        netIns.push_back(netIn);
    }

    std::vector<gridfire::NetOut> results = this->solver->evaluate(*this->solver_ctx, netIns);

    for (size_t zone_idx = 0; zone_idx < this->get_zones(); ++zone_idx) {
        const gridfire::NetOut& netOut = results[zone_idx];
        energy_out[zone_idx] = netOut.energy;
        dEps_dT[zone_idx] = netOut.dEps_dT;;
        dEps_dRho[zone_idx] = netOut.dEps_dRho;
        specific_neutrino_energy_loss[zone_idx] = netOut.specific_neutrino_energy_loss;
        specific_neutrino_flux[zone_idx] = netOut.specific_neutrino_flux;
        instantaneous.write(netOut, zone_idx);

        std::set<fourdst::atomic::Species> seen_species;
        for (size_t i = 0; i < num_species; ++i) {
            fourdst::atomic::Species species = this->speciesList[i];
            Y_out[zone_idx * num_species + i] = netOut.composition.getMolarAbundance(species);
            seen_species.insert(species);
        }
        mass_lost[zone_idx] = 0.0;
        for (const auto& species : netOut.composition.getRegisteredSpecies()) {
            if (!seen_species.contains(species)) {
                mass_lost[zone_idx] += species.mass() * netOut.composition.getMolarAbundance(species);
            }
        }
    }
    return 0;
}

std::unique_ptr<GFContext> make_gf_context(const GFContextType &type) {
    switch (type) {
        case GFContextType::POINT:
            return std::make_unique<GFPointContext>();
        case GFContextType::GRID:
            return std::make_unique<GFGridContext>();
        default:
            throw gridfire::exceptions::GridFireError("Unhandled GFContextType in make_gf_context");
    }
}
