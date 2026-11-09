#include "gridfire/extern/gridfire_context.h"
#include "gridfire/extern/gridfire_extern.h"

#include "fourdst/atomic/species.h"

#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "fourdst/composition/utils/utils.h"
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

namespace {
    /// The context's GraphEngine (the bottom of the policy's engine stack), whose rates the multipliers scale
    gridfire::engine::GraphEngine& graph_engine_of(const GFContext& context) {
        if (context.engine == nullptr || context.policy == nullptr) {
            throw GFConfigError("Rate multipliers need the engine: call gf_construct_engine_from_policy first.");
        }
        auto* graph = dynamic_cast<gridfire::engine::GraphEngine*>(context.policy->get_engine_stack().front().get());
        if (graph == nullptr) {
            throw GFConfigError("Rate multipliers are not supported by this policy's engine.");
        }
        return *graph;
    }
}

void GFContext::set_rate_multiplier(const std::string& reaction_id, const double factor) {
    if (!std::isfinite(factor) || factor < 0.0) {
        throw GFConfigError(std::format("Rate multiplier for {} must be finite and >= 0 (got {}).", reaction_id, factor));
    }
    graph_engine_of(*this).setRateMultiplier(reaction_id, factor);
}

double GFContext::get_rate_multiplier(const std::string& reaction_id) const {
    return graph_engine_of(*this).getRateMultiplier(reaction_id);
}

void GFContext::clear_rate_multipliers() {
    graph_engine_of(*this).clearRateMultipliers();
}

void GFContext::set_passive_species(const std::vector<std::string> &names) {
    if (engine != nullptr) {
        throw GFConfigError("Passive species must be set before the engine is constructed (gf_construct_engine_from_policy).");
    }
    for (const auto& name : names) {
        if (!fourdst::atomic::species.contains(name)) {
            throw GFConfigError(std::format("Passive species {} is not recognized in the atomic species database.", name));
        }
    }
    config.engine.passiveSpecies = names;
}

void GFContext::mass_fractions_to_molar(const double *X, const size_t num_species, const size_t num_zones, double *Y) const {
    if (num_species != speciesList.size()) {
        throw fourdst::composition::exceptions::InvalidCompositionError(std::format(
            "Number of species provided ({}) does not match the registered species count ({}).", num_species, speciesList.size()));
    }
    for (size_t zone = 0; zone < num_zones; ++zone) {
        const double* X_zone = X + zone * num_species;
        fourdst::composition::Composition comp;
        try {
            comp = fourdst::composition::buildCompositionFromMassFractions(
                speciesList, std::vector<double>(X_zone, X_zone + num_species));
        } catch (const fourdst::composition::exceptions::InvalidCompositionError& e) {
            throw GFInvalidStateError(std::format("Zone {}: {}", zone, e.what()));
        }
        for (size_t i = 0; i < num_species; ++i) {
            Y[zone * num_species + i] = comp.getMolarAbundance(speciesList[i]);
        }
    }
}

void GFContext::molar_to_mass_fractions(const double *Y, const size_t num_species, const size_t num_zones, double *X) const {
    if (num_species != speciesList.size()) {
        throw fourdst::composition::exceptions::InvalidCompositionError(std::format(
            "Number of species provided ({}) does not match the registered species count ({}).", num_species, speciesList.size()));
    }
    for (size_t zone = 0; zone < num_zones; ++zone) {
        const double* Y_zone = Y + zone * num_species;
        const fourdst::composition::Composition comp(speciesList, std::vector<double>(Y_zone, Y_zone + num_species));
        const auto fractions = comp.getMassFraction();
        for (size_t i = 0; i < num_species; ++i) {
            X[zone * num_species + i] = fractions.at(speciesList[i]);
        }
    }
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
        MAIN_SEQUENCE_POLICY,
        H_HE_BURNING_POLICY
    };

    static const std::unordered_map<std::string, EnginePolicy> engine_map = {
        {"MAIN_SEQUENCE_POLICY", EnginePolicy::MAIN_SEQUENCE_POLICY},
        {"H_HE_BURNING_POLICY", EnginePolicy::H_HE_BURNING_POLICY}
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
        case EnginePolicy::MAIN_SEQUENCE_POLICY:
            this->policy = std::make_unique<gridfire::policy::MainSequencePolicy>(comp, config);
            break;
        case EnginePolicy::H_HE_BURNING_POLICY:
            this->policy = std::make_unique<gridfire::policy::HydrogenHeliumBurningPolicy>(comp, config);
            break;
        default:
            throw gridfire::exceptions::PolicyError(
                "Unhandled engine policy in GFContext::init_engine_from_policy"
            );
    }

    const auto& [e, ctx] = policy->construct();
    this->engine = &e;
    this->policy_name = policy_name;
    this->engine_ctx = ctx->clone_structure();

    this->network_species_names.clear();
    for (const auto& species : e.getNetworkSpecies(*this->engine_ctx)) {
        this->network_species_names.emplace_back(species.name());
    }
    this->reaction_ids.clear();
    for (const auto& reaction : e.getNetworkReactions(*this->engine_ctx)) {
        this->reaction_ids.emplace_back(reaction->id());
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

void GFContext::validate_zone(const double* Y, const size_t num_species, const double T, const double rho, const size_t zone) const {
    if (num_species != speciesList.size()) {
        throw GFInvalidStateError(std::format(
            "Number of species provided ({}) does not match the registered species count ({}).", num_species, speciesList.size()));
    }
    if (!std::isfinite(T) || T <= 0.0 || !std::isfinite(rho) || rho <= 0.0) {
        throw GFInvalidStateError(std::format(
            "Zone {}: temperature ({}) and density ({}) must be positive and finite.", zone, T, rho));
    }
    const auto& network = network_species();
    for (size_t i = 0; i < num_species; ++i) {
        if (!std::isfinite(Y[i]) || Y[i] < 0.0) {
            throw GFInvalidStateError(std::format(
                "Zone {}: molar abundance of {} ({}) must be finite and non-negative.", zone, speciesList[i].name(), Y[i]));
        }
        if (Y[i] > 0.0 && std::ranges::find(network, speciesList[i]) == network.end()) {
            throw GFInvalidStateError(std::format(
                "Zone {}: {} (molar abundance {:.3e}) is not in the network.", zone, speciesList[i].name(), Y[i]));
        }
    }
}

void GFContext::evaluate(
    const double* Y,
    const size_t num_species,
    const double* T,
    const double* rho,
    const size_t num_zones,
    const GFEvaluationOutputs& out
) const {
    const size_t num_channels = engine->getEnergyChannelNames(*engine_ctx).size();
    for (size_t zone = 0; zone < num_zones; ++zone) {
        const double* Y_zone = Y + zone * num_species;
        // Caller errors are reported as GF_INVALID_STATE before anything reaches the engine
        validate_zone(Y_zone, num_species, T[zone], rho[zone], zone);

        // The same pass-through to the engine that the C++ solvers' evaluateRates uses
        gridfire::NetIn netIn;
        netIn.composition = fourdst::composition::Composition(speciesList, std::vector<double>(Y_zone, Y_zone + num_species));
        netIn.temperature = T[zone];
        netIn.density = rho[zone];
        const gridfire::InstantaneousRates rates = gridfire::solver::evaluate_instantaneous_rates(
            *engine, *engine_ctx, netIn, out.eps_by_channel != nullptr, false);

        out.eps_nuc[zone] = rates.eps_nuc;
        out.dEps_dT[zone] = rates.dEps_dT;
        out.dEps_dRho[zone] = rates.dEps_dRho;
        if (out.neutrino_energy_loss_rate) {
            out.neutrino_energy_loss_rate[zone] = rates.neutrino_energy_loss_rate;
        }
        if (out.eps_by_channel) {
            std::ranges::copy(rates.eps_by_channel, out.eps_by_channel + zone * num_channels);
        }
        if (out.dYdt) {
            for (size_t i = 0; i < num_species; ++i) {
                const auto it = rates.dydt.find(speciesList[i]);
                out.dYdt[zone * num_species + i] = it == rates.dydt.end() ? 0.0 : it->second;
            }
        }
    }
}

void GFContext::evaluate_reactions(
    const double* Y,
    const size_t num_species,
    const double* T,
    const double* rho,
    const size_t num_zones,
    const GFReactionOutputs& out
) const {
    const size_t num_reactions = network_reactions().size();
    for (size_t zone = 0; zone < num_zones; ++zone) {
        const double* Y_zone = Y + zone * num_species;
        validate_zone(Y_zone, num_species, T[zone], rho[zone], zone);

        gridfire::NetIn netIn;
        netIn.composition = fourdst::composition::Composition(speciesList, std::vector<double>(Y_zone, Y_zone + num_species));
        netIn.temperature = T[zone];
        netIn.density = rho[zone];
        const gridfire::InstantaneousRates rates = gridfire::solver::evaluate_instantaneous_rates(
            *engine, *engine_ctx, netIn, false, true);
        if (rates.reaction_energy_rate.size() != num_reactions) {
            throw gridfire::exceptions::GridFireError("This engine does not report per-reaction contributions.");
        }
        const size_t offset = zone * num_reactions;
        if (out.molar_flow) std::ranges::copy(rates.reaction_molar_flow, out.molar_flow + offset);
        if (out.energy_rate) std::ranges::copy(rates.reaction_energy_rate, out.energy_rate + offset);
        if (out.neutrino_energy_loss_rate) std::ranges::copy(rates.reaction_neutrino_energy_loss_rate, out.neutrino_energy_loss_rate + offset);
        if (out.neutrino_rate) std::ranges::copy(rates.reaction_neutrino_rate, out.neutrino_rate + offset);
    }
}

void GFPointContext::init_solver_from_engine() {
    this->solver = std::make_unique<gridfire::solver::PointSolver>(*this->engine, config);
    this->solver_ctx = std::make_unique<gridfire::solver::PointSolverContext>(*engine_ctx);
    this->solver_ctx->set_stdout_logging(false); // the per-step table would go to the host code's stdout
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
    zone_status.assign(1, GF_INVALID_STATE); // until the input has been validated
    zone_time_reached.assign(1, 0.0);
    validate_zone(Y_in, num_species, T, rho, 0);
    zone_status[0] = GF_SUCCESS;

    const std::vector<double> Y_scratch(Y_in, Y_in + num_species);
    const fourdst::composition::Composition comp = init_composition_from_abundance_vector(Y_scratch, num_species);

    gridfire::NetIn netIn;
    netIn.temperature = T;
    netIn.density = rho;
    netIn.dt0 = dt0;
    netIn.tMax = tMax;
    netIn.composition = comp;

    gridfire::NetOut result;
    try {
        result = this->solver->evaluate(*solver_ctx, netIn);
    } catch (...) {
        // A failed zone hands back its input abundances and zero rates, and records how far it got
        std::string message;
        zone_status[0] = gf_error_code(std::current_exception(), message);
        if (const auto* point_ctx = dynamic_cast<const gridfire::solver::PointSolverContext*>(solver_ctx.get())) {
            zone_time_reached[0] = point_ctx->last_accepted_time;
        }
        std::copy_n(Y_in, num_species, Y_out);
        energy_out = dEps_dT = dEps_dRho = specific_neutrino_energy_loss = specific_neutrino_flux = mass_lost = 0.0;
        instantaneous.write_zero(0, engine->getEnergyChannelNames(*engine_ctx).size());
        throw;
    }
    zone_time_reached[0] = result.time_reached;

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

void GFInstantaneousOutputs::write_zero(const size_t index, const size_t num_channels) const {
    if (eps_nuc) eps_nuc[index] = 0.0;
    if (dlnEps_dlnT) dlnEps_dlnT[index] = 0.0;
    if (dlnEps_dlnRho) dlnEps_dlnRho[index] = 0.0;
    if (neutrino_energy_loss_rate) neutrino_energy_loss_rate[index] = 0.0;
    if (eps_by_channel) {
        std::fill_n(eps_by_channel + index * num_channels, num_channels, 0.0);
    }
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
    auto grid_ctx = std::make_unique<gridfire::solver::GridSolverContext>(*engine_ctx);
    grid_ctx->zone_completion_logging = false; // nothing on the host code's stdout
    grid_ctx->set_stdout_logging(false);
    this->solver_ctx = std::move(grid_ctx);
}

void GFGridContext::set_shared_composition_groups(const int* group, const double* mass) {
    if (group == nullptr) {
        shared_composition_groups.clear();
        zone_masses.clear();
        return;
    }
    if (zones == 0) {
        throw GFInvalidStateError("Set the number of zones (gf_set_num_zones) before the shared-composition groups.");
    }
    std::vector<int> groups(group, group + zones);
    std::vector<double> masses(zones, 0.0);
    for (size_t z = 0; z < zones; ++z) {
        if (groups[z] < -1) {
            throw GFInvalidStateError(std::format("Zone {} has group id {}; use -1 for a zone burnt on its own.", z, groups[z]));
        }
        if (groups[z] >= 0) {
            if (mass == nullptr || !std::isfinite(mass[z]) || mass[z] <= 0.0) {
                throw GFInvalidStateError(std::format("Zone {} is in group {} but has no positive, finite mass.", z, groups[z]));
            }
        }
        if (mass != nullptr) {
            masses[z] = mass[z];
        }
    }
    shared_composition_groups = std::move(groups);
    zone_masses = std::move(masses);
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

    zone_status.assign(this->get_zones(), GF_INVALID_STATE); // until the input has been validated
    zone_time_reached.assign(this->get_zones(), 0.0);
    for (size_t i = 0; i < this->get_zones(); ++i) {
        validate_zone(Y_in + i * num_species, num_species, T[i], rho[i], i);
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

    auto& grid_ctx = dynamic_cast<gridfire::solver::GridSolverContext&>(*this->solver_ctx);
    if (shared_composition_groups.empty()) {
        grid_ctx.clear_shared_composition_groups();
    } else {
        grid_ctx.set_shared_composition_groups(shared_composition_groups, zone_masses);
    }
    std::vector<gridfire::NetOut> results = this->solver->evaluate(*this->solver_ctx, netIns);
    const auto& failures = dynamic_cast<const gridfire::solver::GridSolverContext&>(*this->solver_ctx).zone_failures;
    const size_t num_channels = engine->getEnergyChannelNames(*engine_ctx).size();

    zone_status.assign(this->get_zones(), GF_SUCCESS);
    zone_time_reached.assign(this->get_zones(), 0.0);
    size_t num_failed = 0;
    std::string failure_summary;

    for (size_t zone_idx = 0; zone_idx < this->get_zones(); ++zone_idx) {
        const gridfire::NetOut& netOut = results[zone_idx];
        zone_time_reached[zone_idx] = netOut.time_reached;
        if (!netOut.converged) {
            // A failed zone hands back its input abundances and zero rates; the other zones are unaffected
            std::string message;
            zone_status[zone_idx] = zone_idx < failures.size() && failures[zone_idx]
                ? gf_error_code(failures[zone_idx], message)
                : GF_SOLVER_ERROR;
            if (num_failed < 3) {
                failure_summary += std::format("{}zone {}: {}", num_failed ? "; " : "", zone_idx,
                                               message.empty() ? netOut.failure_message : message);
            }
            ++num_failed;
            std::copy_n(Y_in + zone_idx * num_species, num_species, Y_out + zone_idx * num_species);
            energy_out[zone_idx] = dEps_dT[zone_idx] = dEps_dRho[zone_idx] = 0.0;
            specific_neutrino_energy_loss[zone_idx] = specific_neutrino_flux[zone_idx] = mass_lost[zone_idx] = 0.0;
            instantaneous.write_zero(zone_idx, num_channels);
            continue;
        }
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
    if (num_failed > 0) {
        last_error = std::format("{} of {} zones failed ({}{}). See gf_get_zone_status.", num_failed, this->get_zones(),
                                 failure_summary, num_failed > 3 ? "; ..." : "");
        return GF_ZONE_FAILURE;
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
