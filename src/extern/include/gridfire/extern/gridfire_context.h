#ifndef GF_GRIDFIRE_CONTEXT_H
#define GF_GRIDFIRE_CONTEXT_H

#include "gridfire/gridfire.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Raised by the C API for an invalid configuration, or a configuration change after the engine was built.
 * gf_* functions report it as GF_CONFIG_ERROR.
 */
struct GFConfigError final : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * @brief Raised by the C API for an invalid thermodynamic state or composition passed in by the caller (as opposed to a
 * solver failure). gf_* functions report it as GF_INVALID_STATE.
 */
struct GFInvalidStateError final : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// The GF_* / FDSSE_* code for an exception (as every gf_* function reports it), with its message.
int gf_error_code(const std::exception_ptr& error, std::string& message);

/// Destinations for gf_evaluate. Each points to num_zones entries (eps_by_channel and dYdt: per zone, contiguous).
struct GFEvaluationOutputs {
    double* eps_nuc = nullptr;
    double* dEps_dT = nullptr;
    double* dEps_dRho = nullptr;
    double* neutrino_energy_loss_rate = nullptr; ///< May be null
    double* eps_by_channel = nullptr; ///< May be null; num_channels per zone
    double* dYdt = nullptr; ///< May be null; num_species per zone, in registered species order
};

/**
 * @brief Destinations for the instantaneous final-state quantities in NetOut.
 *
 * For a point context each scalar pointer refers to a single double; for a grid context each refers to an array with
 * one entry per zone. eps_by_channel holds the channel energies of each zone contiguously
 * (eps_by_channel[zone * num_channels + channel]).
 */
struct GFInstantaneousOutputs {
    double* eps_nuc = nullptr;
    double* dlnEps_dlnT = nullptr;
    double* dlnEps_dlnRho = nullptr;
    double* neutrino_energy_loss_rate = nullptr;
    double* eps_by_channel = nullptr;

    void write(const gridfire::NetOut& netOut, size_t index) const;
    /// Zero the outputs of zone `index` (for a zone that failed)
    void write_zero(size_t index, size_t num_channels) const;
};

enum class GFContextType {
    POINT,
    GRID
};

struct GFContext {
    virtual ~GFContext() = default;

    std::unique_ptr<gridfire::policy::NetworkPolicy> policy;
    const gridfire::engine::DynamicEngine* engine = nullptr;
    gridfire::config::GridFireConfig config; ///< Used for the engine and solver; set before init_engine_from_policy

    /// Replace the configuration. Throws GFConfigError once the engine has been built, since it would not apply.
    void set_config(const gridfire::config::GridFireConfig& new_config);
    /// Replace the configuration with one read from a TOML file (the format written by save_config).
    void load_config(const std::string& path);
    /// Write the current configuration as TOML, a complete file to edit and pass to load_config.
    void save_config(const std::string& path) const;
    std::unique_ptr<gridfire::engine::scratch::StateBlob> engine_ctx;
    std::vector<fourdst::atomic::Species> speciesList;
    std::vector<std::string> network_species_names; ///< Null-terminated names handed out by gf_get_network_species
    std::vector<std::string> reaction_ids; ///< Null-terminated ids handed out by gf_get_reaction_id

    virtual void init_species_map(const std::vector<std::string>& species_names);
    virtual void init_engine_from_policy(const std::string& policy_name, const double *abundances, size_t num_species);
    virtual void init_solver_from_engine() = 0;

    fourdst::composition::Composition init_composition_from_abundance_vector(const std::vector<double> &abundances, size_t num_species) const;
    /// The engine's network species, in network order. Throws if no engine has been built yet.
    [[nodiscard]] const std::vector<fourdst::atomic::Species>& network_species() const;
    /// Reactions in the engine's network. Throws if no engine has been built yet.
    [[nodiscard]] const gridfire::reaction::ReactionSet& network_reactions() const;
    /// Register the network's own species, in network order, as the species of every later call.
    void register_network_species();

    /**
     * Instantaneous eps_nuc, its derivatives, neutrino losses, channel energies and dY/dt for num_zones states
     * (Y holds num_species abundances per zone, in registered species order). No solver is involved.
     */
    void evaluate(const double* Y, size_t num_species, const double* T, const double* rho, size_t num_zones,
                  const GFEvaluationOutputs& out) const;
    [[nodiscard]] const std::vector<std::string>& get_energy_channel_names() const;
    mutable std::string last_error;

    /// Outcome of each zone of the latest gf_evolve call (GF_SUCCESS or the zone's error code), and the time it reached
    mutable std::vector<int> zone_status;
    mutable std::vector<double> zone_time_reached;

    /// Throws GFInvalidStateError unless T, rho and the abundances of one zone are physical and Y has no positive
    /// abundance outside the network (which the integrator would have to drop).
    void validate_zone(const double* Y, size_t num_species, double T, double rho, size_t zone) const;
};

struct GFPointContext final: GFContext{
    std::unique_ptr<gridfire::solver::SingleZoneDynamicNetworkSolver> solver;
    std::unique_ptr<gridfire::solver::SolverContextBase> solver_ctx;

    void init_solver_from_engine() override;

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
        double& dEps_dRho,
        double& specific_neutrino_energy_loss,
        double& specific_neutrino_flux,
        double& mass_lost,
        const GFInstantaneousOutputs& instantaneous
    ) const;

};

struct GFGridContext final : GFContext {
    std::unique_ptr<gridfire::solver::SingleZoneDynamicNetworkSolver> local_solver;
    std::unique_ptr<gridfire::solver::MultiZoneDynamicNetworkSolver> solver;
    std::unique_ptr<gridfire::solver::SolverContextBase> solver_ctx;

    void init_solver_from_engine() override;

    size_t zones;

    void set_zones(const size_t num_zones) {
        zones = num_zones;
    }

    [[nodiscard]] size_t get_zones() const {
        return zones;
    }

    int evolve(
        const double* Y_in,
        size_t num_species,
        const double* T,
        const double* rho,
        double tMax,
        double dt0,
        double* Y_out,
        double* energy_out,
        double* dEps_dT,
        double* dEps_dRho,
        double* specific_neutrino_energy_loss,
        double* specific_neutrino_flux,
        double* mass_lost,
        const GFInstantaneousOutputs& instantaneous
    ) const;

};

std::unique_ptr<GFContext> make_gf_context(const GFContextType& type);

#endif