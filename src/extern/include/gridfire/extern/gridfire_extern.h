#ifndef GF_GRIDFIRE_EXTERN_H
#define GF_GRIDFIRE_EXTERN_H

#include <stddef.h>

/*
 * Thread safety (OpenMP builds, which include the published wheels)
 *
 *   - Build and free contexts from serial code, while no other GridFire call is running: gf_init through
 *     gf_construct_solver_from_engine, and gf_free. Engine construction and destruction register and remove CppAD
 *     atomic functions in a process-wide table, which CppAD does not allow from several threads.
 *   - Once built, separate contexts may be used concurrently: gf_evolve, gf_evaluate and gf_evaluate_reactions on
 *     different contexts, from any threads. These may be OpenMP threads of any runtime (including GCC's libgomp, as
 *     in a gfortran -fopenmp program, while GridFire carries its own) or plain threads; up to 512 at once
 *     (CPPAD_MAX_NUM_THREADS). Results are bit-identical to running the same contexts serially (test
 *     extern_c_threads, run with both libomp and libgomp hosts).
 *   - A single context must never be used by two threads at once: it holds the solver state and scratch memory.
 *   - Make the first gf_init call outside any parallel region; it prepares CppAD for use from several threads.
 *   - MULTI_ZONE contexts spread their zones over GridFire's own OpenMP threads; call them from serial code (calling
 *     them from inside a host parallel region is untested and oversubscribes the cores).
 *
 * Builds without OpenMP: contexts are independent, but use them from one thread.
 */

#ifdef __cplusplus
extern "C" {
#endif
    enum GF_TYPE {
        SINGLE_ZONE = 1001,
        MULTI_ZONE = 1002
    };


    enum FDSSE_ERROR_CODES {
        FDSSE_NON_4DSTAR_ERROR = -102,
        FDSSE_UNKNOWN_ERROR = -101,
        FDSSE_SUCCESS = 1,
        FDSSE_UNKNOWN_SYMBOL_ERROR = 100,
        FDSSE_SPECIES_ERROR = 101,
        FDSSE_INVALID_COMPOSITION_ERROR = 102,
        FDSSE_COMPOSITION_ERROR = 103
    };

    enum GF_ERROR_CODES {
        GF_NON_GRIDFIRE_ERROR = -2,
        GF_UNKNOWN_ERROR = -1,
        GF_SUCCESS = 0,

        GF_INVALID_QSE_SOLUTION_ERROR = 5,
        GF_FAILED_TO_PARTITION_ENGINE_ERROR = 6,
        GF_NETWORK_RESIZED_ERROR = 7,
        GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR = 8,
        GF_BAD_COLLECTION_ERROR = 9,
        GF_BAD_RHS_ENGINE_ERROR = 10,
        GF_STALE_JACOBIAN_ERROR = 11,
        GF_UNINITIALIZED_JACOBIAN_ERROR = 12,
        GF_UNKNOWN_JACOBIAN_ERROR = 13,
        GF_JACOBIAN_ERROR = 14,
        GF_ENGINE_ERROR = 15,

        GF_MISSING_BASE_REACTION_ERROR = 16,
        GF_MISSING_SEED_SPECIES_ERROR = 17,
        GF_MISSING_KEY_REACTION_ERROR = 18,
        GF_POLICY_ERROR = 19,

        GF_REACTION_PARSING_ERROR = 20,
        GF_REACTION_ERROR = 21,

        GF_SINGULAR_JACOBIAN_ERROR = 22,
        GF_ILL_CONDITIONED_JACOBIAN_ERROR = 23,
        GF_CVODE_SOLVER_FAILURE_ERROR = 24,
        GF_KINSOL_SOLVER_FAILURE_ERROR = 25,
        GF_SUNDIALS_ERROR = 26,
        GF_SOLVER_ERROR = 27,

        GF_HASHING_ERROR = 28,
        GF_UTILITY_ERROR = 29,

        GF_DEBUG_ERROR = 30,

        GF_GRIDFIRE_ERROR = 31,
        GF_UNINITIALIZED_INPUT_MEMORY_ERROR = 32,
        GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR = 33,

        GF_INVALID_NUM_SPECIES = 34,

        GF_INVALID_TIMESTEPS = 35,
        GF_UNKNOWN_FREE_TYPE = 36,

        GF_INVALID_TYPE = 37,

        GF_CONFIG_ERROR = 38,
        GF_INVALID_STATE = 39, /**< Caller error: non-physical T, rho or abundances, or a species outside the network */
        GF_ZONE_FAILURE = 40, /**< MULTI_ZONE gf_evolve: one or more zones failed, the rest are valid (gf_get_zone_status) */
    };

    /** Electron screening model (GridFireConfig::engine.screening). */
    enum GF_SCREENING_TYPE {
        GF_SCREENING_BARE = 0, /**< No screening */
        GF_SCREENING_WEAK = 1, /**< Salpeter (1954) weak screening */
        GF_SCREENING_CD09 = 2, /**< Chugunov & DeWitt (2009), weak through strong coupling */
        GF_SCREENING_BLENDED = 3 /**< Salpeter handing over to CD09 as coupling grows (see engine_screening_blend_*); the default */
    };

    /** How the engine computes the species Jacobian (GridFireConfig::engine.jacobian). */
    enum GF_JACOBIAN_METHOD {
        GF_JACOBIAN_ANALYTIC = 0,                 /**< Analytic, from the precomputed reactions; the default */
        GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION = 1 /**< CppAD automatic differentiation */
    };

    /**
     * Run-time configuration, the C mirror of gridfire::config::GridFireConfig (field names follow its nesting).
     *
     * Start from gf_config_defaults, change what you need, and pass it to gf_set_config before
     * gf_construct_engine_from_policy. Booleans are ints (0 = false). The energy channel fuel list has no field here;
     * it keeps its current value (the default, or what gf_load_config read).
     *
     * Fortran: type(GFConfig) in gridfire_mod has the same layout.
     */
    struct GFConfig {
        /* solver.pointSolver */
        double solver_abs_tol;
        double solver_rel_tol;
        /* solver.pointSolver.trigger */
        double trigger_off_diagonal_threshold;
        double trigger_timestep_collapse_ratio;
        double trigger_max_convergence_failures;
        double trigger_boundary_flux_relative_threshold;
        double trigger_boundary_flux_absolute_threshold;
        /* engine.views.adaptiveEngineView */
        double engine_adaptive_relative_culling_threshold;
        /* engine.screeningBlend: BLENDED screening uses Salpeter below this weak-screening exponent H ... */
        double engine_screening_blend_weak;
        /* ... and CD09 above this one (0 < weak < strong <= 2) */
        double engine_screening_blend_strong;
        /* solver.pointSolver */
        int solver_clamp_negative_iterates;
        int solver_max_step_failure_restarts;
        /* engine */
        int engine_screening;              /**< enum GF_SCREENING_TYPE */
        int engine_jacobian;               /**< enum GF_JACOBIAN_METHOD */
        int engine_split_hydrogen_pp_cno;  /**< engine.energyChannels.splitHydrogenPPCNO */
        int engine_use_reverse_reactions;  /**< engine.useReverseReactions: reverse rates by detailed balance (default 0) */
    };

    char* gf_get_last_error_message(void* ptr);

    char* gf_error_code_to_string(int error_code);

    void* gf_init(const enum GF_TYPE type);

    int gf_free(const enum GF_TYPE type, void *ctx);

    int gf_set_num_zones(const enum GF_TYPE type, void* ptr, const size_t num_zones);

    int gf_register_species(void* ptr, const int num_species, const char** species_names);

    /** Fill config with GridFire's defaults (the same values a C++ GridFireConfig{} has). */
    int gf_config_defaults(struct GFConfig* config);

    /**
     * Use config for the engine and solver this context builds. Must be called before
     * gf_construct_engine_from_policy (GF_CONFIG_ERROR otherwise, or for out-of-range values).
     */
    int gf_set_config(void* ptr, const struct GFConfig* config);

    /** The configuration this context uses (defaults, or whatever gf_set_config / gf_load_config set). */
    int gf_get_config(void* ptr, struct GFConfig* config);

    /**
     * Read the configuration from a TOML file, with the same rules as gf_set_config. The file must be complete, in
     * the format gf_save_config writes, so the easiest way to make one is to save the defaults (gf_save_config, or
     * the gf_generate_config_file tool) and edit them.
     */
    int gf_load_config(void* ptr, const char* path);

    /** Write this context's configuration as TOML (the format gf_load_config reads). */
    int gf_save_config(void* ptr, const char* path);

    /** GridFire's version, e.g. "1.2.0". The string is static. */
    const char* gf_get_version(void);

    /**
     * What this context will compute with, for a host code to record next to its results (e.g. in output headers):
     * one "key: value" line each for the GridFire version, the policy and network (size and a hash of its reactions,
     * once the engine exists), the strong and weak rate libraries with hashes of their compiled data, rate
     * multipliers, screening, Jacobian method, passive species, energy channels and solver tolerances. Two runs with
     * the same description used the same code version, rates, network and physics options. The string is owned by
     * the context and stays valid until the next gf_get_provenance call or gf_free.
     */
    int gf_get_provenance(void* ptr, const char** description);

    /**
     * Species carried in the composition without any reactions, e.g. the metals a stellar code lumps together and
     * represents by one isotope. They count towards Y_e and screening but never change, and the network is not grown
     * from them. They must also be registered (gf_register_species) and none may take part in a reaction of the
     * network (gf_construct_engine_from_policy fails with GF_BAD_COLLECTION_ERROR otherwise). Must be called before
     * gf_construct_engine_from_policy (GF_CONFIG_ERROR otherwise); num_species = 0 clears the list.
     */
    int gf_set_passive_species(void* ptr, size_t num_species, const char** species_names);

    /**
     * Convert mass fractions to the molar abundances every other call takes. GridFire works in molar abundances only;
     * these helpers are for callers that keep mass fractions, and change nothing in the context.
     *
     * X and Y hold num_species values per zone (registered species order, zone-major), for num_zones zones. Each zone
     * goes through the same conversion as the C++ library (fourdst::composition::buildCompositionFromMassFractions):
     * Y_i = X_i / m_i, with m_i the atomic mass [u] reported by gf_get_network_species. Each zone's mass fractions must
     * sum to one to within 1e-10 (GF_INVALID_STATE otherwise, naming the zone). Needs gf_register_species first.
     */
    int gf_mass_fractions_to_molar(void* ptr, const double* X, size_t num_species, size_t num_zones, double* Y);

    /**
     * Convert molar abundances to mass fractions, X_i = Y_i m_i / sum_j Y_j m_j for each zone (the C++ library's
     * Composition::getMassFraction), so each zone's X sums to one. Burning lowers sum_j Y_j m_j (the binding energy
     * released, about 0.7% of the rest mass for complete hydrogen burning); the normalisation absorbs it, so mass
     * fractions taken after gf_evolve are the true fractions of the material's current rest mass. Register the whole
     * network (gf_register_network_species) so that no burnt mass sits in a species the caller does not track.
     */
    int gf_molar_to_mass_fractions(void* ptr, const double* Y, size_t num_species, size_t num_zones, double* X);

    /**
     * Build the engine from a named policy and the starting abundances (registered species order).
     *
     *   "MAIN_SEQUENCE_POLICY" - pp chains and CNO cycle, a network grown from the starting composition; requires
     *                            H-1, He-3, He-4, C-12, N-14, O-16, Ne-20 and Mg-24.
     *   "H_HE_BURNING_POLICY"  - pp chains, CNO cycle and helium burning (triple-alpha, alpha captures to Mg-24,
     *                            N-14 -> Ne-22). The network is the same whatever the starting composition, and no
     *                            species is required; the composition only has to lie within the network.
     */
    int gf_construct_engine_from_policy(void* ptr, const char* policy_name, const double *abundances, size_t num_species);

    int gf_construct_solver_from_engine(void* ptr);

    /**
     * Evolve one zone (SINGLE_ZONE) or every zone (MULTI_ZONE) for tMax seconds.
     *
     * For MULTI_ZONE, T, rho, and every scalar output are arrays with one entry per zone, and Y_in / Y_out hold
     * num_species abundances per zone.
     *
     * Outputs:
     *   energy_out                    - specific energy released over tMax [erg/g]
     *   dEps_dT, dEps_dRho            - partial derivatives of eps_nuc at the final state (fixed composition)
     *   specific_neutrino_energy_loss - neutrino energy lost over tMax [erg/g]
     *   specific_neutrino_flux        - neutrinos emitted over tMax [1/g]
     *   mass_lost                     - mass in species outside the registered species list
     *   eps_nuc                       - instantaneous energy generation rate at the final state, net of neutrinos [erg/g/s]
     *   dlnEps_dlnT, dlnEps_dlnRho    - logarithmic derivatives of eps_nuc (0 if eps_nuc == 0)
     *   neutrino_energy_loss_rate     - instantaneous reaction neutrino loss rate at the final state [erg/g/s]
     *   eps_by_channel                - eps_nuc split by energy channel, num_channels entries per zone
     *                                   (see gf_get_num_energy_channels); sums to eps_nuc
     *
     * eps_nuc, its derivatives, and eps_by_channel are all evaluated at the same final state.
     *
     * Calls share no integrator state: each one integrates from Y_in at t = 0 with a fresh step-size and error
     * history (the solver only keeps its allocated memory between calls), so the result depends on the inputs, the
     * configuration and the rate multipliers alone, never on earlier calls. The network is fixed for the life of the
     * context (see gf_get_num_network_species), so nothing is pruned or repartitioned between calls either.
     */
    int gf_evolve(
        enum GF_TYPE type,
        void* ptr,
        const void* Y_in,
        size_t num_species,
        const void* T,
        const void* rho,
        double tMax,
        double dt0,
        void* Y_out,
        void* energy_out,
        void* dEps_dT,
        void* dEps_dRho,
        void* specific_neutrino_energy_loss,
        void* specific_neutrino_flux,
        void* mass_lost,
        void* eps_nuc,
        void* dlnEps_dlnT,
        void* dlnEps_dlnRho,
        void* neutrino_energy_loss_rate,
        void* eps_by_channel
    );

    /**
     * MULTI_ZONE only: burn groups of zones as one shared composition in every later gf_evolve call (a convective zone
     * under instantaneous mixing).
     *
     *   group - one entry per zone: -1 to burn the zone on its own, or a group id >= 0; zones with the same id share
     *           one composition (the mass-weighted average of their inputs), which burns at the mass-weighted mean of
     *           their rates. NULL clears the groups (every zone on its own, the default).
     *   mass  - one entry per zone: its mass (any units, only ratios within a group matter); must be positive for every
     *           grouped zone. May be NULL only when no zone is grouped.
     *
     * Every member of a group gets the group's Y_out, but its own energy, neutrino losses, eps_nuc, derivatives and
     * channel energies, evaluated at its own T and rho. Members use the zone's tMax and dt0 (one per call anyway).
     * A group that fails fails as a unit: every member gets the group's status (gf_get_zone_status). Groups stay set
     * until changed or cleared; the zone count must not change while they are set.
     *
     * Returns GF_INVALID_STATE for an id below -1 or a grouped zone without a positive mass.
     */
    int gf_set_shared_composition_groups(void* ptr, const int* group, const double* mass);

    /**
     * Group ids for gf_set_shared_composition_groups from a per-zone flag (e.g. 1 = convective), for zones ordered in
     * space: each contiguous run of non-zero flags is one group (ids 0, 1, ...), every other zone gets -1.
     */
    int gf_shared_composition_groups_from_mask(const int* flagged, size_t num_zones, int* group);

    /**
     * Outcome of each zone in the latest gf_evolve call (1 entry for SINGLE_ZONE, gf_set_num_zones entries for
     * MULTI_ZONE). Either array may be NULL.
     *
     *   status       - GF_SUCCESS, or the zone's own error code. GF_INVALID_STATE (a non-physical input) is a caller
     *                  error; solver failures such as GF_CVODE_SOLVER_FAILURE_ERROR are worth retrying with a shorter
     *                  step.
     *   time_reached - how far the zone was integrated [s]: tMax on success, the last accepted step on failure.
     *
     * Failure semantics of gf_evolve: a failed zone's Y_out holds its input abundances and its other outputs are 0;
     * every other zone's outputs are complete. A SINGLE_ZONE call returns the zone's error code; a MULTI_ZONE call
     * returns GF_ZONE_FAILURE if any zone failed (gf_get_last_error_message names the first few). Input errors found
     * before integrating (GF_INVALID_STATE, a wrong species count) fail the whole call and write nothing. Nothing is
     * printed to stdout.
     */
    int gf_get_zone_status(void* ptr, int* status, double* time_reached);

    /**
     * Instantaneous nuclear energy generation at fixed composition, for every zone, without integrating.
     *
     * No solver is set up or touched, so this needs only gf_construct_engine_from_policy and can be called any number
     * of times between (or instead of) gf_evolve calls. For MULTI_ZONE every argument is an array with one entry per
     * zone, Y holds num_species abundances per zone, and dYdt likewise.
     *
     * Inputs:
     *   Y                         - molar abundances in registered species order (gf_register_species, or
     *                               gf_register_network_species); network species not registered are taken as absent
     *   T, rho                    - temperature [K] and density [g/cm^3]
     * Outputs (eps_nuc, dEps_dT, dEps_dRho required; the rest may be NULL):
     *   eps_nuc                   - energy generation rate, net of reaction neutrino losses [erg/g/s]
     *   dEps_dT, dEps_dRho        - partial derivatives of eps_nuc at fixed composition [erg/g/s/K], [erg/g/s / (g/cm^3)],
     *                               always returned as absolute derivatives (also when eps_nuc = 0)
     *   neutrino_energy_loss_rate - reaction neutrino energy loss rate [erg/g/s]
     *   eps_by_channel            - eps_nuc by energy channel, gf_get_num_energy_channels entries per zone
     *   dYdt                      - dY/dt for each registered species [mol/g/s]
     *
     * Returns GF_INVALID_STATE (a caller error) for a non-positive or non-finite T or rho, a negative or non-finite
     * abundance, or a positive abundance of a registered species that is not in the network.
     */
    int gf_evaluate(
        enum GF_TYPE type,
        void* ptr,
        const double* Y,
        size_t num_species,
        const double* T,
        const double* rho,
        double* eps_nuc,
        double* dEps_dT,
        double* dEps_dRho,
        double* neutrino_energy_loss_rate,
        double* eps_by_channel,
        double* dYdt
    );

    /**
     * What each reaction contributes at fixed composition, for every zone, without integrating (like gf_evaluate).
     *
     * Each output holds gf_get_num_reactions entries per zone, in gf_get_reaction_id order (zone-major for MULTI_ZONE);
     * any may be NULL. Over the reactions they sum to gf_evaluate's totals at the same state:
     *   molar_flow                - net molar flow of each reaction [mol/g/s]
     *   energy_rate               - energy released, net of neutrinos [erg/g/s]; sums to eps_nuc
     *   neutrino_energy_loss_rate - energy carried off by reaction neutrinos [erg/g/s]; sums to the neutrino loss rate
     *   neutrino_rate             - reaction neutrinos emitted [1/g/s]; the solar sources are the entries for
     *                               p(p,e+)d (pp), p(e-p,)d (pep), he3(p,e+)he4 (hep), be7(e-,)li7, b8(,e+ a)he4 (B-8),
     *                               n13(,e+)c13, o15(,e+)n15 and f17(,e+)o17
     *
     * Called with gf_evolve's Y_out, T and rho, it gives the per-reaction values at the end of that step (gf_evaluate
     * reproduces gf_evolve's final-state eps exactly). Inputs and errors are as for gf_evaluate.
     */
    int gf_evaluate_reactions(
        enum GF_TYPE type,
        void* ptr,
        const double* Y,
        size_t num_species,
        const double* T,
        const double* rho,
        double* molar_flow,
        double* energy_rate,
        double* neutrino_energy_loss_rate,
        double* neutrino_rate
    );

    /**
     * The network built by gf_construct_engine_from_policy is fixed for the life of the context: its species, their
     * order and its reactions never change, so the indices below stay valid. Network species not registered with
     * gf_register_species are dropped (reported as mass_lost) by gf_evolve; gf_register_network_species avoids that.
     */

    /** Number of species in the network. */
    int gf_get_num_network_species(void* ptr, size_t* num_species);

    /**
     * Network species `index` (0-based, network order): name (e.g. "N-15", owned by the context), mass number A,
     * charge Z and atomic mass [u], the masses GridFire uses. Any output pointer may be NULL.
     */
    int gf_get_network_species(void* ptr, size_t index, const char** name, int* A, int* Z, double* mass);

    /**
     * Register the network's own species, in network order (replacing gf_register_species), so gf_evolve and
     * gf_evaluate take and return the full network composition and no species is dropped between calls.
     */
    int gf_register_network_species(void* ptr);

    /** Number of reactions in the network. */
    int gf_get_num_reactions(void* ptr, size_t* num_reactions);

    /** Identifier of reaction `index` (0-based), e.g. "p(p,e+)d" or "he4(2a,)c12". The string is owned by the context. */
    int gf_get_reaction_id(void* ptr, size_t index, const char** id);

    /**
     * Scale one reaction's rate by a constant factor, e.g. to sample rate uncertainties in Monte Carlo runs.
     *
     * reaction_id is an identifier from gf_get_reaction_id. The factor multiplies that reaction's rate in everything
     * the context computes afterwards (gf_evolve, gf_evaluate, gf_evaluate_reactions, for every zone): energy
     * generation, abundances, neutrino losses, derivatives and both Jacobian methods. A reverse rate is scaled with
     * its forward rate, so detailed balance holds. The network is not rebuilt, so changing factors between runs is
     * cheap. Factors stay until changed; 1 removes one, 0 switches the reaction off. gf_get_provenance lists them.
     *
     * Needs gf_construct_engine_from_policy first (GF_CONFIG_ERROR otherwise). Returns GF_CONFIG_ERROR for a negative
     * or non-finite factor and GF_REACTION_ERROR for a reaction not in the network. Do not call while another thread
     * is evolving or evaluating with the same context.
     */
    int gf_set_rate_multiplier(void* ptr, const char* reaction_id, double factor);

    /** The factor set for a reaction (1 if none). Errors as for gf_set_rate_multiplier. */
    int gf_get_rate_multiplier(void* ptr, const char* reaction_id, double* factor);

    /** Remove every rate multiplier (all rates back to their library values). Needs the engine. */
    int gf_clear_rate_multipliers(void* ptr);

    /**
     * Number of energy channels reported per zone by gf_evolve (valid once the engine has been constructed).
     */
    int gf_get_num_energy_channels(void* ptr, size_t* num_channels);

    /**
     * Name of energy channel `index` (e.g. "H_pp", "H_CNO", "He"). The string is owned by the context.
     */
    int gf_get_energy_channel_name(void* ptr, size_t index, const char** name);

#ifdef __cplusplus
}
#endif

#endif
