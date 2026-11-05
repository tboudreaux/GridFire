/*
 * Run-time configuration through the C API: defaults, building a GFConfig in C, the TOML round trip, rejection of
 * invalid or late configuration, and that the configuration reaches the engine (screening changes the energy).
 *
 * Usage: gf_extern_c_config <scratch directory>. Exits non-zero on the first failed check.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define NUM_SPECIES 8
#define MAX_CHANNELS 16

static const char* species_names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
static const double abundances[NUM_SPECIES] = {
    0.702616602672027, 9.74791583949078e-06, 0.06895512307276903, 0.00025,
    7.855418029399437e-05, 0.0006014411598306529, 8.103062886768109e-05, 2.151340851063217e-05
};

static int failures = 0;

#define CHECK(cond, ...) do { \
    if (!(cond)) { \
        printf("FAIL (line %d): ", __LINE__); \
        printf(__VA_ARGS__); \
        printf("\n"); \
        failures++; \
    } \
} while (0)

#define CHECK_OK(ret, ctx, what) CHECK((ret) == GF_SUCCESS || (ret) == FDSSE_SUCCESS, "%s: %s [%s]", what, \
    gf_get_last_error_message(ctx), gf_error_code_to_string(ret))

static int configs_equal(const struct GFConfig* a, const struct GFConfig* b) {
    return a->solver_abs_tol == b->solver_abs_tol && a->solver_rel_tol == b->solver_rel_tol &&
           a->trigger_off_diagonal_threshold == b->trigger_off_diagonal_threshold &&
           a->trigger_timestep_collapse_ratio == b->trigger_timestep_collapse_ratio &&
           a->trigger_max_convergence_failures == b->trigger_max_convergence_failures &&
           a->trigger_boundary_flux_relative_threshold == b->trigger_boundary_flux_relative_threshold &&
           a->trigger_boundary_flux_absolute_threshold == b->trigger_boundary_flux_absolute_threshold &&
           a->engine_adaptive_relative_culling_threshold == b->engine_adaptive_relative_culling_threshold &&
           a->engine_screening_blend_weak == b->engine_screening_blend_weak &&
           a->engine_screening_blend_strong == b->engine_screening_blend_strong &&
           a->solver_clamp_negative_iterates == b->solver_clamp_negative_iterates &&
           a->solver_max_step_failure_restarts == b->solver_max_step_failure_restarts &&
           a->engine_screening == b->engine_screening && a->engine_jacobian == b->engine_jacobian &&
           a->engine_split_hydrogen_pp_cno == b->engine_split_hydrogen_pp_cno;
}

/* Build a single-zone context with config (NULL: defaults), evolve 1e10 s at solar-core conditions, return eps_nuc. */
static double evolved_eps_nuc(const struct GFConfig* config) {
    void* ctx = gf_init(SINGLE_ZONE);
    int ret;
    if (config) {
        ret = gf_set_config(ctx, config);
        CHECK_OK(ret, ctx, "gf_set_config");
    }
    ret = gf_register_species(ctx, NUM_SPECIES, species_names);
    CHECK_OK(ret, ctx, "gf_register_species");
    ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", abundances, NUM_SPECIES);
    CHECK_OK(ret, ctx, "gf_construct_engine_from_policy");
    ret = gf_construct_solver_from_engine(ctx);
    CHECK_OK(ret, ctx, "gf_construct_solver_from_engine");

    const double T = 1.5e7, rho = 150.0;
    double Y_out[NUM_SPECIES], energy, dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost;
    double eps_nuc = NAN, dlnEps_dlnT, dlnEps_dlnRho, nu_loss_rate, eps_by_channel[MAX_CHANNELS];
    ret = gf_evolve(SINGLE_ZONE, ctx, abundances, NUM_SPECIES, &T, &rho, 1e10, 1e-12, Y_out, &energy, &dEps_dT,
                    &dEps_dRho, &nu_loss, &nu_flux, &mass_lost, &eps_nuc, &dlnEps_dlnT, &dlnEps_dlnRho, &nu_loss_rate,
                    eps_by_channel);
    CHECK_OK(ret, ctx, "gf_evolve");
    gf_free(SINGLE_ZONE, ctx);
    return eps_nuc;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("usage: %s <scratch directory>\n", argv[0]);
        return 2;
    }

    /* --- Defaults --- */
    struct GFConfig defaults;
    CHECK(gf_config_defaults(&defaults) == GF_SUCCESS, "gf_config_defaults");
    CHECK(defaults.engine_screening == GF_SCREENING_BLENDED, "default screening is %d", defaults.engine_screening);
    CHECK(defaults.engine_jacobian == GF_JACOBIAN_ANALYTIC, "default Jacobian is %d", defaults.engine_jacobian);
    CHECK(defaults.solver_abs_tol > 0.0 && defaults.solver_rel_tol > 0.0, "default tolerances");

    void* ctx = gf_init(SINGLE_ZONE);
    struct GFConfig current;
    CHECK_OK(gf_get_config(ctx, &current), ctx, "gf_get_config");
    CHECK(configs_equal(&current, &defaults), "a new context does not have the default configuration");

    /* --- A configuration built in C round-trips --- */
    struct GFConfig custom = defaults;
    custom.engine_screening = GF_SCREENING_BARE;
    custom.engine_jacobian = GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION;
    custom.solver_rel_tol = 1e-7;
    custom.engine_screening_blend_weak = 0.4;
    custom.engine_screening_blend_strong = 1.2;
    custom.solver_max_step_failure_restarts = 5;
    custom.engine_split_hydrogen_pp_cno = 0;
    CHECK_OK(gf_set_config(ctx, &custom), ctx, "gf_set_config");
    CHECK_OK(gf_get_config(ctx, &current), ctx, "gf_get_config");
    CHECK(configs_equal(&current, &custom), "gf_get_config does not return what gf_set_config set");

    /* --- Invalid values are rejected and leave the configuration unchanged --- */
    struct GFConfig invalid = custom;
    invalid.engine_jacobian = 7;
    CHECK(gf_set_config(ctx, &invalid) == GF_CONFIG_ERROR, "an invalid Jacobian method was accepted");
    invalid = custom;
    invalid.engine_screening_blend_strong = 2.5; /* above the weak-screening cap */
    CHECK(gf_set_config(ctx, &invalid) == GF_CONFIG_ERROR, "a blend window above H = 2 was accepted");
    invalid = custom;
    invalid.solver_abs_tol = -1.0;
    CHECK(gf_set_config(ctx, &invalid) == GF_CONFIG_ERROR, "a negative tolerance was accepted");
    CHECK_OK(gf_get_config(ctx, &current), ctx, "gf_get_config");
    CHECK(configs_equal(&current, &custom), "a rejected configuration changed the context");

    /* --- TOML round trip --- */
    char path[1024];
    snprintf(path, sizeof path, "%s/gf_c_config_test.toml", argv[1]);
    CHECK_OK(gf_save_config(ctx, path), ctx, "gf_save_config");
    void* loaded_ctx = gf_init(SINGLE_ZONE);
    CHECK_OK(gf_load_config(loaded_ctx, path), loaded_ctx, "gf_load_config");
    CHECK_OK(gf_get_config(loaded_ctx, &current), loaded_ctx, "gf_get_config");
    CHECK(configs_equal(&current, &custom), "the configuration read back from TOML differs from the one saved");
    snprintf(path, sizeof path, "%s/does_not_exist.toml", argv[1]);
    CHECK(gf_load_config(loaded_ctx, path) == GF_CONFIG_ERROR, "loading a missing file did not fail with GF_CONFIG_ERROR");
    gf_free(SINGLE_ZONE, loaded_ctx);

    /* --- The configuration cannot change once the engine exists --- */
    CHECK_OK(gf_register_species(ctx, NUM_SPECIES, species_names), ctx, "gf_register_species");
    CHECK_OK(gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", abundances, NUM_SPECIES), ctx, "gf_construct_engine_from_policy");
    CHECK(gf_set_config(ctx, &defaults) == GF_CONFIG_ERROR, "gf_set_config after the engine was built was accepted");
    size_t num_channels = 0;
    CHECK_OK(gf_get_num_energy_channels(ctx, &num_channels), ctx, "gf_get_num_energy_channels");
    int has_split_channel = 0;
    for (size_t i = 0; i < num_channels; ++i) {
        const char* name = NULL;
        CHECK_OK(gf_get_energy_channel_name(ctx, i, &name), ctx, "gf_get_energy_channel_name");
        if (name && strcmp(name, "H_pp") == 0) has_split_channel = 1;
    }
    CHECK(!has_split_channel, "engine_split_hydrogen_pp_cno = 0 did not reach the engine");
    gf_free(SINGLE_ZONE, ctx);

    /* --- The configuration reaches the engine: weak screening raises pp burning by a few percent --- */
    struct GFConfig bare = defaults;
    bare.engine_screening = GF_SCREENING_BARE;
    struct GFConfig autodiff = defaults;
    autodiff.engine_jacobian = GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION;
    const double eps_default = evolved_eps_nuc(NULL);
    const double eps_bare = evolved_eps_nuc(&bare);
    const double eps_autodiff = evolved_eps_nuc(&autodiff);
    printf("eps_nuc: default %.10e, bare screening %.10e, CppAD Jacobian %.10e\n", eps_default, eps_bare, eps_autodiff);
    CHECK(eps_default / eps_bare > 1.02, "screening did not change eps_nuc (ratio %.6f)", eps_default / eps_bare);
    /* The Jacobian only steers Newton, so both methods reach the same solution to within the tolerances */
    CHECK(fabs(eps_autodiff / eps_default - 1.0) < 1e-3, "Jacobian methods disagree (ratio %.8f)", eps_autodiff / eps_default);

    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All C configuration checks passed\n");
    return 0;
}
