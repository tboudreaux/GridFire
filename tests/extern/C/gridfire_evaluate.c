/*
 * gf_evaluate and the network listing calls: derivatives against finite differences, channels summing to eps,
 * agreement with gf_evolve's final-state values, multi-zone against single-zone, input validation, and that
 * registering the network's species stops gf_evolve dropping any of them.
 *
 * Exits non-zero if any check fails.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NUM_SEEDS 8
#define NUM_ZONES 3

static const char* seed_names[NUM_SEEDS] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
static const double seed_Y[NUM_SEEDS] = {
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

static double rel_diff(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return scale == 0.0 ? 0.0 : fabs(a - b) / scale;
}

/* A context with the main sequence engine built from the seed composition. */
static void* make_context(enum GF_TYPE type) {
    void* ctx = gf_init(type);
    if (type == MULTI_ZONE) {
        CHECK_OK(gf_set_num_zones(type, ctx, NUM_ZONES), ctx, "gf_set_num_zones");
    }
    CHECK_OK(gf_register_species(ctx, NUM_SEEDS, seed_names), ctx, "gf_register_species");
    CHECK_OK(gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", seed_Y, NUM_SEEDS), ctx, "gf_construct_engine_from_policy");
    return ctx;
}

static double eps_at(void* ctx, const double* Y, size_t n, double T, double rho) {
    double eps, dT, dRho;
    const int ret = gf_evaluate(SINGLE_ZONE, ctx, Y, n, &T, &rho, &eps, &dT, &dRho, NULL, NULL, NULL);
    CHECK_OK(ret, ctx, "gf_evaluate");
    return eps;
}

int main(void) {
    /* --- Listing before the engine exists is an error, not a crash --- */
    void* empty = gf_init(SINGLE_ZONE);
    size_t n = 0;
    CHECK(gf_get_num_network_species(empty, &n) != GF_SUCCESS, "network listing worked before the engine was built");
    gf_free(SINGLE_ZONE, empty);

    void* ctx = make_context(SINGLE_ZONE);

    /* --- Network species and reactions --- */
    size_t num_network = 0, num_reactions = 0;
    CHECK_OK(gf_get_num_network_species(ctx, &num_network), ctx, "gf_get_num_network_species");
    CHECK_OK(gf_get_num_reactions(ctx, &num_reactions), ctx, "gf_get_num_reactions");
    CHECK(num_network > NUM_SEEDS, "the network has only %zu species", num_network);
    printf("network: %zu species, %zu reactions\n", num_network, num_reactions);

    int seeds_found = 0;
    for (size_t i = 0; i < num_network; ++i) {
        const char* name = NULL;
        int A = 0, Z = -1;
        double mass = 0.0;
        CHECK_OK(gf_get_network_species(ctx, i, &name, &A, &Z, &mass), ctx, "gf_get_network_species");
        CHECK(name != NULL && A >= 1 && Z >= 0 && Z <= A, "species %zu: name %s, A %d, Z %d", i, name ? name : "(null)", A, Z);
        CHECK(fabs(mass - A) < 0.1, "species %s: mass %.6f u for A = %d", name, mass, A);
        for (int s = 0; s < NUM_SEEDS; ++s) {
            if (name && strcmp(name, seed_names[s]) == 0) seeds_found++;
        }
    }
    CHECK(seeds_found == NUM_SEEDS, "only %d of the seed species are in the network listing", seeds_found);
    CHECK(gf_get_network_species(ctx, num_network, NULL, NULL, NULL, NULL) != GF_SUCCESS, "an out-of-range index was accepted");

    int has_pp = 0, has_triple_alpha = 0;
    for (size_t r = 0; r < num_reactions; ++r) {
        const char* id = NULL;
        CHECK_OK(gf_get_reaction_id(ctx, r, &id), ctx, "gf_get_reaction_id");
        if (id && strcmp(id, "p(p,e+)d") == 0) has_pp = 1;
        if (id && strcmp(id, "he4(2a,)c12") == 0) has_triple_alpha = 1;
    }
    CHECK(has_pp && has_triple_alpha, "reaction listing is missing p(p,e+)d or he4(2a,)c12");

    /* --- Derivatives against central finite differences, channels summing to eps --- */
    const double T = 1.5e7, rho = 150.0;
    size_t num_channels = 0;
    CHECK_OK(gf_get_num_energy_channels(ctx, &num_channels), ctx, "gf_get_num_energy_channels");
    double eps, dEps_dT, dEps_dRho, nu_loss, channels[32], dYdt[NUM_SEEDS];
    CHECK(num_channels <= 32, "too many channels");
    CHECK_OK(gf_evaluate(SINGLE_ZONE, ctx, seed_Y, NUM_SEEDS, &T, &rho, &eps, &dEps_dT, &dEps_dRho, &nu_loss, channels, dYdt),
             ctx, "gf_evaluate");
    printf("eps_nuc %.10e, dEps/dT %.6e, dEps/drho %.6e, nu loss %.6e\n", eps, dEps_dT, dEps_dRho, nu_loss);
    CHECK(eps > 0.0 && dEps_dT > 0.0 && dEps_dRho > 0.0 && nu_loss > 0.0, "solar-core values have the wrong sign");

    const double hT = 1e-4 * T, hRho = 1e-4 * rho;
    const double fd_T = (eps_at(ctx, seed_Y, NUM_SEEDS, T + hT, rho) - eps_at(ctx, seed_Y, NUM_SEEDS, T - hT, rho)) / (2.0 * hT);
    const double fd_rho = (eps_at(ctx, seed_Y, NUM_SEEDS, T, rho + hRho) - eps_at(ctx, seed_Y, NUM_SEEDS, T, rho - hRho)) / (2.0 * hRho);
    CHECK(rel_diff(dEps_dT, fd_T) < 1e-6, "dEps/dT %.10e vs finite difference %.10e", dEps_dT, fd_T);
    CHECK(rel_diff(dEps_dRho, fd_rho) < 1e-6, "dEps/drho %.10e vs finite difference %.10e", dEps_dRho, fd_rho);

    double channel_sum = 0.0;
    for (size_t c = 0; c < num_channels; ++c) channel_sum += channels[c];
    CHECK(rel_diff(channel_sum, eps) < 1e-12, "channels sum to %.12e, eps is %.12e", channel_sum, eps);
    CHECK(dYdt[0] < 0.0, "hydrogen is not being destroyed (dY/dt = %.3e)", dYdt[0]);

    /* --- Input validation: caller errors are GF_INVALID_STATE --- */
    const double bad_T = -1.0, nan_rho = NAN;
    double Y_bad[NUM_SEEDS];
    memcpy(Y_bad, seed_Y, sizeof Y_bad);
    Y_bad[3] = -1e-3;
    CHECK(gf_evaluate(SINGLE_ZONE, ctx, seed_Y, NUM_SEEDS, &bad_T, &rho, &eps, &dEps_dT, &dEps_dRho, NULL, NULL, NULL) == GF_INVALID_STATE, "negative T accepted");
    CHECK(gf_evaluate(SINGLE_ZONE, ctx, seed_Y, NUM_SEEDS, &T, &nan_rho, &eps, &dEps_dT, &dEps_dRho, NULL, NULL, NULL) == GF_INVALID_STATE, "NaN rho accepted");
    CHECK(gf_evaluate(SINGLE_ZONE, ctx, Y_bad, NUM_SEEDS, &T, &rho, &eps, &dEps_dT, &dEps_dRho, NULL, NULL, NULL) == GF_INVALID_STATE, "negative abundance accepted");
    CHECK(gf_evaluate(SINGLE_ZONE, ctx, seed_Y, NUM_SEEDS - 1, &T, &rho, &eps, &dEps_dT, &dEps_dRho, NULL, NULL, NULL) == GF_INVALID_STATE, "wrong species count accepted");

    /* --- Full-network registration: gf_evolve keeps every species, and gf_evaluate reproduces its final state --- */
    CHECK_OK(gf_register_network_species(ctx), ctx, "gf_register_network_species");
    CHECK_OK(gf_construct_solver_from_engine(ctx), ctx, "gf_construct_solver_from_engine");
    double* Y_full = calloc(num_network, sizeof(double));
    double* Y_out = calloc(num_network, sizeof(double));
    double* channels_out = calloc(num_channels, sizeof(double));
    for (size_t i = 0; i < num_network; ++i) {
        const char* name = NULL;
        gf_get_network_species(ctx, i, &name, NULL, NULL, NULL);
        for (int s = 0; s < NUM_SEEDS; ++s) {
            if (strcmp(name, seed_names[s]) == 0) Y_full[i] = seed_Y[s];
        }
    }
    double energy, evolve_dT, evolve_dRho, sp_nu_loss, sp_nu_flux, mass_lost, evolve_eps, dlnT, dlnRho, evolve_nu_rate;
    CHECK_OK(gf_evolve(SINGLE_ZONE, ctx, Y_full, num_network, &T, &rho, 3e15, 1e-12, Y_out, &energy, &evolve_dT,
                       &evolve_dRho, &sp_nu_loss, &sp_nu_flux, &mass_lost, &evolve_eps, &dlnT, &dlnRho, &evolve_nu_rate,
                       channels_out), ctx, "gf_evolve");
    CHECK(mass_lost == 0.0, "with the network registered gf_evolve still lost %.3e of mass", mass_lost);
    double eval_eps, eval_dT, eval_dRho;
    CHECK_OK(gf_evaluate(SINGLE_ZONE, ctx, Y_out, num_network, &T, &rho, &eval_eps, &eval_dT, &eval_dRho, NULL, NULL, NULL),
             ctx, "gf_evaluate after gf_evolve");
    CHECK(rel_diff(eval_eps, evolve_eps) < 1e-10, "gf_evaluate eps %.12e vs gf_evolve final eps %.12e", eval_eps, evolve_eps);
    CHECK(rel_diff(eval_dT, evolve_dT) < 1e-10, "gf_evaluate dEps/dT %.12e vs gf_evolve %.12e", eval_dT, evolve_dT);
    CHECK(rel_diff(eval_dRho, evolve_dRho) < 1e-10, "gf_evaluate dEps/drho %.12e vs gf_evolve %.12e", eval_dRho, evolve_dRho);

    /* --- Cost per zone --- */
    const int reps = 2000;
    double T_vary = T;
    /* Wall time: clock() would also count GridFire's background logging thread */
    struct timespec start, stop;
    timespec_get(&start, TIME_UTC);
    for (int k = 0; k < reps; ++k) {
        T_vary *= 1.0 + 1e-9;
        gf_evaluate(SINGLE_ZONE, ctx, Y_out, num_network, &T_vary, &rho, &eval_eps, &eval_dT, &eval_dRho, NULL, NULL, NULL);
    }
    timespec_get(&stop, TIME_UTC);
    printf("gf_evaluate (eps + derivatives, %zu species): %.1f us per zone\n", num_network,
           ((double)(stop.tv_sec - start.tv_sec) * 1e6 + (double)(stop.tv_nsec - start.tv_nsec) * 1e-3) / reps);
    free(Y_full);
    free(Y_out);
    free(channels_out);
    gf_free(SINGLE_ZONE, ctx);

    /* --- Multi-zone: one call equals one single-zone call per zone --- */
    void* grid = make_context(MULTI_ZONE);
    void* point = make_context(SINGLE_ZONE);
    double Y_zones[NUM_ZONES * NUM_SEEDS];
    const double T_zones[NUM_ZONES] = {1.2e7, 1.5e7, 2.5e7};
    const double rho_zones[NUM_ZONES] = {50.0, 150.0, 80.0};
    for (int z = 0; z < NUM_ZONES; ++z) memcpy(Y_zones + z * NUM_SEEDS, seed_Y, sizeof seed_Y);
    Y_zones[2 * NUM_SEEDS + 3] *= 2.0; /* a different composition in the last zone */
    double eps_z[NUM_ZONES], dT_z[NUM_ZONES], dRho_z[NUM_ZONES];
    CHECK_OK(gf_evaluate(MULTI_ZONE, grid, Y_zones, NUM_SEEDS, T_zones, rho_zones, eps_z, dT_z, dRho_z, NULL, NULL, NULL),
             grid, "gf_evaluate (multi-zone)");
    for (int z = 0; z < NUM_ZONES; ++z) {
        double e1, d1, r1;
        CHECK_OK(gf_evaluate(SINGLE_ZONE, point, Y_zones + z * NUM_SEEDS, NUM_SEEDS, &T_zones[z], &rho_zones[z], &e1, &d1, &r1, NULL, NULL, NULL),
                 point, "gf_evaluate (single zone)");
        CHECK(eps_z[z] == e1 && dT_z[z] == d1 && dRho_z[z] == r1, "zone %d differs from the single-zone call", z);
    }
    gf_free(MULTI_ZONE, grid);
    gf_free(SINGLE_ZONE, point);

    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All gf_evaluate checks passed\n");
    return 0;
}
