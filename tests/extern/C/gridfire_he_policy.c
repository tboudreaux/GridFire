/*
 * H_HE_BURNING_POLICY through the C API: built from a composition with only H-1 and He-4 (no Ne-20 or Mg-24), the full
 * network is registered, and a helium-burning zone turns He-4 into C-12 and O-16. Exits non-zero on failure.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CHANNELS 16

static int failures = 0;

#define CHECK(cond, ...) do { \
    if (!(cond)) { \
        printf("FAIL (line %d): ", __LINE__); \
        printf(__VA_ARGS__); \
        printf("\n"); \
        failures++; \
    } \
} while (0)

static long index_of(void* ctx, size_t n, const char* wanted) {
    for (size_t i = 0; i < n; ++i) {
        const char* name = NULL;
        gf_get_network_species(ctx, i, &name, NULL, NULL, NULL);
        if (name && strcmp(name, wanted) == 0) return (long)i;
    }
    return -1;
}

int main(void) {
    const char* seeds[2] = {"H-1", "He-4"};
    const double seed_Y[2] = {0.75 / 1.007825, 0.25 / 4.002603};

    void* ctx = gf_init(SINGLE_ZONE);
    int ret = gf_register_species(ctx, 2, seeds);
    CHECK(ret == GF_SUCCESS || ret == FDSSE_SUCCESS, "gf_register_species");
    ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", seed_Y, 2);
    CHECK(ret != GF_SUCCESS, "MAIN_SEQUENCE_POLICY accepted a composition without its seed species");
    ret = gf_construct_engine_from_policy(ctx, "H_HE_BURNING_POLICY", seed_Y, 2);
    CHECK(ret == GF_SUCCESS, "H_HE_BURNING_POLICY: %s [%s]", gf_get_last_error_message(ctx), gf_error_code_to_string(ret));

    size_t n = 0, num_reactions = 0;
    gf_get_num_network_species(ctx, &n);
    gf_get_num_reactions(ctx, &num_reactions);
    printf("H_HE_BURNING_POLICY network: %zu species, %zu reactions\n", n, num_reactions);
    int has_triple_alpha = 0;
    for (size_t r = 0; r < num_reactions; ++r) {
        const char* id = NULL;
        gf_get_reaction_id(ctx, r, &id);
        if (id && strcmp(id, "he4(2a,)c12") == 0) has_triple_alpha = 1;
    }
    CHECK(has_triple_alpha, "no triple-alpha in the network");

    /* Helium burning on the full network */
    CHECK(gf_register_network_species(ctx) == GF_SUCCESS, "gf_register_network_species");
    CHECK(gf_construct_solver_from_engine(ctx) == GF_SUCCESS, "gf_construct_solver_from_engine");
    const long iHe = index_of(ctx, n, "He-4"), iC = index_of(ctx, n, "C-12"), iO = index_of(ctx, n, "O-16");
    CHECK(iHe >= 0 && iC >= 0 && iO >= 0, "He-4, C-12 or O-16 missing from the network");
    double* Y_in = calloc(n, sizeof(double));
    double* Y_out = calloc(n, sizeof(double));
    Y_in[iHe] = 1.0 / 4.002603;
    const double T = 1.5e8, rho = 1e4;
    double energy, dT, dRho, nl, nf, lost, eps, l1, l2, nr, ch[MAX_CHANNELS];
    ret = gf_evolve(SINGLE_ZONE, ctx, Y_in, n, &T, &rho, 3e13, 1e-6, Y_out, &energy, &dT, &dRho, &nl, &nf, &lost, &eps,
                    &l1, &l2, &nr, ch);
    CHECK(ret == GF_SUCCESS, "gf_evolve: %s [%s]", gf_get_last_error_message(ctx), gf_error_code_to_string(ret));
    printf("He-4 %.4e -> %.4e, C-12 -> %.4e, O-16 -> %.4e, eps %.4e erg/g/s\n", Y_in[iHe], Y_out[iHe], Y_out[iC], Y_out[iO], eps);
    CHECK(Y_out[iHe] < Y_in[iHe] && Y_out[iC] > 0.0 && Y_out[iO] > 0.0 && eps > 0.0 && lost == 0.0, "helium was not burnt");

    free(Y_in);
    free(Y_out);
    gf_free(SINGLE_ZONE, ctx);
    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All H_HE_BURNING_POLICY checks passed\n");
    return 0;
}
