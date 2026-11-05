/*
 * Shared-composition groups through the C API: a mask-built group of a hot and a cool zone shares one composition but
 * not its energy, ungrouped zones are unchanged, and invalid groups are rejected. Exits non-zero on failure.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <stdio.h>
#include <string.h>

#define NUM_SPECIES 8
#define NUM_ZONES 4
#define MAX_CHANNELS 16

static const char* names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
static const double Y0[NUM_SPECIES] = {
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

struct Outputs {
    double Y[NUM_ZONES * NUM_SPECIES];
    double energy[NUM_ZONES], dT[NUM_ZONES], dRho[NUM_ZONES], nu_loss[NUM_ZONES], nu_flux[NUM_ZONES], lost[NUM_ZONES];
    double eps[NUM_ZONES], dlnT[NUM_ZONES], dlnRho[NUM_ZONES], nu_rate[NUM_ZONES], channels[NUM_ZONES * MAX_CHANNELS];
};

static int evolve(void* ctx, const double* T, const double* rho, struct Outputs* out) {
    double Y_in[NUM_ZONES * NUM_SPECIES];
    for (int z = 0; z < NUM_ZONES; ++z) memcpy(Y_in + z * NUM_SPECIES, Y0, sizeof Y0);
    return gf_evolve(MULTI_ZONE, ctx, Y_in, NUM_SPECIES, T, rho, 1e15, 1e-12, out->Y, out->energy, out->dT, out->dRho,
                     out->nu_loss, out->nu_flux, out->lost, out->eps, out->dlnT, out->dlnRho, out->nu_rate, out->channels);
}

int main(void) {
    void* ctx = gf_init(MULTI_ZONE);
    gf_set_num_zones(MULTI_ZONE, ctx, NUM_ZONES);
    gf_register_species(ctx, NUM_SPECIES, names);
    gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", Y0, NUM_SPECIES);
    gf_construct_solver_from_engine(ctx);

    const double T[NUM_ZONES] = {2.5e7, 1.2e7, 1.5e7, 2e7};
    const double rho[NUM_ZONES] = {120.0, 60.0, 150.0, 100.0};
    const double mass[NUM_ZONES] = {1e31, 3e31, 1e31, 1e31};

    struct Outputs separate, grouped;
    int ret = evolve(ctx, T, rho, &separate);
    CHECK(ret == GF_SUCCESS, "ungrouped evolve: %s", gf_get_last_error_message(ctx));

    /* Zones 0 and 1 form a convective zone */
    const int convective[NUM_ZONES] = {1, 1, 0, 0};
    int group[NUM_ZONES];
    CHECK(gf_shared_composition_groups_from_mask(convective, NUM_ZONES, group) == GF_SUCCESS, "groups from mask");
    CHECK(group[0] == 0 && group[1] == 0 && group[2] == -1 && group[3] == -1, "mask gave %d %d %d %d", group[0], group[1], group[2], group[3]);
    CHECK(gf_set_shared_composition_groups(ctx, group, mass) == GF_SUCCESS, "set groups: %s", gf_get_last_error_message(ctx));
    ret = evolve(ctx, T, rho, &grouped);
    CHECK(ret == GF_SUCCESS, "grouped evolve: %s", gf_get_last_error_message(ctx));

    CHECK(memcmp(grouped.Y, grouped.Y + NUM_SPECIES, sizeof(double) * NUM_SPECIES) == 0, "group members differ in composition");
    CHECK(grouped.energy[0] > 10.0 * grouped.energy[1], "hot member energy %.4e vs cool %.4e", grouped.energy[0], grouped.energy[1]);
    CHECK(grouped.eps[0] > grouped.eps[1], "members' eps not evaluated at their own T");
    const double weighted = (mass[0] * grouped.energy[0] + mass[1] * grouped.energy[1]) / (mass[0] + mass[1]);
    printf("group: H-1 %.6f, energies %.4e / %.4e erg/g (mass-weighted %.4e); separate: H-1 %.6f / %.6f\n",
           grouped.Y[0], grouped.energy[0], grouped.energy[1], weighted, separate.Y[0], separate.Y[NUM_SPECIES]);
    CHECK(grouped.Y[0] < separate.Y[NUM_SPECIES] && grouped.Y[0] > separate.Y[0], "group H-1 not between its members' separate burns");
    for (int z = 2; z < NUM_ZONES; ++z) {
        CHECK(grouped.energy[z] == separate.energy[z], "ungrouped zone %d changed (%.10e vs %.10e)", z, grouped.energy[z], separate.energy[z]);
    }

    /* Clearing restores independent zones */
    CHECK(gf_set_shared_composition_groups(ctx, NULL, NULL) == GF_SUCCESS, "clear groups");
    struct Outputs cleared;
    evolve(ctx, T, rho, &cleared);
    CHECK(cleared.energy[0] == separate.energy[0] && cleared.energy[1] == separate.energy[1], "clearing did not restore independent zones");

    /* Invalid groups */
    const int bad[NUM_ZONES] = {0, -2, -1, -1};
    CHECK(gf_set_shared_composition_groups(ctx, bad, mass) == GF_INVALID_STATE, "group id -2 accepted");
    const double no_mass[NUM_ZONES] = {1.0, 0.0, 1.0, 1.0};
    CHECK(gf_set_shared_composition_groups(ctx, group, no_mass) == GF_INVALID_STATE, "a grouped zone without mass accepted");
    gf_free(MULTI_ZONE, ctx);

    void* point = gf_init(SINGLE_ZONE);
    CHECK(gf_set_shared_composition_groups(point, group, mass) == GF_INVALID_TYPE, "groups accepted on a single-zone context");
    gf_free(SINGLE_ZONE, point);

    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All shared-composition checks passed\n");
    return 0;
}
