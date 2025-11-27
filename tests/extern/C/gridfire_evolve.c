#include "gridfire/extern/gridfire_extern.h"
#include <stdio.h>
#define NUM_SPECIES 8

// Define a macro to check return codes
#define CHECK_RET_CODE(ret, ctx, msg) \
    if (ret != 0 && ret != 1) { \
        printf("Error %s: %s\n", msg, gf_get_last_error_message(ctx)); \
        gf_free(ctx); \
        return 1; \
    }

int main() {
    void* ctx = gf_init();

    const char* species_names[NUM_SPECIES];
    species_names[0] = "H-1";
    species_names[1] = "He-3";
    species_names[2] = "He-4";
    species_names[3] = "C-12";
    species_names[4] = "N-14";
    species_names[5] = "O-16";
    species_names[6] = "Ne-20";
    species_names[7] = "Mg-24";
    const double abundances[NUM_SPECIES] = {
        0.702616602672027,
        9.74791583949078e-06,
        0.06895512307276903,
        0.00025,
        7.855418029399437e-05,
        0.0006014411598306529,
        8.103062886768109e-05,
        2.151340851063217e-05
    };

    int ret = gf_register_species(ctx, NUM_SPECIES, species_names);
    CHECK_RET_CODE(ret, ctx, "SPECIES");

    ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", abundances, NUM_SPECIES);
    CHECK_RET_CODE(ret, ctx, "MAIN_SEQUENCE_POLICY");

    ret = gf_construct_solver_from_engine(ctx, "CVODE");
    CHECK_RET_CODE(ret, ctx, "CVODE");

    double Y_out[NUM_SPECIES];
    double energy_out;
    double dEps_dT;
    double dEps_dRho;
    double mass_lost;

    ret = gf_evolve(
        ctx,
        abundances,
        NUM_SPECIES,
        1.5e7,    // Temperature in K
        1.5e2,    // Density in g/cm^3
        3e17,      // Time step in seconds
        1e-12, // Initial time step in seconds
        Y_out,
        &energy_out,
        &dEps_dT,
        &dEps_dRho, &mass_lost
    );

    CHECK_RET_CODE(ret, ctx, "EVOLUTION");


    printf("Evolved abundances:\n");
    for (size_t i = 0; i < NUM_SPECIES; i++) {
        printf("Species %s: %e\n", species_names[i], Y_out[i]);
    }
    printf("Energy output: %e\n", energy_out);
    printf("dEps/dT: %e\n", dEps_dT);
    printf("dEps/dRho: %e\n", dEps_dRho);
    printf("Mass lost: %e\n", mass_lost);

    gf_free(ctx);

    return 0;
}