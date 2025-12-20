#include "gridfire/extern/gridfire_extern.h"
#include <stdio.h>
#define NUM_SPECIES 8

// Define a macro to check return codes
#define CHECK_RET_CODE(ret, ctx, msg) \
    if ((ret) != 0 && (ret) != 1) { \
        printf("Error %s: %s [%s]\n", msg, gf_get_last_error_message(ctx), gf_error_code_to_string(ret)); \
        gf_free(SINGLE_ZONE, ctx); \
        return 1; \
    }

int main() {
    void* ctx = gf_init(SINGLE_ZONE);

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

    ret = gf_construct_solver_from_engine(ctx);
    CHECK_RET_CODE(ret, ctx, "CVODE");

    double Y_out[NUM_SPECIES];
    double energy_out;
    double dEps_dT;
    double dEps_dRho;
    double specific_neutrino_energy_loss;
    double specific_neutrino_flux;
    double mass_lost;

    const double T_in = 1e7; // Temperature in K
    const double rho_in = 1.5e2; // Density in g/cm^3

    ret = gf_evolve(
        SINGLE_ZONE,
        ctx,
        abundances,
        NUM_SPECIES,
        &T_in,    // Temperature in K
        &rho_in,    // Density in g/cm^3
        3e17,      // Time step in seconds
        1e-12,
        Y_out,
        &energy_out,
        &dEps_dT,
        &dEps_dRho,
        &specific_neutrino_energy_loss,
        &specific_neutrino_flux, &mass_lost
    );

    CHECK_RET_CODE(ret, ctx, "EVOLUTION");


    printf("Evolved abundances:\n");
    for (size_t i = 0; i < NUM_SPECIES; i++) {
        printf("Species %s: %e\n", species_names[i], Y_out[i]);
    }
    printf("Energy output: %e\n", energy_out);
    printf("dEps/dT: %e\n", dEps_dT);
    printf("dEps/dRho: %e\n", dEps_dRho);
    printf("Specific neutrino energy loss: %e\n", specific_neutrino_energy_loss);
    printf("Specific neutrino flux: %e\n", specific_neutrino_flux);
    printf("Mass lost: %e\n", mass_lost);

    gf_free(SINGLE_ZONE, ctx);

    return 0;
}