#include "gridfire/extern/gridfire_extern.h"
#include <stdio.h>

#define NUM_SPECIES 8
#define ZONES 100
#define ZONE_POLICY MULTI_ZONE

// Define a macro to check return codes
#define CHECK_RET_CODE(ret, ctx, msg) \
    if ((ret) != 0 && (ret) != 1) { \
        printf("Error %s: %s (No. %d) [%s]\n", msg, gf_error_code_to_string(ret), ret, gf_get_last_error_message(ctx)); \
        gf_free(ZONE_POLICY, ctx); \
        return 1; \
    }

int main() {
    printf("Testing GridFireEvolve Multi\n");
    printf(" Number of zones: %d\n", ZONES);
    printf(" Number of species: %d\n", NUM_SPECIES);
    printf(" Initializing..");
    void* ctx = gf_init(ZONE_POLICY);
    printf(" Done.\n");
    printf(" Setting number of zones to %d..", ZONES);
    gf_set_num_zones(ZONE_POLICY, ctx, ZONES);
    printf(" Done.\n");

    const char* species_names[NUM_SPECIES];
    species_names[0] = "H-1";
    species_names[1] = "He-3";
    species_names[2] = "He-4";
    species_names[3] = "C-12";
    species_names[4] = "N-14";
    species_names[5] = "O-16";
    species_names[6] = "Ne-20";
    species_names[7] = "Mg-24";
    const double abundance_root[NUM_SPECIES] = {
        0.702616602672027,
        9.74791583949078e-06,
        0.06895512307276903,
        0.00025,
        7.855418029399437e-05,
        0.0006014411598306529,
        8.103062886768109e-05,
        2.151340851063217e-05
    };

    double abundances[ZONES][NUM_SPECIES];
    for (size_t zone = 0; zone < ZONES; zone++) {
        for (size_t i = 0; i < NUM_SPECIES; i++) {
            abundances[zone][i] = abundance_root[i];
        }
    }

    double Temps[ZONES];
    double Rhos[ZONES];

    for (size_t zone = 0; zone < ZONES; zone++) {
        Temps[zone] = 1.0e7 + (double)zone * 1.0e5;  // From 10 million K to 20 million K
        Rhos[zone] = 1.5e2;
        printf("Zone %zu - Temp: %e K, Rho: %e g/cm^3\n", zone, Temps[zone], Rhos[zone]);
    }
    return 0;

    printf(" Registering species...");
    int ret = gf_register_species(ctx, NUM_SPECIES, species_names);
    CHECK_RET_CODE(ret, ctx, "SPECIES");
    printf(" Done.\n");

    printf(" Constructing engine from policy...");
    ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", abundance_root, NUM_SPECIES);
    CHECK_RET_CODE(ret, ctx, "MAIN_SEQUENCE_POLICY");
    printf(" Done.\n");

    printf(" Constructing solver from engine...");
    ret = gf_construct_solver_from_engine(ctx);
    CHECK_RET_CODE(ret, ctx, "CVODE");
    printf(" Done.\n");

    double Y_out[ZONES][NUM_SPECIES];
    double energy_out[ZONES];
    double dEps_dT[ZONES];
    double dEps_dRho[ZONES];
    double specific_neutrino_energy_loss[ZONES];
    double specific_neutrino_flux[ZONES];
    double mass_lost[ZONES];

    printf(" Evolving...\n");
    ret = gf_evolve(
        ZONE_POLICY,
        ctx,
        abundances,
        NUM_SPECIES,
        Temps,    // Temperature in K
        Rhos,    // Density in g/cm^3
        3e17,      // Time step in seconds
        1e-12,
        Y_out,
        energy_out,
        dEps_dT,
        dEps_dRho,
        specific_neutrino_energy_loss,
        specific_neutrino_flux, &mass_lost
    );
    CHECK_RET_CODE(ret, ctx, "EVOLUTION");
    printf(" Done.\n");


    printf("Evolved abundances:\n");
    for (size_t zone = 0; zone < ZONES; zone++) {
        printf("=== Zone %zu ===\n", zone);
        for (size_t i = 0; i < NUM_SPECIES; i++) {
            printf("  Species %s: %e\n", species_names[i], Y_out[zone][i]);
        }
        printf("  Energy output: %e\n", energy_out[zone]);
        printf("  dEps/dT: %e\n", dEps_dT[zone]);
        printf("  dEps/dRho: %e\n", dEps_dRho[zone]);
        printf("  Specific neutrino energy loss: %e\n", specific_neutrino_energy_loss[zone]);
        printf("  Specific neutrino flux: %e\n", specific_neutrino_flux[zone]);
        printf("  Mass lost: %e\n", mass_lost[zone]);
    }

    gf_free(ZONE_POLICY, ctx);

    return 0;
}