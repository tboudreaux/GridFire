/*
 * A downstream C program built only from what an installed GridFire reports (gf-config, pkg-config, CMake): solar
 * core burning for one zone through the C API.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
#include <stdio.h>

#define NUM_SPECIES 8

int main(void) {
    const char* names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
    const double Y[NUM_SPECIES] = {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05};
    const double T = 1.5e7, rho = 150.0;
    double Y_out[NUM_SPECIES], energy, dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost;
    double eps_nuc, dlnT, dlnRho, nu_rate, channels[16];

    void* ctx = gf_init(SINGLE_ZONE);
    int ret = gf_register_species(ctx, NUM_SPECIES, names);
    if (ret == GF_SUCCESS || ret == FDSSE_SUCCESS) ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", Y, NUM_SPECIES);
    if (ret == GF_SUCCESS) ret = gf_construct_solver_from_engine(ctx);
    if (ret == GF_SUCCESS) {
        ret = gf_evolve(SINGLE_ZONE, ctx, Y, NUM_SPECIES, &T, &rho, 1e10, 1e-12, Y_out, &energy, &dEps_dT, &dEps_dRho,
                        &nu_loss, &nu_flux, &mass_lost, &eps_nuc, &dlnT, &dlnRho, &nu_rate, channels);
    }
    if (ret != GF_SUCCESS) {
        printf("FAIL: %s [%s]\n", gf_get_last_error_message(ctx), gf_error_code_to_string(ret));
        return 1;
    }
    gf_free(SINGLE_ZONE, ctx);
    if (!(eps_nuc > 1.0 && eps_nuc < 1e3)) {
        printf("FAIL: implausible solar-core eps_nuc %g erg/g/s\n", eps_nuc);
        return 1;
    }
    printf("C consumer OK: eps_nuc = %.6e erg/g/s\n", eps_nuc);
    return 0;
}
