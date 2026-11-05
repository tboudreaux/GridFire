/*
 * Per-reaction contributions, passive species and mass-fraction units through the C API (DSEP requirements 7, 8):
 *   - gf_evaluate_reactions sums to gf_evaluate's totals, and pp and pep are separate solar neutrino sources;
 *   - a passive species (Fe-56) survives gf_evolve unchanged, strengthens screening, and is rejected when it would
 *     take part in a reaction or is set after the engine exists;
 *   - the stateless mass-fraction helpers round-trip, sum to one, reject bad input, and after burning the drop in
 *     sum_i Y_i m_i is the rest mass released as energy and neutrinos.
 *
 * Exits non-zero if any check fails.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_SPECIES 9
#define MAX_CHANNELS 16
#define I_FE 8

static const char* species_names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24", "Fe-56"};
static const double solar_Y[NUM_SPECIES] = {
    0.702616602672027, 9.74791583949078e-06, 0.06895512307276903, 0.00025,
    7.855418029399437e-05, 0.0006014411598306529, 8.103062886768109e-05, 2.151340851063217e-05, 2.0e-5
};
static const char* passive[] = {"Fe-56"};
static const double T_sun = 1.57e7, rho_sun = 150.0;

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

/* Single-zone context, Fe-56 passive, engine and solver built. */
static void* make_context(const double* initial) {
    void* ctx = gf_init(SINGLE_ZONE);
    CHECK_OK(gf_register_species(ctx, NUM_SPECIES, species_names), ctx, "gf_register_species");
    CHECK_OK(gf_set_passive_species(ctx, 1, passive), ctx, "gf_set_passive_species");
    CHECK_OK(gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", initial, NUM_SPECIES), ctx, "gf_construct_engine_from_policy");
    CHECK_OK(gf_construct_solver_from_engine(ctx), ctx, "gf_construct_solver_from_engine");
    return ctx;
}

static int evolve(void* ctx, const double* Y_in, double* Y_out, double tMax) {
    double energy, dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost, eps, dlnT, dlnRho, nu_rate, channels[MAX_CHANNELS];
    return gf_evolve(SINGLE_ZONE, ctx, Y_in, NUM_SPECIES, &T_sun, &rho_sun, tMax, 1e-12, Y_out, &energy, &dEps_dT,
                     &dEps_dRho, &nu_loss, &nu_flux, &mass_lost, &eps, &dlnT, &dlnRho, &nu_rate, channels);
}

static double eps_nuc(void* ctx, const double* Y, double* nu_loss) {
    double eps = NAN, dT, dRho;
    CHECK_OK(gf_evaluate(SINGLE_ZONE, ctx, Y, NUM_SPECIES, &T_sun, &rho_sun, &eps, &dT, &dRho, nu_loss, NULL, NULL),
             ctx, "gf_evaluate");
    return eps;
}

static void check_reactions(void* ctx) {
    size_t num_reactions = 0;
    CHECK_OK(gf_get_num_reactions(ctx, &num_reactions), ctx, "gf_get_num_reactions");
    double* flow = calloc(num_reactions, sizeof(double));
    double* energy = calloc(num_reactions, sizeof(double));
    double* nu_energy = calloc(num_reactions, sizeof(double));
    double* nu_rate = calloc(num_reactions, sizeof(double));
    CHECK_OK(gf_evaluate_reactions(SINGLE_ZONE, ctx, solar_Y, NUM_SPECIES, &T_sun, &rho_sun, flow, energy, nu_energy, nu_rate),
             ctx, "gf_evaluate_reactions");

    double nu_loss = NAN;
    const double eps = eps_nuc(ctx, solar_Y, &nu_loss);
    double sum_energy = 0.0, sum_nu = 0.0;
    double pp = NAN, pep = NAN;
    for (size_t j = 0; j < num_reactions; ++j) {
        sum_energy += energy[j];
        sum_nu += nu_energy[j];
        const char* id = NULL;
        CHECK_OK(gf_get_reaction_id(ctx, j, &id), ctx, "gf_get_reaction_id");
        if (id && strcmp(id, "p(p,e+)d") == 0) pp = nu_rate[j];
        if (id && strcmp(id, "p(e-p,)d") == 0) pep = nu_rate[j];
    }
    printf("eps_nuc %.6e (sum over reactions %.6e), neutrino loss %.6e (sum %.6e), pep/pp %.4e\n",
           eps, sum_energy, nu_loss, sum_nu, pep / pp);
    CHECK(rel_diff(sum_energy, eps) < 1e-10, "reaction energies sum to %.10e, eps_nuc is %.10e", sum_energy, eps);
    CHECK(rel_diff(sum_nu, nu_loss) < 1e-10, "reaction neutrino losses sum to %.10e, the total is %.10e", sum_nu, nu_loss);
    CHECK(pp > 0.0 && pep > 0.0, "pp (%g) and pep (%g) must both emit neutrinos", pp, pep);
    CHECK(pep / pp > 1e-3 && pep / pp < 1e-2, "pep/pp = %g, expected a few tenths of a percent", pep / pp);

    /* A NULL output is allowed */
    CHECK_OK(gf_evaluate_reactions(SINGLE_ZONE, ctx, solar_Y, NUM_SPECIES, &T_sun, &rho_sun, NULL, NULL, NULL, nu_rate),
             ctx, "gf_evaluate_reactions with NULL outputs");
    free(flow);
    free(energy);
    free(nu_energy);
    free(nu_rate);
}

static void check_passive(void* ctx) {
    double Y_out[NUM_SPECIES];
    CHECK_OK(evolve(ctx, solar_Y, Y_out, 3.15e13), ctx, "gf_evolve");
    CHECK(Y_out[I_FE] == solar_Y[I_FE], "passive Fe-56 changed: %.17g -> %.17g", solar_Y[I_FE], Y_out[I_FE]);
    CHECK(Y_out[0] < solar_Y[0], "hydrogen did not burn");

    /* Iron adds charge to the plasma, so it raises the screened rates */
    double Y_more_iron[NUM_SPECIES];
    memcpy(Y_more_iron, solar_Y, sizeof solar_Y);
    Y_more_iron[I_FE] = 2e-3;
    const double ratio = eps_nuc(ctx, Y_more_iron, NULL) / eps_nuc(ctx, solar_Y, NULL);
    CHECK(ratio > 1.001 && ratio < 1.2, "more passive iron changed eps_nuc by %.6f", ratio);

    /* Too late: the engine exists */
    CHECK(gf_set_passive_species(ctx, 1, passive) == GF_CONFIG_ERROR, "gf_set_passive_species after the engine was accepted");

    /* A species that burns cannot be passive */
    void* bad = gf_init(SINGLE_ZONE);
    const char* burning[] = {"C-12"};
    CHECK_OK(gf_register_species(bad, NUM_SPECIES, species_names), bad, "gf_register_species");
    CHECK_OK(gf_set_passive_species(bad, 1, burning), bad, "gf_set_passive_species");
    const int ret = gf_construct_engine_from_policy(bad, "MAIN_SEQUENCE_POLICY", solar_Y, NUM_SPECIES);
    CHECK(ret != GF_SUCCESS, "C-12 was accepted as a passive species");
    printf("passive C-12 rejected: %s\n", gf_get_last_error_message(bad));
    const char* unknown[] = {"Qq-12"};
    CHECK(gf_set_passive_species(bad, 1, unknown) == GF_CONFIG_ERROR, "an unknown passive species was accepted");
    gf_free(SINGLE_ZONE, bad);
}

#define MAX_NETWORK 256

static void check_mass_fraction_helpers(void* ctx) {
    /* Molar -> mass fractions: each zone sums to one */
    double X[NUM_SPECIES], Y_back[NUM_SPECIES], X_back[NUM_SPECIES];
    CHECK_OK(gf_molar_to_mass_fractions(ctx, solar_Y, NUM_SPECIES, 1, X), ctx, "gf_molar_to_mass_fractions");
    double sum_X = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) sum_X += X[i];
    CHECK(fabs(sum_X - 1.0) < 1e-14, "mass fractions sum to %.16f", sum_X);

    /* Mass fractions -> molar -> mass fractions is exact, and the molar abundances keep their ratios */
    CHECK_OK(gf_mass_fractions_to_molar(ctx, X, NUM_SPECIES, 1, Y_back), ctx, "gf_mass_fractions_to_molar");
    CHECK_OK(gf_molar_to_mass_fractions(ctx, Y_back, NUM_SPECIES, 1, X_back), ctx, "gf_molar_to_mass_fractions");
    for (int i = 0; i < NUM_SPECIES; ++i) {
        CHECK(rel_diff(X_back[i], X[i]) < 1e-14, "%s: X round trip %.17g -> %.17g", species_names[i], X[i], X_back[i]);
        CHECK(rel_diff(Y_back[i] / solar_Y[i], Y_back[0] / solar_Y[0]) < 1e-14, "%s: molar ratios changed", species_names[i]);
    }

    /* Several zones at once, and a zone whose mass fractions do not sum to one is reported by number */
    double X2[2 * NUM_SPECIES], Y2[2 * NUM_SPECIES];
    memcpy(X2, X, sizeof X);
    memcpy(X2 + NUM_SPECIES, X, sizeof X);
    CHECK_OK(gf_mass_fractions_to_molar(ctx, X2, NUM_SPECIES, 2, Y2), ctx, "gf_mass_fractions_to_molar, 2 zones");
    CHECK(memcmp(Y2, Y2 + NUM_SPECIES, sizeof X) == 0 && memcmp(Y2, Y_back, sizeof X) == 0, "zones convert independently");
    X2[NUM_SPECIES] *= 1.01;
    const int ret = gf_mass_fractions_to_molar(ctx, X2, NUM_SPECIES, 2, Y2);
    CHECK(ret == GF_INVALID_STATE, "mass fractions summing to %.4f accepted [%s]", 1.0 + 0.01 * X[0], gf_error_code_to_string(ret));
    CHECK(strstr(gf_get_last_error_message(ctx), "Zone 1") != NULL, "error does not name the zone: %s", gf_get_last_error_message(ctx));
    CHECK(gf_mass_fractions_to_molar(ctx, X, NUM_SPECIES - 1, 1, Y_back) != GF_SUCCESS, "wrong species count accepted");
}

/* Burning lowers sum_i Y_i m_i by the rest mass turned into energy: (drop) * u * N_A * c^2 = energy + neutrino losses */
static void check_mass_defect(void* ctx) {
    CHECK_OK(gf_register_network_species(ctx), ctx, "gf_register_network_species");
    static double Y[MAX_NETWORK], Y_out[MAX_NETWORK], X_out[MAX_NETWORK], mass[MAX_NETWORK];
    static const char* names[MAX_NETWORK];
    size_t n = 0;
    CHECK_OK(gf_get_num_network_species(ctx, &n), ctx, "gf_get_num_network_species");
    CHECK(n <= MAX_NETWORK, "network too large for the test (%zu)", n);
    for (size_t i = 0; i < n; ++i) {
        CHECK_OK(gf_get_network_species(ctx, i, &names[i], NULL, NULL, &mass[i]), ctx, "gf_get_network_species");
        Y[i] = 0.0;
        for (int r = 0; r < NUM_SPECIES; ++r) {
            if (strcmp(names[i], species_names[r]) == 0) Y[i] = solar_Y[r];
        }
    }
    double energy, dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost, eps, dlnT, dlnRho, nu_rate, channels[MAX_CHANNELS];
    CHECK_OK(gf_evolve(SINGLE_ZONE, ctx, Y, n, &T_sun, &rho_sun, 3.15e15, 1e-12, Y_out, &energy, &dEps_dT, &dEps_dRho,
                       &nu_loss, &nu_flux, &mass_lost, &eps, &dlnT, &dlnRho, &nu_rate, channels), ctx, "gf_evolve");
    double M_in = 0.0, M_out = 0.0, sum_X = 0.0;
    for (size_t i = 0; i < n; ++i) {
        M_in += Y[i] * mass[i];
        M_out += Y_out[i] * mass[i];
    }
    CHECK_OK(gf_molar_to_mass_fractions(ctx, Y_out, n, 1, X_out), ctx, "gf_molar_to_mass_fractions");
    for (size_t i = 0; i < n; ++i) sum_X += X_out[i];
    const double u_NA_c2 = 1.66053906660e-24 * 6.02214076e23 * 2.99792458e10 * 2.99792458e10; /* erg per (g/mol)/(mol/g) */
    const double released = energy + nu_loss;
    const double defect = (M_in - M_out) * u_NA_c2;
    printf("rest mass per gram %.12f -> %.12f; defect c^2 %.6e erg/g, energy + neutrinos %.6e erg/g (ratio %.6f); sum X %.15f\n",
           M_in, M_out, defect, released, defect / released, sum_X);
    CHECK(M_out < M_in, "burning did not lower the rest mass");
    CHECK(fabs(defect / released - 1.0) < 1e-3, "mass defect and released energy disagree (ratio %.6f)", defect / released);
    CHECK(fabs(sum_X - 1.0) < 1e-14, "mass fractions after burning sum to %.16f", sum_X);
}

int main(void) {
    void* ctx = make_context(solar_Y);
    check_reactions(ctx);
    check_passive(ctx);
    check_mass_fraction_helpers(ctx);
    check_mass_defect(ctx);
    gf_free(SINGLE_ZONE, ctx);

    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All C per-reaction, passive species and mass fraction checks passed\n");
    return 0;
}
