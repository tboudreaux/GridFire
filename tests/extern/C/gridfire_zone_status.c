/*
 * Failure reporting from gf_evolve: one zone that cannot be integrated (T = 1e10 K, no step-failure restarts) among
 * zones that can, per-zone status and time reached, the documented outputs of a failed zone, the input-error code,
 * and that nothing is printed to stdout or stderr.
 *
 * Usage: gf_extern_c_zone_status <scratch directory>. Exits non-zero if any check fails.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NUM_SPECIES 8
#define NUM_ZONES 3
#define MAX_CHANNELS 16

static const char* names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
static const double Y0[NUM_SPECIES] = {
    0.702616602672027, 9.74791583949078e-06, 0.06895512307276903, 0.00025,
    7.855418029399437e-05, 0.0006014411598306529, 8.103062886768109e-05, 2.151340851063217e-05
};

static int failures = 0;

#define CHECK(cond, ...) do { \
    if (!(cond)) { \
        fprintf(report, "FAIL (line %d): ", __LINE__); \
        fprintf(report, __VA_ARGS__); \
        fprintf(report, "\n"); \
        failures++; \
    } \
} while (0)

static FILE* report; /* checks are reported on a copy of the original stdout, which is redirected while GridFire runs */

static void* make_context(enum GF_TYPE type, size_t zones) {
    void* ctx = gf_init(type);
    if (type == MULTI_ZONE) gf_set_num_zones(type, ctx, zones);
    struct GFConfig config;
    gf_config_defaults(&config);
    config.solver_max_step_failure_restarts = 0; /* fail on the first CVODE failure */
    gf_set_config(ctx, &config);
    gf_register_species(ctx, NUM_SPECIES, names);
    gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", Y0, NUM_SPECIES);
    gf_construct_solver_from_engine(ctx);
    return ctx;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printf("usage: %s <scratch directory>\n", argv[0]);
        return 2;
    }

    /* Send GridFire's stdout and stderr to files so we can check it writes nothing there */
    report = fdopen(dup(STDOUT_FILENO), "w");
    char out_path[1024], err_path[1024];
    snprintf(out_path, sizeof out_path, "%s/zone_status_stdout.txt", argv[1]);
    snprintf(err_path, sizeof err_path, "%s/zone_status_stderr.txt", argv[1]);
    fflush(stdout);
    fflush(stderr);
    if (!freopen(out_path, "w", stdout) || !freopen(err_path, "w", stderr)) {
        fprintf(report, "could not redirect stdout/stderr\n");
        return 2;
    }

    const double tMax = 3e15, dt0 = 1e-12;

    /* --- MULTI_ZONE: the middle zone fails, the others complete --- */
    void* grid = make_context(MULTI_ZONE, NUM_ZONES);
    double Y_in[NUM_ZONES * NUM_SPECIES], Y_out[NUM_ZONES * NUM_SPECIES];
    for (int z = 0; z < NUM_ZONES; ++z) memcpy(Y_in + z * NUM_SPECIES, Y0, sizeof Y0);
    const double T[NUM_ZONES] = {1.5e7, 1e10, 2e7};
    const double rho[NUM_ZONES] = {150.0, 1e4, 100.0};
    double energy[NUM_ZONES], dT[NUM_ZONES], dRho[NUM_ZONES], nu_loss[NUM_ZONES], nu_flux[NUM_ZONES], lost[NUM_ZONES];
    double eps[NUM_ZONES], dlnT[NUM_ZONES], dlnRho[NUM_ZONES], nu_rate[NUM_ZONES], channels[NUM_ZONES * MAX_CHANNELS];
    for (int i = 0; i < NUM_ZONES * NUM_SPECIES; ++i) Y_out[i] = -1.0;

    int ret = gf_evolve(MULTI_ZONE, grid, Y_in, NUM_SPECIES, T, rho, tMax, dt0, Y_out, energy, dT, dRho, nu_loss,
                        nu_flux, lost, eps, dlnT, dlnRho, nu_rate, channels);
    CHECK(ret == GF_ZONE_FAILURE, "multi-zone returned %d (%s), expected GF_ZONE_FAILURE", ret, gf_error_code_to_string(ret));
    CHECK(strstr(gf_get_last_error_message(grid), "zone 1") != NULL, "message does not name zone 1: %s", gf_get_last_error_message(grid));

    int status[NUM_ZONES];
    double reached[NUM_ZONES];
    CHECK(gf_get_zone_status(grid, status, reached) == GF_SUCCESS, "gf_get_zone_status");
    for (int z = 0; z < NUM_ZONES; ++z) {
        fprintf(report, "zone %d: status %d (%s), time reached %.3e s, eps %.4e\n", z, status[z],
                gf_error_code_to_string(status[z]), reached[z], eps[z]);
    }
    CHECK(status[0] == GF_SUCCESS && status[2] == GF_SUCCESS, "the good zones did not succeed");
    CHECK(status[1] == GF_CVODE_SOLVER_FAILURE_ERROR, "zone 1 status %d, expected a CVODE failure", status[1]);
    CHECK(reached[0] == tMax && reached[2] == tMax, "the good zones did not reach tMax");
    CHECK(reached[1] >= 0.0 && reached[1] < tMax, "zone 1 reached %.3e", reached[1]);
    CHECK(memcmp(Y_out + NUM_SPECIES, Y_in + NUM_SPECIES, sizeof(double) * NUM_SPECIES) == 0, "failed zone's Y_out is not its input");
    CHECK(energy[1] == 0.0 && eps[1] == 0.0 && dT[1] == 0.0 && lost[1] == 0.0, "failed zone's outputs are not zero");
    CHECK(energy[0] > 0.0 && eps[0] > 0.0 && Y_out[0] < Y_in[0], "zone 0 was not burnt");
    CHECK(energy[2] > 0.0 && eps[2] > 0.0 && Y_out[2 * NUM_SPECIES] < Y_in[2 * NUM_SPECIES], "zone 2 was not burnt");

    /* --- An input error is reported before integrating, for the whole call --- */
    const double bad_T[NUM_ZONES] = {1.5e7, -5.0, 2e7};
    ret = gf_evolve(MULTI_ZONE, grid, Y_in, NUM_SPECIES, bad_T, rho, tMax, dt0, Y_out, energy, dT, dRho, nu_loss,
                    nu_flux, lost, eps, dlnT, dlnRho, nu_rate, channels);
    CHECK(ret == GF_INVALID_STATE, "a negative temperature returned %d, expected GF_INVALID_STATE", ret);
    gf_get_zone_status(grid, status, NULL);
    CHECK(status[1] == GF_INVALID_STATE, "zone status after an input error is %d", status[1]);
    gf_free(MULTI_ZONE, grid);

    /* --- SINGLE_ZONE: the zone's own code, and the same failed-zone outputs --- */
    void* point = make_context(SINGLE_ZONE, 1);
    double Y_point[NUM_SPECIES], e1, d1, r1, nl1, nf1, m1, eps1, l1, l2, nr1, ch1[MAX_CHANNELS];
    ret = gf_evolve(SINGLE_ZONE, point, Y0, NUM_SPECIES, &T[1], &rho[1], tMax, dt0, Y_point, &e1, &d1, &r1, &nl1, &nf1,
                    &m1, &eps1, &l1, &l2, &nr1, ch1);
    CHECK(ret == GF_CVODE_SOLVER_FAILURE_ERROR, "single zone returned %d (%s)", ret, gf_error_code_to_string(ret));
    int point_status;
    double point_reached;
    gf_get_zone_status(point, &point_status, &point_reached);
    CHECK(point_status == ret, "single-zone status %d differs from the return code %d", point_status, ret);
    CHECK(memcmp(Y_point, Y0, sizeof Y0) == 0 && e1 == 0.0 && eps1 == 0.0, "failed single zone's outputs");

    ret = gf_evolve(SINGLE_ZONE, point, Y0, NUM_SPECIES, &T[0], &rho[0], tMax, dt0, Y_point, &e1, &d1, &r1, &nl1, &nf1,
                    &m1, &eps1, &l1, &l2, &nr1, ch1);
    CHECK(ret == GF_SUCCESS, "the same context did not recover after a failure (%d)", ret);
    gf_get_zone_status(point, &point_status, &point_reached);
    CHECK(point_status == GF_SUCCESS && point_reached == tMax, "status %d, time reached %.3e after success", point_status, point_reached);
    gf_free(SINGLE_ZONE, point);

    /* --- Nothing on stdout or stderr --- */
    fflush(stdout);
    fflush(stderr);
    FILE* out = fopen(out_path, "r");
    FILE* err = fopen(err_path, "r");
    fseek(out, 0, SEEK_END);
    fseek(err, 0, SEEK_END);
    CHECK(ftell(out) == 0, "GridFire wrote %ld bytes to stdout", ftell(out));
    CHECK(ftell(err) == 0, "GridFire wrote %ld bytes to stderr", ftell(err));
    fclose(out);
    fclose(err);

    if (failures) {
        fprintf(report, "%d check(s) failed\n", failures);
        return 1;
    }
    fprintf(report, "All zone status checks passed\n");
    return 0;
}
