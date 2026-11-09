/*
 * Thread safety of the C API: separate contexts, built serially, used concurrently from the caller's own threads.
 *
 * Every case runs the same zones twice on fresh contexts, once serially and once concurrently, and requires
 * bit-identical results (each zone's history is identical, so any difference is interference between threads):
 *   1. contexts evolved from an OpenMP parallel loop (analytic Jacobian, the default);
 *   2. the same with the CppAD Jacobian (CppAD's per-thread memory and tapes in use);
 *   3. gf_evaluate and gf_evaluate_reactions from an OpenMP parallel loop;
 *   4. contexts evolved from plain POSIX threads (not OpenMP), with the CppAD Jacobian.
 *
 * Build it with any OpenMP runtime: with GCC's libgomp the caller's threads are invisible to GridFire's own runtime,
 * which is the situation of a gfortran -fopenmp host code.
 *
 * Exits non-zero if any check fails.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_SPECIES 8
#define NUM_ZONES 16
#define MAX_CHANNELS 16

static const char* species_names[NUM_SPECIES] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
static const double Y0[NUM_SPECIES] = {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05};

static int failures = 0;

#define CHECK(cond, ...) do { \
    if (!(cond)) { \
        _Pragma("omp critical") { printf("FAIL (line %d): ", __LINE__); printf(__VA_ARGS__); printf("\n"); failures++; } \
    } \
} while (0)

struct ZoneResult {
    int status;
    double Y_out[NUM_SPECIES];
    double eps_nuc;
    double dlnT;
    double nu_loss;
};

static double zone_T(int zone) { return 1.2e7 + 1.0e6 * zone; }

static void* build_context(int jacobian) {
    void* ctx = gf_init(SINGLE_ZONE);
    struct GFConfig config;
    gf_config_defaults(&config);
    config.engine_jacobian = jacobian;
    int ret = gf_set_config(ctx, &config);
    if (ret == GF_SUCCESS) ret = gf_register_species(ctx, NUM_SPECIES, species_names);
    if (ret == GF_SUCCESS || ret == FDSSE_SUCCESS) ret = gf_construct_engine_from_policy(ctx, "MAIN_SEQUENCE_POLICY", Y0, NUM_SPECIES);
    if (ret == GF_SUCCESS) ret = gf_construct_solver_from_engine(ctx);
    CHECK(ret == GF_SUCCESS, "building a context: %s", gf_get_last_error_message(ctx));
    return ctx;
}

static void evolve_zone(void* ctx, int zone, struct ZoneResult* r) {
    const double T = zone_T(zone), rho = 100.0 + 5.0 * zone;
    double energy, dEps_dT, dEps_dRho, nu_flux, mass_lost, dlnRho, nu_rate, channels[MAX_CHANNELS];
    r->status = gf_evolve(SINGLE_ZONE, ctx, Y0, NUM_SPECIES, &T, &rho, 3.15e13, 1e-12, r->Y_out, &energy, &dEps_dT,
                          &dEps_dRho, &r->nu_loss, &nu_flux, &mass_lost, &r->eps_nuc, &r->dlnT, &dlnRho, &nu_rate, channels);
}

static void compare(const char* label, const struct ZoneResult* serial, const struct ZoneResult* parallel) {
    int differing = 0;
    for (int z = 0; z < NUM_ZONES; ++z) {
        CHECK(serial[z].status == GF_SUCCESS && parallel[z].status == GF_SUCCESS, "%s: zone %d failed (%d, %d)", label, z,
              serial[z].status, parallel[z].status);
        if (memcmp(&serial[z], &parallel[z], sizeof(struct ZoneResult)) != 0) ++differing;
    }
    CHECK(differing == 0, "%s: %d of %d zones differ between serial and concurrent runs", label, differing, NUM_ZONES);
    printf("%-55s %s (%d threads)\n", label, differing ? "DIFFERS" : "identical", omp_get_max_threads());
}

static void prebuilt_contexts(const char* label, int jacobian) {
    struct ZoneResult serial[NUM_ZONES], parallel[NUM_ZONES];
    memset(serial, 0, sizeof serial);
    memset(parallel, 0, sizeof parallel);
    void* a[NUM_ZONES];
    void* b[NUM_ZONES];
    for (int z = 0; z < NUM_ZONES; ++z) {
        a[z] = build_context(jacobian);
        b[z] = build_context(jacobian);
    }
    for (int z = 0; z < NUM_ZONES; ++z) evolve_zone(a[z], z, &serial[z]);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int z = 0; z < NUM_ZONES; ++z) evolve_zone(b[z], z, &parallel[z]);
    compare(label, serial, parallel);
    for (int z = 0; z < NUM_ZONES; ++z) {
        gf_free(SINGLE_ZONE, a[z]);
        gf_free(SINGLE_ZONE, b[z]);
    }
}

struct PthreadJob {
    void* ctx;
    int zone;
    struct ZoneResult* result;
};

static void* pthread_evolve(void* arg) {
    struct PthreadJob* job = arg;
    evolve_zone(job->ctx, job->zone, job->result);
    return NULL;
}

/* Plain POSIX threads, created per batch so that thread numbers are released and reused */
static void posix_threads(void) {
    struct ZoneResult serial[NUM_ZONES], parallel[NUM_ZONES];
    memset(serial, 0, sizeof serial);
    memset(parallel, 0, sizeof parallel);
    void* a[NUM_ZONES];
    void* b[NUM_ZONES];
    for (int z = 0; z < NUM_ZONES; ++z) {
        a[z] = build_context(GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION);
        b[z] = build_context(GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION);
    }
    for (int z = 0; z < NUM_ZONES; ++z) evolve_zone(a[z], z, &serial[z]);
    const int batch = 8;
    for (int first = 0; first < NUM_ZONES; first += batch) {
        pthread_t threads[8];
        struct PthreadJob jobs[8];
        for (int k = 0; k < batch && first + k < NUM_ZONES; ++k) {
            jobs[k] = (struct PthreadJob){b[first + k], first + k, &parallel[first + k]};
            pthread_create(&threads[k], NULL, pthread_evolve, &jobs[k]);
        }
        for (int k = 0; k < batch && first + k < NUM_ZONES; ++k) pthread_join(threads[k], NULL);
    }
    compare("contexts evolved from POSIX threads (CppAD Jacobian)", serial, parallel);
    for (int z = 0; z < NUM_ZONES; ++z) {
        gf_free(SINGLE_ZONE, a[z]);
        gf_free(SINGLE_ZONE, b[z]);
    }
}

static void concurrent_evaluation(void) {
    void* ctx_serial = build_context(GF_JACOBIAN_ANALYTIC);
    size_t num_reactions = 0;
    gf_get_num_reactions(ctx_serial, &num_reactions);
    double* serial = calloc(NUM_ZONES * (3 + num_reactions), sizeof(double));
    double* parallel = calloc(NUM_ZONES * (3 + num_reactions), sizeof(double));
    void* ctxs[NUM_ZONES];
    for (int z = 0; z < NUM_ZONES; ++z) ctxs[z] = build_context(GF_JACOBIAN_ANALYTIC);

    for (int pass = 0; pass < 2; ++pass) {
        double* out = pass ? parallel : serial;
        #pragma omp parallel for schedule(dynamic, 1) if (pass)
        for (int z = 0; z < NUM_ZONES; ++z) {
            void* ctx = pass ? ctxs[z] : ctx_serial;
            const double T = zone_T(z), rho = 100.0 + 5.0 * z;
            double* row = out + z * (3 + num_reactions);
            int ret = gf_evaluate(SINGLE_ZONE, ctx, Y0, NUM_SPECIES, &T, &rho, &row[0], &row[1], &row[2], NULL, NULL, NULL);
            if (ret == GF_SUCCESS) {
                ret = gf_evaluate_reactions(SINGLE_ZONE, ctx, Y0, NUM_SPECIES, &T, &rho, NULL, row + 3, NULL, NULL);
            }
            CHECK(ret == GF_SUCCESS, "evaluation in zone %d failed", z);
        }
    }
    const int same = memcmp(serial, parallel, NUM_ZONES * (3 + num_reactions) * sizeof(double)) == 0;
    CHECK(same, "concurrent gf_evaluate / gf_evaluate_reactions differ from serial");
    printf("%-55s %s (%d threads)\n", "gf_evaluate and gf_evaluate_reactions concurrently", same ? "identical" : "DIFFERS",
           omp_get_max_threads());
    for (int z = 0; z < NUM_ZONES; ++z) gf_free(SINGLE_ZONE, ctxs[z]);
    gf_free(SINGLE_ZONE, ctx_serial);
    free(serial);
    free(parallel);
}

int main(void) {
    setvbuf(stdout, NULL, _IONBF, 0);
    if (omp_get_max_threads() < 2) omp_set_num_threads(4);
    prebuilt_contexts("prebuilt contexts evolved concurrently (analytic)", GF_JACOBIAN_ANALYTIC);
    prebuilt_contexts("prebuilt contexts evolved concurrently (CppAD Jacobian)", GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION);
    concurrent_evaluation();
    posix_threads();
    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All thread-safety checks passed\n");
    return 0;
}
