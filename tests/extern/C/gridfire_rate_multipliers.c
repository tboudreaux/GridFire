/*
 * Per-reaction rate multipliers through the C API (DSEP requirement 5, the Monte Carlo pipeline), and the reverse
 * reaction switch (requirement 4):
 *   - a factor scales exactly that reaction's flow, energy and neutrino emission, and nothing else; eps moves by the
 *     reaction's own share; 0 switches it off; clearing restores the library rates bit for bit;
 *   - it reaches gf_evolve with both Jacobian methods, including on a context that has already integrated (cached
 *     results and the CppAD tape copies are refreshed), and a burn with a factor matches a fresh context with it;
 *   - with reverse reactions on, a factor scales forward and reverse rate together (detailed balance), so the net
 *     flow of a reversible reaction scales by exactly the factor;
 *   - errors: unknown reaction, negative or non-finite factor, no engine yet; provenance lists the factors;
 *   - consecutive gf_evolve calls share no integrator state (requirement 6).
 *
 * Exits non-zero if any check fails.
 */
#include "gridfire/extern/gridfire_extern.h"

#include <math.h>
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
#define CHECK_OK(ret, ctx, what) CHECK((ret) == GF_SUCCESS || (ret) == FDSSE_SUCCESS, "%s: %s [%s]", what, \
    gf_get_last_error_message(ctx), gf_error_code_to_string(ret))

static double rel_diff(double a, double b) {
    const double scale = fmax(fabs(a), fabs(b));
    return scale == 0.0 ? 0.0 : fabs(a - b) / scale;
}

/* H+He network, every network species registered, in network order */
struct Net {
    void* ctx;
    size_t ns, nr;
    double* Y; /* a solar-like composition over the whole network */
};

static size_t species_index(void* ctx, size_t ns, const char* name) {
    for (size_t i = 0; i < ns; ++i) {
        const char* n = NULL;
        gf_get_network_species(ctx, i, &n, NULL, NULL, NULL);
        if (strcmp(n, name) == 0) return i;
    }
    printf("species %s not in the network\n", name);
    exit(2);
}

static size_t reaction_index(void* ctx, size_t nr, const char* id) {
    for (size_t i = 0; i < nr; ++i) {
        const char* r = NULL;
        gf_get_reaction_id(ctx, i, &r);
        if (strcmp(r, id) == 0) return i;
    }
    printf("reaction %s not in the network\n", id);
    exit(2);
}

/* tight = solver tolerances well below the default ones, so two burns that take different steps agree closely */
static struct Net make_net_tol(int jacobian, int reverse, int tight) {
    struct Net net = {0};
    net.ctx = gf_init(SINGLE_ZONE);
    struct GFConfig config;
    gf_config_defaults(&config);
    config.engine_jacobian = jacobian;
    config.engine_use_reverse_reactions = reverse;
    if (tight) {
        config.solver_abs_tol = 1e-14;
        config.solver_rel_tol = 1e-10;
    }
    CHECK_OK(gf_set_config(net.ctx, &config), net.ctx, "gf_set_config");
    const char* seed[] = {"H-1", "He-4"};
    const double seed_Y[] = {0.7, 0.07};
    CHECK_OK(gf_register_species(net.ctx, 2, seed), net.ctx, "gf_register_species");
    CHECK_OK(gf_construct_engine_from_policy(net.ctx, "H_HE_BURNING_POLICY", seed_Y, 2), net.ctx, "construct engine");
    CHECK_OK(gf_register_network_species(net.ctx), net.ctx, "gf_register_network_species");
    CHECK_OK(gf_construct_solver_from_engine(net.ctx), net.ctx, "construct solver");
    gf_get_num_network_species(net.ctx, &net.ns);
    gf_get_num_reactions(net.ctx, &net.nr);
    net.Y = calloc(net.ns, sizeof(double));
    const char* names[] = {"H-1", "He-3", "He-4", "C-12", "N-14", "O-16", "Ne-20", "Mg-24"};
    const double values[] = {0.7026166, 9.7479e-06, 0.0689551, 0.00025, 7.8554e-05, 0.00060144, 8.103e-05, 2.1513e-05};
    for (size_t k = 0; k < 8; ++k) net.Y[species_index(net.ctx, net.ns, names[k])] = values[k];
    return net;
}

static struct Net make_net(int jacobian, int reverse) { return make_net_tol(jacobian, reverse, 0); }

static void free_net(struct Net* net) {
    gf_free(SINGLE_ZONE, net->ctx);
    free(net->Y);
}

struct Rates {
    double eps, dT, dRho;
    double *flow, *energy, *nu_rate;
};

static struct Rates evaluate(const struct Net* net, const double* Y, double T, double rho) {
    struct Rates r = {0};
    r.flow = calloc(net->nr, sizeof(double));
    r.energy = calloc(net->nr, sizeof(double));
    r.nu_rate = calloc(net->nr, sizeof(double));
    CHECK_OK(gf_evaluate(SINGLE_ZONE, net->ctx, Y, net->ns, &T, &rho, &r.eps, &r.dT, &r.dRho, NULL, NULL, NULL),
             net->ctx, "gf_evaluate");
    CHECK_OK(gf_evaluate_reactions(SINGLE_ZONE, net->ctx, Y, net->ns, &T, &rho, r.flow, r.energy, NULL, r.nu_rate),
             net->ctx, "gf_evaluate_reactions");
    return r;
}

static void free_rates(struct Rates* r) {
    free(r->flow);
    free(r->energy);
    free(r->nu_rate);
}

struct Burn {
    int status;
    double* Y_out;
    double energy, eps;
};

static struct Burn burn(const struct Net* net, double tMax) {
    struct Burn b = {0};
    b.Y_out = calloc(net->ns, sizeof(double));
    const double T = 1.57e7, rho = 150.0;
    double dEps_dT, dEps_dRho, nu_loss, nu_flux, mass_lost, dlnT, dlnRho, nu_rate, channels[MAX_CHANNELS];
    b.status = gf_evolve(SINGLE_ZONE, net->ctx, net->Y, net->ns, &T, &rho, tMax, 1e-12, b.Y_out, &b.energy, &dEps_dT,
                         &dEps_dRho, &nu_loss, &nu_flux, &mass_lost, &b.eps, &dlnT, &dlnRho, &nu_rate, channels);
    CHECK_OK(b.status, net->ctx, "gf_evolve");
    return b;
}

static double burn_diff(const struct Net* net, const struct Burn* a, const struct Burn* b) {
    double worst = rel_diff(a->energy, b->energy);
    for (size_t i = 0; i < net->ns; ++i) {
        if (fabs(a->Y_out[i]) > 1e-20 || fabs(b->Y_out[i]) > 1e-20) worst = fmax(worst, rel_diff(a->Y_out[i], b->Y_out[i]));
    }
    return worst;
}

/* A factor scales exactly its reaction; eps moves by that reaction's share; 0 switches it off; clear restores */
static void check_evaluation(void) {
    struct Net net = make_net(GF_JACOBIAN_ANALYTIC, 0);
    const double T = 1.57e7, rho = 150.0;
    const size_t pp = reaction_index(net.ctx, net.nr, "p(p,e+)d");
    const size_t cno = reaction_index(net.ctx, net.nr, "n14(p,g)o15");

    struct Rates base = evaluate(&net, net.Y, T, rho);
    CHECK(base.flow[pp] > 0.0 && base.flow[cno] > 0.0, "pp and n14(p,g) should both run in the solar core");

    CHECK_OK(gf_set_rate_multiplier(net.ctx, "p(p,e+)d", 2.0), net.ctx, "set pp x2");
    CHECK_OK(gf_set_rate_multiplier(net.ctx, "n14(p,g)o15", 0.5), net.ctx, "set n14(p,g) x0.5");
    double factor = 0.0;
    CHECK_OK(gf_get_rate_multiplier(net.ctx, "p(p,e+)d", &factor), net.ctx, "get");
    CHECK(factor == 2.0, "gf_get_rate_multiplier returned %g, expected 2", factor);
    CHECK_OK(gf_get_rate_multiplier(net.ctx, "he4(2a,)c12", &factor), net.ctx, "get unset");
    CHECK(factor == 1.0, "an unset reaction should report 1, got %g", factor);

    struct Rates scaled = evaluate(&net, net.Y, T, rho);
    CHECK(rel_diff(scaled.flow[pp], 2.0 * base.flow[pp]) < 1e-13, "pp flow: %.17g vs 2 x %.17g", scaled.flow[pp], base.flow[pp]);
    CHECK(rel_diff(scaled.nu_rate[pp], 2.0 * base.nu_rate[pp]) < 1e-13, "pp neutrino rate not doubled");
    CHECK(rel_diff(scaled.flow[cno], 0.5 * base.flow[cno]) < 1e-13, "n14(p,g) flow not halved");
    size_t untouched_changed = 0;
    for (size_t r = 0; r < net.nr; ++r) {
        if (r != pp && r != cno && scaled.flow[r] != base.flow[r]) ++untouched_changed;
    }
    CHECK(untouched_changed == 0, "%zu reactions without a multiplier changed", untouched_changed);
    const double expected_eps = base.eps + base.energy[pp] - 0.5 * base.energy[cno];
    CHECK(rel_diff(scaled.eps, expected_eps) < 1e-12, "eps %.17g, expected %.17g", scaled.eps, expected_eps);

    /* dEps/dT at fixed composition: rates are scaled by constants, so it must match a finite difference */
    const double h = 1e-5 * T, Tp = T + h, Tm = T - h;
    double ep, em, d1, d2;
    gf_evaluate(SINGLE_ZONE, net.ctx, net.Y, net.ns, &Tp, &rho, &ep, &d1, &d2, NULL, NULL, NULL);
    gf_evaluate(SINGLE_ZONE, net.ctx, net.Y, net.ns, &Tm, &rho, &em, &d1, &d2, NULL, NULL, NULL);
    CHECK(rel_diff(scaled.dT, (ep - em) / (2.0 * h)) < 1e-5, "dEps/dT %.10g vs finite difference %.10g", scaled.dT,
          (ep - em) / (2.0 * h));

    CHECK_OK(gf_set_rate_multiplier(net.ctx, "p(p,e+)d", 0.0), net.ctx, "switch pp off");
    struct Rates off = evaluate(&net, net.Y, T, rho);
    CHECK(off.flow[pp] == 0.0, "pp flow with factor 0 is %g", off.flow[pp]);

    CHECK_OK(gf_clear_rate_multipliers(net.ctx), net.ctx, "clear");
    struct Rates cleared = evaluate(&net, net.Y, T, rho);
    CHECK(cleared.eps == base.eps && cleared.dT == base.dT && cleared.dRho == base.dRho,
          "after clearing, eps and derivatives should equal the library-rate values exactly");
    CHECK(memcmp(cleared.flow, base.flow, net.nr * sizeof(double)) == 0, "after clearing, every flow should be restored");

    printf("evaluation: pp x2 -> flow x%.15f, n14(p,g) x0.5 -> flow x%.15f, eps %.6e -> %.6e\n",
           scaled.flow[pp] / base.flow[pp], scaled.flow[cno] / base.flow[cno], base.eps, scaled.eps);
    free_rates(&base);
    free_rates(&scaled);
    free_rates(&off);
    free_rates(&cleared);
    free_net(&net);
}

/* gf_evolve sees the factor with both Jacobian methods, also after a burn without it on the same context */
static void check_evolve(void) {
    const double tMax = 3.15e13; /* 1 Myr */
    struct Burn plain[2], fresh[2], reused[2];
    for (int method = 0; method < 2; ++method) {
        const int jacobian = method ? GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION : GF_JACOBIAN_ANALYTIC;

        /* the same context: burn with library rates, then set the factor and burn again */
        struct Net a = make_net_tol(jacobian, 0, 1);
        plain[method] = burn(&a, tMax);
        CHECK_OK(gf_set_rate_multiplier(a.ctx, "p(p,e+)d", 1.5), a.ctx, "set on a used context");
        reused[method] = burn(&a, tMax);

        /* a context that has the factor from the start */
        struct Net b = make_net_tol(jacobian, 0, 1);
        CHECK_OK(gf_set_rate_multiplier(b.ctx, "p(p,e+)d", 1.5), b.ctx, "set on a fresh context");
        fresh[method] = burn(&b, tMax);

        CHECK(burn_diff(&a, &reused[method], &fresh[method]) < 1e-12,
              "%s: setting the factor on a used context differs from a fresh one by %.3e",
              method ? "CppAD" : "analytic", burn_diff(&a, &reused[method], &fresh[method]));
        CHECK(fresh[method].energy > 1.1 * plain[method].energy,
              "%s: pp x1.5 should release clearly more energy (%.6e vs %.6e)", method ? "CppAD" : "analytic",
              fresh[method].energy, plain[method].energy);
        free_net(&a);
        free_net(&b);
    }
    /* Both Jacobian methods integrate the same scaled network: with tight tolerances they differ only as much as they
       do without a multiplier (different Jacobians take different steps). A tape that missed the factor would put
       CppAD's burn on the unscaled network, ~18% away in energy. */
    struct Net ref = make_net(GF_JACOBIAN_ANALYTIC, 0);
    const double between = burn_diff(&ref, &fresh[0], &fresh[1]);
    const double baseline = burn_diff(&ref, &plain[0], &plain[1]);
    CHECK(between < fmax(10.0 * baseline, 1e-7), "analytic and CppAD disagree by %.3e with a multiplier (%.3e without)",
          between, baseline);
    CHECK(rel_diff(fresh[1].energy, fresh[0].energy) < 1e-7, "CppAD energy %.10e vs analytic %.10e", fresh[1].energy,
          fresh[0].energy);
    printf("evolve (1 Myr, pp x1.5): energy %.6e -> %.6e erg/g; analytic vs CppAD %.2e (%.2e without); "
           "used vs fresh context %.2e / %.2e\n", plain[0].energy, fresh[0].energy, between, baseline,
           burn_diff(&ref, &reused[0], &fresh[0]), burn_diff(&ref, &reused[1], &fresh[1]));
    free_net(&ref);
    for (int m = 0; m < 2; ++m) {
        free(plain[m].Y_out);
        free(fresh[m].Y_out);
        free(reused[m].Y_out);
    }
}

/* With reverse reactions on, a factor scales forward and reverse together: the net flow scales by the factor */
static void check_reverse(void) {
    struct Net off = make_net(GF_JACOBIAN_ANALYTIC, 0);
    struct Net on = make_net(GF_JACOBIAN_ANALYTIC, 1);
    /* Hot, with every species present, so reverse rates are a sizeable part of some net flows */
    const double T = 3.0e9, rho = 1.0e6;
    double* flat = calloc(on.ns, sizeof(double));
    for (size_t i = 0; i < on.ns; ++i) flat[i] = 1.0e-3;
    struct Rates r_off = evaluate(&off, flat, T, rho);
    struct Rates r_on = evaluate(&on, flat, T, rho);

    /* the reversible reaction whose net flow the reverse rate changes most */
    size_t best = (size_t)-1;
    double best_change = 0.0;
    for (size_t r = 0; r < on.nr; ++r) {
        if (r_on.flow[r] == 0.0 || r_off.flow[r] == 0.0) continue;
        const double change = rel_diff(r_on.flow[r], r_off.flow[r]);
        if (change > best_change && change < 0.9) {
            best_change = change;
            best = r;
        }
    }
    CHECK(best != (size_t)-1 && best_change > 1e-3, "no reaction's net flow changed with reverse reactions on");
    if (best != (size_t)-1) {
        const char* id = NULL;
        gf_get_reaction_id(on.ctx, best, &id);
        CHECK_OK(gf_set_rate_multiplier(on.ctx, id, 3.0), on.ctx, "set on a reversible reaction");
        struct Rates scaled = evaluate(&on, flat, T, rho);
        CHECK(rel_diff(scaled.flow[best], 3.0 * r_on.flow[best]) < 1e-12,
              "%s: net flow %.17g, expected 3 x %.17g (forward and reverse scaled together)", id, scaled.flow[best],
              r_on.flow[best]);
        printf("reverse on: %s net flow changes by %.3f%% at T = %.1e K; x3 multiplier -> net flow x%.15f\n", id,
               100.0 * best_change, T, scaled.flow[best] / r_on.flow[best]);
        free_rates(&scaled);
    }

    /* a solar-core burn with reverse reactions on runs and matches the default to well within tolerances */
    struct Burn b_off = burn(&off, 3.15e13);
    struct Burn b_on = burn(&on, 3.15e13);
    CHECK(rel_diff(b_on.energy, b_off.energy) < 1e-6, "solar-core energy with reverse on %.10e vs off %.10e",
          b_on.energy, b_off.energy);
    printf("reverse on: solar-core burn energy differs by %.2e from reverse off\n", rel_diff(b_on.energy, b_off.energy));
    free(b_off.Y_out);
    free(b_on.Y_out);
    free(flat);
    free_rates(&r_off);
    free_rates(&r_on);
    free_net(&off);
    free_net(&on);
}

/* Consecutive gf_evolve calls share no integrator state (DSEP requirement 6): burning the previous call's output on
   a used context (the solver reuses its memory) matches the same burn on a fresh context bit for bit */
static void check_no_carried_state(void) {
    struct Net used = make_net(GF_JACOBIAN_ANALYTIC, 0);
    struct Burn first = burn(&used, 3.15e13);
    memcpy(used.Y, first.Y_out, used.ns * sizeof(double));
    struct Burn second = burn(&used, 3.15e13);

    struct Net fresh = make_net(GF_JACOBIAN_ANALYTIC, 0);
    memcpy(fresh.Y, first.Y_out, fresh.ns * sizeof(double));
    struct Burn alone = burn(&fresh, 3.15e13);

    CHECK(memcmp(second.Y_out, alone.Y_out, used.ns * sizeof(double)) == 0 && second.energy == alone.energy,
          "a burn on a used context differs from the same burn on a fresh one (max rel %.3e)",
          burn_diff(&used, &second, &alone));
    printf("no carried state: chained burn on a used context == fresh context, bit for bit: %s\n",
           burn_diff(&used, &second, &alone) == 0.0 ? "yes" : "no");
    free(first.Y_out);
    free(second.Y_out);
    free(alone.Y_out);
    free_net(&used);
    free_net(&fresh);
}

static void check_errors_and_provenance(void) {
    void* bare = gf_init(SINGLE_ZONE);
    int ret = gf_set_rate_multiplier(bare, "p(p,e+)d", 2.0);
    CHECK(ret == GF_CONFIG_ERROR, "setting a factor before the engine exists returned %s", gf_error_code_to_string(ret));
    gf_free(SINGLE_ZONE, bare);

    struct Net net = make_net(GF_JACOBIAN_ANALYTIC, 0);
    ret = gf_set_rate_multiplier(net.ctx, "x(y,z)w", 2.0);
    CHECK(ret == GF_REACTION_ERROR, "an unknown reaction returned %s", gf_error_code_to_string(ret));
    ret = gf_set_rate_multiplier(net.ctx, "p(p,e+)d", -1.0);
    CHECK(ret == GF_CONFIG_ERROR, "a negative factor returned %s", gf_error_code_to_string(ret));
    ret = gf_set_rate_multiplier(net.ctx, "p(p,e+)d", NAN);
    CHECK(ret == GF_CONFIG_ERROR, "a NaN factor returned %s", gf_error_code_to_string(ret));
    double factor = 0.0;
    ret = gf_get_rate_multiplier(net.ctx, "x(y,z)w", &factor);
    CHECK(ret == GF_REACTION_ERROR, "get for an unknown reaction returned %s", gf_error_code_to_string(ret));

    const char* text = NULL;
    CHECK_OK(gf_get_provenance(net.ctx, &text), net.ctx, "provenance");
    CHECK(text && strstr(text, "rate_multipliers: none"), "provenance without factors:\n%s", text ? text : "(null)");
    CHECK(text && strstr(text, "reverse_reactions: off"), "provenance should report reverse reactions off");
    gf_set_rate_multiplier(net.ctx, "p(p,e+)d", 1.1);
    gf_set_rate_multiplier(net.ctx, "c12(a,g)o16", 0.8);
    CHECK_OK(gf_get_provenance(net.ctx, &text), net.ctx, "provenance");
    CHECK(text && strstr(text, "p(p,e+)d=1.1") && strstr(text, "c12(a,g)o16=0.8"), "provenance with factors:\n%s",
          text ? text : "(null)");
    gf_set_rate_multiplier(net.ctx, "p(p,e+)d", 1.0); /* 1 removes the factor */
    CHECK_OK(gf_get_provenance(net.ctx, &text), net.ctx, "provenance");
    CHECK(text && !strstr(text, "p(p,e+)d=") && strstr(text, "c12(a,g)o16=0.8"), "factor 1 should drop out:\n%s",
          text ? text : "(null)");
    free_net(&net);
}

int main(void) {
    setvbuf(stdout, NULL, _IONBF, 0);
    check_evaluation();
    check_evolve();
    check_reverse();
    check_no_carried_state();
    check_errors_and_provenance();
    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("All rate multiplier checks passed\n");
    return 0;
}
