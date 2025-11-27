#ifndef GF_GRIDFIRE_EXTERN_H
#define GF_GRIDFIRE_EXTERN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
    enum FDSSE_ERROR_CODES {
        FDSSE_NON_4DSTAR_ERROR = -102,
        FDSSE_UNKNOWN_ERROR = -101,
        FDSSE_SUCCESS = 1,
        FDSSE_UNKNOWN_SYMBOL_ERROR = 100,
        FDSSE_SPECIES_ERROR = 101,
        FDSSE_INVALID_COMPOSITION_ERROR = 102,
        FDSSE_COMPOSITION_ERROR = 103
    };

    enum GF_ERROR_CODES {
        GF_NON_GRIDFIRE_ERROR = -2,
        GF_UNKNOWN_ERROR = -1,
        GF_SUCCESS = 0,

        GF_INVALID_QSE_SOLUTION_ERROR = 5,
        GF_FAILED_TO_PARTITION_ENGINE_ERROR = 6,
        GF_NETWORK_RESIZED_ERROR = 7,
        GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR = 8,
        GF_BAD_COLLECTION_ERROR = 9,
        GF_BAD_RHS_ENGINE_ERROR = 10,
        GF_STALE_JACOBIAN_ERROR = 11,
        GF_UNINITIALIZED_JACOBIAN_ERROR = 12,
        GF_UNKNOWN_JACOBIAN_ERROR = 13,
        GF_JACOBIAN_ERROR = 14,
        GF_ENGINE_ERROR = 15,

        GF_MISSING_BASE_REACTION_ERROR = 16,
        GF_MISSING_SEED_SPECIES_ERROR = 17,
        GF_MISSING_KEY_REACTION_ERROR = 18,
        GF_POLICY_ERROR = 19,

        GF_REACTION_PARSING_ERROR = 20,
        GF_REACTION_ERROR = 21,

        GF_SINGULAR_JACOBIAN_ERROR = 22,
        GF_ILL_CONDITIONED_JACOBIAN_ERROR = 23,
        GF_CVODE_SOLVER_FAILURE_ERROR = 24,
        GF_KINSOL_SOLVER_FAILURE_ERROR = 25,
        GF_SUNDIALS_ERROR = 26,
        GF_SOLVER_ERROR = 27,

        GF_HASHING_ERROR = 28,
        GF_UTILITY_ERROR = 29,

        GF_DEBUG_ERROR = 30,

        GF_GRIDFIRE_ERROR = 31,
    };

    char* gf_get_last_error_message(void* ptr);

    char* gf_error_code_to_string(int error_code);

    void* gf_init();

    void gf_free(void* ctx);

    int gf_register_species(void* ptr, const int num_species, const char** species_names);

    int gf_construct_engine_from_policy(void* ptr, const char* policy_name, const double *abundances, size_t num_species);

    int gf_construct_solver_from_engine(void* ptr, const char* solver_name);

    int gf_evolve(
        void* ptr,
        const double* Y_in,
        size_t num_species,
        double T,
        double rho,
        double tMax,
        double dt0,
        double* Y_out,
        double* energy_out,
        double* dEps_dT,
        double* dEps_dRho,
        double* specific_neutrino_energy_loss,
        double* specific_neutrino_flux,
        double* mass_lost
    );

#ifdef __cplusplus
}
#endif

#endif
