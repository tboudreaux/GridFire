#include "gridfire/gridfire.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/extern/gridfire_context.h"
#include "gridfire/extern/gridfire_extern.h"

extern "C" {
    void* gf_init() {
        return new GridFireContext();
    }

    void gf_free(void* ctx) {
        delete static_cast<GridFireContext*>(ctx);
    }

    int gf_register_species(void* ptr, const int num_species, const char** species_names) {
        auto* ctx = static_cast<GridFireContext*>(ptr);
        try {
            std::vector<std::string> names;
            for(int i=0; i<num_species; ++i) {
                names.emplace_back(species_names[i]);
            }
            ctx->init_species_map(names);
            return FDSSE_SUCCESS;
        } catch (const fourdst::composition::exceptions::UnknownSymbolError& e) {
            ctx->last_error = e.what();
            return FDSSE_UNKNOWN_SYMBOL_ERROR;
        } catch (const fourdst::composition::exceptions::SpeciesError& e) {
            ctx->last_error = e.what();
            return FDSSE_SPECIES_ERROR;
        } catch (const std::exception& e) {
            ctx->last_error = e.what();
            return FDSSE_NON_4DSTAR_ERROR;
        } catch (...) {
            ctx->last_error = "Unknown error occurred during species registration.";
            return FDSSE_UNKNOWN_ERROR;
        }
    }

    int gf_construct_engine_from_policy(
        void* ptr,
        const char* policy_name,
        const double *abundances,
        const size_t num_species
    ) {
        auto* ctx = static_cast<GridFireContext*>(ptr);
        try {
            ctx->init_engine_from_policy(std::string(policy_name), abundances, num_species);
            return GF_SUCCESS;
        } catch (const gridfire::exceptions::MissingBaseReactionError& e) {
            ctx->last_error = e.what();
            return GF_MISSING_BASE_REACTION_ERROR;
        } catch (const gridfire::exceptions::MissingSeedSpeciesError& e) {
            ctx->last_error = e.what();
            return GF_MISSING_SEED_SPECIES_ERROR;
        }  catch (const gridfire::exceptions::MissingKeyReactionError& e) {
            ctx->last_error = e.what();
            return GF_MISSING_KEY_REACTION_ERROR;
        } catch (const gridfire::exceptions::PolicyError& e) {
            ctx->last_error = e.what();
            return GF_POLICY_ERROR;
        } catch (std::exception& e) {
            ctx->last_error = e.what();
            return GF_NON_GRIDFIRE_ERROR;
        } catch (...) {
            ctx->last_error = "Unknown error occurred during engine construction.";
            return GF_UNKNOWN_ERROR;
        }
    }

    int gf_construct_solver_from_engine(
        void* ptr,
        const char* solver_name
    ) {
        auto* ctx = static_cast<GridFireContext*>(ptr);
        try {
            ctx->init_solver_from_engine(std::string(solver_name));
            return GF_SUCCESS;
        } catch (std::exception& e) {
            ctx->last_error = e.what();
            return GF_NON_GRIDFIRE_ERROR;
        } catch (...) {
            ctx->last_error = "Unknown error occurred during solver construction.";
            return GF_UNKNOWN_ERROR;
        }
    }

    int gf_evolve(
        void* ptr,
        const double* Y,
        const size_t num_species,
        const double T,
        const double rho,
        const double tMax,
        const double dt0,
        double* Y_out,
        double* energy_out,
        double* dEps_dT,
        double* dEps_dRho, double* mass_lost
    ) {
        auto* ctx = static_cast<GridFireContext*>(ptr);
        try {
            const int result = ctx->evolve(
                Y,
                num_species,
                T,
                rho,
                tMax,
                dt0,
                Y_out,
                *energy_out,
                *dEps_dT,
                *dEps_dRho, *mass_lost
            );
            if (result != 0) {
                return result;
            }
            return GF_SUCCESS;
        } catch (fourdst::composition::exceptions::UnknownSymbolError& e) {
            ctx->last_error = e.what();
            return FDSSE_UNKNOWN_SYMBOL_ERROR;
        } catch (const fourdst::composition::exceptions::SpeciesError& e) {
            ctx->last_error = e.what();
            return FDSSE_SPECIES_ERROR;
        } catch (const fourdst::composition::exceptions::InvalidCompositionError& e) {
            ctx->last_error = e.what();
            return FDSSE_INVALID_COMPOSITION_ERROR;
        } catch (const fourdst::composition::exceptions::CompositionError& e) {
            ctx->last_error = e.what();
            return FDSSE_COMPOSITION_ERROR;
        } catch (const gridfire::exceptions::InvalidQSESolutionError& e) {
            ctx->last_error = e.what();
            return GF_INVALID_QSE_SOLUTION_ERROR;
        } catch (const gridfire::exceptions::FailedToPartitionEngineError& e) {
            ctx->last_error = e.what();
            return GF_FAILED_TO_PARTITION_ENGINE_ERROR;
        } catch (const gridfire::exceptions::NetworkResizedError& e) {
            ctx->last_error = e.what();
            return GF_NETWORK_RESIZED_ERROR;
        } catch (const gridfire::exceptions::UnableToSetNetworkReactionsError& e) {
            ctx->last_error = e.what();
            return GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR;
        } catch (const gridfire::exceptions::BadCollectionError& e) {
            ctx->last_error = e.what();
            return GF_BAD_COLLECTION_ERROR;
        } catch (const gridfire::exceptions::BadRHSEngineError& e) {
            ctx->last_error = e.what();
            return GF_BAD_RHS_ENGINE_ERROR;
        } catch (const gridfire::exceptions::StaleJacobianError& e) {
            ctx->last_error = e.what();
            return GF_STALE_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::UninitializedJacobianError& e) {
            ctx->last_error = e.what();
            return GF_UNINITIALIZED_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::UnknownJacobianError& e) {
            ctx->last_error = e.what();
            return GF_UNKNOWN_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::JacobianError& e) {
            ctx->last_error = e.what();
            return GF_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::EngineError& e) {
            ctx->last_error = e.what();
            return GF_ENGINE_ERROR;
        } catch (const gridfire::exceptions::ReactionParsingError& e) {
            ctx->last_error = e.what();
            return GF_REACTION_PARSING_ERROR;
        } catch (const gridfire::exceptions::ReactionError& e) {
            ctx->last_error = e.what();
            return GF_REACTION_ERROR;
        } catch (const gridfire::exceptions::SingularJacobianError& e) {
            ctx->last_error = e.what();
            return GF_SINGULAR_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::IllConditionedJacobianError& e) {
            ctx->last_error = e.what();
            return GF_ILL_CONDITIONED_JACOBIAN_ERROR;
        } catch (const gridfire::exceptions::CVODESolverFailureError& e) {
            ctx->last_error = e.what();
            return GF_CVODE_SOLVER_FAILURE_ERROR;
        } catch (const gridfire::exceptions::KINSolSolverFailureError& e) {
            ctx->last_error = e.what();
            return GF_KINSOL_SOLVER_FAILURE_ERROR;
        } catch (const gridfire::exceptions::SUNDIALSError& e) {
            ctx->last_error = e.what();
            return GF_SUNDIALS_ERROR;
        } catch (const gridfire::exceptions::SolverError& e) {
            ctx->last_error = e.what();
            return GF_SOLVER_ERROR;
        } catch (const gridfire::exceptions::HashingError& e) {
            ctx->last_error = e.what();
            return GF_HASHING_ERROR;
        } catch (const gridfire::exceptions::UtilityError& e) {
            ctx->last_error = e.what();
            return GF_UTILITY_ERROR;
        } catch (const gridfire::exceptions::DebugException& e) {
            ctx->last_error = e.what();
            return GF_DEBUG_ERROR;
        } catch (const gridfire::exceptions::GridFireError& e) {
            ctx->last_error = e.what();
            return GF_GRIDFIRE_ERROR;
        } catch (std::exception& e) {
            ctx->last_error = e.what();
            return GF_NON_GRIDFIRE_ERROR;
        } catch (...) {
            ctx->last_error = "Unknown error occurred during evolution.";
            return GF_UNKNOWN_ERROR;
        }
    }

    char* gf_get_last_error_message(void* ptr) {
        const auto* ctx = static_cast<GridFireContext*>(ptr);
        return const_cast<char*>(ctx->last_error.c_str());
    }

    char* gf_error_code_to_string(const int error_code) {
        switch (error_code) {
            case GF_SUCCESS:
                return const_cast<char*>("GF_SUCCESS");
            case GF_UNKNOWN_ERROR:
                return const_cast<char*>("GF_UNKNOWN_ERROR");
            case GF_NON_GRIDFIRE_ERROR:
                return const_cast<char*>("GF_NON_GRIDFIRE_ERROR");
            case GF_INVALID_QSE_SOLUTION_ERROR:
                return const_cast<char*>("GF_INVALID_QSE_SOLUTION_ERROR");
            case GF_FAILED_TO_PARTITION_ENGINE_ERROR:
                return const_cast<char*>("GF_FAILED_TO_PARTITION_ENGINE_ERROR");
            case GF_NETWORK_RESIZED_ERROR:
                return const_cast<char*>("GF_NETWORK_RESIZED_ERROR");
            case GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR:
                return const_cast<char*>("GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR");
            case GF_BAD_COLLECTION_ERROR:
                return const_cast<char*>("GF_BAD_COLLECTION_ERROR");
            case GF_BAD_RHS_ENGINE_ERROR:
                return const_cast<char*>("GF_BAD_RHS_ENGINE_ERROR");
            case GF_STALE_JACOBIAN_ERROR:
                return const_cast<char*>("GF_STALE_JACOBIAN_ERROR");
            case GF_UNINITIALIZED_JACOBIAN_ERROR:
                return const_cast<char*>("GF_UNINITIALIZED_JACOBIAN_ERROR");
            case GF_UNKNOWN_JACOBIAN_ERROR:
                return const_cast<char*>("GF_UNKNOWN_JACOBIAN_ERROR");
            case GF_JACOBIAN_ERROR:
                return const_cast<char*>("GF_JACOBIAN_ERROR");
            case GF_ENGINE_ERROR:
                return const_cast<char*>("GF_ENGINE_ERROR");
            case GF_MISSING_BASE_REACTION_ERROR:
                return const_cast<char*>("GF_MISSING_BASE_REACTION_ERROR");
            case GF_MISSING_SEED_SPECIES_ERROR:
                return const_cast<char*>("GF_MISSING_SEED_SPECIES_ERROR");
            case GF_MISSING_KEY_REACTION_ERROR:
                return const_cast<char*>("GF_MISSING_KEY_REACTION_ERROR");
            case GF_POLICY_ERROR:
                return const_cast<char*>("GF_POLICY_ERROR");
            case GF_REACTION_PARSING_ERROR:
                return const_cast<char*>("GF_REACTION_PARSING_ERROR");
            case GF_REACTION_ERROR:
                return const_cast<char*>("GF_REACTION_ERROR");
            case GF_SINGULAR_JACOBIAN_ERROR:
                return const_cast<char*>("GF_SINGULAR_JACOBIAN_ERROR");
            case GF_ILL_CONDITIONED_JACOBIAN_ERROR:
                return const_cast<char*>("GF_ILL_CONDITIONED_JACOBIAN_ERROR");
            case GF_CVODE_SOLVER_FAILURE_ERROR:
                return const_cast<char*>("GF_CVODE_SOLVER_FAILURE_ERROR");
            case GF_KINSOL_SOLVER_FAILURE_ERROR:
                return const_cast<char*>("GF_KINSOL_SOLVER_FAILURE_ERROR");
            case GF_SUNDIALS_ERROR:
                return const_cast<char*>("GF_SUNDIALS_ERROR");
            case GF_SOLVER_ERROR:
                return const_cast<char*>("GF_SOLVER_ERROR");
            case GF_HASHING_ERROR:
                return const_cast<char*>("GF_HASHING_ERROR");
            case GF_UTILITY_ERROR:
                return const_cast<char*>("GF_UTILITY_ERROR");
            case GF_DEBUG_ERROR:
                return const_cast<char*>("GF_DEBUG_ERROR");
            case GF_GRIDFIRE_ERROR:
                return const_cast<char*>("GF_GRIDFIRE_ERROR");
            case FDSSE_NON_4DSTAR_ERROR:
                return const_cast<char*>("FDSSE_NON_4DSTAR_ERROR");
            case FDSSE_UNKNOWN_ERROR:
                return const_cast<char*>("FDSSE_UNKNOWN_ERROR");
            case FDSSE_SUCCESS:
                return const_cast<char*>("FDSSE_SUCCESS");
            case FDSSE_UNKNOWN_SYMBOL_ERROR:
                return const_cast<char*>("FDSSE_UNKNOWN_SYMBOL_ERROR");
            case FDSSE_SPECIES_ERROR:
                return const_cast<char*>("FDSSE_SPECIES_ERROR");
            case FDSSE_INVALID_COMPOSITION_ERROR:
                return const_cast<char*>("FDSSE_INVALID_COMPOSITION_ERROR");
            case FDSSE_COMPOSITION_ERROR:
                return const_cast<char*>("FDSSE_COMPOSITION_ERROR");
            default:
                return const_cast<char*>("GF_UNRECOGNIZED_ERROR_CODE");
        }
    }
}