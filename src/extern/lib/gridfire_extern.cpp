#include "gridfire/gridfire.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "gridfire/extern/gridfire_context.h"
#include "gridfire/extern/gridfire_extern.h"

namespace {

    template<typename T>
    concept ErrorTrackable = requires(T a) {
        { a.last_error } -> std::convertible_to<std::string>;
    };

    template <ErrorTrackable Context, typename Func>
    int execute_guarded(Context* ctx, Func&& action) {
        try {
            const int result = action();

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

}

extern "C" {


    void* gf_init(const enum GF_TYPE type) {
        if (type == MULTI_ZONE) {
            return new GFGridContext();
        }
        if (type == SINGLE_ZONE) {
            return new GFPointContext();
        }
        return nullptr;
    }

    int gf_free(const enum GF_TYPE type, void *ctx) {
        if (!ctx) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }
        if (type == MULTI_ZONE) {
            delete static_cast<GFGridContext*>(ctx);
            return GF_SUCCESS;
        }
        if (type == SINGLE_ZONE) {
            delete static_cast<GFPointContext*>(ctx);
            return GF_SUCCESS;
        }
        return GF_UNKNOWN_FREE_TYPE;
    }

    int gf_set_num_zones(const enum GF_TYPE type, void* ptr, const size_t num_zones) {
        if (type != MULTI_ZONE) {
            return GF_INVALID_TYPE;
        }

        if (!ptr) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }

        auto* ctx = static_cast<GFGridContext*>(ptr);
        return execute_guarded(ctx, [&]() {
            ctx->set_zones(num_zones);
            return GF_SUCCESS;
        });
    }

    int gf_register_species(void* ptr, const int num_species, const char** species_names) {
        if (num_species < 0) return GF_INVALID_NUM_SPECIES;

        if (num_species == 0) return GF_SUCCESS;

        if (!ptr || !species_names) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }

        for (int i=0; i < num_species; ++i) {
            if (!species_names[i]) {
                return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
            }
        }

        auto* ctx = static_cast<GFContext*>(ptr);
        return execute_guarded(ctx, [&]() {
           std::vector<std::string> names;
            for (int i=0; i<num_species; ++i) {
                names.emplace_back(species_names[i]);
            }
            ctx->init_species_map(names);
            return FDSSE_SUCCESS;
        });
    }

    int gf_construct_engine_from_policy(
        void* ptr,
        const char* policy_name,
        const double *abundances,
        const size_t num_species
    ) {
        auto* ctx = static_cast<GFContext*>(ptr);
        return execute_guarded(ctx, [&]() {
            ctx->init_engine_from_policy(std::string(policy_name), abundances, num_species);
            return GF_SUCCESS;
        });
    }

    int gf_construct_solver_from_engine(
        void* ptr
    ) {
        auto* ctx = static_cast<GFContext*>(ptr);
        return execute_guarded(ctx, [&]() {
            ctx->init_solver_from_engine();
            return GF_SUCCESS;
        });
    }

    int gf_evolve(
        const enum GF_TYPE type,
        void* ptr,
        const void* Y_in,
        const size_t num_species,
        const void* T,
        const void* rho,
        const double tMax,
        const double dt0,
        void* Y_out,
        void* energy_out,
        void* dEps_dT,
        void* dEps_dRho,
        void* specific_neutrino_energy_loss,
        void* specific_neutrino_flux,
        void* mass_lost,
        void* eps_nuc,
        void* dlnEps_dlnT,
        void* dlnEps_dlnRho,
        void* neutrino_energy_loss_rate,
        void* eps_by_channel
    ) {
        if (!ptr || !Y_in || !T || !rho) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }

        if (!Y_out || !energy_out || !dEps_dT || !dEps_dRho || !specific_neutrino_energy_loss || !specific_neutrino_flux || !mass_lost ||
            !eps_nuc || !dlnEps_dlnT || !dlnEps_dlnRho || !neutrino_energy_loss_rate || !eps_by_channel) {
            return GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR;
        }

        if (tMax <= 0 || dt0 <= 0) {
            return GF_INVALID_TIMESTEPS;
        }

        if (num_species <= 0) {
            return GF_INVALID_NUM_SPECIES;
        }

        const GFInstantaneousOutputs instantaneous{
            .eps_nuc = static_cast<double*>(eps_nuc),
            .dlnEps_dlnT = static_cast<double*>(dlnEps_dlnT),
            .dlnEps_dlnRho = static_cast<double*>(dlnEps_dlnRho),
            .neutrino_energy_loss_rate = static_cast<double*>(neutrino_energy_loss_rate),
            .eps_by_channel = static_cast<double*>(eps_by_channel)
        };

        auto* Y_out_local = static_cast<double*>(Y_out);
        auto* energy_out_local = static_cast<double*>(energy_out);
        auto* dEps_dT_local = static_cast<double*>(dEps_dT);
        auto* dEps_dRho_local = static_cast<double*>(dEps_dRho);
        auto* specific_neutrino_energy_loss_local = static_cast<double*>(specific_neutrino_energy_loss);
        auto* specific_neutrino_flux_local = static_cast<double*>(specific_neutrino_flux);
        auto* mass_lost_local = static_cast<double*>(mass_lost);

        switch (type) {
            case SINGLE_ZONE : {
                auto* ctx = static_cast<GFPointContext*>(ptr);
                return execute_guarded(ctx, [&]() {
                    return ctx->evolve(
                        static_cast<const double*>(Y_in),
                        num_species,
                        *static_cast<const double*>(T),
                        *static_cast<const double*>(rho),
                        tMax,
                        dt0,
                        Y_out_local,
                        *energy_out_local,
                        *dEps_dT_local,
                        *dEps_dRho_local,
                        *specific_neutrino_energy_loss_local,
                        *specific_neutrino_flux_local,
                        *mass_lost_local,
                        instantaneous
                    );
                });
            }
            case MULTI_ZONE : {
                auto* ctx = static_cast<GFGridContext*>(ptr);
                return execute_guarded(ctx, [&]() {
                    return ctx->evolve(
                        static_cast<const double*>(Y_in),
                        num_species,
                        static_cast<const double*>(T),
                        static_cast<const double*>(rho),
                        tMax,
                        dt0,
                        Y_out_local,
                        energy_out_local,
                        dEps_dT_local,
                        dEps_dRho_local,
                        specific_neutrino_energy_loss_local,
                        specific_neutrino_flux_local,
                        mass_lost_local,
                        instantaneous
                    );
                });
            }
            default :
                return GF_UNKNOWN_ERROR;
        }
    }

    int gf_get_num_energy_channels(void* ptr, size_t* num_channels) {
        if (!ptr) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }
        if (!num_channels) {
            return GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR;
        }
        auto* ctx = static_cast<GFContext*>(ptr);
        return execute_guarded(ctx, [&]() {
            *num_channels = ctx->get_energy_channel_names().size();
            return GF_SUCCESS;
        });
    }

    int gf_get_energy_channel_name(void* ptr, const size_t index, const char** name) {
        if (!ptr) {
            return GF_UNINITIALIZED_INPUT_MEMORY_ERROR;
        }
        if (!name) {
            return GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR;
        }
        auto* ctx = static_cast<GFContext*>(ptr);
        return execute_guarded(ctx, [&]() {
            *name = ctx->get_energy_channel_names().at(index).c_str();
            return GF_SUCCESS;
        });
    }

    char* gf_get_last_error_message(void* ptr) {
        if (!ptr) {
            return const_cast<char*>("GF_UNINITIALIZED_INPUT_MEMORY_ERROR");
        }
        const auto* ctx = static_cast<GFContext*>(ptr);
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
            case GF_UNINITIALIZED_INPUT_MEMORY_ERROR:
                return const_cast<char*>("GF_UNINITIALIZED_INPUT_MEMORY_ERROR");
            case GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR:
                return const_cast<char*>("GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR");
            case GF_INVALID_NUM_SPECIES:
                return const_cast<char*>("GF_INVALID_NUM_SPECIES");
            case GF_INVALID_TIMESTEPS:
                return const_cast<char*>("GF_INVALID_TIMESTEPS");
            case GF_UNKNOWN_FREE_TYPE:
                return const_cast<char*>("GF_UNKNOWN_FREE_TYPE");
            case GF_INVALID_TYPE:
                return const_cast<char*>("GF_INVALID_TYPE");
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