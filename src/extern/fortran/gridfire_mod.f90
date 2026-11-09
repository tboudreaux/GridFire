! Fortran interface to GridFire's C API (gridfire/extern/gridfire_extern.h). Compile this file with the compiler
! that builds your program (gf-config --fortran-source gives its installed path) and link against gridfire_extern.
!
! Thread safety, as for the C API: set GridFire objects up (gff_init through gff_setup_solver) and free them from
! serial code; once set up, separate objects may evolve and evaluate concurrently from OpenMP threads (one object per
! thread at a time). MULTI_ZONE objects only from serial code.
module gridfire_mod
    use iso_c_binding
    implicit none

    type, public :: GF_TYPE
        integer(c_int) :: value
    end type GF_TYPE

    type(GF_TYPE), parameter, public :: &
            SINGLE_ZONE   = GF_TYPE(1001), &
            MULTI_ZONE    = GF_TYPE(1002)

    enum, bind (C)
        enumerator :: FDSSE_NON_4DSTAR_ERROR = -102
        enumerator :: FDSSE_UNKNOWN_ERROR = -101

        enumerator :: FDSSE_SUCCESS = 1
        enumerator :: FDSSE_UNKNOWN_SYMBOL_ERROR = 100

        enumerator :: FDSSE_SPECIES_ERROR = 101
        enumerator :: FDSSE_INVALID_COMPOSITION_ERROR = 102
        enumerator :: FDSSE_COMPOSITION_ERROR = 103

        enumerator :: GF_NON_GRIDFIRE_ERROR = -2
        enumerator :: GF_UNKNOWN_ERROR = -1
        enumerator :: GF_SUCCESS = 0

        enumerator :: GF_INVALID_QSE_SOLUTION_ERROR = 5
        enumerator :: GF_FAILED_TO_PARTITION_ERROR = 6
        enumerator :: GF_NETWORK_RESIZED_ERROR = 7
        enumerator :: GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR = 8
        enumerator :: GF_BAD_COLLECTION_ERROR = 9
        enumerator :: GF_BAD_RHS_ENIGNE_ERROR = 10
        enumerator :: GF_STALE_JACOBIAN_ERROR = 11
        enumerator :: GF_UNINITIALIZED_JACOBIAN_ERROR = 12
        enumerator :: GF_UNKNONWN_JACOBIAN_ERROR = 13
        enumerator :: GF_JACOBIAN_ERROR = 14
        enumerator :: GF_ENGINE_ERROR = 15

        enumerator :: GF_MISSING_BASE_REACTION_ERROR = 16
        enumerator :: GF_MISSING_SEED_SPECIES_ERROR = 17
        enumerator :: GF_MISSING_KEY_REACTION_ERROR = 18
        enumerator :: GF_POLICY_ERROR = 19

        enumerator :: GF_REACTION_PARSING_ERROR = 20
        enumerator :: GF_REACTOION_ERROR = 21

        enumerator :: GF_SINGULAR_JACOBIAN_ERROR = 22
        enumerator :: GF_ILL_CONDITIONED_JACOBIAN_ERROR = 23
        enumerator :: GF_CVODE_SOLVER_FAILURE_ERROR = 24
        enumerator :: GF_KINSOL_SOLVER_FAILURE_ERROR = 25
        enumerator :: GF_SUNDIALS_ERROR = 26
        enumerator :: GF_SOLVER_ERROR = 27

        enumerator :: GF_HASHING_ERROR = 28
        enumerator :: GF_UTILITY_ERROR = 29

        enumerator :: GF_DEBUG_ERRROR = 30

        enumerator :: GF_GRIDFIRE_ERROR = 31
        enumerator :: GF_UNINITIALIZED_INPUT_MEMORY_ERROR = 32
        enumerator :: GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR = 33

        enumerator :: GF_INVALD_NUM_SPECIES = 34
        enumerator :: GF_INVALID_TIMESTEPS = 35
        enumerator :: GF_UNKNONWN_FREE_TYPE = 36

        enumerator :: GF_INVALID_TYPE = 37

        enumerator :: GF_CONFIG_ERROR = 38
        enumerator :: GF_INVALID_STATE = 39
        enumerator :: GF_ZONE_FAILURE = 40

        enumerator :: GF_SINGLE_ZONE = 1001
        enumerator :: GF_MULTI_ZONE = 1002
        ! Correctly spelled names for codes whose original names above are misspelled (both are kept, same values)
        enumerator :: GF_FAILED_TO_PARTITION_ENGINE_ERROR = 6
        enumerator :: GF_BAD_RHS_ENGINE_ERROR = 10
        enumerator :: GF_UNKNOWN_JACOBIAN_ERROR = 13
        enumerator :: GF_REACTION_ERROR = 21
        enumerator :: GF_DEBUG_ERROR = 30
        enumerator :: GF_INVALID_NUM_SPECIES = 34
        enumerator :: GF_UNKNOWN_FREE_TYPE = 36
    end enum

    ! Values for GFConfig%engine_screening (enum GF_SCREENING_TYPE)
    enum, bind(C)
        enumerator :: GF_SCREENING_BARE = 0
        enumerator :: GF_SCREENING_WEAK = 1
        enumerator :: GF_SCREENING_CD09 = 2
        enumerator :: GF_SCREENING_BLENDED = 3
    end enum

    ! Values for GFConfig%engine_jacobian (enum GF_JACOBIAN_METHOD)
    enum, bind(C)
        enumerator :: GF_JACOBIAN_ANALYTIC = 0
        enumerator :: GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION = 1
    end enum

    ! Run-time configuration; same layout as struct GFConfig in gridfire_extern.h. Start from gff_default_config(),
    ! change fields, and pass it to gff_set_config before gff_setup_policy. Logical options are integers (0 = false).
    type, bind(C), public :: GFConfig
        ! solver.pointSolver
        real(c_double) :: solver_abs_tol
        real(c_double) :: solver_rel_tol
        ! solver.pointSolver.trigger
        real(c_double) :: trigger_off_diagonal_threshold
        real(c_double) :: trigger_timestep_collapse_ratio
        real(c_double) :: trigger_max_convergence_failures
        real(c_double) :: trigger_boundary_flux_relative_threshold
        real(c_double) :: trigger_boundary_flux_absolute_threshold
        ! engine.views.adaptiveEngineView
        real(c_double) :: engine_adaptive_relative_culling_threshold
        ! engine.screeningBlend: BLENDED screening uses Salpeter below engine_screening_blend_weak (weak-screening
        ! exponent H) and CD09 above engine_screening_blend_strong (0 < weak < strong <= 2)
        real(c_double) :: engine_screening_blend_weak
        real(c_double) :: engine_screening_blend_strong
        ! solver.pointSolver
        integer(c_int) :: solver_clamp_negative_iterates
        integer(c_int) :: solver_max_step_failure_restarts
        ! engine
        integer(c_int) :: engine_screening        ! GF_SCREENING_BARE, _WEAK, _CD09 or _BLENDED
        integer(c_int) :: engine_jacobian         ! GF_JACOBIAN_ANALYTIC or GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION
        integer(c_int) :: engine_split_hydrogen_pp_cno
        integer(c_int) :: engine_use_reverse_reactions  ! reverse rates by detailed balance (default 0)
    end type GFConfig

    interface
        ! void* gf_init()
        function gf_init(ctx_type) bind(C, name="gf_init")
            import :: c_ptr, c_int
            type(c_ptr) :: gf_init
            integer(c_int), value :: ctx_type
        end function gf_init

        ! int gf_free(void* gf)
        function gf_free(ctx_type, ptr) result(c_res) bind(C, name="gf_free")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            integer(c_int), value :: ctx_type
            integer(c_int) :: c_res
        end function gf_free

        function gf_set_num_zones(ctx_type, ptr, num_zones) result(c_res) bind(C, name="gf_set_num_zones")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_int), value :: ctx_type
            integer(c_size_t), value :: num_zones
            integer(c_int) :: c_res
        end function gf_set_num_zones

        ! char* gf_get_last_error_message(void* ptr);
        function gf_get_last_error_message(ptr) result(c_msg) bind(C, name="gf_get_last_error_message")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            type(c_ptr) :: c_msg
        end function

        ! const char* gf_get_version(void);
        function gf_get_version() result(c_str) bind(C, name="gf_get_version")
            import :: c_ptr
            type(c_ptr) :: c_str
        end function

        ! int gf_get_provenance(void* ptr, const char** description);
        function gf_get_provenance(ptr, description) result(ierr) bind(C, name="gf_get_provenance")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            type(c_ptr), intent(out) :: description
            integer(c_int) :: ierr
        end function

        ! size_t strlen(const char*), to copy C strings of any length
        function c_strlen(str) result(n) bind(C, name="strlen")
            import :: c_ptr, c_size_t
            type(c_ptr), value :: str
            integer(c_size_t) :: n
        end function

        ! int gf_config_defaults(struct GFConfig* config);
        function gf_config_defaults(config) result(ierr) bind(C, name="gf_config_defaults")
            import :: c_int, GFConfig
            type(GFConfig), intent(out) :: config
            integer(c_int) :: ierr
        end function

        ! int gf_set_config(void* ptr, const struct GFConfig* config);
        function gf_set_config(ptr, config) result(ierr) bind(C, name="gf_set_config")
            import :: c_ptr, c_int, GFConfig
            type(c_ptr), value :: ptr
            type(GFConfig), intent(in) :: config
            integer(c_int) :: ierr
        end function

        ! int gf_get_config(void* ptr, struct GFConfig* config);
        function gf_get_config(ptr, config) result(ierr) bind(C, name="gf_get_config")
            import :: c_ptr, c_int, GFConfig
            type(c_ptr), value :: ptr
            type(GFConfig), intent(out) :: config
            integer(c_int) :: ierr
        end function

        ! int gf_load_config(void* ptr, const char* path);
        function gf_load_config(ptr, path) result(ierr) bind(C, name="gf_load_config")
            import :: c_ptr, c_int, c_char
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int) :: ierr
        end function

        ! int gf_save_config(void* ptr, const char* path);
        function gf_save_config(ptr, path) result(ierr) bind(C, name="gf_save_config")
            import :: c_ptr, c_int, c_char
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int) :: ierr
        end function

        ! int gf_set_passive_species(void* ptr, size_t num_species, const char** species_names);
        function gf_set_passive_species(ptr, num_species, species_names) result(ierr) bind(C, name="gf_set_passive_species")
            import
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: num_species
            type(c_ptr), dimension(*), intent(in) :: species_names
            integer(c_int) :: ierr
        end function

        ! int gf_mass_fractions_to_molar(void* ptr, const double* X, size_t num_species, size_t num_zones, double* Y);
        function gf_mass_fractions_to_molar(ptr, X, num_species, num_zones, Y) result(ierr) &
                bind(C, name="gf_mass_fractions_to_molar")
            import
            type(c_ptr), value :: ptr
            real(c_double), dimension(*), intent(in) :: X
            integer(c_size_t), value :: num_species, num_zones
            real(c_double), dimension(*), intent(out) :: Y
            integer(c_int) :: ierr
        end function

        ! int gf_molar_to_mass_fractions(void* ptr, const double* Y, size_t num_species, size_t num_zones, double* X);
        function gf_molar_to_mass_fractions(ptr, Y, num_species, num_zones, X) result(ierr) &
                bind(C, name="gf_molar_to_mass_fractions")
            import
            type(c_ptr), value :: ptr
            real(c_double), dimension(*), intent(in) :: Y
            integer(c_size_t), value :: num_species, num_zones
            real(c_double), dimension(*), intent(out) :: X
            integer(c_int) :: ierr
        end function

        ! int gf_register_species(void* ptr, const int num_species, const char** species_names);
        function gf_register_species(ptr, num_species, species_names) result(ierr) bind(C, name="gf_register_species")
            import
            type(c_ptr), value :: ptr
            integer(c_int), value :: num_species
            type(c_ptr), dimension(*), intent(in) :: species_names ! Array of C pointers
            integer(c_int) :: ierr
        end function

        ! int gf_construct_engine_from_policy(void* ptr, const char* policy_name, const double *abundances, size_t num_species);
        function gf_construct_engine_from_policy(ptr, policy_name, abundances, num_species) result(ierr) &
                bind(C, name="gf_construct_engine_from_policy")
            import
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: policy_name
            real(c_double), dimension(*), intent(in) :: abundances
            integer(c_size_t), value :: num_species
            integer(c_int) :: ierr
        end function

        ! int gf_construct_solver_from_engine(void* ptr);
        function gf_construct_solver_from_engine(ptr) result(ierr) &
                bind(C, name="gf_construct_solver_from_engine")
            import
            type(c_ptr), value :: ptr
            integer(c_int) :: ierr
        end function

        ! int gf_evolve(...)
        function gf_evolve_c_scalar(ctx_type, ptr, Y_in, num_species, T, rho, tMax, dt0, &
                Y_out, energy, dedt, dedrho, &
                nue_loss, nu_flux, mass_lost, &
                eps_nuc, dlneps_dlnt, dlneps_dlnrho, nue_loss_rate, eps_by_channel) result(ierr) &
                bind(C, name="gf_evolve")
            import :: c_ptr, c_int, c_double, c_size_t
            type(c_ptr), value :: ptr
            integer(c_int), value :: ctx_type
            integer(c_size_t), value :: num_species

            ! Arrays
            real(c_double), dimension(*), intent(in) :: Y_in
            real(c_double), dimension(*), intent(out) :: Y_out
            real(c_double), dimension(*), intent(out) :: eps_by_channel

            ! Scalars (Passed by Reference -> matches void*)
            real(c_double), intent(in) :: T, rho
            real(c_double), intent(out) :: energy, dedt, dedrho, nue_loss, nu_flux, mass_lost
            real(c_double), intent(out) :: eps_nuc, dlneps_dlnt, dlneps_dlnrho, nue_loss_rate

            ! Scalars (Passed by Value)
            real(c_double), value :: tMax, dt0

            integer(c_int) :: ierr
        end function

        ! 2. Interface for Multi Zone (Arrays)
        function gf_evolve_c_array(ctx_type, ptr, Y_in, num_species, T, rho, tMax, dt0, &
                Y_out, energy, dedt, dedrho, &
                nue_loss, nu_flux, mass_lost, &
                eps_nuc, dlneps_dlnt, dlneps_dlnrho, nue_loss_rate, eps_by_channel) result(ierr) &
                bind(C, name="gf_evolve")
            import :: c_ptr, c_int, c_double, c_size_t
            type(c_ptr), value :: ptr
            integer(c_int), value :: ctx_type
            integer(c_size_t), value :: num_species

            ! All Arrays (dimension(*))
            real(c_double), dimension(*), intent(in) :: Y_in
            real(c_double), dimension(*), intent(in) :: T, rho

            real(c_double), dimension(*), intent(out) :: Y_out
            real(c_double), dimension(*), intent(out) :: energy, dedt, dedrho, nue_loss, nu_flux, mass_lost
            real(c_double), dimension(*), intent(out) :: eps_nuc, dlneps_dlnt, dlneps_dlnrho, nue_loss_rate
            real(c_double), dimension(*), intent(out) :: eps_by_channel

            ! Scalars (Passed by Value)
            real(c_double), value :: tMax, dt0

            integer(c_int) :: ierr
        end function

        ! int gf_evaluate(...): T and rho are scalars for SINGLE_ZONE and arrays for MULTI_ZONE; both are passed by
        ! reference, so one interface serves both. Optional outputs are passed as C pointers (c_null_ptr to skip).
        function gf_evaluate_c(ctx_type, ptr, Y, num_species, T, rho, eps_nuc, dedt, dedrho, &
                nue_loss_rate, eps_by_channel, dydt) result(ierr) bind(C, name="gf_evaluate")
            import :: c_ptr, c_int, c_double, c_size_t
            integer(c_int), value :: ctx_type
            type(c_ptr), value :: ptr
            real(c_double), dimension(*), intent(in) :: Y
            integer(c_size_t), value :: num_species
            real(c_double), dimension(*), intent(in) :: T, rho
            real(c_double), dimension(*), intent(out) :: eps_nuc, dedt, dedrho
            type(c_ptr), value :: nue_loss_rate, eps_by_channel, dydt
            integer(c_int) :: ierr
        end function

        ! int gf_set_shared_composition_groups(void* ptr, const int* group, const double* mass);
        function gf_set_shared_composition_groups(ptr, group, mass) result(ierr) bind(C, name="gf_set_shared_composition_groups")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            type(c_ptr), value :: group, mass
            integer(c_int) :: ierr
        end function

        ! int gf_shared_composition_groups_from_mask(const int* flagged, size_t num_zones, int* group);
        function gf_shared_composition_groups_from_mask(flagged, num_zones, group) result(ierr) &
                bind(C, name="gf_shared_composition_groups_from_mask")
            import :: c_int, c_size_t
            integer(c_int), dimension(*), intent(in) :: flagged
            integer(c_size_t), value :: num_zones
            integer(c_int), dimension(*), intent(out) :: group
            integer(c_int) :: ierr
        end function

        ! int gf_get_zone_status(void* ptr, int* status, double* time_reached);
        function gf_get_zone_status(ptr, status, time_reached) result(ierr) bind(C, name="gf_get_zone_status")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            type(c_ptr), value :: status, time_reached
            integer(c_int) :: ierr
        end function

        ! int gf_evaluate_reactions(...): T and rho by reference (scalar or array); outputs as C pointers (c_null_ptr to skip)
        function gf_evaluate_reactions_c(ctx_type, ptr, Y, num_species, T, rho, molar_flow, energy_rate, &
                nue_loss_rate, nu_rate) result(ierr) bind(C, name="gf_evaluate_reactions")
            import :: c_ptr, c_int, c_double, c_size_t
            integer(c_int), value :: ctx_type
            type(c_ptr), value :: ptr
            real(c_double), dimension(*), intent(in) :: Y
            integer(c_size_t), value :: num_species
            real(c_double), dimension(*), intent(in) :: T, rho
            type(c_ptr), value :: molar_flow, energy_rate, nue_loss_rate, nu_rate
            integer(c_int) :: ierr
        end function

        ! int gf_get_num_network_species(void* ptr, size_t* num_species);
        function gf_get_num_network_species(ptr, num_species) result(ierr) bind(C, name="gf_get_num_network_species")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_size_t), intent(out) :: num_species
            integer(c_int) :: ierr
        end function

        ! int gf_get_network_species(void* ptr, size_t index, const char** name, int* A, int* Z, double* mass);
        function gf_get_network_species(ptr, index, name, A, Z, mass) result(ierr) bind(C, name="gf_get_network_species")
            import :: c_ptr, c_int, c_size_t, c_double
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: index
            type(c_ptr), intent(out) :: name
            integer(c_int), intent(out) :: A, Z
            real(c_double), intent(out) :: mass
            integer(c_int) :: ierr
        end function

        ! int gf_register_network_species(void* ptr);
        function gf_register_network_species(ptr) result(ierr) bind(C, name="gf_register_network_species")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            integer(c_int) :: ierr
        end function

        ! int gf_get_num_reactions(void* ptr, size_t* num_reactions);
        function gf_get_num_reactions(ptr, num_reactions) result(ierr) bind(C, name="gf_get_num_reactions")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_size_t), intent(out) :: num_reactions
            integer(c_int) :: ierr
        end function

        ! int gf_get_reaction_id(void* ptr, size_t index, const char** id);
        function gf_get_reaction_id(ptr, index, id) result(ierr) bind(C, name="gf_get_reaction_id")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: index
            type(c_ptr), intent(out) :: id
            integer(c_int) :: ierr
        end function

        ! int gf_set_rate_multiplier(void* ptr, const char* reaction_id, double factor);
        function gf_set_rate_multiplier(ptr, reaction_id, factor) result(ierr) bind(C, name="gf_set_rate_multiplier")
            import :: c_ptr, c_int, c_char, c_double
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: reaction_id
            real(c_double), value :: factor
            integer(c_int) :: ierr
        end function

        ! int gf_get_rate_multiplier(void* ptr, const char* reaction_id, double* factor);
        function gf_get_rate_multiplier(ptr, reaction_id, factor) result(ierr) bind(C, name="gf_get_rate_multiplier")
            import :: c_ptr, c_int, c_char, c_double
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: reaction_id
            real(c_double), intent(out) :: factor
            integer(c_int) :: ierr
        end function

        ! int gf_clear_rate_multipliers(void* ptr);
        function gf_clear_rate_multipliers(ptr) result(ierr) bind(C, name="gf_clear_rate_multipliers")
            import :: c_ptr, c_int
            type(c_ptr), value :: ptr
            integer(c_int) :: ierr
        end function

        ! int gf_get_num_energy_channels(void* ptr, size_t* num_channels);
        function gf_get_num_energy_channels(ptr, num_channels) result(ierr) bind(C, name="gf_get_num_energy_channels")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_size_t), intent(out) :: num_channels
            integer(c_int) :: ierr
        end function

        ! int gf_get_energy_channel_name(void* ptr, size_t index, const char** name);
        function gf_get_energy_channel_name(ptr, index, name) result(ierr) bind(C, name="gf_get_energy_channel_name")
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: index
            type(c_ptr), intent(out) :: name
            integer(c_int) :: ierr
        end function
    end interface

    type :: GridFire
        type(c_ptr) :: ctx = c_null_ptr
        integer(c_int) :: ctx_type = SINGLE_ZONE%value
        integer(c_size_t) :: num_species = 0
        integer(c_size_t) :: num_zones = 1
    contains
        procedure :: gff_init
        procedure :: gff_free
        procedure :: gff_register_species
        procedure :: gff_set_config
        procedure :: gff_get_config
        procedure :: gff_load_config
        procedure :: gff_save_config
        procedure :: gff_set_passive_species
        procedure :: gff_mass_fractions_to_molar_single
        procedure :: gff_mass_fractions_to_molar_multi
        generic :: gff_mass_fractions_to_molar => gff_mass_fractions_to_molar_single, gff_mass_fractions_to_molar_multi
        procedure :: gff_molar_to_mass_fractions_single
        procedure :: gff_molar_to_mass_fractions_multi
        generic :: gff_molar_to_mass_fractions => gff_molar_to_mass_fractions_single, gff_molar_to_mass_fractions_multi
        procedure :: gff_setup_policy
        procedure :: gff_setup_solver
        procedure :: gff_get_last_error
        procedure :: gff_get_provenance

        procedure :: gff_evolve_single
        procedure :: gff_evolve_multi

        generic :: gff_evolve => gff_evolve_single, gff_evolve_multi

        procedure :: gff_evaluate_single
        procedure :: gff_evaluate_multi
        generic :: gff_evaluate => gff_evaluate_single, gff_evaluate_multi

        procedure :: gff_evaluate_reactions_single
        procedure :: gff_evaluate_reactions_multi
        generic :: gff_evaluate_reactions => gff_evaluate_reactions_single, gff_evaluate_reactions_multi

        procedure :: gff_get_zone_status
        procedure :: gff_set_shared_composition_groups
        procedure :: gff_clear_shared_composition_groups
        procedure :: gff_get_num_network_species
        procedure :: gff_get_network_species
        procedure :: gff_register_network_species
        procedure :: gff_get_num_reactions
        procedure :: gff_get_reaction_ids
        procedure :: gff_set_rate_multiplier
        procedure :: gff_get_rate_multiplier
        procedure :: gff_clear_rate_multipliers

        procedure :: gff_get_num_energy_channels
        procedure :: gff_get_energy_channel_names
    end type GridFire

    contains
        subroutine gff_init(self, type, zones)
            class(GridFire), intent(out) :: self
            type(GF_TYPE), intent(in)    :: type
            integer(c_size_t), intent(in), optional :: zones
            integer(c_int) :: ierr

            if (type%value==1002) then
                if (.not. present(zones)) then
                    print *, "GridFire Error: Multi-zone type requires number of zones to be specficied in the GridFire init method (i.e. GridFire(MULTI_ZONE, 10) for 10 zones)."
                    error stop
                end if

                self%num_zones = zones
            end if

            self%ctx_type = type%value

            self%ctx = gf_init(self%ctx_type)

            if (type%value==1002) then
                ierr = gf_set_num_zones(self%ctx_type, self%ctx, self%num_zones)
                if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                    print *, "GridFire Multi-Zone Error: ", self%gff_get_last_error()
                    error stop
                end if
            end if
        end subroutine gff_init

        subroutine gff_free(self)
            class(GridFire), intent(inout) :: self
            integer(c_int) :: ierr

            if (c_associated(self%ctx)) then
                ierr = gf_free(self%ctx_type, self%ctx)
                if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                    print *, "GridFire Free Error: ", self%gff_get_last_error()
                    error stop
                end if
                self%ctx = c_null_ptr
            end if
        end subroutine gff_free

        function gff_get_last_error(self) result(msg)
            class(GridFire), intent(in) :: self
            character(len=:), allocatable :: msg
            type(c_ptr) :: c_msg_ptr
            character(kind=c_char), pointer :: char_ptr(:)
            integer :: i, len_str

            c_msg_ptr = gf_get_last_error_message(self%ctx)
            if (.not. c_associated(c_msg_ptr)) then
                msg = "GridFire: Unknown Error (Null Pointer returned)"
                return
            end if

            call c_f_pointer(c_msg_ptr, char_ptr, [1024])
            len_str = 0
            do i = 1, 1024
                if (char_ptr(i) == c_null_char) exit
                len_str = len_str + 1
            end do

            msg = repeat(' ', len_str+10)
            msg(1:10) = "GridFire: "
            do i = 1, len_str
                msg(i+10:i+10) = char_ptr(i)
            end do
        end function gff_get_last_error

        ! A C string as a Fortran string ("" for a null pointer)
        function c_to_fortran_string(c_str) result(str)
            type(c_ptr), intent(in) :: c_str
            character(len=:), allocatable :: str
            character(kind=c_char), pointer :: chars(:)
            integer :: i, n

            if (.not. c_associated(c_str)) then
                str = ""
                return
            end if
            n = int(c_strlen(c_str))
            call c_f_pointer(c_str, chars, [max(n, 1)])
            allocate(character(len=n) :: str)
            do i = 1, n
                str(i:i) = chars(i)
            end do
        end function c_to_fortran_string

        ! GridFire's version, e.g. "1.2.0"
        function gff_get_version() result(version)
            character(len=:), allocatable :: version
            version = c_to_fortran_string(gf_get_version())
        end function gff_get_version

        ! What this context computes with, one "key: value" line each (see gf_get_provenance in gridfire_extern.h);
        ! lines are separated by new_line('a'), ready to write into an output header.
        function gff_get_provenance(self) result(description)
            class(GridFire), intent(in) :: self
            character(len=:), allocatable :: description
            type(c_ptr) :: c_str
            integer(c_int) :: ierr

            ierr = gf_get_provenance(self%ctx, c_str)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Error: ", self%gff_get_last_error()
                error stop
            end if
            description = c_to_fortran_string(c_str)
        end function gff_get_provenance

        ! GridFire's default configuration, to modify and pass to gff_set_config.
        function gff_default_config() result(config)
            type(GFConfig) :: config
            integer(c_int) :: ierr

            ierr = gf_config_defaults(config)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: could not read the default configuration."
                error stop
            end if
        end function gff_default_config

        ! Use config for the engine and solver; call before gff_setup_policy.
        subroutine gff_set_config(self, config)
            class(GridFire), intent(in) :: self
            type(GFConfig), intent(in) :: config
            integer(c_int) :: ierr

            ierr = gf_set_config(self%ctx, config)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_set_config

        function gff_get_config(self) result(config)
            class(GridFire), intent(in) :: self
            type(GFConfig) :: config
            integer(c_int) :: ierr

            ierr = gf_get_config(self%ctx, config)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: ", self%gff_get_last_error()
                error stop
            end if
        end function gff_get_config

        ! Read the configuration from a TOML file (as written by gff_save_config); call before gff_setup_policy.
        subroutine gff_load_config(self, path)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: path
            integer(c_int) :: ierr

            ierr = gf_load_config(self%ctx, trim(path) // c_null_char)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_load_config

        subroutine gff_save_config(self, path)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: path
            integer(c_int) :: ierr

            ierr = gf_save_config(self%ctx, trim(path) // c_null_char)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_save_config

        ! Species carried without reactions (count towards Y_e and screening, never change). They must also be passed
        ! to gff_register_species, and none may take part in a reaction. Call before gff_setup_policy.
        subroutine gff_set_passive_species(self, species_list)
            class(GridFire), intent(in) :: self
            character(len=*), dimension(:), intent(in) :: species_list

            type(c_ptr), allocatable, dimension(:) :: c_ptrs
            character(kind=c_char, len=:), allocatable, target :: temp_strs(:)
            integer :: i, n
            integer(c_int) :: ierr

            n = size(species_list)
            allocate(c_ptrs(max(n, 1)))
            allocate(character(len=len(species_list)+1) :: temp_strs(max(n, 1)))
            do i = 1, n
                temp_strs(i) = trim(species_list(i)) // c_null_char
                c_ptrs(i) = c_loc(temp_strs(i))
            end do

            ierr = gf_set_passive_species(self%ctx, int(n, c_size_t), c_ptrs)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Config Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_set_passive_species

        ! Scale one reaction's rate (an id from gff_get_reaction_ids) by factor >= 0 in everything computed afterwards;
        ! 1 removes it, 0 switches the reaction off. Call after gff_setup_policy. With ierr present, failures (e.g.
        ! GF_REACTION_ERROR for a reaction not in the network) are returned; without it they stop the program.
        subroutine gff_set_rate_multiplier(self, reaction_id, factor, ierr)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: reaction_id
            real(c_double), intent(in) :: factor
            integer, intent(out), optional :: ierr
            integer(c_int) :: c_ierr

            c_ierr = gf_set_rate_multiplier(self%ctx, trim(reaction_id) // c_null_char, factor)
            call gff_report_or_stop(self, c_ierr, "GridFire Rate Multiplier Error: ", ierr)
        end subroutine gff_set_rate_multiplier

        ! The factor set for a reaction (1 if none); errors as for gff_set_rate_multiplier.
        subroutine gff_get_rate_multiplier(self, reaction_id, factor, ierr)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: reaction_id
            real(c_double), intent(out) :: factor
            integer, intent(out), optional :: ierr
            integer(c_int) :: c_ierr

            c_ierr = gf_get_rate_multiplier(self%ctx, trim(reaction_id) // c_null_char, factor)
            call gff_report_or_stop(self, c_ierr, "GridFire Rate Multiplier Error: ", ierr)
        end subroutine gff_get_rate_multiplier

        ! Remove every rate multiplier.
        subroutine gff_clear_rate_multipliers(self, ierr)
            class(GridFire), intent(in) :: self
            integer, intent(out), optional :: ierr
            integer(c_int) :: c_ierr

            c_ierr = gf_clear_rate_multipliers(self%ctx)
            call gff_report_or_stop(self, c_ierr, "GridFire Rate Multiplier Error: ", ierr)
        end subroutine gff_clear_rate_multipliers

        ! Hand an error code back through ierr when the caller passed one; otherwise stop on failure with the message
        subroutine gff_report_or_stop(self, c_ierr, prefix, ierr)
            class(GridFire), intent(in) :: self
            integer(c_int), intent(in) :: c_ierr
            character(len=*), intent(in) :: prefix
            integer, intent(out), optional :: ierr

            if (present(ierr)) then
                ierr = int(c_ierr)
            else if (c_ierr /= GF_SUCCESS) then
                print *, prefix, self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_report_or_stop

        ! Mass fractions to molar abundances (registered species order), Y_i = X_i / m_i; X must sum to one to 1e-10.
        ! GridFire itself takes molar abundances only; ierr is GF_INVALID_STATE for a bad X.
        subroutine gff_mass_fractions_to_molar_single(self, X, Y, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:), intent(in) :: X
            real(c_double), dimension(:), intent(out) :: Y
            integer, intent(out) :: ierr
            ierr = int(gf_mass_fractions_to_molar(self%ctx, X, int(size(X), c_size_t), 1_c_size_t, Y))
        end subroutine gff_mass_fractions_to_molar_single

        ! As above for X(num_species, num_zones)
        subroutine gff_mass_fractions_to_molar_multi(self, X, Y, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:,:), intent(in) :: X
            real(c_double), dimension(:,:), intent(out) :: Y
            integer, intent(out) :: ierr
            ierr = int(gf_mass_fractions_to_molar(self%ctx, X, int(size(X, 1), c_size_t), int(size(X, 2), c_size_t), Y))
        end subroutine gff_mass_fractions_to_molar_multi

        ! Molar abundances to mass fractions, X_i = Y_i m_i / sum_j Y_j m_j (each zone sums to one; see gridfire_extern.h)
        subroutine gff_molar_to_mass_fractions_single(self, Y, X, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:), intent(in) :: Y
            real(c_double), dimension(:), intent(out) :: X
            integer, intent(out) :: ierr
            ierr = int(gf_molar_to_mass_fractions(self%ctx, Y, int(size(Y), c_size_t), 1_c_size_t, X))
        end subroutine gff_molar_to_mass_fractions_single

        subroutine gff_molar_to_mass_fractions_multi(self, Y, X, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:,:), intent(in) :: Y
            real(c_double), dimension(:,:), intent(out) :: X
            integer, intent(out) :: ierr
            ierr = int(gf_molar_to_mass_fractions(self%ctx, Y, int(size(Y, 1), c_size_t), int(size(Y, 2), c_size_t), X))
        end subroutine gff_molar_to_mass_fractions_multi

        subroutine gff_register_species(self, species_list)
            class(GridFire), intent(inout) :: self
            character(len=*), dimension(:), intent(in) :: species_list

            type(c_ptr), allocatable, dimension(:) :: c_ptrs
            character(kind=c_char, len=:), allocatable, target :: temp_strs(:)
            integer :: i, n, ierr

            n = size(species_list)
            self%num_species = int(n, c_size_t)

            allocate(c_ptrs(n))
            allocate(character(len=len(species_list(1))+1) :: temp_strs(n)) ! +1 for null terminator

            do i = 1, n
                temp_strs(i) = trim(species_list(i)) // c_null_char
                c_ptrs(i) = c_loc(temp_strs(i))
            end do

            ierr = gf_register_species(self%ctx, int(n, c_int), c_ptrs)
            if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                print *, "GridFire: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_register_species

        subroutine gff_setup_policy(self, policy_name, abundances)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: policy_name
            real(c_double), dimension(:), intent(in) :: abundances
            integer(c_int) :: ierr

            if (size(abundances) /= self%num_species) then
                print *, "GridFire Error: Abundance array size mismatch."
                error stop
            end if

            ierr = gf_construct_engine_from_policy(self%ctx, &
                    trim(policy_name) // c_null_char, &
                    abundances, &
                    self%num_species)

            if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                print *, "GridFire Policy Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_setup_policy

        ! Build the solver (CVODE; for MULTI_ZONE a grid of CVODE solvers). solver_name is optional and kept for
        ! existing callers: GridFire has one solver, so only "CVODE" (any case) is accepted.
        subroutine gff_setup_solver(self, solver_name)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in), optional :: solver_name
            integer(c_int) :: ierr
            character(len=:), allocatable :: name
            integer :: i

            if (present(solver_name)) then
                name = trim(adjustl(solver_name))
                do i = 1, len(name)
                    if (name(i:i) >= "a" .and. name(i:i) <= "z") name(i:i) = achar(iachar(name(i:i)) - 32)
                end do
                if (name /= "CVODE") then
                    print *, "GridFire Solver Error: unknown solver '", trim(solver_name), "' (the only solver is CVODE)"
                    error stop
                end if
            end if

            ierr = gf_construct_solver_from_engine(self%ctx)
            if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                print *, "GridFire Solver Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_setup_solver

        subroutine gff_evolve_single(self, Y_in, T, rho, tMax, dt0, Y_out, energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost, &
                eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate, eps_by_channel, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:), intent(in) :: Y_in
            real(c_double), intent(in) :: T, rho
            real(c_double), value :: tMax, dt0

            real(c_double), dimension(:), intent(out) :: Y_out
            real(c_double), intent(out) :: energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost
            real(c_double), intent(out) :: eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate
            real(c_double), dimension(:), intent(out) :: eps_by_channel ! size gff_get_num_energy_channels()
            integer, intent(out) :: ierr
            integer(c_int) :: c_ierr

            c_ierr = gf_evolve_c_scalar(self%ctx_type, self%ctx, &
                    Y_in, self%num_species, &
                    T, rho, tMax, dt0, &
                    Y_out, &
                    energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost, &
                    eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate, eps_by_channel)
            ierr = int(c_ierr)
        end subroutine gff_evolve_single

        subroutine gff_evolve_multi(self, Y_in, T, rho, tMax, dt0, Y_out, energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost, &
                eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate, eps_by_channel, ierr)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:,:), intent(in) :: Y_in
            real(c_double), dimension(:), intent(in) :: T, rho
            real(c_double), value :: tMax, dt0

            real(c_double), dimension(:,:), intent(out) :: Y_out
            real(c_double), dimension(:), intent(out) :: energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost
            real(c_double), dimension(:), intent(out) :: eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate
            real(c_double), dimension(:,:), intent(out) :: eps_by_channel ! (num_channels, num_zones)
            integer, intent(out) :: ierr
            integer(c_int) :: c_ierr

            c_ierr = gf_evolve_c_array(self%ctx_type, self%ctx, &
                    Y_in, self%num_species, &
                    T, rho, tMax, dt0, &
                    Y_out, &
                    energy, dedt, dedrho, nu_e_loss, nu_flux, mass_lost, &
                    eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_e_loss_rate, eps_by_channel)
            ierr = int(c_ierr)
        end subroutine gff_evolve_multi

        ! Copy a C string owned by the context into a blank-padded Fortran string.
        subroutine copy_c_string(c_str, f_str)
            type(c_ptr), intent(in) :: c_str
            character(len=*), intent(out) :: f_str
            character(kind=c_char), pointer :: chars(:)
            integer :: k

            f_str = ' '
            call c_f_pointer(c_str, chars, [len(f_str) + 1])
            do k = 1, len(f_str)
                if (chars(k) == c_null_char) exit
                f_str(k:k) = chars(k)
            end do
        end subroutine copy_c_string

        ! Instantaneous eps_nuc and its derivatives at fixed composition (no integration). Y in registered species order.
        subroutine gff_evaluate_single(self, Y, T, rho, eps_nuc, dedt, dedrho, ierr, nu_e_loss_rate, eps_by_channel, dydt)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:), intent(in) :: Y
            real(c_double), intent(in) :: T, rho
            real(c_double), intent(out) :: eps_nuc, dedt, dedrho
            integer, intent(out) :: ierr
            real(c_double), intent(out), optional, target :: nu_e_loss_rate
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: eps_by_channel ! num_channels
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: dydt ! num_species
            real(c_double) :: T_a(1), rho_a(1), eps_a(1), dedt_a(1), dedrho_a(1)
            type(c_ptr) :: p_nu, p_channels, p_dydt

            p_nu = c_null_ptr
            p_channels = c_null_ptr
            p_dydt = c_null_ptr
            if (present(nu_e_loss_rate)) p_nu = c_loc(nu_e_loss_rate)
            if (present(eps_by_channel)) p_channels = c_loc(eps_by_channel)
            if (present(dydt)) p_dydt = c_loc(dydt)

            T_a(1) = T
            rho_a(1) = rho
            ierr = int(gf_evaluate_c(self%ctx_type, self%ctx, Y, self%num_species, T_a, rho_a, &
                    eps_a, dedt_a, dedrho_a, p_nu, p_channels, p_dydt))
            eps_nuc = eps_a(1)
            dedt = dedt_a(1)
            dedrho = dedrho_a(1)
        end subroutine gff_evaluate_single

        ! As gff_evaluate_single for every zone: Y(num_species, num_zones), T(num_zones), ...
        subroutine gff_evaluate_multi(self, Y, T, rho, eps_nuc, dedt, dedrho, ierr, nu_e_loss_rate, eps_by_channel, dydt)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:,:), intent(in) :: Y
            real(c_double), dimension(:), intent(in) :: T, rho
            real(c_double), dimension(:), intent(out) :: eps_nuc, dedt, dedrho
            integer, intent(out) :: ierr
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: nu_e_loss_rate
            real(c_double), dimension(:,:), intent(out), optional, target, contiguous :: eps_by_channel ! (num_channels, num_zones)
            real(c_double), dimension(:,:), intent(out), optional, target, contiguous :: dydt ! (num_species, num_zones)
            type(c_ptr) :: p_nu, p_channels, p_dydt

            p_nu = c_null_ptr
            p_channels = c_null_ptr
            p_dydt = c_null_ptr
            if (present(nu_e_loss_rate)) p_nu = c_loc(nu_e_loss_rate)
            if (present(eps_by_channel)) p_channels = c_loc(eps_by_channel)
            if (present(dydt)) p_dydt = c_loc(dydt)

            ierr = int(gf_evaluate_c(self%ctx_type, self%ctx, Y, self%num_species, T, rho, &
                    eps_nuc, dedt, dedrho, p_nu, p_channels, p_dydt))
        end subroutine gff_evaluate_multi

        ! MULTI_ZONE: burn zones with the same group id (>= 0) as one shared composition in every later gff_evolve
        ! (a convective zone under instantaneous mixing); -1 burns a zone on its own. mass(i) is zone i's mass (any
        ! units), required for grouped zones. Every member gets the group's Y_out but its own energy and rates.
        subroutine gff_set_shared_composition_groups(self, group, mass)
            class(GridFire), intent(in) :: self
            integer(c_int), dimension(:), intent(in), target, contiguous :: group
            real(c_double), dimension(:), intent(in), target, contiguous :: mass
            integer(c_int) :: ierr

            if (size(group) /= self%num_zones .or. size(mass) /= self%num_zones) then
                print *, "GridFire Error: shared-composition group and mass arrays need one entry per zone."
                error stop
            end if
            ierr = gf_set_shared_composition_groups(self%ctx, c_loc(group), c_loc(mass))
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Shared Composition Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_set_shared_composition_groups

        ! Burn every zone on its own again.
        subroutine gff_clear_shared_composition_groups(self)
            class(GridFire), intent(in) :: self
            integer(c_int) :: ierr

            ierr = gf_set_shared_composition_groups(self%ctx, c_null_ptr, c_null_ptr)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Shared Composition Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_clear_shared_composition_groups

        ! Group ids from a per-zone flag (e.g. 1 = convective) for zones ordered in space: each contiguous run of
        ! non-zero flags is one group, every other zone -1.
        function gff_shared_composition_groups_from_mask(flagged) result(group)
            integer(c_int), dimension(:), intent(in) :: flagged
            integer(c_int), dimension(size(flagged)) :: group
            integer(c_int) :: ierr

            ierr = gf_shared_composition_groups_from_mask(flagged, int(size(flagged), c_size_t), group)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Shared Composition Error: could not build groups from the mask."
                error stop
            end if
        end function gff_shared_composition_groups_from_mask

        ! Outcome of each zone in the latest gff_evolve call: status(i) is GF_SUCCESS or the zone's error code and
        ! time_reached(i) how far it was integrated [s]. Arrays have num_zones entries (1 for SINGLE_ZONE); both optional.
        ! A failed zone's Y_out holds its input abundances and its other outputs are 0. A MULTI_ZONE gff_evolve returns
        ! GF_ZONE_FAILURE in ierr if any zone failed.
        subroutine gff_get_zone_status(self, status, time_reached)
            class(GridFire), intent(in) :: self
            integer(c_int), dimension(:), intent(out), optional, target, contiguous :: status
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: time_reached
            type(c_ptr) :: p_status, p_time
            integer(c_int) :: ierr

            p_status = c_null_ptr
            p_time = c_null_ptr
            if (present(status)) p_status = c_loc(status)
            if (present(time_reached)) p_time = c_loc(time_reached)
            ierr = gf_get_zone_status(self%ctx, p_status, p_time)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Zone Status Error: ", self%gff_get_last_error()
                error stop
            end if
        end subroutine gff_get_zone_status

        ! Each reaction's molar flow [mol/g/s], energy rate net of neutrinos [erg/g/s], neutrino energy loss rate
        ! [erg/g/s] and neutrino emission rate [1/g/s] at fixed composition (no integration). Arrays have
        ! gff_get_num_reactions() entries, in gff_get_reaction_ids order; all four are optional.
        subroutine gff_evaluate_reactions_single(self, Y, T, rho, ierr, molar_flow, energy_rate, nu_e_loss_rate, nu_rate)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:), intent(in) :: Y
            real(c_double), intent(in) :: T, rho
            integer, intent(out) :: ierr
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: molar_flow, energy_rate
            real(c_double), dimension(:), intent(out), optional, target, contiguous :: nu_e_loss_rate, nu_rate
            real(c_double) :: T_a(1), rho_a(1)
            type(c_ptr) :: p_flow, p_energy, p_loss, p_rate

            p_flow = c_null_ptr
            p_energy = c_null_ptr
            p_loss = c_null_ptr
            p_rate = c_null_ptr
            if (present(molar_flow)) p_flow = c_loc(molar_flow)
            if (present(energy_rate)) p_energy = c_loc(energy_rate)
            if (present(nu_e_loss_rate)) p_loss = c_loc(nu_e_loss_rate)
            if (present(nu_rate)) p_rate = c_loc(nu_rate)
            T_a(1) = T
            rho_a(1) = rho
            ierr = int(gf_evaluate_reactions_c(self%ctx_type, self%ctx, Y, self%num_species, T_a, rho_a, &
                    p_flow, p_energy, p_loss, p_rate))
        end subroutine gff_evaluate_reactions_single

        ! As gff_evaluate_reactions_single for every zone: Y(num_species, num_zones), outputs (num_reactions, num_zones)
        subroutine gff_evaluate_reactions_multi(self, Y, T, rho, ierr, molar_flow, energy_rate, nu_e_loss_rate, nu_rate)
            class(GridFire), intent(in) :: self
            real(c_double), dimension(:,:), intent(in) :: Y
            real(c_double), dimension(:), intent(in) :: T, rho
            integer, intent(out) :: ierr
            real(c_double), dimension(:,:), intent(out), optional, target, contiguous :: molar_flow, energy_rate
            real(c_double), dimension(:,:), intent(out), optional, target, contiguous :: nu_e_loss_rate, nu_rate
            type(c_ptr) :: p_flow, p_energy, p_loss, p_rate

            p_flow = c_null_ptr
            p_energy = c_null_ptr
            p_loss = c_null_ptr
            p_rate = c_null_ptr
            if (present(molar_flow)) p_flow = c_loc(molar_flow)
            if (present(energy_rate)) p_energy = c_loc(energy_rate)
            if (present(nu_e_loss_rate)) p_loss = c_loc(nu_e_loss_rate)
            if (present(nu_rate)) p_rate = c_loc(nu_rate)
            ierr = int(gf_evaluate_reactions_c(self%ctx_type, self%ctx, Y, self%num_species, T, rho, &
                    p_flow, p_energy, p_loss, p_rate))
        end subroutine gff_evaluate_reactions_multi

        function gff_get_num_network_species(self) result(num_species)
            class(GridFire), intent(in) :: self
            integer(c_size_t) :: num_species
            integer(c_int) :: ierr

            ierr = gf_get_num_network_species(self%ctx, num_species)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Network Error: ", self%gff_get_last_error()
                error stop
            end if
        end function gff_get_num_network_species

        ! Names, A, Z and atomic masses [u] of the network species, in network order. Each array has
        ! gff_get_num_network_species() entries; A, Z and mass are optional.
        subroutine gff_get_network_species(self, names, A, Z, mass)
            class(GridFire), intent(in) :: self
            character(len=*), dimension(:), intent(out) :: names
            integer(c_int), dimension(:), intent(out), optional :: A, Z
            real(c_double), dimension(:), intent(out), optional :: mass
            type(c_ptr) :: c_name
            integer(c_int) :: ierr, A_i, Z_i
            real(c_double) :: mass_i
            integer :: i

            do i = 1, size(names)
                ierr = gf_get_network_species(self%ctx, int(i - 1, c_size_t), c_name, A_i, Z_i, mass_i)
                if (ierr /= GF_SUCCESS) then
                    print *, "GridFire Network Error: ", self%gff_get_last_error()
                    error stop
                end if
                call copy_c_string(c_name, names(i))
                if (present(A)) A(i) = A_i
                if (present(Z)) Z(i) = Z_i
                if (present(mass)) mass(i) = mass_i
            end do
        end subroutine gff_get_network_species

        ! Register the network's species (network order) so gff_evolve / gff_evaluate carry the full network.
        subroutine gff_register_network_species(self)
            class(GridFire), intent(inout) :: self
            integer(c_int) :: ierr

            ierr = gf_register_network_species(self%ctx)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Network Error: ", self%gff_get_last_error()
                error stop
            end if
            self%num_species = self%gff_get_num_network_species()
        end subroutine gff_register_network_species

        function gff_get_num_reactions(self) result(num_reactions)
            class(GridFire), intent(in) :: self
            integer(c_size_t) :: num_reactions
            integer(c_int) :: ierr

            ierr = gf_get_num_reactions(self%ctx, num_reactions)
            if (ierr /= GF_SUCCESS) then
                print *, "GridFire Network Error: ", self%gff_get_last_error()
                error stop
            end if
        end function gff_get_num_reactions

        ! Reaction identifiers (e.g. "p(p,e+)d"); ids has gff_get_num_reactions() entries.
        subroutine gff_get_reaction_ids(self, ids)
            class(GridFire), intent(in) :: self
            character(len=*), dimension(:), intent(out) :: ids
            type(c_ptr) :: c_id
            integer(c_int) :: ierr
            integer :: i

            do i = 1, size(ids)
                ierr = gf_get_reaction_id(self%ctx, int(i - 1, c_size_t), c_id)
                if (ierr /= GF_SUCCESS) then
                    print *, "GridFire Network Error: ", self%gff_get_last_error()
                    error stop
                end if
                call copy_c_string(c_id, ids(i))
            end do
        end subroutine gff_get_reaction_ids

        function gff_get_num_energy_channels(self) result(num_channels)
            class(GridFire), intent(in) :: self
            integer(c_size_t) :: num_channels
            integer(c_int) :: ierr

            ierr = gf_get_num_energy_channels(self%ctx, num_channels)
            if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                print *, "GridFire Energy Channel Error: ", self%gff_get_last_error()
                error stop
            end if
        end function gff_get_num_energy_channels

        subroutine gff_get_energy_channel_names(self, names)
            class(GridFire), intent(in) :: self
            character(len=*), dimension(:), intent(out) :: names ! size gff_get_num_energy_channels()
            type(c_ptr) :: c_name
            character(kind=c_char), pointer :: char_ptr(:)
            integer(c_int) :: ierr
            integer :: i, k

            do i = 1, size(names)
                ierr = gf_get_energy_channel_name(self%ctx, int(i - 1, c_size_t), c_name)
                if (ierr /= GF_SUCCESS .AND. ierr /= FDSSE_SUCCESS) then
                    print *, "GridFire Energy Channel Error: ", self%gff_get_last_error()
                    error stop
                end if

                call c_f_pointer(c_name, char_ptr, [len(names(i)) + 1])
                names(i) = ' '
                do k = 1, len(names(i))
                    if (char_ptr(k) == c_null_char) exit
                    names(i)(k:k) = char_ptr(k)
                end do
            end do
        end subroutine gff_get_energy_channel_names

end module gridfire_mod
