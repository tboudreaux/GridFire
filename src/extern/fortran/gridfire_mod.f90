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

        enumerator :: GF_SINGLE_ZONE = 1001
        enumerator :: GF_MULTI_ZONE = 1002
    end enum

    ! Values for GFConfig%engine_screening (enum GF_SCREENING_TYPE)
    enum, bind(C)
        enumerator :: GF_SCREENING_BARE = 0
        enumerator :: GF_SCREENING_WEAK = 1
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
        ! solver.pointSolver
        integer(c_int) :: solver_clamp_negative_iterates
        integer(c_int) :: solver_max_step_failure_restarts
        ! engine
        integer(c_int) :: engine_screening        ! GF_SCREENING_BARE or GF_SCREENING_WEAK
        integer(c_int) :: engine_jacobian         ! GF_JACOBIAN_ANALYTIC or GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION
        integer(c_int) :: engine_split_hydrogen_pp_cno
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

        ! int gf_construct_solver_from_engine(void* ptr, const char* solver_name);
        function gf_construct_solver_from_engine(ptr, solver_name) result(ierr) &
                bind(C, name="gf_construct_solver_from_engine")
            import
            type(c_ptr), value :: ptr
            character(kind=c_char), dimension(*), intent(in) :: solver_name
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
        procedure :: gff_setup_policy
        procedure :: gff_setup_solver
        procedure :: gff_get_last_error

        procedure :: gff_evolve_single
        procedure :: gff_evolve_multi

        generic :: gff_evolve => gff_evolve_single, gff_evolve_multi

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

        subroutine gff_setup_solver(self, solver_name)
            class(GridFire), intent(in) :: self
            character(len=*), intent(in) :: solver_name
            integer(c_int) :: ierr

            ierr = gf_construct_solver_from_engine(self%ctx, trim(solver_name) // c_null_char)
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
