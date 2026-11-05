! Run-time configuration through the Fortran module: a GFConfig built in Fortran, the TOML round trip, and that the
! configuration reaches the engine (bare screening lowers eps_nuc; the Jacobian method does not change the answer).
!
! Usage: gf_fortran_config_test <scratch directory>. Stops with a non-zero code on the first failed check.
program gridfire_config_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    type(GFConfig) :: defaults, custom, readback
    type(GridFire) :: net
    character(len=1024) :: scratch_dir, toml_path
    real(c_double) :: eps_default, eps_bare, eps_autodiff

    call get_command_argument(1, scratch_dir)
    if (len_trim(scratch_dir) == 0) then
        print *, "usage: gf_fortran_config_test <scratch directory>"
        error stop 2
    end if

    ! --- Defaults ---
    defaults = gff_default_config()
    call check(defaults%engine_screening == GF_SCREENING_WEAK, "default screening is weak")
    call check(defaults%engine_jacobian == GF_JACOBIAN_ANALYTIC, "default Jacobian is analytic")

    ! --- A configuration built in Fortran round-trips through the context and through TOML ---
    custom = defaults
    custom%engine_screening = GF_SCREENING_BARE
    custom%engine_jacobian = GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION
    custom%solver_rel_tol = 1.0d-7
    custom%solver_max_step_failure_restarts = 5

    call net%gff_init(SINGLE_ZONE)
    call net%gff_set_config(custom)
    readback = net%gff_get_config()
    call check(same_config(readback, custom), "gff_get_config returns what gff_set_config set")

    toml_path = trim(scratch_dir) // "/gf_fortran_config_test.toml"
    call net%gff_save_config(toml_path)
    call net%gff_free()

    call net%gff_init(SINGLE_ZONE)
    call net%gff_load_config(toml_path)
    readback = net%gff_get_config()
    call check(same_config(readback, custom), "the configuration read back from TOML matches the one saved")
    call net%gff_free()

    ! --- The configuration reaches the engine ---
    eps_default = evolved_eps_nuc(defaults)
    custom = defaults
    custom%engine_screening = GF_SCREENING_BARE
    eps_bare = evolved_eps_nuc(custom)
    custom = defaults
    custom%engine_jacobian = GF_JACOBIAN_AUTOMATIC_DIFFERENTIATION
    eps_autodiff = evolved_eps_nuc(custom)
    print '(a, 3es18.10)', " eps_nuc (default, bare screening, CppAD Jacobian):", eps_default, eps_bare, eps_autodiff
    call check(eps_default / eps_bare > 1.02d0, "weak screening raises eps_nuc")
    call check(abs(eps_autodiff / eps_default - 1.0d0) < 1.0d-3, "both Jacobian methods reach the same solution")

    print *, "All Fortran configuration checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

    logical function same_config(a, b)
        type(GFConfig), intent(in) :: a, b
        same_config = a%solver_abs_tol == b%solver_abs_tol .and. a%solver_rel_tol == b%solver_rel_tol .and. &
                a%trigger_off_diagonal_threshold == b%trigger_off_diagonal_threshold .and. &
                a%trigger_timestep_collapse_ratio == b%trigger_timestep_collapse_ratio .and. &
                a%trigger_max_convergence_failures == b%trigger_max_convergence_failures .and. &
                a%trigger_boundary_flux_relative_threshold == b%trigger_boundary_flux_relative_threshold .and. &
                a%trigger_boundary_flux_absolute_threshold == b%trigger_boundary_flux_absolute_threshold .and. &
                a%engine_adaptive_relative_culling_threshold == b%engine_adaptive_relative_culling_threshold .and. &
                a%solver_clamp_negative_iterates == b%solver_clamp_negative_iterates .and. &
                a%solver_max_step_failure_restarts == b%solver_max_step_failure_restarts .and. &
                a%engine_screening == b%engine_screening .and. a%engine_jacobian == b%engine_jacobian .and. &
                a%engine_split_hydrogen_pp_cno == b%engine_split_hydrogen_pp_cno
    end function same_config

    ! Evolve 1e10 s at solar-core conditions with the given configuration and return eps_nuc at the end.
    real(c_double) function evolved_eps_nuc(config)
        type(GFConfig), intent(in) :: config
        type(GridFire) :: run
        character(len=5), dimension(8) :: species_names = [ &
                "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
        real(c_double), dimension(8) :: Y_in = [ &
                0.702616602672027d0, 9.74791583949078d-06, 0.06895512307276903d0, 0.00025d0, &
                7.855418029399437d-05, 0.0006014411598306529d0, 8.103062886768109d-05, 2.151340851063217d-05 ]
        real(c_double), dimension(8) :: Y_out
        real(c_double) :: energy, dedt, dedrho, nu_loss, nu_flux, mass_lost
        real(c_double) :: eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_loss_rate
        real(c_double), allocatable :: eps_by_channel(:)
        integer :: ierr

        call run%gff_init(SINGLE_ZONE)
        call run%gff_set_config(config)
        call run%gff_register_species(species_names)
        call run%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y_in)
        call run%gff_setup_solver("CVODE")
        allocate(eps_by_channel(run%gff_get_num_energy_channels()))
        call run%gff_evolve(Y_in, 1.5d7, 150.0d0, 1.0d10, 1.0d-12, Y_out, energy, dedt, dedrho, nu_loss, nu_flux, &
                mass_lost, eps_nuc, dlneps_dlnt, dlneps_dlnrho, nu_loss_rate, eps_by_channel, ierr)
        call check(ierr == GF_SUCCESS, "gff_evolve succeeds")
        call run%gff_free()
        evolved_eps_nuc = eps_nuc
    end function evolved_eps_nuc

end program gridfire_config_test
