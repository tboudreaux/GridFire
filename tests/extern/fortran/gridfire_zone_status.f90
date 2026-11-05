! Per-zone failure reporting from Fortran: one failing zone (T = 1e10 K, no restarts) among good ones.
program gridfire_zone_status_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(8) :: names = [ &
            "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
    real(c_double), dimension(8) :: Y0 = [ &
            0.702616602672027d0, 9.74791583949078d-06, 0.06895512307276903d0, 0.00025d0, &
            7.855418029399437d-05, 0.0006014411598306529d0, 8.103062886768109d-05, 2.151340851063217d-05 ]
    integer, parameter :: nz = 3
    real(c_double), parameter :: tMax = 3.0d15
    type(GridFire) :: net
    type(GFConfig) :: config
    real(c_double) :: Y_in(8, nz), Y_out(8, nz), T(nz), rho(nz)
    real(c_double), dimension(nz) :: energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, eps, dlnt, dlnrho, nu_rate
    real(c_double), allocatable :: channels(:,:)
    integer(c_int) :: status(nz)
    real(c_double) :: reached(nz)
    integer :: i, ierr

    call net%gff_init(MULTI_ZONE, int(nz, c_size_t))
    config = gff_default_config()
    config%solver_max_step_failure_restarts = 0
    call net%gff_set_config(config)
    call net%gff_register_species(names)
    call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y0)
    call net%gff_setup_solver("CVODE")
    allocate(channels(net%gff_get_num_energy_channels(), nz))

    do i = 1, nz
        Y_in(:, i) = Y0
    end do
    T = [1.5d7, 1.0d10, 2.0d7]
    rho = [150.0d0, 1.0d4, 100.0d0]
    call net%gff_evolve(Y_in, T, rho, tMax, 1.0d-12, Y_out, energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, &
            eps, dlnt, dlnrho, nu_rate, channels, ierr)
    call check(ierr == GF_ZONE_FAILURE, "gff_evolve returns GF_ZONE_FAILURE")

    call net%gff_get_zone_status(status, reached)
    print '(a, 3i4)', " zone status:", status
    call check(status(1) == GF_SUCCESS .and. status(3) == GF_SUCCESS, "the good zones succeed")
    call check(status(2) == GF_CVODE_SOLVER_FAILURE_ERROR, "the hot zone reports a CVODE failure")
    call check(reached(1) == tMax .and. reached(3) == tMax .and. reached(2) < tMax, "time reached")
    call check(all(Y_out(:, 2) == Y_in(:, 2)) .and. energy(2) == 0.0d0, "a failed zone hands back its input")
    call check(energy(1) > 0.0d0 .and. energy(3) > 0.0d0, "the good zones burn")

    call net%gff_free()
    print *, "All Fortran zone status checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

end program gridfire_zone_status_test
