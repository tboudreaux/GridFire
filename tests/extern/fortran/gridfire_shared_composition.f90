! Shared-composition groups from Fortran: a mask-built group shares one composition but not its energy.
program gridfire_shared_composition_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(8) :: names = [ &
            "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
    real(c_double), dimension(8) :: Y0 = [ &
            0.702616602672027d0, 9.74791583949078d-06, 0.06895512307276903d0, 0.00025d0, &
            7.855418029399437d-05, 0.0006014411598306529d0, 8.103062886768109d-05, 2.151340851063217d-05 ]
    integer, parameter :: nz = 3
    type(GridFire) :: net
    real(c_double) :: Y_in(8, nz), Y_out(8, nz), T(nz), rho(nz), mass(nz)
    real(c_double), dimension(nz) :: energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, eps, dlnt, dlnrho, nu_rate
    real(c_double), allocatable :: channels(:,:)
    integer(c_int) :: group(nz)
    integer :: i, ierr

    call net%gff_init(MULTI_ZONE, int(nz, c_size_t))
    call net%gff_register_species(names)
    call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y0)
    call net%gff_setup_solver("CVODE")
    allocate(channels(net%gff_get_num_energy_channels(), nz))

    do i = 1, nz
        Y_in(:, i) = Y0
    end do
    T = [2.5d7, 1.2d7, 1.5d7]
    rho = [120.0d0, 60.0d0, 150.0d0]
    mass = [1.0d0, 3.0d0, 1.0d0]

    group = gff_shared_composition_groups_from_mask([1_c_int, 1_c_int, 0_c_int])
    call check(all(group == [0, 0, -1]), "groups from mask")
    call net%gff_set_shared_composition_groups(group, mass)
    call net%gff_evolve(Y_in, T, rho, 1.0d15, 1.0d-12, Y_out, energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, &
            eps, dlnt, dlnrho, nu_rate, channels, ierr)
    call check(ierr == GF_SUCCESS, "grouped evolve succeeds")
    call check(all(Y_out(:, 1) == Y_out(:, 2)), "group members share one composition")
    call check(energy(1) > 10.0d0 * energy(2), "each member keeps its own energy")
    print '(a, 3es14.5)', " energies:", energy

    call net%gff_clear_shared_composition_groups()
    call net%gff_free()
    print *, "All Fortran shared-composition checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

end program gridfire_shared_composition_test
