! A downstream Fortran program built only from what an installed GridFire reports: the module is compiled from the
! installed gridfire_mod.f90 with this compiler, and linked against the C API library.
program gridfire_consumer
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(8) :: names = [ "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
    real(c_double), dimension(8) :: Y = [ 0.7026166d0, 9.7479d-06, 0.0689551d0, 0.00025d0, 7.8554d-05, 0.00060144d0, &
                                          8.103d-05, 2.1513d-05 ]
    real(c_double), dimension(8) :: Y_out
    real(c_double), allocatable :: channels(:)
    real(c_double) :: energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, eps_nuc, dlnT, dlnRho, nu_rate
    type(GridFire) :: net
    integer :: ierr

    call net%gff_init(SINGLE_ZONE)
    call net%gff_register_species(names)
    call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", Y)
    call net%gff_setup_solver("CVODE")
    allocate(channels(net%gff_get_num_energy_channels()))
    call net%gff_evolve(Y, 1.5d7, 150.0d0, 1d10, 1d-12, Y_out, energy, dedt, dedrho, nu_loss, nu_flux, mass_lost, &
            eps_nuc, dlnT, dlnRho, nu_rate, channels, ierr)
    call net%gff_free()
    if (ierr /= GF_SUCCESS .or. .not. (eps_nuc > 1d0 .and. eps_nuc < 1d3)) then
        print *, "FAIL: ierr", ierr, "eps_nuc", eps_nuc
        error stop 1
    end if
    print '(a, es13.6, a)', " Fortran consumer OK: eps_nuc = ", eps_nuc, " erg/g/s"
end program gridfire_consumer
