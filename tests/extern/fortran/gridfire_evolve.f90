program main
    use iso_c_binding
    use gridfire_mod
    implicit none

    type(GridFire) :: net
    integer(c_int) :: ierr
    integer :: i

    ! --- 1. Define Species and Initial Conditions ---
    ! Note: String lengths must match or exceed the longest name.
    ! We pad with spaces, which 'trim' handles inside the module.
    character(len=5), dimension(8) :: species_names = [ &
            "H-1  ", &
            "He-3 ", &
            "He-4 ", &
            "C-12 ", &
            "N-14 ", &
            "O-16 ", &
            "Ne-20", &
            "Mg-24"  &
        ]

    ! Initial Mass Fractions (converted to Molar Abundances Y = X/A)
    ! Standard solar-ish composition
    real(c_double), dimension(8) :: Y_in = [ &
        0.702616602672027,     &
        9.74791583949078e-06,  &
        0.06895512307276903,   &
        0.00025,               &
        7.855418029399437e-05, &
        0.0006014411598306529, &
        8.103062886768109e-05, &
        2.151340851063217e-05  &
        ]

    ! Output buffers
    real(c_double), dimension(8) :: Y_out
    real(c_double) :: energy_out, dedt, dedrho, snu_e_loss, snu_flux, dmass

    ! Thermodynamic Conditions (Solar Core-ish)
    real(c_double) :: T = 1.5e7      ! 15 Million K
    real(c_double) :: rho = 150.0e0  ! 150 g/cm^3
    real(c_double) :: dt = 3.0e17     ! 1 second timestep

    ! --- 2. Initialize GridFire ---
    print *, "Initializing GridFire..."
    call net%gff_init()

    ! --- 3. Register Species ---
    print *, "Registering species..."
    call net%register_species(species_names)

    ! --- 4. Configure Engine & Solver ---
    print *, "Setting up Main Sequence Policy..."
    call net%setup_policy("MAIN_SEQUENCE_POLICY", Y_in)

    print *, "Setting up CVODE Solver..."
    call net%setup_solver("CVODE")

    ! --- 5. Evolve ---
    print *, "Evolving system (dt =", dt, "s)..."
    call net%evolve(Y_in, T, rho, dt, Y_out, energy_out, dedt, dedrho, snu_e_loss, snu_flux, dmass, ierr)

    if (ierr /= 0) then
        print *, "Evolution Failed with error code: ", ierr
        print *, "Error Message: ", net%get_last_error()
        call net%gff_free() ! Always cleanup
        stop
    end if

    ! --- 6. Report Results ---
    print *, ""
    print *, "--- Results ---"
    print '(A, ES12.5, A)', "Energy Generation: ", energy_out, " erg/g/s"
    print '(A, ES12.5)',    "dEps/dT:           ", dedt
    print '(A, ES12.5)',    "dEps/drho:         ", dedrho
    print '(A, ES12.5)',    "Neutrino Energy Loss:      ", snu_e_loss
    print '(A, ES12.5)',    "Neutrino Flux:          ", snu_flux
    print '(A, ES12.5)',    "Mass Change:       ", dmass

    print *, ""
    print *, "Abundances:"
    do i = 1, size(species_names)
        print '(A, " : ", ES12.5, " -> ", ES12.5)', &
                trim(species_names(i)), Y_in(i), Y_out(i)
    end do

    ! --- 7. Cleanup ---
    call net%gff_free()

end program main