program main_multi
    use iso_c_binding
    use gridfire_mod
    implicit none

    ! --- Constants ---
    integer, parameter :: NUM_SPECIES = 8
    integer, parameter :: ZONES = 100

    type(GridFire) :: net
    integer(c_int) :: ierr
    integer :: i, z

    ! --- 1. Define Species ---
    character(len=5), dimension(NUM_SPECIES) :: species_names = [ &
            "H-1  ", &
            "He-3 ", &
            "He-4 ", &
            "C-12 ", &
            "N-14 ", &
            "O-16 ", &
            "Ne-20", &
            "Mg-24"  &
    ]

    ! Initial Mass Fractions
    ! Standard solar-ish composition template
    real(c_double), dimension(NUM_SPECIES) :: abundance_root = [ &
            0.702616602672027d0,     &
                    9.74791583949078d-06,  &
                    0.06895512307276903d0,   &
                    0.00025d0,               &
                    7.855418029399437d-05, &
                    0.0006014411598306529d0, &
                    8.103062886768109d-05, &
                    2.151340851063217d-05  &
            ]

    ! --- Multi-Zone Arrays ---
    ! Fortran is Column-Major. To match C's Row-Major [ZONES][SPECIES],
    ! we dimension as (SPECIES, ZONES) so that Species are contiguous for each Zone.
    real(c_double), dimension(NUM_SPECIES, ZONES) :: Y_in, Y_out

    ! Thermodynamic Conditions Arrays
    real(c_double), dimension(ZONES) :: T_arr
    real(c_double), dimension(ZONES) :: rho_arr

    ! Output Arrays
    real(c_double), dimension(ZONES) :: energy_out
    real(c_double), dimension(ZONES) :: dedt
    real(c_double), dimension(ZONES) :: dedrho
    real(c_double), dimension(ZONES) :: snu_e_loss
    real(c_double), dimension(ZONES) :: snu_flux
    real(c_double), dimension(ZONES) :: dmass

    ! Time settings
    real(c_double) :: tMax = 3.0e17  ! 10 Gyr total time
    real(c_double) :: dt0 = 1e-12    ! Starting Timestep

    ! --- 2. Setup Data ---
    print *, "Testing GridFireEvolve Multi (Fortran)"
    print *, " Number of zones: ", ZONES
    print *, " Number of species: ", NUM_SPECIES

    do z = 1, ZONES
        ! Initialize Abundances (Copy root to every zone)
        Y_in(:, z) = abundance_root

        ! Initialize T and Rho gradient
        ! T: 1.0e7 -> 2.0e7 in steps of 1.0e5
        T_arr(z) = 1.0d7 + dble(z-1) * 1.0d5
        rho_arr(z) = 1.5d2

        if (z <= 3) then
            print '(A, I0, A, ES12.5, A, ES12.5, A)', &
                    " Zone ", z-1, " - Temp: ", T_arr(z), " K, Rho: ", rho_arr(z), " g/cm^3"
        end if
    end do

    ! --- 3. Initialize GridFire Multi-Zone ---
    print *, "Initializing GridFire..."
    ! Note: Pass integer(c_size_t) for zone count
    call net%gff_init(MULTI_ZONE, int(ZONES, c_size_t))

    ! --- 4. Register Species ---
    print *, "Registering species..."
    call net%gff_register_species(species_names)

    ! --- 5. Configure Engine & Solver ---
    print *, "Setting up Main Sequence Policy..."
    call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", abundance_root)

    print *, "Setting up CVODE Solver..."
    call net%gff_setup_solver("CVODE")

    ! --- 6. Evolve ---
    print *, "Evolving system..."

    ! Note: We pass the arrays T_arr and rho_arr.
    call net%gff_evolve(Y_in, T_arr, rho_arr, tMax, dt0, &
            Y_out, energy_out, dedt, dedrho, &
            snu_e_loss, snu_flux, dmass, ierr)

    if (ierr /= 0) then
        print *, "Evolution Failed with error code: ", ierr
        print *, "Error Message: ", net%gff_get_last_error()
        call net%gff_free()
        stop
    end if

    ! --- 7. Report Results ---
    print *, ""
    print *, "--- Results (First 3 Zones) ---"

    do z = 1, 3
        print *, "=== Zone ", z-1, " ==="
        print '(A, ES12.5, A)', " Energy Gen:  ", energy_out(z), " erg/g/s"
        print '(A, ES12.5)',    " Mass Lost:   ", dmass(z)
        print '(A, ES12.5)',    " T:           ", T_arr(z)

        print *, " Key Abundances (H-1, He-4, C-12):"
        print '(3(ES12.5, 1X))', Y_out(1, z), Y_out(3, z), Y_out(4, z)
        print *, ""
    end do

    print *, "... (Zones 4-", ZONES, " omitted) ..."

    ! --- 8. Cleanup ---
    call net%gff_free()

end program main_multi