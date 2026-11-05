! Per-reaction contributions, passive species and the mass-fraction helpers from Fortran (DSEP requirements 7, 8):
! reaction energies sum to eps, pep is a separate source, a passive Fe-56 survives gff_evolve unchanged, and mass
! fractions converted to molar abundances and back round-trip and sum to one before and after burning.
program gridfire_dsep_composition_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(9) :: seed_names = [ &
            "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24", "Fe-56" ]
    real(c_double), dimension(9) :: seed_Y = [ &
            0.702616602672027d0, 9.74791583949078d-06, 0.06895512307276903d0, 0.00025d0, &
            7.855418029399437d-05, 0.0006014411598306529d0, 8.103062886768109d-05, 2.151340851063217d-05, 2.0d-5 ]
    type(GridFire) :: net
    character(len=8), allocatable :: names(:)
    character(len=32), allocatable :: ids(:)
    integer(c_int), allocatable :: A(:), Z(:)
    real(c_double), allocatable :: mass(:), X(:), X_out(:), Y(:), Y_out(:), X_back(:), flow(:), energy_rate(:), nu_loss(:), nu_rate(:), channels(:)
    real(c_double) :: T, rho, eps, dedt, dedrho, nu_total
    real(c_double) :: energy, dE_dT, dE_drho, nu_e_loss, nu_flux, mass_lost, eps_nuc, dlnT, dlnRho, nu_e_rate
    integer(c_size_t) :: n
    integer :: i, s, ierr, i_pp, i_pep, i_fe, i_h

    T = 1.57d7
    rho = 150.0d0

    call net%gff_init(SINGLE_ZONE)
    call net%gff_register_species(seed_names)
    call net%gff_set_passive_species(["Fe-56"])
    call net%gff_setup_policy("MAIN_SEQUENCE_POLICY", seed_Y)
    call net%gff_setup_solver("CVODE")

    ! Full network; the caller keeps mass fractions and converts at the boundary
    n = net%gff_get_num_network_species()
    allocate(names(n), A(n), Z(n), mass(n), X(n), X_out(n), Y(n), Y_out(n), X_back(n))
    call net%gff_get_network_species(names, A, Z, mass)
    call net%gff_register_network_species()
    X = 0.0d0
    i_fe = 0
    i_h = 0
    do i = 1, int(n)
        do s = 1, size(seed_names)
            if (trim(names(i)) == trim(seed_names(s))) X(i) = seed_Y(s) * mass(i)
        end do
        if (trim(names(i)) == "Fe-56") i_fe = i
        if (trim(names(i)) == "H-1") i_h = i
    end do
    call check(i_fe > 0 .and. i_h > 0, "Fe-56 (passive) and H-1 are network species")
    X = X / sum(X)
    call net%gff_mass_fractions_to_molar(X, Y, ierr)
    call check(ierr == GF_SUCCESS, "gff_mass_fractions_to_molar succeeds")
    call net%gff_molar_to_mass_fractions(Y, X_back, ierr)
    call check(ierr == GF_SUCCESS, "gff_molar_to_mass_fractions succeeds")
    call check(maxval(abs(X_back - X)) < 1d-15, "mass fractions round-trip")
    call net%gff_mass_fractions_to_molar(X * 1.01d0, Y_out, ierr)
    call check(ierr == GF_INVALID_STATE, "mass fractions summing to 1.01 are rejected")

    ! --- Per-reaction contributions ---
    allocate(ids(net%gff_get_num_reactions()))
    call net%gff_get_reaction_ids(ids)
    allocate(flow(size(ids)), energy_rate(size(ids)), nu_loss(size(ids)), nu_rate(size(ids)))
    call net%gff_evaluate_reactions(Y, T, rho, ierr, molar_flow=flow, energy_rate=energy_rate, &
            nu_e_loss_rate=nu_loss, nu_rate=nu_rate)
    call check(ierr == GF_SUCCESS, "gff_evaluate_reactions succeeds")
    call net%gff_evaluate(Y, T, rho, eps, dedt, dedrho, ierr, nu_e_loss_rate=nu_total)
    call check(ierr == GF_SUCCESS, "gff_evaluate succeeds")
    print '(a, 2es16.7)', " eps and sum of reaction energies:", eps, sum(energy_rate)
    call check(abs(sum(energy_rate) - eps) <= 1.0d-10 * abs(eps), "reaction energies sum to eps")
    call check(abs(sum(nu_loss) - nu_total) <= 1.0d-10 * abs(nu_total), "reaction neutrino losses sum to the total")
    i_pp = findloc(ids, "p(p,e+)d", dim=1)
    i_pep = findloc(ids, "p(e-p,)d", dim=1)
    call check(i_pp > 0 .and. i_pep > 0, "pp and pep are separate reactions")
    print '(a, es12.4)', " pep/pp neutrino rate:", nu_rate(i_pep) / nu_rate(i_pp)
    call check(nu_rate(i_pep) / nu_rate(i_pp) > 1d-3 .and. nu_rate(i_pep) / nu_rate(i_pp) < 1d-2, "pep/pp is a few 0.1%")

    ! --- Evolve one Myr: passive iron unchanged, hydrogen burns, mass fractions still sum to one ---
    allocate(channels(net%gff_get_num_energy_channels()))
    call net%gff_evolve(Y, T, rho, 3.15d13, 1d-12, Y_out, energy, dE_dT, dE_drho, nu_e_loss, nu_flux, mass_lost, &
            eps_nuc, dlnT, dlnRho, nu_e_rate, channels, ierr)
    call check(ierr == GF_SUCCESS, "gff_evolve succeeds")
    call check(Y_out(i_fe) == Y(i_fe), "passive Fe-56 is unchanged")
    call net%gff_molar_to_mass_fractions(Y_out, X_out, ierr)
    call check(ierr == GF_SUCCESS, "gff_molar_to_mass_fractions after burning succeeds")
    print '(a, 2f20.15, a, 2f12.8)', " sum X in, out:", sum(X), sum(X_out), "; X(H-1) in, out:", X(i_h), X_out(i_h)
    call check(abs(sum(X_out) - 1d0) < 1d-14, "mass fractions after burning sum to one")
    call check(X_out(i_h) < X(i_h), "hydrogen burns")

    call net%gff_free()
    print *, "All Fortran per-reaction, passive species and mass fraction checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

end program gridfire_dsep_composition_test
