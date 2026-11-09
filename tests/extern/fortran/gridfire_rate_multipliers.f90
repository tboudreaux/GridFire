! Rate multipliers from Fortran (DSEP requirement 5, the Monte Carlo pipeline): a factor scales exactly its reaction,
! gff_get_rate_multiplier reads it back, an unknown reaction comes back as GF_REACTION_ERROR through ierr (so a caller
! can check its list at start-up and stop cleanly), provenance lists the factors, and clearing restores the rates.
program gridfire_rate_multipliers_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(2) :: seed_names = ["H-1  ", "He-4 "]
    real(c_double), dimension(2) :: seed_Y = [0.7d0, 0.07d0]
    character(len=5), dimension(8) :: solar_names = [ &
            "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
    real(c_double), dimension(8) :: solar_Y = [ &
            0.7026166d0, 9.7479d-06, 0.0689551d0, 0.00025d0, 7.8554d-05, 0.00060144d0, 8.103d-05, 2.1513d-05 ]
    type(GridFire) :: net
    character(len=8), allocatable :: names(:)
    character(len=32), allocatable :: ids(:)
    real(c_double), allocatable :: Y(:), flow(:), flow_scaled(:), flow_cleared(:)
    real(c_double) :: T, rho, factor
    character(len=:), allocatable :: provenance
    integer(c_size_t) :: n, nr
    integer :: i, s, ierr, i_pp

    T = 1.57d7
    rho = 150.0d0

    call net%gff_init(SINGLE_ZONE)
    call net%gff_register_species(seed_names)
    call net%gff_setup_policy("H_HE_BURNING_POLICY", seed_Y)
    call net%gff_register_network_species()
    call net%gff_setup_solver()

    n = net%gff_get_num_network_species()
    nr = net%gff_get_num_reactions()
    allocate(names(n), Y(n), ids(nr), flow(nr), flow_scaled(nr), flow_cleared(nr))
    call net%gff_get_network_species(names)
    call net%gff_get_reaction_ids(ids)
    Y = 0.0d0
    do i = 1, int(n)
        do s = 1, size(solar_names)
            if (trim(names(i)) == trim(solar_names(s))) Y(i) = solar_Y(s)
        end do
    end do
    i_pp = 0
    do i = 1, int(nr)
        if (trim(ids(i)) == "p(p,e+)d") i_pp = i
    end do
    call check(i_pp > 0, "p(p,e+)d is in the network")

    call net%gff_evaluate_reactions_single(Y, T, rho, ierr, molar_flow=flow)
    call check(ierr == GF_SUCCESS, "gff_evaluate_reactions")

    ! Start-up check of a perturbation list: an unknown reaction is reported, not fatal, when ierr is passed
    call net%gff_set_rate_multiplier("x(y,z)w", 2.0d0, ierr)
    call check(ierr == GF_REACTION_ERROR, "an unknown reaction gives GF_REACTION_ERROR")
    call net%gff_set_rate_multiplier("p(p,e+)d", -1.0d0, ierr)
    call check(ierr == GF_CONFIG_ERROR, "a negative factor gives GF_CONFIG_ERROR")

    call net%gff_set_rate_multiplier("p(p,e+)d", 1.3d0)
    call net%gff_get_rate_multiplier("p(p,e+)d", factor)
    call check(factor == 1.3d0, "gff_get_rate_multiplier reads the factor back")
    call net%gff_evaluate_reactions_single(Y, T, rho, ierr, molar_flow=flow_scaled)
    call check(ierr == GF_SUCCESS, "gff_evaluate_reactions with a factor")
    call check(abs(flow_scaled(i_pp) / (1.3d0 * flow(i_pp)) - 1d0) < 1d-13, "the pp flow is scaled by 1.3")
    call check(count(flow_scaled /= flow) == 1, "only the multiplied reaction changes")

    provenance = net%gff_get_provenance()
    call check(index(provenance, "p(p,e+)d=1.3") > 0, "provenance lists the factor")

    call net%gff_clear_rate_multipliers()
    call net%gff_evaluate_reactions_single(Y, T, rho, ierr, molar_flow=flow_cleared)
    call check(all(flow_cleared == flow), "clearing restores every flow exactly")

    print '(a, f18.15)', " pp flow ratio with factor 1.3:", flow_scaled(i_pp) / flow(i_pp)
    call net%gff_free()
    print *, "All Fortran rate multiplier checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

end program gridfire_rate_multipliers_test
