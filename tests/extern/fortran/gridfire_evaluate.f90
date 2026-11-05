! gff_evaluate and the network listing from Fortran: species metadata, full-network registration, finite-difference
! derivatives, and a multi-zone call agreeing with single-zone calls.
program gridfire_evaluate_test
    use iso_c_binding
    use gridfire_mod
    implicit none

    character(len=5), dimension(8) :: seed_names = [ &
            "H-1  ", "He-3 ", "He-4 ", "C-12 ", "N-14 ", "O-16 ", "Ne-20", "Mg-24" ]
    real(c_double), dimension(8) :: seed_Y = [ &
            0.702616602672027d0, 9.74791583949078d-06, 0.06895512307276903d0, 0.00025d0, &
            7.855418029399437d-05, 0.0006014411598306529d0, 8.103062886768109d-05, 2.151340851063217d-05 ]
    integer, parameter :: num_zones = 3
    type(GridFire) :: point, grid
    character(len=8), allocatable :: names(:)
    character(len=32), allocatable :: reaction_ids(:)
    integer(c_int), allocatable :: A(:), Z(:)
    real(c_double), allocatable :: mass(:), Y(:), Y_zones(:,:), dydt(:), channels(:)
    real(c_double) :: T, rho, eps, dedt, dedrho, eps_up, eps_down, d1, d2, nu_loss
    real(c_double) :: T_zones(num_zones), rho_zones(num_zones), eps_z(num_zones), dedt_z(num_zones), dedrho_z(num_zones)
    integer(c_size_t) :: n
    integer :: i, s, ierr

    ! --- Network listing and full-network registration ---
    call point%gff_init(SINGLE_ZONE)
    call point%gff_register_species(seed_names)
    call point%gff_setup_policy("MAIN_SEQUENCE_POLICY", seed_Y)
    n = point%gff_get_num_network_species()
    allocate(names(n), A(n), Z(n), mass(n), Y(n), dydt(n))
    call point%gff_get_network_species(names, A, Z, mass)
    call check(all(A >= 1) .and. all(Z >= 0) .and. all(abs(mass - A) < 0.1d0), "species metadata is sane")
    allocate(reaction_ids(point%gff_get_num_reactions()))
    call point%gff_get_reaction_ids(reaction_ids)
    call check(any(reaction_ids == "p(p,e+)d"), "p(p,e+)d is listed")

    call point%gff_register_network_species()
    call check(point%num_species == n, "registering the network sets num_species")
    Y = 0.0d0
    do i = 1, int(n)
        do s = 1, size(seed_names)
            if (trim(names(i)) == trim(seed_names(s))) Y(i) = seed_Y(s)
        end do
    end do

    ! --- Single zone: derivatives against finite differences ---
    T = 1.5d7
    rho = 150.0d0
    allocate(channels(point%gff_get_num_energy_channels()))
    call point%gff_evaluate(Y, T, rho, eps, dedt, dedrho, ierr, nu_e_loss_rate=nu_loss, eps_by_channel=channels, dydt=dydt)
    call check(ierr == GF_SUCCESS, "gff_evaluate succeeds")
    print '(a, 4es16.7)', " eps, dEps/dT, dEps/drho, nu loss:", eps, dedt, dedrho, nu_loss
    call check(abs(sum(channels) - eps) <= 1.0d-12 * abs(eps), "channels sum to eps")
    call check(dydt(1) < 0.0d0 .or. trim(names(1)) /= "H-1", "hydrogen is destroyed")

    call point%gff_evaluate(Y, T * (1 + 1d-4), rho, eps_up, d1, d2, ierr)
    call point%gff_evaluate(Y, T * (1 - 1d-4), rho, eps_down, d1, d2, ierr)
    call check(abs((eps_up - eps_down) / (2d-4 * T) / dedt - 1.0d0) < 1.0d-6, "dEps/dT matches finite differences")
    call point%gff_evaluate(Y, T, rho * (1 + 1d-4), eps_up, d1, d2, ierr)
    call point%gff_evaluate(Y, T, rho * (1 - 1d-4), eps_down, d1, d2, ierr)
    call check(abs((eps_up - eps_down) / (2d-4 * rho) / dedrho - 1.0d0) < 1.0d-6, "dEps/drho matches finite differences")

    call point%gff_evaluate(Y, -1.0d0, rho, eps, dedt, dedrho, ierr)
    call check(ierr == GF_INVALID_STATE, "a negative temperature is GF_INVALID_STATE")

    ! --- Multi-zone equals single-zone, zone by zone ---
    call grid%gff_init(MULTI_ZONE, int(num_zones, c_size_t))
    call grid%gff_register_species(seed_names)
    call grid%gff_setup_policy("MAIN_SEQUENCE_POLICY", seed_Y)
    call grid%gff_register_network_species()
    allocate(Y_zones(n, num_zones))
    do i = 1, num_zones
        Y_zones(:, i) = Y
    end do
    T_zones = [1.2d7, 1.5d7, 2.5d7]
    rho_zones = [50.0d0, 150.0d0, 80.0d0]
    call grid%gff_evaluate(Y_zones, T_zones, rho_zones, eps_z, dedt_z, dedrho_z, ierr)
    call check(ierr == GF_SUCCESS, "multi-zone gff_evaluate succeeds")
    do i = 1, num_zones
        call point%gff_evaluate(Y, T_zones(i), rho_zones(i), eps, dedt, dedrho, ierr)
        call check(eps == eps_z(i) .and. dedt == dedt_z(i) .and. dedrho == dedrho_z(i), "zone matches single-zone call")
    end do

    call grid%gff_free()
    call point%gff_free()
    print *, "All Fortran gff_evaluate checks passed"

contains

    subroutine check(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        if (.not. condition) then
            print *, "FAIL: ", description
            error stop 1
        end if
    end subroutine check

end program gridfire_evaluate_test
