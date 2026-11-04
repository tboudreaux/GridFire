! Per-reaction energy generation of a MESA net at one state, for diagnosing GridFire/MESA differences.
!
!   mesa_rates <network.net> <composition file> <T [K]> <rho [g/cm^3]> [screening: none|extended]
!
! The composition file has one "<MESA iso name> <mass fraction>" per line. Output (CSV on stdout):
!   total,<eps_nuc>,<eps_neu_total>
!   reaction,<name>,<rate_raw>,<Q [MeV]>,<neuQ [MeV]>,<eps_nuc with only this reaction on>,<eps_neu with only it on>,<dX(h1)/dt with only it on>
! Per-reaction values come from calling net_get with every other rate_factor set to zero, so they need no knowledge of
! MESA's rate normalisation conventions. eta = 0, as in the bbq runs (EOS skipped).
program mesa_rates
   use const_lib, only: const_init
   use math_lib
   use chem_def
   use chem_lib
   use rates_def
   use rates_lib, only: rates_init, rates_warning_init, screening_option
   use net_def
   use net_lib
   use utils_lib, only: mesa_error
   implicit none

   character(len=512) :: net_file, comp_file, arg, my_mesa_dir
   character(len=32) :: iso
   integer :: ierr, handle, species, num_reactions, i, j, cid, u, ios, screening_mode
   integer, pointer :: chem_id(:), net_iso(:)
   real(dp) :: T, rho, x_in
   real(dp), allocatable, target :: rate_factors(:)
   real(dp), allocatable :: xa(:), dxdt(:), d_dxdt_dRho(:), d_dxdt_dT(:), d_dxdt_dx(:,:), &
      d_eps_nuc_dx(:), dabar_dx(:), dzbar_dx(:), dmc_dx(:)
   real(dp) :: eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, eps_neu_total, eps_nuc_categories(num_categories), &
      xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, xsum, eta, d_eta_dlnT, d_eta_dlnRho, weak_rate_factor
   type(Net_Info) :: n
   type(Net_General_Info), pointer :: g

   call get_command_argument(1, net_file)
   call get_command_argument(2, comp_file)
   call get_command_argument(3, arg); read(arg, *) T
   call get_command_argument(4, arg); read(arg, *) rho
   arg = ''
   if (command_argument_count() >= 5) call get_command_argument(5, arg)
   if (trim(arg) == 'none') arg = ''
   screening_mode = screening_option(arg, ierr)  ! same mapping bbq uses ('' = no screening)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'unknown screening option')
   call get_environment_variable('MESA_DIR', my_mesa_dir)

   ierr = 0
   call math_init()
   call const_init(my_mesa_dir, ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call chem_init('isotopes.data', ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call rates_init('reactions.list', '', 'rate_tables', .false., .false., '', '', '', ierr)
   if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call rates_warning_init(.false., 10d0)
   call net_init(ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   handle = alloc_net_handle(ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call net_start_def(handle, ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call read_net_file(net_file, handle, ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call net_finish_def(handle, ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   call net_setup_tables(handle, '', ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)
   species = net_num_isos(handle, ierr)
   num_reactions = net_num_reactions(handle, ierr)
   call get_chem_id_table_ptr(handle, chem_id, ierr)
   call get_net_iso_table_ptr(handle, net_iso, ierr)
   call get_net_ptr(handle, g, ierr); if (ierr /= 0) call mesa_error(__FILE__, __LINE__)

   allocate(xa(species), rate_factors(num_reactions), dxdt(species), d_dxdt_dRho(species), d_dxdt_dT(species), &
      d_dxdt_dx(species, species), d_eps_nuc_dx(species), dabar_dx(species), dzbar_dx(species), dmc_dx(species))

   xa = 0
   open(newunit=u, file=comp_file, status='old')
   do
      read(u, *, iostat=ios) iso, x_in
      if (ios /= 0) exit
      cid = get_nuclide_index(iso)
      if (cid <= 0) then
         write(0, *) 'unknown iso ', trim(iso); stop 1
      end if
      if (net_iso(cid) <= 0) then
         write(0, *) 'iso not in net ', trim(iso); stop 1
      end if
      xa(net_iso(cid)) = x_in
   end do
   close(u)

   call composition_info(species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, xsum, &
      dabar_dx, dzbar_dx, dmc_dx)
   eta = 0; d_eta_dlnT = 0; d_eta_dlnRho = 0; weak_rate_factor = 1

   rate_factors = 1
   call eval()
   write(*, '(a,2(",",es24.16))') 'total', eps_nuc, eps_neu_total
   if (command_argument_count() >= 6) then
      do j = 1, num_reactions
         i = g% reaction_id(j)
         if (trim(reaction_Name(i)) == 'r_h1_h1_wk_h2' .or. trim(reaction_Name(i)) == 'r_h2_pg_he3') then
            write(0, *) trim(reaction_Name(i)), ' raw ', n% rate_raw(j), ' screened ', n% rate_screened(j)
         end if
      end do
      write(0, *) 'y(h1) ', n% y(net_iso(ih1)), ' x(h1) ', xa(net_iso(ih1)), ' dxdt(h1) ', dxdt(net_iso(ih1)), ' ye ', ye, ' abar ', abar
      write(0, *) 'logTcut ', g% logTcut_lo, g% logTcut_lim, ' T ', T, ' rho ', rho
   end if

   do j = 1, num_reactions
      rate_factors = 0
      rate_factors(j) = 1
      call eval()
      i = g% reaction_id(j)
      write(*, '(a,",",a,5(",",es24.16),",",es24.16)') 'reaction', trim(reaction_Name(i)), n% rate_raw(j), &
         std_reaction_Qs(i), std_reaction_neuQs(i), eps_nuc, eps_neu_total, dxdt(net_iso(ih1))
   end do

contains

   subroutine eval()
      call net_get(handle, .false., n, species, num_reactions, xa, T, log10(T), rho, log10(rho), &
         abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, rate_factors, weak_rate_factor, &
         std_reaction_Qs, std_reaction_neuQs, eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
         dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, screening_mode, eps_nuc_categories, eps_neu_total, ierr)
      if (ierr /= 0) call mesa_error(__FILE__, __LINE__, 'net_get')
   end subroutine eval

end program mesa_rates
