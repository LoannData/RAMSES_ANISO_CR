module feedback_module
  use amr_parameters,only:dp,ndim
  implicit none

  public


  logical::sn_feedback_sink = .false. !SN feedback emanates from the sink

  !mass, energy and momentum of supernova for forcing by sinks
  ! sn_e in erg, sn_p in g cm/s , sn_mass in Ms
  real(dp):: sn_e=1.d51 , sn_p=4.d43 , sn_mass=1.
  real(dp):: sn_e_ref , sn_p_ref , sn_mass_ref

  !limit speed and temperatue in supernova remnants, minimum radius for the SN remnant
  real(dp):: Tsat=1.d99 , Vsat=1.d99, sn_r_sat=0.

  !dispersion velocity of the stellar objects in km/s
  real(dp):: Vdisp=1. 


  ! Stellar object related arrays, those parameters are read in  read_stellar_params 
  logical:: stellar = .false.
  logical:: sn_direct = .false.
  logical:: make_stellar_glob = .false. !if used, the objects are created when the total mass in sinks exceeds stellar_msink_th
  integer:: nstellarmax ! maximum number of stellar objects
  integer:: nstellar = 0 ! current number of stellar objects
  integer:: nstellar_tot = 0 ! total number of stellar objects
  real(dp):: imf_index, imf_low, imf_high ! power-law IMF model: PDF index, lower and higher mass bounds (Msun)
  real(dp):: lt_t0, lt_m0, lt_a, lt_b ! Stellar lifetime model: t(M) = lt_t0 * exp(lt_a * (log(lt_m0 / M))**lt_b)

!  real(dp):: stf_K, stf_m0, stf_a, stf_b, stf_c 
  ! Stellar ionizing flux model: S(M) = stf_K * (M / stf_m0)**stf_a / (1 + (M / stf_m0)**stf_b)**stf_c
  ! This is a fit from Vacca et al. 1996
  ! Corresponding routine : vaccafit
  real(dp)::stf_K=9.634642584812752d48 ! s**(-1) then normalised in code units in read_stellar
  real(dp)::stf_m0=2.728098824280431d1 ! Msun then normalised in code units in read_stellar
  real(dp)::stf_a=6.840015602892084d0
  real(dp)::stf_b=4.353614230584390d0
  real(dp)::stf_c=1.142166657042991d0 

  !     hii_w: density profile exponent (n = n_0 * (r / r_0)**(-hii_w))
  !     hii_alpha: recombination constant in code units
  !     hii_c: HII region sound speed
  !     hii_t: fiducial HII region lifetime, it is normalised in code units in read_stellar 
  !     hii_T2: HII region temperature
  real(dp):: hii_w, hii_alpha, hii_c, hii_t, hii_T2
  real(dp):: mH_code ! mH in code units
  real(dp):: stellar_msink_th ! sink mass threshold for stellar object creation (Msun)
  real(dp), allocatable, dimension(:, :):: xstellar ! stellar object position
  real(dp), allocatable, dimension(:):: mstellar, tstellar, ltstellar ! stellar object mass, birth time, life time
  integer, allocatable, dimension(:):: id_stellar !the id  of the sink to which it belongs



  ! Proto-stellar jet feedback
!  real(dp),allocatable,dimension(:)::vol_gas_jet,vol_gas_jet_all
  logical::jet_feedback_sink=.false.         ! Put jet on sink 
  real(dp)::mass_jet_sink=0.                 ! Mass above which a jets is included
  real(dp)::frac_acc_ej=0.333333333333d0     ! Fraction of mass ejected into jet
  real(dp)::cone_jet=90.                     ! Outflow cone opening angle of jet; in degree
  real(dp)::v_jet = 1.6d6                    ! 16 km/s from Wang et al. 2010
  real(dp)::expo_jet=5.d-1                   ! Powerlaw exponent of jet velocity on mass
  logical::verbose_jet=.false.

  ! Use the supernova module?
  logical::SN_on = .false.

  !series of supernovae specified by "hand"
  ! Number of supernovae (max limit and number active in namelist)
  integer,parameter::NSNMAX=1000
  integer::SN_nsource

  ! Supernova times and whether they're done or not
  real(dp),dimension(1:NSNMAX)::SN_time = 1d10
  logical,dimension(1:NSNMAX)::SN_done = .false.
  
  ! SN position in units from 0 to 1
  real(dp),dimension(1:NSNMAX)::SN_pos_x = 0.5d0
  real(dp),dimension(1:NSNMAX)::SN_pos_y = 0.5d0
  real(dp),dimension(1:NSNMAX)::SN_pos_z = 0.5d0

  ! Ejecta mass in solar masses
  real(dp),dimension(1:NSNMAX)::SN_mejecta = 1.d0

  ! Energy of supernova in ergs
  real(dp),dimension(1:NSNMAX)::SN_energy = 1.d51

  ! Momentum of supernova in g cm s^-1 (needed if momentum is added instead of energy)
  real(dp),dimension(1:NSNMAX)::SN_momentum = 4.d43

  ! Use a thermal blast? (Otherwise add kinetic energy)
  ! Note that if SN_radius is 0 it's thermal anyway
  logical,dimension(1:NSNMAX)::SN_thermal = .false.

  ! Radius to deposit energy inside in number of cells (at highest level)
  real(dp),dimension(1:NSNMAX)::SN_radius = 12d0

  ! Radius in number of cells at highest level to refine fully around SN
  integer,dimension(1:NSNMAX)::SN_r_refine = 10


CONTAINS



!################################################################
!################################################################
!################################################################
!################################################################
!*************************************************************************
SUBROUTINE read_feedback_params(nml_ok)

! Read FEEDBACK_PARAMS namelist
!-------------------------------------------------------------------------
  use amr_parameters,only:dp
  use amr_commons
  implicit none
  logical::nml_ok
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
!-------------------------------------------------------------------------
  namelist/FEEDBACK_PARAMS/  Tsat, Vsat, sn_r_sat, sn_p, sn_e, sn_mass, &
       & SN_nsource, SN_on, SN_time, &
       & SN_pos_x, SN_pos_y, SN_pos_z, &
       & SN_mejecta, SN_energy, SN_thermal, SN_momentum, &
       & SN_radius, SN_r_refine, Vdisp, &
       & jet_feedback_sink, mass_jet_sink, frac_acc_ej, cone_jet, & 
       & v_jet, expo_jet, verbose_jet
  rewind(1)
  read(1,NML=FEEDBACK_PARAMS,END=101)
101 continue


      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

      !normalise the supernova quantities
      sn_p_ref = sn_p / (scale_d * scale_v * scale_l**3)
      sn_e_ref = sn_e / (scale_d * scale_v**2 * scale_l**3)
      sn_mass_ref = sn_mass / (scale_d * scale_l**3) !10 solar mass ejected

      !normalise Vsat which is assumed to be in KM/S
      Vsat = Vsat * 1.e5 / scale_v


      !normalise Vdisp which is assumed to be in KM/S
      Vdisp = Vdisp * 1.e5 / scale_v

END SUBROUTINE read_feedback_params
!*************************************************************************
!*************************************************************************
!*************************************************************************
!*************************************************************************
SUBROUTINE vaccafit(M,S)
  ! This routine is called in sink_RT_feedback
  ! perform a fit of the Vacca et al. 96 ionising flux
  ! M - stellar mass / solar masses
  ! S - photon emission rate in / s

!  use amr_parameters

  real(dp),intent(in)::M
  real(dp),intent(out)::S
  
  S = stf_K * (M / stf_m0)**stf_a / (1. + (M / stf_m0)**stf_b)**stf_c

END SUBROUTINE
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sn_stellar
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx
  real(dp):: scale, dx_min, dx_loc, vol_loc
  real(dp), dimension(1:3):: xbound, skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_p, sn_e, sn_vol, sn_d, sn_ed, sn_rp
  real(dp):: rr, pi,dens_max,pgas,dgas,ekin,mass_sn_tot,dens_max_all,mass_sn_tot_all
  integer:: n_sn,n_sn_all,info
  integer ,dimension(1:nvector)::cc
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:ndim):: xshift, x_sn

  logical, save:: first = .true.
  real(dp), save:: xseed

  real(dp) ::dens_max_loc,dens_max_loc_all

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m, pc

  logical, dimension(1:nstellarmax):: mark_del
  integer:: istellar

  real(dp)::T_sn,sn_ed_lim,pnorm_sn,pnorm_sn_all,vol_sn,vol_sn_all,vol_rap, mass_sn, mass_sn_all,dens_moy

  real(dp)::pgas_check,pgas_check_all

  integer, parameter:: navg = 3
  integer, parameter:: nsph = 1
  real(dp), dimension(1:nsph, 1:ndim):: avg_center
  real(dp), dimension(1:nsph):: avg_radius
  real(dp), dimension(1:navg):: avg_rpow
  real(dp), dimension(1:navg, 1:nvar+3):: avg_upow
  real(dp), dimension(1:navg, 1:nsph):: avg

  real(dp):: norm, rad_sn

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn_stellar'

  pi = acos(-1.0)

  if (first) then
     xseed = 0.5
     call random_number(xseed)
     first = .false.
  endif

  ! Mesh spacing in that level
  xbound(1:3) = (/ dble(nx), dble(ny), dble(nz) /)
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
  skip_loc(1) = dble(icoarse_min)
  skip_loc(2) = dble(jcoarse_min)
  skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*(scale_l**3)
  pc = 3.08d18 / scale_l

  sn_r = 3.0d0*(0.5d0**levelmin)*scale
  if(sn_r_sat .ne. 0) sn_r = max(sn_r, sn_r_sat * pc) !impose a minimum size of 12 pc for the radius
!  sn_r = 2.*(0.5**levelmin)*scale
  sn_m = sn_mass_ref
  sn_p = sn_p_ref
  sn_e = sn_e_ref
  sn_rp = 0.

!  if(sn_r /= 0.0) then
    sn_vol = 4. / 3. * pi * sn_r**3
!    sn_d = sn_m / sn_vol
!    sn_ed = sn_e / sn_vol
!  end if

  !we loop over stellar objets to determine whether one is turning supernovae
  !after it happens, the object is removed from the list
  mark_del = .false.
  do istellar = 1, nstellar
    if(t - tstellar(istellar) < ltstellar(istellar)) cycle
    mark_del(istellar) = .true.

    
    !the velocity dispersion times the life time of the object
    rad_sn = ltstellar(istellar)*Vdisp

    !find a random point within a sphere of radius < 1
    norm=2.
    do while (norm .gt. 1)
       norm=0.
       do idim = 1, ndim
          call random_number(xseed)
          xshift(idim) = (xseed-0.5)*2.
          norm = norm + xshift(idim)**2
       end do
    end do
    do idim = 1, ndim - 1 
        xshift(idim) = xshift(idim) * rad_sn
    end do
    !special treatment for the z-coordinates to maintain it at low altitude
    xshift(3) = xshift(3) * min(rad_sn,100.)

!    x_sn(:) = xstellar(istellar, :) + xshift(:)
    !place the supernovae around sink particles
    x_sn(:) = xsink(id_stellar(istellar), :) + xshift(:)

    !apply periodic boundary conditions (only along x and y)
    if( x_sn(1) .lt. 0) x_sn(1) = boxlen - x_sn(1) 
    if( x_sn(2) .lt. 0) x_sn(2) = boxlen - x_sn(2) 
    if( x_sn(1) .gt. boxlen) x_sn(1) = - boxlen + x_sn(1) 
    if( x_sn(2) .gt. boxlen) x_sn(2) = - boxlen + x_sn(2) 



    if(.true.) then
      avg_center(1, :) = x_sn(:)
      avg_radius = sn_r
      avg_rpow = 0.0d0
      avg_upow = 0.0d0
      ! avg_rpow(1) = 0 ; avg_upow(1, :) = 0 -> integrand(1) = 1
      ! avg_rpow(2) = 0 ; avg_upow(2, 1) = 1 -> integrand(2) = density
      avg_upow(2, 1) = 1.0d0
      ! avg_rpow(3) = 1 ; avg_upow(3, :) = 0 -> integrand(3) = radius
      avg_rpow(3) = 1.0d0
      call sphere_average(navg, nsph, avg_center, avg_radius, avg_rpow, avg_upow, avg)
      vol_sn = avg(1, 1)
      mass_sn = avg(2, 1) + vol_sn * sn_d ! region average + ejecta
      pnorm_sn = avg(3, 1)
    else
    !do a first path to compute the volume of the cells that are enclosed in the supernovae radius
    !this is to correct for the grid effects
    vol_sn = 0. ; vol_sn_all=0.
    mass_sn = 0. ; mass_sn_all=0.
    pnorm_sn = 0. ; pnorm_sn_all=0.

    do ilevel = levelmin, nlevelmax
      ! Computing local volume (important for averaging hydro quantities)
      dx = 0.5d0**ilevel
      dx_loc = dx * scale
      vol_loc = dx_loc**ndim

      ! Cell center position relative to grid center position
      do ind=1,twotondim
        iz = (ind - 1) / 4
        iy = (ind - 1 - 4 * iz) / 2
        ix = (ind - 1 - 2 * iy - 4 * iz)
        xc(ind,1) = (dble(ix) - 0.5d0) * dx
        xc(ind,2) = (dble(iy) - 0.5d0) * dx
        xc(ind,3) = (dble(iz) - 0.5d0) * dx
      end do

      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
        ngrid = min(nvector, ncache - igrid + 1)
        do i = 1, ngrid
          ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        end do

        ! Loop over cells
        do ind = 1, twotondim
          ! Gather cell indices
          iskip = ncoarse + (ind - 1) * ngridmax
          do i = 1, ngrid
            ind_cell(i) = iskip + ind_grid(i)
          end do

          ! Gather cell center positions
          do i = 1, ngrid
            xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
          end do
          ! Rescale position from coarse grid units to code units
          do i = 1, ngrid
             xx(i, :) = (xx(i, :) - skip_loc(:)) * scale
          end do

          ! Flag leaf cells
          do i = 1, ngrid
            ok(i) = (son(ind_cell(i)) == 0)
          end do

          do i = 1, ngrid
            if(ok(i)) then
                rr = sum(((xx(i, :) - x_sn(:)) / sn_r)**2)

                if(rr < 1.) then
                  vol_sn = vol_sn + vol_loc
                  mass_sn = mass_sn + vol_loc * (uold(ind_cell(i), 1) + sn_d)
                  pnorm_sn = pnorm_sn + vol_loc * sqrt(rr)
                endif
            endif
          end do
          !  End loop over sublist of cells
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(vol_sn,vol_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mass_sn,mass_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pnorm_sn,pnorm_sn_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

    vol_sn  = vol_sn_all
    mass_sn = mass_sn_all
    pnorm_sn = pnorm_sn_all
#endif
    end if

    !compute energy and mass density
    sn_d = sn_m / vol_sn
    sn_ed = sn_e / vol_sn

    dens_moy = mass_sn / vol_sn

    dens_max_loc = 0.
    mass_sn_tot = 0.
    dens_max_loc_all = 0.
    mass_sn_tot_all = 0.
    n_sn=0.

    pgas_check=0.

    !now loop over cells again and damp energies, mass and momentum
    !loop over levels 
    do ilevel = levelmin, nlevelmax
      ! Computing local volume (important for averaging hydro quantities)
      dx = 0.5d0**ilevel
      dx_loc = dx * scale
      vol_loc = dx_loc**ndim


      ! Cell center position relative to grid center position
      do ind=1,twotondim
        iz = (ind - 1) / 4
        iy = (ind - 1 - 4 * iz) / 2
        ix = (ind - 1 - 2 * iy - 4 * iz)
        xc(ind,1) = (dble(ix) - 0.5d0) * dx
        xc(ind,2) = (dble(iy) - 0.5d0) * dx
        xc(ind,3) = (dble(iz) - 0.5d0) * dx
      end do

      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid = 1, ncache, nvector
        ngrid = min(nvector, ncache - igrid + 1)
        do i = 1, ngrid
          ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
        end do

        ! Loop over cells
        do ind = 1, twotondim
          ! Gather cell indices
          iskip = ncoarse + (ind - 1) * ngridmax
          do i = 1, ngrid
            ind_cell(i) = iskip + ind_grid(i)
          end do

          ! Gather cell center positions
          do i = 1, ngrid
            xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
          end do
          ! Rescale position from coarse grid units to code units
          do i = 1, ngrid
             xx(i, :) = (xx(i, :) - skip_loc(:)) * scale
          end do

          ! Flag leaf cells
          do i = 1, ngrid
            ok(i) = (son(ind_cell(i)) == 0)
          end do

          do i = 1, ngrid
            if(ok(i)) then
                rr = sqrt(sum(((xx(i, :) - x_sn(:)) / sn_r)**2))

                if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  dgas = uold(ind_cell(i), 1)

                  if(dgas .gt. dens_max_loc) dens_max_loc = dgas

                  mass_sn_tot = mass_sn_tot + dgas*vol_loc


                  !compute velocity of the gas within this cell assuming 
                  !energy equipartition
                  !pgas = sqrt(eff_sn*sn_ed / max(dens_moy,dgas) ) * dgas 
                  !for cells where dgas < dens_moy, take the same velocity as if
                  !the density were equal to dens_moy
                  pgas = min(sn_p / pnorm_sn * rr /  dgas , Vsat) * dgas
                  pgas_check = pgas_check + pgas * vol_loc

                  ekin = ((uold(ind_cell(i),2))**2 + (uold(ind_cell(i),3))**2 + (uold(ind_cell(i),4))**2) / dgas / 2.
                  uold(ind_cell(i), 5) = uold(ind_cell(i), 5) - ekin

                  uold(ind_cell(i), 2:4) = uold(ind_cell(i), 2:4) + pgas * (xx(i, 1:3) - x_sn(1:3)) / (rr * sn_r)

                  ekin = ((uold(ind_cell(i),2))**2 + (uold(ind_cell(i),3))**2 + (uold(ind_cell(i),4))**2) / dgas / 2.

                  !before adding thermal energy make sure the temperature is not too high (too small timesteps otherwise)
                  T_sn = (sn_ed / dgas * (gamma-1.) ) * scale_T2
                  T_sn = min( T_sn , Tsat) / scale_T2
                  sn_ed_lim = T_sn * dgas / (gamma-1.)

                  uold(ind_cell(i), 5) = uold(ind_cell(i), 5) + ekin + sn_ed_lim

                  n_sn = n_sn + 1

                  !write(*,*) 'put SN, myid , n_sn ',myid, n_sn, 'x,y,z: ',x_sn(1),x_sn(2),x_sn(3), 'density ',dgas

                end if

            end if
          end do
          ! End loop over sublist of cells
        end do
        ! End loop over cells
      end do
      ! End loop over grids
    end do
    ! End loop over levels


#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(n_sn,n_sn_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(mass_sn_tot,mass_sn_tot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(dens_max_loc,dens_max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
    call MPI_ALLREDUCE(pgas_check,pgas_check_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
    n_sn_all = n_sn
    mass_sn_tot_all = mass_sn_tot
    dens_max_loc_all = dens_max_loc
    pgas_check_all = pgas_check
#endif

    if(myid == 1) write(*, *) "SN momentum (injected, expected):", pgas_check_all, sn_p
    if(myid == 1) write(*, *) "Physical units:", pgas_check_all * scale_d * scale_l**3 * scale_v, sn_p * scale_d * scale_l**3 * scale_v

  !write(*,*) '4 n_sn ', n_sn_all

  !  if(myid .eq. cc(1)) write(*,*) '4 n_sn ', n_sn_all

    !calculate grid effect
    vol_rap = vol_sn / sn_vol

    if(myid .eq. 1) then 
       open(103,file='supernovae2.txt',form='formatted',status='unknown',access='append')
         write(103,112) t,x_sn(1),x_sn(2),x_sn(3),dens_max_loc_all,mass_sn_tot_all,vol_rap,pgas_check_all,sn_p
       close(103)
    endif

112 format(9e12.4)

  end do ! end of the loop over stellar objects

  call delete_stellar(mark_del)

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo
end subroutine make_sn_stellar
!################################################################
!################################################################
!################################################################
!################################################################
subroutine sphere_average(navg, nsph, center, radius, rpow, upow, avg)
    use amr_parameters, only: boxlen, dp, hydro, icoarse_max, icoarse_min &
        & , jcoarse_min, kcoarse_min, levelmin, ndim, ngridmax, nlevelmax &
        & , nvector, twotondim, verbose
    use amr_commons, only: active, ncoarse, son, xg, myid
    use hydro_parameters, only: nvar
    use hydro_commons, only: uold
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif

    ! Integrate quantities over spheres
    ! The integrand is (r / radius(isph))**rpow(iavg) * product(u(ivar)**upow(iavg, ivar), ivar=1:nvar)

    integer, intent(in):: navg                               ! Number of quantities
    integer, intent(in):: nsph                               ! Number of spheres
    real(dp), dimension(1:nsph, 1:ndim), intent(in):: center ! Sphere centers
    real(dp), dimension(1:nsph), intent(in):: radius         ! Sphere radii
    real(dp), dimension(1:navg), intent(in):: rpow           ! Power of radius in the integral
#ifdef SOLVERmhd
    real(dp), dimension(1:navg, 1:nvar+3), intent(in):: upow ! Power of hydro variables in the integral
#else
    real(dp), dimension(1:navg, 1:nvar), intent(in):: upow   ! Power of hydro variables in the integral
#endif
    real(dp), dimension(1:navg, 1:nsph), intent(out):: avg   ! Averages

    integer:: i, ivar, ilevel, igrid, ind, ix, iy, iz, iskip, isph
    integer:: nx_loc, ncache, ngrid
    integer, dimension(1:nvector):: ind_grid, ind_cell
    logical, dimension(1:nvector):: ok

    real(dp):: scale, dx, dx_loc, vol_loc, rr
    real(dp), dimension(1:ndim):: skip_loc
    real(dp), dimension(1:twotondim, 1:ndim):: xc
    real(dp), dimension(1:nvector, 1:ndim):: xx

    integer:: info
    real(dp), dimension(1:navg, 1:nsph):: avg_loc
    real(dp), dimension(1:navg):: integrand
    real(dp), dimension(1:navg):: utemp

    if(.not. hydro)return
    if(ndim .ne. 3)return

    if(verbose .and. myid == 1) write(*, *) 'Entering sphere_average'

    ! Mesh spacing in that level
    nx_loc = icoarse_max - icoarse_min + 1
    skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
    skip_loc(1) = dble(icoarse_min)
    skip_loc(2) = dble(jcoarse_min)
    skip_loc(3) = dble(kcoarse_min)
    scale = boxlen / dble(nx_loc)

    avg_loc = 0.0d0

    do ilevel = levelmin, nlevelmax
        ! Computing local volume (important for averaging hydro quantities)
        dx = 0.5d0**ilevel
        dx_loc = dx * scale
        vol_loc = dx_loc**ndim

        ! Cell center position relative to grid center position
        do ind = 1, twotondim
            iz = (ind - 1) / 4
            iy = (ind - 1 - 4 * iz) / 2
            ix = (ind - 1 - 2 * iy - 4 * iz)
            xc(ind, 1) = (dble(ix) - 0.5d0) * dx
            xc(ind, 2) = (dble(iy) - 0.5d0) * dx
            xc(ind, 3) = (dble(iz) - 0.5d0) * dx
        end do

        ! Loop over grids
        ncache = active(ilevel)%ngrid
        do igrid = 1, ncache, nvector
            ngrid = min(nvector, ncache - igrid + 1)
            do i = 1, ngrid
                ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
            end do

            ! Loop over cells
            do ind = 1, twotondim
                ! Gather cell indices
                iskip = ncoarse + (ind - 1) * ngridmax
                do i = 1, ngrid
                    ind_cell(i) = iskip + ind_grid(i)
                end do

                ! Gather cell center positions
                do i = 1, ngrid
                    xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
                end do

                ! Rescale position from coarse grid units to code units
                do i = 1, ngrid
                    xx(i, :) = (xx(i, :) - skip_loc(:)) * scale
                end do

                ! Flag leaf cells
                do i = 1, ngrid
                    ok(i) = (son(ind_cell(i)) == 0)
                end do

                do i = 1, ngrid
                    if(ok(i)) then
                        do isph = 1, nsph
                            rr = sqrt(sum(((xx(i, :) - center(isph, :)) / radius(isph))**2))

                            if(rr < 1.) then
                                integrand = rr**rpow
                                where(abs(rpow) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                    integrand = 1.0d0
                                end where
#ifdef SOLVERmhd
                                do ivar = 1, nvar + 3
#else
                                do ivar = 1, nvar
#endif
                                    utemp(:) = uold(ind_cell(i), ivar)
                                    where(abs(upow(:, ivar)) < 1.0d-10) ! Avoid NaNs of the form 0**0
                                        utemp = 1.0d0
                                    end where
                                    integrand = integrand * utemp**upow(:, ivar)
                                end do
                                avg_loc(:, isph) = avg_loc(:, isph) + vol_loc * integrand
                            endif
                            ! End test on radius
                        end do
                        ! End loop over spheres
                    endif
                    ! End test on leaf cells
                end do
                ! End loop over sublist of cells
            end do
            ! End loop over cells
        end do
        ! End loop over grids
    end do
    ! End loop over levels

#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(avg_loc, avg, navg * nsph, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
    avg = avg_loc
#endif
end subroutine sphere_average
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine feedback_sn(ilevel)
  use amr_commons
  use hydro_parameters
  implicit none
  integer::ilevel,isn
!---------------------
  ! Supernova explosions
  !---------------------
  if(SN_on .and. ilevel == levelmin) then
     do isn=1,SN_nsource
        if(t >= sn_time(isn) .and. .not. sn_done(isn)) then
           if(myid == 1) write (*,*) 'Supernova ',isn,' @ t = ', t, &
                & 'sn_time =', sn_time(isn)
           call make_sn_blast(isn)
           sn_done(isn) = .true.
        endif
     end do
  endif
end subroutine feedback_sn

subroutine make_sn_blast(isn)
  ! Adapted from O. Iffrig
  use amr_commons
  use hydro_commons
  implicit none

  integer, intent(in) :: isn

  integer:: ivar
  integer:: ilevel, ind, ix, iy, iz, ngrid, iskip, idim
  integer:: i, nx_loc, igrid, ncache
  integer, dimension(1:nvector), save:: ind_grid, ind_cell
  real(dp):: dx
  real(dp):: scale, dx_min, dx_loc, vol_loc

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun, scale_ecgs

  real(dp), dimension(1:3):: xbound, skip_loc
  real(dp), dimension(1:twotondim, 1:3):: xc
  logical, dimension(1:nvector), save:: ok

  real(dp),dimension(1:ndim):: sn_cent
  real(dp), dimension(1:nvector, 1:ndim), save:: xx
  real(dp):: sn_r, sn_m, sn_e, sn_vol, sn_d, sn_ed, dx_sel, sn_p, sn_v
  real(dp):: rr, pi
  real(dp), dimension(1:ndim)::rvec
  logical:: sel = .false.
  real(dp),parameter::m_sun=1.9891d33  ! Solar mass [g]

  !logical, save:: first = .true.
  !logical,dimension(1:n_source), save:: sn_done = .false.

  if(.not. hydro)return
  if(ndim .ne. 3)return

  if(verbose)write(*,*)'Entering make_sn'

  pi = acos(-1.0)

  ! Mesh spacing in that level
  xbound(1:3) = (/ dble(nx), dble(ny), dble(nz) /)
  nx_loc = icoarse_max - icoarse_min + 1
  skip_loc = (/ 0.0d0, 0.0d0, 0.0d0 /)
  skip_loc(1) = dble(icoarse_min)
  skip_loc(2) = dble(jcoarse_min)
  skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)
  dx_min = scale * 0.5d0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_d * scale_l**ndim / m_sun
  scale_ecgs = scale_d * scale_v**2 * scale_l**ndim

  ! Hard-code sn properties to centre of box
  sn_r = sn_radius(isn)*(0.5**nlevelmax)*scale

  !sn_m = 80d0 / scale_msun ! I dunno?
  sn_m = SN_mejecta(isn) / scale_msun ! Put in 10 solar masses
  sn_e = SN_energy(isn) / scale_ecgs
  sn_v = 0d0
  sn_v = sqrt(2.0*(sn_e*scale_ecgs)/(sn_m*scale_msun*m_sun))
  sn_v = sn_v / scale_v
  sn_cent(1)= SN_pos_x(isn)*boxlen
  sn_cent(2)= SN_pos_y(isn)*boxlen
  sn_cent(3)= SN_pos_z(isn)*boxlen
     
  ! HACK !!! - KINETIC BLAST ONLY WORKS FOR sn_r > 0.0 !!!
  if(sn_r /= 0.0) then
     sn_vol = 4. / 3. * pi * sn_r**3
     sn_d = sn_m / sn_vol
     sn_ed = sn_e / sn_vol
     sn_p = sn_d*sn_v ! uniform momentum of blast ejecta
  end if
     
  if(myid .eq. 1) then
     write(*,*) 'Supernova blast! Wow!'
     write(*,*) 'x_sn, y_sn, z_sn, ',sn_cent(1),sn_cent(2),sn_cent(3)
  endif

  ! Loop over levels
  do ilevel = levelmin, nlevelmax
    ! Computing local volume (important for averaging hydro quantities)
    dx = 0.5d0**ilevel
    dx_loc = dx * scale
    vol_loc = dx_loc**ndim
    !if(.not. sel) then
      ! dx_sel will contain the size of the biggest leaf cell around the center
      !dx_sel = dx_loc
      !sn_vol = vol_loc
    !end if

    ! Cell center position relative to grid center position
    do ind=1,twotondim
      iz = (ind - 1) / 4
      iy = (ind - 1 - 4 * iz) / 2
      ix = (ind - 1 - 2 * iy - 4 * iz)
      xc(ind,1) = (dble(ix) - 0.5d0) * dx
      xc(ind,2) = (dble(iy) - 0.5d0) * dx
      xc(ind,3) = (dble(iz) - 0.5d0) * dx
    end do

    ! Loop over grids
    ncache=active(ilevel)%ngrid
    do igrid = 1, ncache, nvector
      ngrid = min(nvector, ncache - igrid + 1)
      do i = 1, ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid + i - 1)
      end do

      ! Loop over cells
      do ind = 1, twotondim
        ! Gather cell indices
        iskip = ncoarse + (ind - 1) * ngridmax
        do i = 1, ngrid
          ind_cell(i) = iskip + ind_grid(i)
        end do

        ! Gather cell center positions
        do i = 1, ngrid
          xx(i, :) = xg(ind_grid(i), :) + xc(ind, :)
        end do
        ! Rescale position from coarse grid units to code units
        do i = 1, ngrid
           xx(i, :) = (xx(i, :) - skip_loc(:)) * scale
        end do

        ! Flag leaf cells
        do i = 1, ngrid
          ok(i) = (son(ind_cell(i)) == 0)
        end do

        do i = 1, ngrid
          if(ok(i)) then
            if(sn_r == 0.0) then
              sn_d = sn_m / vol_loc ! XXX
              sn_ed = sn_e / vol_loc ! XXX
              rr = 1.0
              do idim = 1, 3
                !rr = rr * max(1.0 - abs(xx(i, idim) - sn_center(sn_i, idim)) / dx_sel, 0.0)
                rr = rr * max(1.0 - abs(xx(i, idim) - sn_cent(idim)) / dx_loc, 0.0)
              end do
              !if(rr > 0.0) then
                !if(.not. sel) then
                  !! We found a leaf cell near the supernova center
                  !sel = .true.
                  !sn_d = sn_m / sn_vol
                  !sn_ed = sn_e / sn_vol
                !end if
                uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d * rr
                uold(ind_cell(i), 5) = uold(ind_cell(i), 5) + sn_ed * rr
              !end if
            else
               ! Get direction to point the explosion in
               rvec = (xx(i, :) - sn_cent(:)) / sn_r
               rr = sqrt(sum(rvec**2))
               rvec = rvec/rr
   
               if(rr < 1.) then
                  uold(ind_cell(i), 1) = uold(ind_cell(i), 1) + sn_d
                  ! If not entirely thermal injection, add some velocity
                  if (.not.SN_thermal(isn)) then
                     uold(ind_cell(i),2:4) = uold(ind_cell(i),2:4) + &
                          & sn_p * rvec
                  end if
                  uold(ind_cell(i), 5) = uold(ind_cell(i), 5) + sn_ed
               endif
            endif
          endif
        end do
      end do
      ! End loop over cells
    end do
    ! End loop over grids
  end do
  ! End loop over levels

  ! Update hydro quantities for split cells
  do ilevel = nlevelmax, levelmin, -1
    call upload_fine(ilevel)
    do ivar = 1, nvar
      call make_virtual_fine_dp(uold(1, ivar), ilevel)
    enddo
  enddo
end subroutine make_sn_blast

!################################################################
!################################################################
!################################################################
!################################################################

SUBROUTINE feedback_refine(xx,ok,ncell,ilevel)

! This routine flags cells immediately around SN sources to the finest
! level of refinement. The criteria for refinement at a point are:
! a) The point is less than one ilevel cell width from an SN source.
! b) The point is within SN_r_wind finest level cell widths from
!    the SN source.
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ncell,ilevel,i,k,nx_loc,isn
  real(dp),dimension(1:nvector,1:ndim)::xx
  logical ,dimension(1:nvector)::ok
  real(dp)::dx_loc,rvec(ndim),w,rmag,rFB
!-------------------------------------------------------------------------
  nx_loc=(icoarse_max-icoarse_min+1)
  dx_loc = boxlen*0.5D0**ilevel/dble(nx_loc)
  ! Loop over regions
  do isn=1,SN_nsource
     do i=1,ncell
        rFB = SN_r_refine(isn)*boxlen*0.5D0**nlevelmax
        rvec(1)=xx(i,1)-SN_pos_x(isn)*boxlen
        rvec(2)=xx(i,2)-SN_pos_y(isn)*boxlen
        rvec(3)=xx(i,3)-SN_pos_z(isn)*boxlen
        rmag=sqrt(sum(rvec**2))
        if(rmag .le. 2*rFB+dx_loc) then
           ok(i)=.true.
        endif
     end do
  end do
  
END SUBROUTINE feedback_refine





END MODULE feedback_module

