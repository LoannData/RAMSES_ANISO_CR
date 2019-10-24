module cloud_module
  use amr_parameters

  ! TODO - CLEAN THIS OUT

  !initial temperature used for the isothermal run
  real(dp)::temper
  real(dp)::temper_iso

  !feedback from jet
  logical:: jet = .false., rad_jet=.false. 
  real(dp)::Ucoef=1.
  real(dp):: mass_jet_sink=0. !mass above which a jets is included

  !Initial conditions parameter for the dense core
  real(dp)::bx_bound=0.
  real(dp)::by_bound=0.
  real(dp)::bz_bound=0.
  real(dp)::turb=0.
  real(dp)::dens0=0.
  real(dp)::V0=0.
  real(dp)::Height0=0.

  real(dp):: switch_solv=1.d20

  !Initial conditions parameter for the dense core
  real(dp)::mass_c=0d0   !cloud mass in solar mass
  real(dp)::rap=1.      !axis ratio
  real(dp)::cont=1.     !density contras
  real(dp)::ff_sct=1.   !freefall time / sound crossing time
  real(dp)::ff_rt=1.    !freefall time / rotation time
  real(dp)::ff_act=1.   !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.   !freefall time / Vrms crossing time
  real(dp)::thet_mag=0. !angle between magnetic field and rotation axis
  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor


end module cloud_module

!================================================================
!================================================================
!================================================================
!================================================================

subroutine calc_dmin(d_c)
  use amr_commons
  use hydro_commons
  use cloud_module
  implicit none

  ! NOTE!! - IS THIS REALLY NECESSARY? - SAM GEEN OCTOBER 2015

  real(dp):: d_c, cont_ic, dmin

  cont_ic = 10.
  dmin = d_c / cont / cont_ic

  if (myid == 1) then
    write(*,*) "dmin = ", dmin
  endif
end subroutine calc_dmin
!================================================================
!================================================================
!================================================================
!================================================================
subroutine calc_boxlen
  use amr_commons
  use amr_parameters
  use hydro_commons
  use poisson_parameters
  use cloud_module
!  use const
  implicit none
  !================================================================
  !this routine calculate boxlen
  !================================================================
  integer :: i
  real(dp):: pi
  real(dp):: d_c,zeta
  real(dp):: res_int,r_0,C_s
  integer::  np
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),save:: first
  real(dp):: mu=1.4d0 ! NOTE - MUST BE THE SAME AS IN units.f90!!
!  real(dp)::myid

!   myid=1

    if (first .eq. 0.) then

    pi=acos(-1.0d0)

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    !calculate the mass in code units (Msolar / Mparticle / pc^3
    mass_c = mass_c * (2.d33 / (scale_d * scale_l**3) )

    !calculate the sound speed
    C_s = sqrt( T2_star / scale_T2)

    !calculate  zeta=r_ext/r_0
    zeta = sqrt(cont - 1.)

    !calculate an integral used to compute the cloud radius
    np=1000
    res_int=0.
    do i=1,np
     res_int = res_int + log(1.+(zeta/np*i)**2) * zeta/np
    enddo
    res_int = zeta*log(1.+zeta**2) - res_int

    !now we determine the central density and the external cloud radius
    !we have mass = 2 pi rho_c r_0^2 z_0 * res_int
    !which results from the integration of rho = dc/(1.+(x^2+y^2)/r_O^2+z^2/z_0^2)
    !for (x^2+y^2)/r_O^2+z^2/z_0^2 < zeta
    !we also have ff_sct = sqrt(3. pi / 32 / G / d_c) C_s / (r_0)
    !which just state the ratio of freefall time over sound crossing time
    !from these 2 formula, rho_c and r_0 are found to be:



    r_0 = mass_c / (2.*pi*rap*res_int) * (ff_sct)**2 / (3.*pi/32.) / C_s**2

    d_c = mass_c / (2.*pi*rap*res_int) / r_0**3

    !it is equal to twice the length of the major axis
    boxlen = r_0 * zeta * max(rap,1.) * 4.

    ! Multiply boxlen by an extra factor
    boxlen = bl_fac * boxlen

    if (myid == 1) then
    write(*,*) '** Cloud parameters estimated in calc-boxlen **'
    write(*,*) 'inner radius (pc) ', r_0
    write(*,*) 'peak density (cc) ', d_c
    write(*,*) 'total box length (pc) ', boxlen
    write(*,*) 'cloud mass (code units) ', mass_c
    write(*,*) 'boxlen (code units) ',boxlen
    write(*,*)
    endif



    first=1.
    endif

    call calc_dmin(d_c)

end subroutine calc_boxlen

subroutine read_cloud_params(nml_ok)

  use amr_parameters
  use feedback_module
  use clfind_commons
  use cloud_module

  implicit none
  logical::nml_ok
  real(dp)::cellsize
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),parameter::pcincm=3.086d18

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cloud_params/mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,thet_mag &
       & ,bl_fac,switch_solv

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Calculate boxlen
  if (mass_c .gt. 0) then
     call calc_boxlen
  end if

  !since boxlen is not known initialy we must multiply the
  !refining parameters by boxlen here
  x_refine = x_refine*boxlen
  y_refine = y_refine*boxlen
  z_refine = z_refine*boxlen
  r_refine = r_refine*boxlen

  ! Set the sink formation threshold based on the Jeans criterion
  cellsize = boxlen * 0.5**nlevelmax * pcincm / scale_l
  n_sink = 881.0 / cellsize**2 ! Scaled to give 1e6 for 30pc/1024
  n_clfind = 0.1 * n_sink
  if(myid==1) write(*,*) "SETTING n_sink, n_clfind TO", n_sink, n_clfind

  ! Feedback parameters
  call read_feedback_params(nml_ok)


end subroutine read_cloud_params

