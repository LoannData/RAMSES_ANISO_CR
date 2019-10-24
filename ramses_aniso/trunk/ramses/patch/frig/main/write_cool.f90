!-------------------------------------------------------------------!
! Write the cooling table being used out to file for debugging      !
! Sam Geen, October-November 2015                                   !
!-------------------------------------------------------------------!

SUBROUTINE output_testcool(filename,Ds,Ts,out_rate)
  ! Write a single table to file
  implicit none

  character (len=*) :: filename

  integer,parameter :: nD=4
  integer,parameter :: nT=200 
  real(kind=8), dimension(nD) :: Ds
  real(kind=8), dimension(nT) :: Ts
  real(kind=8), dimension(nD,nT) :: out_rate

  ! Write no flux file
  write(*,*) "Writing file ", filename
  open(unit=8075,file=filename,form='unformatted')
  write(8075) nD,nT
  write(8075) Ds
  write(8075) Ts
  write(8075) out_rate
  write(*,*) "outrate sample: ", out_rate(1,1:10)
  close(8075)

END SUBROUTINE output_testcool

SUBROUTINE write_cool
  ! Write various tables
  ! Function called for reference
  ! rt_metal_cool(dT2,nH(icell),flux,dXion(1),mu,ne,metal_tot,metal_prime)

  ! This is hard-coded for STG's metal cooling function
  !use rt_metal_cooling_module
  
  use cooling_module,only:X, Y,solve_cooling
  use rt_parameters
  use rt_cooling_module
  use rt_metal_cooling_module
  use coolrates_module,only:compCoolrate,inp_coolrates_table

  implicit none

  integer :: iD,iT
  integer,parameter :: nD=4
  integer,parameter :: nT=200 
  real(kind=8) :: ne, nH, T, mu, afact, rate, dtcool
  real(kind=8) :: nHe, xH, xHeI, xHeII, xHeIII
  real(kind=8) :: metal_tot,metal_prime
  real(kind=8), dimension(nD) :: Ds
  real(kind=8), dimension(nT) :: Ts
  real(kind=8), dimension(nD,nT) :: out_rate
  real(kind=8), dimension(1) :: nHa, Ta, boost, zsolar, rateout

  real(kind=8), dimension(3) :: nN,nI

  write(*,*) "Writing cooling table..."

  ! Set temperature, density arrays
  !Ds = 0d0
  Ds(1) = -2d0
  Ds(2) = 0d0
  Ds(3) = 2d0
  Ds(4) = 4d0
  Ts = 0d0
  do iT=1,nT
     ! 0 to 7 in log space
     Ts(iT) = real(iT-1,kind=8)/real(nT,kind=8)*7d0
  end do
  ! Cast to linear space
  Ds = 10d0**Ds
  Ts = 10d0**Ts

  ! Some values
  afact = 1d0
  mu = 1./(0.76*2.0) ! Try ionised hydrogen (~ ionised helium too)

  ! Vanilla RT metals
  do iD=1,nD
     nH = Ds(iD)
     do iT=1,nT
        T = Ts(iT)
        call rt_cmp_metals(T,nH,mu,metal_tot,metal_prime,afact)
        out_rate(iD,iT) = metal_tot / (nH*nH)
     end do
  end do
  call output_testcool("./cool_vanillart.out",Ds,Ts,out_rate)

  ! FRIG + Osterbrock RT metals
  do iD=1,nD
     nH = Ds(iD)
     do iT=1,nT
        T = Ts(iT)
        if (T < 2e3) then
           ne = 2.4E-3*((T/100.)**0.25)/0.5 !formule C15 Wolfire et al. 2003
        else
           ne = nH
        endif
        call rt_metal_cool(T,nH,1d0,mu,metal_tot,metal_prime)
        out_rate(iD,iT) = metal_tot / (nH*nH)
     end do
  end do
  call output_testcool("./cool_frigosterbrock.out",Ds,Ts,out_rate)

  if (rt) then
  ! Hydrogen + helium RT cooling
  do iD=1,nD
     ! Set the neutral and ionised densities
     ! ASSUME COMPLETELY IONISED HYDROGEN + HELIUM !!!
     nH = Ds(iD)
     nHe = 0.25*nH*Y/X
     xH  = 1d0-1d-6
     xHeI = 1d-6
     xHeII = 1d-6
     xHeIII = 1d0-xHeI-xHeII
     ! ASSUME COMPLETELY NEUTRAL
     !xH = 0.1
     !xHeI = 1d0-2d6
     !xHeII = 1d-6
     !xHeIII = 1d-6
     nN(1) = nH*(1d0-xH)
     nN(2) = nHe*xHeI
     nN(3) = nHe*xHeII
     nI(1) = nH*xH
     nI(2) = nN(3) ! idk
     nI(3) = nHe*xHeIII
     ne = nH*xH + nHe*(xHeII + 2d0*xHeIII)
     do iT=1,nT
        T = Ts(iT)
        rate = compCoolrate(T, ne, nN(1), nI(1), nN(2), nN(3), nI(3)    &
             ,afact, metal_prime, RT_OTSA) 
        out_rate(iD,iT) = rate / (nH*nH)
     end do
  end do
  call output_testcool("./cool_HHeRT.out",Ds,Ts,out_rate)
  endif

  if (.not.rt) then
  ! Vanilla solve_cooling (no rt)
  do iD=1,nD
     nH = Ds(iD)
     nHa(1) = nH
     boost(1) = 1d0
     dtcool = 1d-7 ! idk
     zsolar(1) = 1d0
     do iT=1,nT
        T = Ts(iT)
        Ta(1) = T
        call solve_cooling(nHa,Ta,zsolar,boost,dtcool,rateout,1)
        out_rate(iD,iT) = rateout(1)! / (nH*nH)
     end do
  end do
  call output_testcool("./cool_vanillanort.out",Ds,Ts,out_rate)
  endif

  write(*,*) "Cooling tables written, stopping..."

  call clean_stop

END SUBROUTINE write_cool
