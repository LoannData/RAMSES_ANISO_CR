! Non-equlibrium (in H and He) cooling module for radiation-hydrodynamics.
! For details, see Rosdahl et al. 2013, and Rosdahl & Teyssier 2015.
! Joki Rosdahl, Andreas Bleuler, and Romain Teyssier, September 2015.

module rt_cooling_module
  use amr_commons,only:myid  
  use cooling_module,only:X, Y
  use rt_parameters
  use coolrates_module
  implicit none

  private   ! default

  public rt_set_model, rt_solve_cooling, update_UVrates, cmp_chem_eq     &
         , isHe, is_mu_H2, X, Y, rhoc, kB, mH, T2_min_fix, twopi         &
         , signc, sigec, PHrate, UVrates, rt_isIR, kappaAbs, kappaSc     &
         , is_kIR_T, iIR, rt_isIRtrap, iIRtrapVar, rt_pressBoost         &
         , rt_isoPress, rt_T_rad, rt_vc, a_r

  ! NOTE: T2=T/mu
  ! Np = photon density, Fp = photon flux, 

  real(dp),parameter::rhoc      = 1.88000d-29    !  Crit. density [g cm-3]
  real(dp),parameter::mH        = 1.66000d-24    !         H atom mass [g]
  real(dp),parameter::kB        = 1.38062d-16    ! Boltzm.const. [erg K-1]
  real(dp),parameter::a_r       = 7.5657d-15   ! Rad.const. [erg cm-3 K-4]
  real(dp),parameter::mu_mol    = 1.2195D0
  real(dp),parameter::T2_min_fix=1.d-2           !     Min temperature [K]
  real(dp),parameter::twopi     = 6.2831853d0    !            Two times pi

  real(dp)::T_min, T_frac, x_min, x_frac, Np_min, Np_frac, Fp_min, Fp_frac

  integer,parameter::iIR=1                       !          IR group index
  integer::iIRtrapVar=1                          ! IRtrap passScalar index
  ! Namelist parameters:
  logical::isHe=.true.
  logical::is_mu_H2=.false.
  logical::rt_isoPress=.false.         ! Use cE, not F, for rad. pressure
  real(dp)::rt_pressBoost=1d0          ! Boost on RT pressure            
  logical::rt_isIR=.false.             ! Using IR scattering on dust?    
  logical::rt_isIRtrap=.false.         ! IR trapping in passive scalar?  
  logical::is_kIR_T=.false.            ! k_IR propto T^2?               
  logical::rt_T_rad=.false.            ! Use T_gas = T_rad
  logical::rt_vc=.false.               ! (semi-) relativistic RT
  real(dp)::Tmu_dissoc=1d3             ! Dissociation temperature [K]
  real(dp),dimension(nGroups)::kappaAbs=0! Dust absorption opacity    
  real(dp),dimension(nGroups)::kappaSc=0 ! Dust scattering opacity    
  
  ! Cooling constants, updated on SED and c-change [cm3 s-1],[erg cm3 s-1]
  real(dp),dimension(nGroups,nIons)::signc,sigec,PHrate

  real(dp),dimension(nIons, 2)::UVrates     !UV backgr. heating/ion. rates

CONTAINS 

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_set_model(Nmodel, J0in_in, J0min_in, alpha_in              &
     ,normfacJ0_in, zreioniz_in, correct_cooling, realistic_ne, h        &
     ,omegab, omega0, omegaL, astart_sim, T2_sim)
! Initialize cooling. All these parameters are unused at the moment and
! are only there for the original cooling-module.
! Nmodel(integer)     =>     Model for UV background and metals
! J0in_in  (dble)     => Default UV intensity
! J0min_in (dble)     => Minimum UV intensity
! alpha_in (dble)     => Slope of the UV spectrum
! zreioniz_in (dble)  => Reionization redshift
! normfacJ0_in (dble) => Normalization factor fot a Harrdt&Madau UV model
! correct_cooling (integer) => Cooling correction
! realistic_ne (integer) => Use realistic electron density at high z?
! h (dble)            => H0/100
! omegab (dble)       => Omega baryons
! omega0 (dble)       => Omega materal total
! omegaL (dble)       => Omega Lambda
! astart_sim (dble)   => Redshift at which we start the simulation
! T2_sim (dble)      <=  Starting temperature in simulation?
!-------------------------------------------------------------------------
  use UV_module
  use coolrates_module,only: init_coolrates_tables
  real(kind=8) :: J0in_in, zreioniz_in, J0min_in, alpha_in, normfacJ0_in
  real(kind=8) :: astart_sim, T2_sim, h, omegab, omega0, omegaL
  integer  :: Nmodel, correct_cooling, realistic_ne, ig
  real(kind=8) :: astart=0.0001, aend, dasura, T2end=T2_min_fix, mu, ne
  character(128) :: metalcoolfile = '.'
!-------------------------------------------------------------------------
  if(myid==1) write(*,*) &
       '==================RT momentum pressure is turned ON=============='
  if(myid==1 .and. rt_isIR) &
       write(*,*) 'There is an IR group, with index ',iIR        
  if(myid==1 .and. rt_isIRtrap) write(*,*) &
       '=========IR trapping is turned ON=============='
  ! do initialization
  isHe=.true. ; if(Y .le. 0.) isHe=.false.
  T_MIN           = 0.1                  !                      Minimum T2
  T_FRAC          = 0.1            

  x_MIN           = 1.d-6                !    Minimum ionization fractions
  x_FRAC          = 0.1    

  Np_MIN = 1.d-13                        !            Photon density floor
  Np_FRAC = 0.2    

  Fp_MIN  = 1D-13*rt_c_cgs               !           Minimum photon fluxes
  Fp_FRAC = 0.5

  ! Calculate initial temperature
  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif
  aend=astart_sim
  dasura=0.02d0

  call update_rt_c
  call init_UV_background
  if(cosmo) then
     call update_UVrates(aexp)
     call init_coolrates_tables(aexp)
  else
     call update_UVrates(astart_sim)
     call init_coolrates_tables(astart_sim)
  endif

  if(nrestart==0 .and. cosmo)                                            &
       call rt_evol_single_cell(astart,aend,dasura,h,omegab,omega0       &
                               ,omegaL,-1.0d0,T2end,mu,ne,.false.)
  T2_sim=T2end

END SUBROUTINE rt_set_model

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_UVrates(aexp)
! Set the UV ionization and heating rates according to the given a_exp.
!-------------------------------------------------------------------------
  use UV_module
  use amr_parameters,only:haardt_madau
  integer::i
  real(dp)::aexp
!------------------------------------------------------------------------
  UVrates=0.
  if(.not. haardt_madau) RETURN
  
  call inp_UV_rates_table(1./aexp - 1., UVrates, .true.)

  !if(myid==1) then
  !   write(*,*) 'The UV rates have changed to:'
  !   do i=1,nIons
  !      write(*,910) UVrates(i,:)
  !   enddo
  !endif
910 format (1pe21.6, ' s-1', 1pe21.6,' erg s-1')

END SUBROUTINE update_UVrates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_solve_cooling(T2, xion, Np, Fp, p_gas, dNpdt, dFpdt        &
                           ,nH, c_switch, Zsolar, dt, a_exp, nCell)
! Semi-implicitly solve for new temperature, ionization states, 
! photon density/flux, and gas velocity in a number of cells.
! Parameters: 
! T2     <=> T/mu [K] 
! xion   <=> NION ionization fractions 
! Np     <=> NGROUPS photon number densities [cm-3]
! Fp     <=> NGROUPS * ndim photon number fluxes [cm-2 s-1]
! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
! dNpdt   =>  Op split increment in photon densities during dt
! dFpdt   =>  Op split increment in photon flux magnitudes during dt
! nH      =>  Hydrogen number densities [cm-3]
! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
! Zsolar  =>  Cell metallicities [solar fraction]
! dt      =>  Timestep size             [s]
! a_exp   =>  Cosmic expansion
! nCell   =>  Number of cells (length of all the above vectors)
!
! We use a slightly modified method of Anninos et al. (1997).
!-------------------------------------------------------------------------
  use amr_commons
  implicit none  
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:nIons, 1:nvector):: xion
  real(dp),dimension(1:nGroups, 1:nvector):: Np, dNpdt
  real(dp),dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
  real(dp),dimension(1:ndim, 1:nvector):: p_gas
  real(dp),dimension(1:nvector):: nH, Zsolar
  logical,dimension(1:nvector):: c_switch
  real(dp)::dt, a_exp
  integer::ncell !--------------------------------------------------------
  real(dp),dimension(1:nvector):: tLeft, ddt
  logical:: dt_ok
  real(dp)::dt_rec
  real(dp):: dT2
  real(dp),dimension(nIons):: dXion
  real(dp),dimension(nGroups):: dNp
  real(dp),dimension(1:ndim, 1:nGroups):: dFp
  real(dp),dimension(1:ndim):: dp_gas
  integer::i, ia, ig,  nAct, nAct_next, loopcnt, code
  integer,dimension(1:nvector):: indAct              ! Active cell indexes
  real(dp)::one_over_rt_c_cgs, one_over_egy_IR_erg, one_over_x_FRAC
  real(dp)::one_over_Np_FRAC, one_over_Fp_FRAC, one_over_T_FRAC
  real(dp),dimension(1:nGroups) :: group_egy_ratio, group_egy_erg
  real(dp):: mu

  ! Store some temporary variables reduce computations
  one_over_rt_c_cgs = 1d0 / rt_c_cgs
  one_over_Np_FRAC = 1d0 / Np_FRAC
  one_over_Fp_FRAC = 1d0 / Fp_FRAC
  one_over_T_FRAC = 1d0 / T_FRAC
  one_over_x_FRAC = 1d0 / x_FRAC
#if NGROUPS>0 
  if(rt .and. nGroups .gt. 0) then 
     group_egy_erg(1:nGroups) = group_egy(1:nGroups) * ev_to_erg
     if(rt_isIR) then
        group_egy_ratio(1:nGroups) = group_egy(1:nGroups) / group_egy(iIR)
        one_over_egy_IR_erg = 1.d0 / group_egy_erg(iIR)
     endif
  endif
#endif
  !-----------------------------------------------------------------------
  tleft(1:ncell) = dt                !       Time left in dt for each cell
  ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths

  do i=1,ncell
     indact(i) = i                   !      Set up indexes of active cells
     ! Ensure all state vars are legal:
     T2(i) = MAX(T2(i), T2_min_fix)
     ! HACK - if a high temperature, set fully ionised
     if (T2(i) .gt. 1e5) then
        xion(1,i) = 1d0
        xion(2,i) = 1d-6
        xion(3,i) = 1d0
     endif
     ! END HACK
     xion(1:nIons,i) = MIN(MAX(xion(1:nIons,i), x_MIN),1.d0)
     if(xion(2,i)+xion(3,i) .gt. 1.d0) then
        if(xion(2,i) .gt. xion(3,i)) then
           xion(2,i)=1.d0-xion(3,i)
        else
           xion(3,i)=1.d0-xion(2,i)
        endif
     endif
     if(rt) then
        do ig=1,ngroups
           Np(ig,i) = MAX(smallNp, Np(ig,i))
           call reduce_flux(Fp(:,ig,i),Np(ig,i)*rt_c_cgs)
        end do
     endif
  end do

  ! Loop until all cells have tleft=0
  ! **********************************************
  nAct=nCell                                      ! Currently active cells
  loopcnt=0 ; n_cool_cells=n_cool_cells+nCell     !             Statistics
  do while (nAct .gt. 0)      ! Iterate while there are still active cells
     loopcnt=loopcnt+1   ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct 
     nAct_next=0                     ! Active cells for the next iteration
     do ia=1,nAct                             ! Loop over the active cells
        i = indAct(ia)                        !                 Cell index
        call cool_step(i)
        if(loopcnt .gt. 100000) then
           call display_coolinfo(.true., loopcnt, i, dt-tleft(i), dt     &
                            ,ddt(i), nH(i), T2(i),  xion(:,i),  Np(:,i)  &
                            ,Fp(:,:,i),  p_gas(:,i)                      &
                            ,dT2, dXion, dNp, dFp, dp_gas, code)
        endif
        if(.not. dt_ok) then
           ddt(i)=ddt(i)/2.                    ! Try again with smaller dt 
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
           loopCodes(code) = loopCodes(code)+1
           cycle 
        endif
        ! Update the cell state (advance the time by ddt):
        T2(i)     = T2(i) + dT2
        xion(:,i) = xion(:,i) + dXion(:)
        if(nGroups .gt. 0) then 
           Np(:,i)   = Np(:,i) + dNp(:)
           Fp(:,:,i) = Fp(:,:,i) + dFp(:,:)
        endif
        p_gas(:,i)   = p_gas(:,i) + dp_gas(:)

        tleft(i)=tleft(i)-ddt(i)
        if(tleft(i) .gt. 0.) then           ! Not finished with this cell
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
        else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
           print*,'In rt_solve_cooling: tleft < 0  !!'
           stop
        endif
        ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from cool_step
     end do ! end loop over active cells
     nAct=nAct_next
  end do ! end iterative loop

  ! loop statistics
  max_cool_loopcnt=max(max_cool_loopcnt,loopcnt)
contains
  
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE cool_step(icell)
  ! Compute change in cell state in timestep ddt(icell), or set in dt_rec
  ! a recommendation for new timestep if ddt(icell) proves too large.
  ! T2      => T/mu [K]                                -- dT2 is new value
  ! xion    => NION ionization fractions               --     dXion is new
  ! Np      => NGROUPS photon number densities [cm-3]  -- dNp is new value
  ! Fp      => NGROUPS * ndim photon fluxes [cm-2 s-1] -- dFp is new value 
  ! p_gas   => ndim gas momenta [cm s-1 g cm-3]        --    dp_gas is new
  ! dNpdt   =>  Op split increment in photon densities during dt
  ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
  ! nH      =>  Hydrogen number densities [cm-3]  
  ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat)
  ! Zsolar  =>  Cell metallicities [solar fraction]
  ! dt      =>  Timestep size [s]
  ! a_exp   =>  Cosmic expansion
  ! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
  ! dt_rec  <=  Recommended timesteps for next iteration
  ! code    <= Error code in cool step, if dt_ok=.f.
  !
  ! The original values, T2, xion etc, must stay unchanged, while dT2,
  ! dxion etc contain the new values (the difference at the end of the
  ! routine).
  !-----------------------------------------------------------------------
    use amr_commons
    use rt_metal_cooling_module
    use const
    implicit none  
    integer, intent(in)::icell
    real(dp),dimension(nDim),save:: dmom
    real(dp),dimension(nDim), save:: u_gas ! Gas velocity
    real(dp),dimension(nIons),save:: alpha, beta, nN, nI
    real(dp),save:: dUU, fracMax
    real(dp),save:: xHeI, mu, TK, nHe, ne, neInit, Hrate, dAlpha, dBeta
    real(dp),save:: s, jac, q, Crate, dCdT2, X_nHkb, rate, dRate, cr, de
    real(dp),save:: photoRate, metal_tot,metal_prime, ss_factor
    integer,save:: iion,igroup,idim
    real(dp),dimension(nGroups),save:: recRad, phAbs, phSc, dustAbs
    real(dp),dimension(nGroups),save:: dustSc, kAbs_loc,kSc_loc
    real(dp),save:: rho, TR, one_over_C_v, E_rad, dE_T, fluxMag, mom_fact
    !---------------------------------------------------------------------
    dt_ok=.false.
    nHe=0.25*nH(icell)*Y/X  !         Helium number density
    ! U contains the original values, dU the updated ones:
    dT2=T2(icell) ; dXion(:)=xion(:,icell) ; dNp(:)=Np(:,icell)
    dFp(:,:)=Fp(:,:,icell) ; dp_gas(:)=p_gas(:,icell)
    ! xHI = MAX(1.-dXion(1),0.) ; xHII = dXion(1)
    ! xHeII=dXion(2) ; xHeIII=dXion(3)
    xHeI=MAX(1.-dXion(2)-dXion(3),0.d0)
    ! nN='neutral' species (pre-ionized), nI=their ionized counterparts
    nN(1)  = nH(icell) * (1.d0-dXion(1))                         !     nHI
    nN(2)  = nHe*xHeI                                            !    nHeI
    nN(3)  = nHe*dXion(2)                                        !   nHeII
    nI(1)  = nH(icell) *dXion(1)                                 !    nHII
    nI(2)  = nN(3)                                               !   nHeII
    nI(3)  = nHe*dXion(3)                                        !  nHeIII
    mu = getMu(dXion(1), dXion(2), dXion(3), dT2)
    TK = dT2 * mu                                           !  Temperature
    if(rt_isTconst) TK=rt_Tconst                       !  Force constant T
    ne= nH(icell)*dXion(1)+nHE*(dXion(2)+2.*dXion(3))  !  Electron density
    neInit=ne
    fracMax=0d0   ! Max fractional update, to check if dt can be increased
    ss_factor=1d0                    ! UV background self_shielding factor
    if(self_shielding) ss_factor = exp(-nH(icell)/1d-2)

    rho = nH(icell) / X * mH
#if NGROUPS>0 
    ! Set dust opacities--------------------------------------------------
    if(rt .and. nGroups .gt. 0) then
       kAbs_loc = kappaAbs
       kSc_loc  = kappaSc
       if(is_kIR_T) then ! k_IR depends on T
          ! Special stuff for Krumholz/Davis experiment
          if(rt_T_rad) then  ! Use radiation temperature for kappa
             E_rad = group_egy_erg(iIR) * dNp(iIR)
             TR = max(T2_min_fix,(E_rad*rt_c_fraction/a_r)**0.25)
             dT2 = TR/mu ;   TK = TR
          endif
          kAbs_loc(iIR) = kappaAbs(iIR) * (TK/10d0)**2
          kSc_loc(iIR)  = kappaSc(iIR)  * (TK/10d0)**2
       endif
       ! Set dust absorption and scattering rates [s-1]:
       dustAbs(:)  = kAbs_loc(:) *rho*Zsolar(icell)*rt_c_cgs
       dustSc(iIR) = kSc_loc(iIR)*rho*Zsolar(icell)*rt_c_cgs
    endif

    !(i) UPDATE PHOTON DENSITY AND FLUX **********************************
    if(rt .and. rt_advect) then 
       recRad(1:nGroups)=0. ; phAbs(1:nGroups)=0.              
       ! Scattering rate; reduce the photon flux, but not photon density:
       phSc(1:nGroups)=0.

       ! EMISSION FROM GAS
       if(.not. rt_OTSA .and. rt_advect) then ! ----------- Rec. radiation
          alpha(1) = inp_coolrates_table(tbl_alphaA_HII, TK) &
                   - inp_coolrates_table(tbl_alphaB_HII, TK)
          ! alpha(2) A-B becomes negative around 1K, hence the max
          alpha(2) = MAX(0.d0,  inp_coolrates_table(tbl_alphaA_HeII, TK) &
                              - inp_coolrates_table(tbl_alphaB_HeII, TK))
          alpha(3) = inp_coolrates_table(tbl_alphaA_HeIII, TK) &
                   - inp_coolrates_table(tbl_alphaB_HeIII, TK)
          do iion=1,nIons
             if(spec2group(iion) .gt. 0) &  ! Contribution of ion -> group
                  recRad(spec2group(iion)) = &
                  recRad(spec2group(iion)) + alpha(iion) * nI(iion) * ne
          enddo
       endif

       ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
       do igroup=1,nGroups      ! -------------------Ionization absorbtion
          phAbs(igroup) = SUM(nN(:)*signc(igroup,:)) ! s-1
       end do
       ! IR, optical and UV depletion by dust absorption: ----------------
       if(rt_isIR) & !IR scattering/abs on dust (abs after T update)        
            phSc(iIR)  = phSc(iIR) + dustSc(iIR)                        
       do igroup=1,nGroups        ! Deplete photons, since they go into IR
          if( .not. (rt_isIR .and. igroup.eq.iIR) ) &  ! IR done elsewhere
               phAbs(igroup) = phAbs(igroup) + dustAbs(igroup)
       end do

       dmom(1:nDim)=0d0
       do igroup=1,nGroups  ! ------------------- Do the update of N and F
          dNp(igroup)= MAX(smallNp,                                      &
                        (ddt(icell)*(recRad(igroup)+dNpdt(igroup,icell)) &
                                    +dNp(igroup))                        &
                        / (1.d0+ddt(icell)*phAbs(igroup)))

          dUU = ABS(dNp(igroup)-Np(igroup,icell))                        &
                /(Np(igroup,icell)+Np_MIN) * one_over_Np_FRAC
          if(dUU .gt. 1.d0) then
             code=1 ;   RETURN                        ! ddt(icell) too big
          endif
          fracMax=MAX(fracMax,dUU)      ! To check if ddt can be increased

          do idim=1,nDim
             dFp(idim,igroup) = &
                  (ddt(icell)*dFpdt(idim,igroup,icell)+dFp(idim,igroup)) &
                  /(1d0+ddt(icell)*(phAbs(igroup)+phSc(igroup)))
          end do
          call reduce_flux(dFp(:,igroup),dNp(igroup)*rt_c_cgs)

          do idim=1,nDim
             dUU = ABS(dFp(idim,igroup)-Fp(idim,igroup,icell))           &
                  / (ABS(Fp(idim,igroup,icell))+Fp_MIN) * one_over_Fp_FRAC
             if(dUU .gt. 1.d0) then
                code=2 ;   RETURN                     ! ddt(icell) too big
             endif
             fracMax=MAX(fracMax,dUU)   ! To check if ddt can be increased
          end do

       end do

       do igroup=1,nGroups ! -------Momentum transfer from photons to gas:
          mom_fact = ddt(icell) * (phAbs(igroup) + phSc(igroup)) &
               * group_egy_erg(igroup) * one_over_c_cgs

          if(rt_isoPress .and. .not. (rt_isIR .and. igroup==iIR)) then 
             ! rt_isoPress: assume f=1, where f is reduced flux.
             fluxMag=sqrt(sum((dFp(:,igroup))**2))
             if(fluxMag .gt. 0d0) then
                mom_fact = mom_fact * dNp(igroup) / fluxMag
             else
                mom_fact = 0d0
             endif
          else
             mom_fact = mom_fact * one_over_rt_c_cgs 
          end if

          do idim = 1, nDim
             dmom(idim) = dmom(idim) + dFp(idim,igroup) * mom_fact
          end do
       end do
       dp_gas = dp_gas + dmom * rt_pressBoost        ! update gas momentum

       ! Add absorbed UV/optical energy to IR:----------------------------  
       if(rt_isIR) then   
          do igroup=iIR+1,nGroups
             dNp(iIR) = dNp(iIR) + dustAbs(igroup) * ddt(icell)          &
                  * dNp(igroup) * group_egy_ratio(igroup)
          end do
       endif
       ! -----------------------------------------------------------------
    endif !if(rt)
#endif
    !(ii) UPDATE TEMPERATURE *********************************************
    if(c_switch(icell) .and. cooling .and. .not. rt_T_rad) then
       Hrate=0.                             !  Heating rate [erg cm-3 s-1]
       if(rt .and. rt_advect) then
          do igroup=1,nGroups                              !  Photoheating
             Hrate = Hrate + dNp(igroup) * SUM(nN(:) * PHrate(igroup,:))
          end do
       endif
       if(haardt_madau) Hrate= Hrate + SUM(nN(:)*UVrates(:,2)) * ss_factor
       Crate = compCoolrate(TK, ne, nN(1), nI(1), nN(2), nN(3), nI(3)    &
            ,a_exp, dCdT2, RT_OTSA)                  ! Cooling
       dCdT2 = dCdT2 * mu                            ! dC/dT2 = mu * dC/dT
       metal_tot=0.d0 ; metal_prime=0.d0             ! Metal cooling
       Zsolar(icell) = 1d0 ! HACK
       !if(Zsolar(icell) .gt. 0d0) &
       !     call rt_cmp_metals(T2(icell),nH(icell),mu,metal_tot          &
       !                       ,metal_prime,a_exp)
       call rt_metal_cool(T2(icell),nH(icell),dXion(1),mu,metal_tot,metal_prime)
       X_nHkb= X/(1.5 * nH(icell) * kB)            ! Multiplication factor   
       rate  = X_nHkb*(Hrate - Crate - Zsolar(icell)*metal_tot)
       dRate = -X_nHkb*(dCdT2 + Zsolar(icell)*metal_prime)     ! dRate/dT2
       ! 1st order dt constr
       dUU   = ABS(MAX(T2_min_fix, T2(icell)+rate*ddt(icell))-T2(icell))
       ! New T2 value 
       dT2   = MAX(T2_min_fix &
                  ,T2(icell)+rate*ddt(icell)/(1.-dRate*ddt(icell)))
       dUU   = MAX(dUU, ABS(dT2-T2(icell))) / (T2(icell)+T_MIN) &
                        *one_over_T_FRAC
       if(dUU .gt. 1.) then                                     ! 10% rule
          code=3 ; RETURN
       endif
       fracMax=MAX(fracMax,dUU)
       TK=dT2*mu
    endif

#if NGROUPS>0 
    if(rt_isIR) then
       if(kAbs_loc(iIR) .gt. 0d0 .and. .not. rt_T_rad) then
          ! Evolve IR-Dust equilibrium temperature------------------------
          ! Delta (Cv T)= ( c_red/lambda E - c/lambda a T^4) 
          !           / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
          one_over_C_v = mh*mu*(gamma-1d0) / (rho*kb)
          E_rad = group_egy_erg(iIR) * dNp(iIR)
          dE_T = (rt_c_cgs * E_rad - c_cgs*a_r*TK**4)                    &
               /(1d0/(kAbs_loc(iIR) * Zsolar(icell) * rho * ddt(icell))  &
               +4d0*c_cgs * one_over_C_v *a_r*TK**3+rt_c_cgs)
          dT2 = dT2 + 1d0/mu * one_over_C_v * dE_T
          dNp(iIR) = dNp(iIR) - dE_T * one_over_egy_IR_erg

          dT2 = max(T2_min_fix,dT2)                                   
          dNp(iIR) = max(dNp(iIR), smallNp)
          ! 10% rule for photon density:
          dUU = ABS(dNp(iIR)-Np(iIR,icell)) / (Np(iIR,icell)+Np_MIN)     &
                                            * one_over_Np_FRAC
          if(dUU .gt. 1.) then                 
             code=4 ;   RETURN                          
          endif
          fracMax=MAX(fracMax,dUU)                                           

          dUU   = ABS(dT2-T2(icell)) / (T2(icell)+T_MIN) * one_over_T_FRAC
          if(dUU .gt. 1.) then                           ! 10% rule for T2
             code=5 ; RETURN                                                  
          endif
          fracMax=MAX(fracMax,dUU)
          TK=dT2*mu
          call reduce_flux(dFp(:,iIR),dNp(iIR)*rt_c_cgs)           
       endif
    endif
#endif
    !(iii) UPDATE xHII****************************************************
    ! First recompute interaction rates since T is updated
    if(rt_OTSA .or. .not. rt_advect) then           !  Recombination rates
       alpha(1) = inp_coolrates_table(tbl_alphaB_HII, TK, dalpha)
    else                               
       alpha(1) = inp_coolrates_table(tbl_alphaA_HII, TK, dalpha)
    endif
    beta(1) = inp_coolrates_table(tbl_beta_HI, TK, dBeta) !  Coll-ion rate
    cr = beta(1) * ne                             !               Creation
    if(rt) cr = cr + SUM(signc(:,1)*dNp)          !                  [s-1]
    if(haardt_madau) cr = cr + UVrates(1,1) * ss_factor
    de = alpha(1) * ne                            !            Destruction

    ! Not Anninos, but more stable (this IS neccessary, as the one-cell  !
    ! tests oscillate wildly in the Anninos method):                     !
    S  = cr*(1.-dXion(1))-de*dXion(1)
    dUU= ABS(MIN(MAX(dXion(1)+ddt(icell)*S, x_MIN), 1.)-dXion(1))
    jac=(1.-dXion(1))*(beta(1)*nH(icell)-ne*TK*mu*X*dBeta) & !jac=dS/dxHII
         - cr - de - dXion(1) * (alpha(1)*nH(icell)-ne*TK*mu*X*dAlpha)
    dXion(1) = xion(1,icell)                                             &
             + ddt(icell)*(cr*(1.-xion(1,icell))-de*xion(1,icell))       &
             / (1.-ddt(icell)*jac)
    dXion(1) = MIN(MAX(dXion(1), x_MIN),1.d0)
    dUU = MAX(dUU, ABS(dXion(1)-xion(1,icell))) / (xion(1,icell)+x_MIN)  &
                                                * one_over_x_FRAC
    if(dUU .gt. 1.) then
       code=6 ; RETURN
    endif
    fracMax=MAX(fracMax,dUU)
    !End a more stable and accurate integration---------------------------
    if(isHe) then
       ne= nH(icell)*dXion(1)+nHE*(dXion(2)+2.*dXion(3)) ! Bc changed xhii
       mu = getMu(dXion(1), dXion(2), dXion(3), dT2)
       if(.not. rt_isTconst) TK=dT2*mu !  Update TK because of changed  mu

       !(iv) UPDATE xHeI *************************************************
       if(rt_OTSA .or. .not. rt_advect) then
          alpha(2) = inp_coolrates_table(tbl_alphaB_HeII, TK)
          alpha(3) = inp_coolrates_table(tbl_alphaB_HeIII, TK)
       else                               
          alpha(2) = inp_coolrates_table(tbl_alphaA_HeII, TK)
          alpha(3) = inp_coolrates_table(tbl_alphaA_HeIII, TK)
       endif
       beta(2) =  inp_coolrates_table(tbl_beta_HeI, TK)
       beta(3) = inp_coolrates_table(tbl_beta_HeII, TK)
       ! Creation = recombination of HeII and electrons
       cr = alpha(2) * ne * dXion(2)
       ! Destruction = collisional ionization+photoionization of HeI
       de = beta(2) * ne
       if(rt) de = de + SUM(signc(:,2)*dNp)
       if(haardt_madau) de = de + UVrates(2,1) * ss_factor
       xHeI = (cr*ddt(icell)+xHeI)/(1.+de*ddt(icell))        !  The update
       xHeI = MIN(MAX(xHeI, 0.),1.)

       !(v) UPDATE xHeII *************************************************
       ! Creation = coll.- and photo-ionization of HI + rec. of HeIII
       cr = de * xHeI + alpha(3) * ne * dXion(3)
       ! Destruction = rec. of HeII + coll.- and photo-ionization of HeII
       photoRate=0.
       if(rt) photoRate = SUM(signc(:,3)*dNp)
       if(haardt_madau) photoRate = photoRate + UVrates(3,1) * ss_factor
       de = (alpha(2) + beta(3)) * ne + photoRate
       dXion(2) = (cr*ddt(icell)+dXion(2))/(1.+de*ddt(icell)) ! The update
       dXion(2) = MIN(MAX(dXion(2), x_MIN),1.)

       !(vii) UPDATE xHeIII **********************************************
       ! Creation = coll.- and photo-ionization of HeII
       cr = (beta(3) * ne + photoRate) * dXion(2)          !  xHeII is new
       ! Destruction = rec. of HeIII and e
       de = alpha(3) * ne
       dXion(3) = (cr*ddt(icell)+dXion(3))/(1.+de*ddt(icell)) ! The update
       dXion(3) = MIN(MAX(dXion(3), x_MIN),1.)

       !(viii) ATOMIC CONSERVATION OF He *********************************
       if(xHeI .ge. dXion(3)) then   ! Either HeI or HeII is most abundant 
          if(xHeI .le. dXion(2)) dXion(2) = 1.-xHeI-dXion(3) !HeII most ab
       else                        ! Either HeII or HeIII is most abundant 
          if(dXion(2) .le. dXion(3)) then
             dXion(3) = 1. - xHeI-dXion(2)                         ! HeIII
          else
             dXion(2) = 1. - xHeI-dXion(3)                         !  HeII
          endif
       endif
    endif

    ne = nH(icell)*dXion(1)+nHe*(dXion(2)+2.*dXion(3))
    dUU=ABS((ne-neInit)) / (neInit+x_MIN) * one_over_x_FRAC
    if(dUU .gt. 1.) then
       code=7 ; RETURN
    endif
    fracMax=MAX(fracMax,dUU)

    if(rt_isTconst) dT2=rt_Tconst/mu

    dT2 = dT2-T2(icell) ; dXion(:) = dXion(:)-xion(:,icell)
    dNp(:) = dNp(:)-Np(:,icell) ; dFp(:,:) = dFp(:,:)-Fp(:,:,icell)
    dp_gas(:)= dp_gas(:)-p_gas(:,icell)
    ! Now the dUs are really changes, not new values
    !(ix) Check if we are safe to use a bigger timestep in next iteration:
    if(fracMax .lt. 0.5) then
       dt_rec=ddt(icell)*2.
    else
       dt_rec=ddt(icell)
    endif
    dt_ok=.true.
    code=0

  END SUBROUTINE cool_step

END SUBROUTINE rt_solve_cooling

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE display_coolinfo(stopRun, loopcnt, i, dtDone, dt, ddt, nH    &
                            ,T2,  xion,  Np,  Fp,  p_gas                &
                            ,dT2, dXion, dNp, dFp, dp_gas, code)
! Print cooling information to standard output, and maybe stop execution.
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  real(dp),dimension(nIons):: xion, dXion
  real(dp),dimension(nGroups):: Np, dNp
  real(dp),dimension(nDim, nGroups):: Fp, dFp
  real(dp),dimension(nDim):: p_gas, dp_gas
  real(dp)::T2, dT2, dtDone, dt, ddt, nH
  logical::stopRun
  integer::loopcnt,i, code
!------------------------------------------------------------------------
  if(stopRun) write(*, 111) loopcnt
  if(.true.) then
     write(*,900) loopcnt, myid, code, i, dtDone, dt, ddt, rt_c_cgs, nH
     write(*,901) T2,      xion,      Np,      Fp,      p_gas
     write(*,902) dT2,     dXion,     dNp,     dFp,     dp_gas
     write(*,903) dT2/ddt, dXion/ddt, dNp/ddt, dFp/ddt, dp_gas/ddt
     write(*,904) abs(dT2)/(T2+T_MIN), abs(dxion)/(xion+x_MIN),          &
                  abs(dNp)/(Np+Np_MIN), abs(dFp)/(Fp+Fp_MIN)
  endif
  print*,loopcodes
  print*,group_egy(:)
  if(stopRun) then
     print *,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
     STOP
  endif

111 format(' Stopping because of large number of timestesps in', &
           ' rt_solve_cooling (', I6, ')')
900 format (I3, '  myid=', I2, ' code=', I2, ' i=', I5, ' t=', 1pe12.3,xs&
            '/', 1pe12.3, ' ddt=', 1pe12.3, ' c=', 1pe12.3, &
            ' nH=', 1pe12.3)
901 format ('  U      =', 20(1pe12.3))
902 format ('  dU     =', 20(1pe12.3))
903 format ('  dU/dt  =', 20(1pe12.3))
904 format ('  dU/U % =', 20(1pe12.3))
END SUBROUTINE display_coolinfo

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_chem_eq(TK, nH, t_rad_spec, nSpec, nTot, mu)

! Compute chemical equilibrium abundances of e, HI, HII, HeI, HeII, HeIII
! r_rad_spec => photoionization rates [s-1] for HI, HeI, HeII
!------------------------------------------------------------------------
  implicit none
  real(dp),intent(in)::TK, nH
  real(dp),intent(out)::nTot, mu
  real(dp),dimension(1:3),intent(in)::t_rad_spec
  real(dp),dimension(1:6),intent(out)::nSpec!------------------------
  real(dp)::xx, yy
  real(dp)::n_HI, n_HII, n_HEI, n_HEII, n_HEIII, n_E
  real(dp)::t_rad_HI,  t_rad_HEI,  t_rad_HEII
  real(dp)::t_rec_HI,  t_rec_HEI,  t_rec_HEII
  real(dp)::t_ion_HI,  t_ion_HEI,  t_ion_HEII
  real(dp)::t_ion2_HI, t_ion2_HEI, t_ion2_HEII
  real(dp)::x1, err_nE
  integer,parameter::HI=1, HeI=2, HeII=3
!------------------------------------------------------------------------
  xx=(1.-Y)
  yy=Y/(1.-Y)/4.
  
  t_rad_HI   = t_rad_spec(HI)                !      Photoionization [s-1]
  t_rad_HEI  = t_rad_spec(HeI)
  t_rad_HEII = t_rad_spec(HeII)

  if(rt_OTSA) then                           !    Recombination [cm3 s-1]
     t_rec_HI   = inp_coolrates_table(tbl_alphaB_HII, TK)
     t_rec_HEI  = inp_coolrates_table(tbl_alphaB_HeII, TK)
     t_rec_HEII = inp_coolrates_table(tbl_alphaB_HeIII, TK)
  else 
     t_rec_HI   = inp_coolrates_table(tbl_alphaA_HII, TK)
     t_rec_HEI  = inp_coolrates_table(tbl_alphaA_HeII, TK)
     t_rec_HEII = inp_coolrates_table(tbl_alphaA_HeIII, TK)
  endif

  t_ion_HI   = inp_coolrates_table(tbl_beta_HI, TK) ! Coll. ion. [cm3 s-1]
  t_ion_HEI  = inp_coolrates_table(tbl_beta_HeI, TK)
  t_ion_HEII = inp_coolrates_table(tbl_beta_HeII, TK)
  
  n_E = nH        
  err_nE = 1.
  
  do while(err_nE > 1.d-8)
     t_ion2_HI   = t_ion_HI   + t_rad_HI  /MAX(n_E,1e-15*nH)  ! [cm3 s-1]
     t_ion2_HEI  = t_ion_HEI  + t_rad_HEI /MAX(n_E,1e-15*nH)
     t_ion2_HEII = t_ion_HEII + t_rad_HEII/MAX(n_E,1e-15*nH)
     
     n_HI  = t_rec_HI/(t_ion2_HI+t_rec_HI)*nH
     n_HII = t_ion2_HI/(t_ion2_HI+t_rec_HI)*nH
     
     x1=(                                                                &
          t_rec_HEII*t_rec_HEI                                           &
          +t_ion2_HEI*t_rec_HEII                                         &
          +t_ion2_HEII*t_ion2_HEI)                               ! cm6 s-2
     
     n_HEIII = yy*t_ion2_HEII*t_ion2_HEI/x1*nH
     n_HEII  = yy*t_ion2_HEI *t_rec_HEII/x1*nH
     n_HEI   = yy*t_rec_HEII *t_rec_HEI /x1*nH
     
     err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/nH)
     n_E = 0.5*n_E+0.5*(n_HII + n_HEII + 2.*n_HEIII)
     
  end do
    
  nTOT     = n_E+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       = nH/xx/nTOT
  nSpec(1) = n_E
  nSpec(2) = n_HI
  nSpec(3) = n_HII
  nSpec(4) = n_HEI
  nSpec(5) = n_HEII
  nSpec(6) = n_HEIII
  
END SUBROUTINE cmp_chem_eq

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE rt_evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL &
                           ,J0min_in,T2end,mu,ne,if_write_result)
!-------------------------------------------------------------------------
! Used for initialization of thermal state in cosmological simulations.
!
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100 
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!-------------------------------------------------------------------------
  use amr_commons,only:myid
  use UV_module
  implicit none
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,J0min_in,ne,dasura
  logical :: if_write_result
  real(dp)::aexp,daexp,dt_cool,coeff,T2_com, nH_com  
  real(dp),dimension(nIons)::pHI_rates=0., h_rad_spec=0.
  real(kind=8) ::mu
  real(dp) ::cool_tot,heat_tot, mu_dp,diff
  integer::niter
  real(dp) :: n_spec(1:6)
  real(dp),dimension(1:nvector):: T2
  real(dp),dimension(1:nIons, 1:nvector):: xion
  real(dp),dimension(1:nGroups, 1:nvector):: Np, dNpdt
  real(dp),dimension(1:ndim, 1:nGroups, 1:nvector):: Fp, dFpdt
  real(dp),dimension(1:ndim, 1:nvector):: p_gas
  real(dp),dimension(1:nvector)::nH=0., Zsolar=0.
  logical,dimension(1:nvector)::c_switch=.true.
!-------------------------------------------------------------------------
  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2 / mu_mol
  nH_com = omegab*rhoc*h**2*X/mH

  mu_dp=mu
  call cmp_Equilibrium_Abundances(                                       &
                 T2_com/aexp**2, nH_com/aexp**3, pHI_rates, mu_dp, n_Spec)
  ! Initialize cell state
  T2(1)=T2_com                                          !      Temperature
  xion(1,1)=n_Spec(3)/(nH_com/aexp**3)                  !   HII   fraction
  xion(2,1)=n_Spec(5)/(nH_com/aexp**3)                  !   HeII  fraction
  xion(3,1)=n_Spec(6)/(nH_com/aexp**3)                  !   HeIII fraction
  p_gas(:,1)=0.
  Np(:,1)=0. ; Fp(:,:,1)=0.                  ! Photon densities and fluxes
  dNpdt(:,1)=0. ; dFpdt(:,:,1)=0.                              
  do while (aexp < aend)
     call update_UVrates(aexp)
     call update_coolrates_tables(aexp)

     daexp = dasura*aexp
     dt_cool = daexp                                                     &
             / (aexp*100.*h*3.2408608e-20)                               &
             / HsurH0(1.0/dble(aexp)-1.,omega0,omegaL,1.-omega0-omegaL)
     
     nH(1) = nH_com/aexp**3
     T2(1) = T2(1)/aexp**2
     call rt_solve_cooling(T2,xion,Np,Fp,p_gas,dNpdt,dFpdt,nH,c_switch   &
                           ,Zsolar,dt_cool,aexp,1)
     T2(1)=T2(1)*aexp**2
     aexp = aexp + daexp
     if (if_write_result) write(*,'(4(1pe10.3))')                        &
                              aexp,nH(1),T2_com*mu/aexp**2,n_spec(1)/nH(1)
  end do
  T2end=T2(1)/(aexp-daexp)**2
  ne=(n_spec(3)+(n_spec(5)+2.*n_spec(6))*0.25*Y/X)
end subroutine rt_evol_single_cell

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION HsurH0(z,omega0,omegaL,OmegaR)
!-------------------------------------------------------------------------
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
!-------------------------------------------------------------------------
  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
END FUNCTION HsurH0

! NOTE - rt_cmp_metals HACKED FROM HERE TO rt_metal_cooling_module
!        TO PREVENT CIRCULAR MODULE USE

!*************************************************************************
FUNCTION getMu(xHII, xHeII, xHeIII, Tmu)
! Returns the mean particle mass, in units of the proton mass.
! xHII, xHeII, xHeIII => Hydrogen and helium ionisation fractions
! Tmu => T/mu in Kelvin  
!-------------------------------------------------------------------------
  implicit none
  real(kind=8),intent(in) :: xHII, xHeII, xHeIII, Tmu
  real(kind=8) :: mu
  real(kind=8) :: getMu
!-------------------------------------------------------------------------
  getMu = 1./(X*(1.+xHII) + 0.25*Y*(1.+xHeII+2.*xHeIII))   
  if(is_kIR_T .or. is_mu_H2) &
       getMu = getMu + exp(-1.d0*(Tmu/Tmu_dissoc)**2) * (2.33-getMu)
END FUNCTION getMu


END MODULE rt_cooling_module

!************************************************************************
SUBROUTINE updateRTGroups_CoolConstants()
! Update photon group cooling and heating constants, to reflect an update
! in rt_c_cgs and in the cross-sections and energies in the groups.
!------------------------------------------------------------------------
  use rt_cooling_module
  use rt_parameters
  implicit none
  integer::iP, iI
!------------------------------------------------------------------------
  signc=group_csn*rt_c_cgs                                    ! [cm3 s-1]
  sigec=group_cse*rt_c_cgs                                    ! [cm3 s-1]
  do iP=1,nGroups
     do iI=1,nIons               ! Photoheating rates for photons on ions
        PHrate(iP,iI) =  ev_to_erg * &        ! See eq (19) in Aubert(08)
             (sigec(iP,iI) * group_egy(iP) - signc(iP,iI)*ionEvs(iI))
        PHrate(iP,iI) = max(PHrate(iP,iI),0d0) !      No negative heating
     end do
  end do
END SUBROUTINE updateRTGroups_CoolConstants

!************************************************************************
SUBROUTINE reduce_flux(Fp, cNp)
! Make sure the reduced photon flux is less than one
!------------------------------------------------------------------------
  use rt_parameters
  implicit none
  real(dp),dimension(ndim):: Fp
  real(dp):: cNp, fred
!------------------------------------------------------------------------
  fred = sqrt(sum(Fp**2))/cNp
  if(fred .gt. 1.d0) Fp = Fp/fred
END SUBROUTINE reduce_flux


