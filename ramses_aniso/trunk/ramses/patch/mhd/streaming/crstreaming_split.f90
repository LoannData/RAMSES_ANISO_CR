! ---------------------------------------------------------------
!  COND_SPLIT  This routine solves the heat flux following the
!              anistropic heat conduction.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho, \rho u, \rho v, \rho w, Etot, A, B, C)
!  the hydro variable are cell-centered
!  whereas the magnetic field B=(A,B,C) are face-centered.
!  Note that here we have 3 components for v and B whatever ndim.
!
!  This routine was written by Yohan Dubois & Benoit Commerçon
! ----------------------------------------------------------------
subroutine crstreaming_split(uin,flux,dx,dy,dz,dt,ngrid,fdx)
  use amr_parameters
  use const
  use hydro_parameters
!!$  use pm_commons, ONLY: localseed,iseed
!!$  use amr_commons, ONLY:ncpu,myid
  implicit none

  integer ::ngrid,icr
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::fdx

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr)::Xflux,Yflux,Zflux

  ! Primitive variables

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::jsat

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

!!$
!!$! ------ Importation des variables depuis uloc == uin ici ------
!!$  do k = ku1, ku2
!!$  do j = ju1, ju2
!!$  do i = iu1, iu2
!!$     do l = 1, ngrid
!!$
!!$     end do
!!$  end do
!!$  end do
!!$  end do
! ---------------------------------------------------------------
 flux = 0.0d0
 ! Compute the heat flux in X direction
  call driftXcrflx(uin,Xflux,dx,dy,dz,dt,ngrid,fdx)
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        do icr=1,ncr
           flux(l,i,j,k,icr,1)=Xflux(l, i, j, k,icr)
        end do
    enddo
  enddo
  enddo
  enddo
#if NDIM>1
  ! Compute the heat flux in Y direction
  call driftYcrflx(uin,Yflux,dx,dy,dz,dt,ngrid)
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        do icr=1,ncr
           flux(l,i,j,k,icr,2)=Yflux(l, i, j, k, icr)
        end do
     enddo
  enddo
  enddo
  enddo
#endif
#if NDIM>2
  ! Compute the heat flux in Z direction
  call driftZcrflx(uin,Zflux,dx,dy,dz,dt,ngrid,fdx)
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        do icr=1,ncr
           flux(l,i,j,k,icr,3)=Zflux(l,i,j,k,icr)
        enddo
     enddo
  enddo
  enddo
  enddo
#endif
end subroutine crstreaming_split
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine driftXcrflx(uin,myflux,dx,dy,dz,dt,ngrid,fdx)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid,icr
  integer,dimension(1:ncr)::idcr
  real(dp)::dx,dy,dz,dt

  
  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr)::myflux
  real(dp),dimension(1:ncr)::fluxh_L
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::fdx

  ! Primitive variables
  real(dp)::Bx,Bx_L,Bxh_L,By
  real(dp)::BxL_L,BxL_R,BxL
  real(dp)::BxR_L,BxR_R,BxR
  real(dp)::Bx_xL_L,Bx_xL_R,Bx_xL
  real(dp)::Bx_xR_L,Bx_xR_R,Bx_xR
  real(dp)::By_xL_L,By_xL_R,By_xL
  real(dp)::By_xR_L,By_xR_R,By_xR
  real(dp)::Bx_xhL,By_xhL
  real(dp)::rho_xL,rho_xhL
  real(dp)::Bx_yL_L,Bx_yL_R,Bx_yL
  real(dp)::By_yL_L,By_yL_R,By_yL
  real(dp)::Bx_yhL
  real(dp)::rho,rho_L,rho_R,rhoh_L,rhoh_R
  real(dp)::rho_yL,rho_yhL,facdx
  real(dp),dimension(1:ncr)::ecr,ecr_L,ecr_R,ecr_RR,sigma,sigma_R
  real(dp),dimension(1:ncr)::vstxh_L,vstxh_R,vstx_xhL,vstx_yhL
  real(dp),dimension(1:ncr)::ecr_xL,ecr_xR,ecr_yL,Bgradecr_xhL
  real(dp),dimension(1:ncr)::ecr_xL_yL,ecr_yR,ecr_xL_yR
  real(dp),dimension(1:ncr)::ecr_xhL,ecr_xhL_yR,ecr_xhL_yL
  real(dp),dimension(1:ncr)::Bgradecrh_L
  real(dp),dimension(1:ncr)::ecr_xR_yhL,ecr_xL_yhL,ecr_yhL
  real(dp),dimension(1:ncr)::Bgradecr_yhL,vsty_yhL,By_yhL


  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

 do k=klo,khi
 do j=jlo,jhi
 do i=if1,if2
    do l=1,ngrid
#if NDIM==1
! ---- On recupere les champs magnetiques --------------------------------------
       ! Bx au centre
       Bx   = 0.5d0*(uin(l, i  , j, k, 6) + uin(l, i  , j, k, nvar+1))
       Bx_L = 0.5d0*(uin(l, i-1, j, k, 6) + uin(l, i-1, j, k, nvar+1))
       ! Bx sur les faces
       Bxh_L = 0.5d0*(Bx_L + Bx)
       
! ---- On calcule les valeurs d energie + gradient aux faces -------------------
       do icr=1,ncr
          ecr(icr)   = uin(l, i  , j, k, 8+icr)
          ecr_L(icr) = uin(l, i-1, j, k, 8+icr)
          ecr_R(icr) = uin(l, i+1, j, k, 8+icr)
          !ecr_RR(icr) = uin(l,i+2, j, k, 8+icr)
       end do

       facdx=dx

!!$       facdx=max(fdx(l,i,j,k),fdx(l,i-1,j,k))

!!$       do icr=1,ncr
!!$          call minmod_cr((ecr(icr)-ecr_L(icr))/facdx,(ecr_R(icr)-ecr(icr))/facdx,sigma(icr))
!!$          call minmod_cr((ecr(icr)-ecr_L(icr))/facdx,(ecr_R(icr)-ecr(icr))/facdx,sigma_R(icr))
!!$       end do

       Bgradecrh_L(1:ncr) = Bxh_L*(ecr(1:ncr) - ecr_L(1:ncr))/(dx*facdx)

! ---- On calcule la vitesse de streaming dot u --------------------------------
       do icr=1,ncr
          if(Bgradecrh_L(icr).gt. 0.0d0)then
             rhoh_L = 0.5d0*(uin(l, i  , j, k, 1) + uin(l, i-1, j, k, 1))
             
!!$             write(*,*)'rho',rhoh_L
!!$             write(*,'(A,4e11.3)')'Bgradecrh_L(1)',Bgradecrh_L(1),Bxh_L,ecr(1),ecr_L(1)
             vstxh_L(1:ncr) = - Bxh_L*Bgradecrh_L(1:ncr)/(sqrt(rhoh_L)*abs(Bgradecrh_L(1:ncr)))
!!$       vstxh_L(1:ncr) = 1.0d0 ! YD: test
             
! ---- First order terms --------------------------------
! Streaming velocity at i-1/2 (left cell interface)
! Energy from left or right cell, i-1 or i+1
             if (vstxh_L(icr) > 0.0d0) fluxh_L(icr) = vstxh_L(icr)*ecr_L(icr)*gamma_rad(icr)
             if (vstxh_L(icr) < 0.0d0) fluxh_L(icr) = vstxh_L(icr)*ecr_R(icr)*gamma_rad(icr)
! ---- Second order terms -------------------------------
! with slope limiter
             if (vstxh_L(icr) > 0.0d0) idcr(icr) = i-1
             if (vstxh_L(icr) < 0.0d0) idcr(icr) = i  
             call minmod_cr( (uin(l,idcr(icr),j,k,8+icr)-uin(l,idcr(icr)-1,j,k,8+icr))/facdx &
                  &, (uin(l,idcr(icr)+1,j,k,8+icr)-uin(l,idcr(icr),j,k,8+icr))/facdx, sigma(icr))
             fluxh_L(icr) = fluxh_L(icr) + 0.5d0*abs(vstxh_L(icr))*(facdx-abs(vstxh_L(icr))*dt)*sigma(icr)
             
             ! fluxh_L(1:ncr) = ecrh_L(1:ncr)*vstxh_L(1:ncr)
             
             myflux(l, i, j, k, 1:ncr) = fluxh_L(1:ncr)*dt/dx
          else
             myflux(l, i, j, k, 1:ncr) = 0.0d0
          endif
       enddo
#endif
#if NDIM==2
! ---- 1. Calcul des Flux_x ----------------------------------------------------
! ------------------------------------------------------------------------------
! ---- On recupere les champs magnetiques --------------------------------------
          Bx = 0.5d0*(uin(l, i, j, k, 6)+uin(l, i, j, k, nvar+1))
          By = 0.5d0*(uin(l, i, j, k, 7)+uin(l, i, j, k, nvar+2))

          Bx_xL = 0.5d0*(uin(l, i-1, j, k, 6)+uin(l, i-1, j, k, nvar+1))
          By_xL = 0.5d0*(uin(l, i-1, j, k, 7)+uin(l, i-1, j, k, nvar+2))

          Bx_xhL = 0.5d0*(Bx + Bx_xL)
          By_xhL = 0.5d0*(By + By_xL)

! ---- On calcule les valeurs d energie + gradient aux faces -------------------
          do icr=1,ncr
             ecr(icr)    = uin(l, i  , j, k, 8+icr)
             ecr_xL(icr) = uin(l, i-1, j, k, 8+icr)

             ecr_yL(icr)    = uin(l, i  , j-1, k, 8+icr)
             ecr_xL_yL(icr) = uin(l, i-1, j-1, k, 8+icr)
             ecr_yR(icr)    = uin(l, i  , j+1, k, 8+icr)
             ecr_xL_yR(icr) = uin(l, i-1, j+1, k, 8+icr)

          end do

          ecr_xhL = 0.5d0*(ecr + ecr_xL)

          ecr_xhL_yR = 0.5d0*(ecr_yR + ecr_xL_yR)
          ecr_xhL_yL = 0.5d0*(ecr_yL + ecr_xL_yL)

          Bgradecr_xhL(1:ncr) = Bx_xhL*(ecr(1:ncr) - ecr_xL(1:ncr))/dx &
               & + By_xhL*(ecr_xhL_yR(1:ncr) - ecr_xhL_yL(1:ncr))/(2.0d0*dy)

! ---- On calcule la vitesse de streaming dot u --------------------------------
          rho_xhL = 0.5d0*(uin(l, i, j, k, 1) + uin(l, i-1, j, k, 1))
          
          vstx_xhL(1:ncr) =-vstreaming ! - Bx_xhL*Bgradecr_xhL(1:ncr)/(sqrt(rho_xhL)*abs(Bgradecr_xhL(1:ncr)))

! ---- On calcule les flux xh inf --------------------------------------------

          fluxh_L(1:ncr) = gamma_rad(1:ncr)*ecr_xhL(1:ncr)*vstx_xhL(1:ncr)
          myflux(l, i, j, k, 1:ncr) = -fluxh_L(1:ncr)*dt/dx
#endif

#if NDIM==3

#endif

       enddo
    enddo
 enddo
enddo

end subroutine driftXcrflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine driftYcrflx(uin,myflux,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid,icr
  real(dp)::dx,dy,dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr)::myflux
  real(dp),dimension(1:ncr)::fluxh_L
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin

  ! Primitive variables
  real(dp)::Bx,Bx_L,By
  real(dp)::BxL_L,BxL_R,BxL
  real(dp)::BxR_L,BxR_R,BxR
  real(dp)::Bxh_L,Bxh_R
  real(dp)::Bx_xL_L,Bx_xL_R,Bx_xL
  real(dp)::Bx_xR_L,Bx_xR_R,Bx_xR
  real(dp)::By_xL_L,By_xL_R,By_xL
  real(dp)::By_xR_L,By_xR_R,By_xR
  real(dp)::By_xhL
  real(dp)::rho_xL,rho_xhL
  real(dp)::Bx_yL_L,Bx_yL_R,Bx_yL
  real(dp)::By_yL_L,By_yL_R,By_yL
  real(dp)::Bx_yhL,By_yhL
  real(dp)::rho,rho_L,rho_R,rhoh_L,rhoh_R
  real(dp)::rho_yL,rho_yhL
  real(dp),dimension(1:ncr)::ecr,ecr_L,ecr_R,ecrh_L,ecrh_R
  real(dp),dimension(1:ncr)::Bgrad_ecrh_R,Bgrad_ecrh_L,ecr_yhL_xL
  real(dp),dimension(1:ncr)::vstxh_L,vstxh_R,vstx_xhL,vstx_yhL
  real(dp),dimension(1:ncr)::ecr_xL,ecr_xR,ecr_yL
  real(dp),dimension(1:ncr)::ecr_xL_yL,ecr_yR,ecr_xL_yR,ecr_yL_xL
  real(dp),dimension(1:ncr)::ecr_xhL,ecr_xhL_yR,ecr_xhL_yL,ecr_yhL_xR
  real(dp),dimension(1:ncr)::Bgradecr_xhL,ecr_yL_xR
  real(dp),dimension(1:ncr)::ecr_xR_yhL,ecr_xL_yhL,ecr_yhL
  real(dp),dimension(1:ncr)::Bgradecr_yhL,vsty_yhL


  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ! Conversion factor from user units to cgs units

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
#if NDIM==2

! ---- 2. Calcul des Flux_y ----------------------------------------------------
! ------------------------------------------------------------------------------
! ---- On recupere les champs magnetiques --------------------------------------
        Bx = 0.5d0*(uin(l, i, j, k, 6)+uin(l, i, j, k, nvar+1))
        By = 0.5d0*(uin(l, i, j, k, 7)+uin(l, i, j, k, nvar+2))

          Bx_yL = 0.5d0*(uin(l, i  , j-1, k, 6)+uin(l, i  , j-1, k, nvar+1))
          By_yL = 0.5d0*(uin(l, i  , j-1, k, 7)+uin(l, i  , j-1, k, nvar+2))

          Bx_yhL = 0.5d0*(Bx + Bx_yL)
          By_yhL = 0.5d0*(By + By_yL)

! ---- On calcule les valeurs d energie + gradient aux faces -------------------
          do icr=1,ncr
             ecr(icr)    = uin(l, i, j  , k, 8+icr)
             ecr_yL(icr) = uin(l, i, j-1, k, 8+icr)

             ecr_xL(icr)    = uin(l, i-1, j  , k, 8+icr)
             ecr_yL_xL(icr) = uin(l, i-1, j-1, k, 8+icr)
             ecr_xR(icr)    = uin(l, i+1, j  , k, 8+icr)
             ecr_yL_xR(icr) = uin(l, i+1, j-1, k, 8+icr)

          end do

          ecr_yhL = 0.5d0*(ecr + ecr_yL)

          ecr_yhL_xR = 0.5d0*(ecr_xR + ecr_yL_xR)
          ecr_yhL_xL = 0.5d0*(ecr_xL + ecr_yL_xL)

          Bgradecr_yhL(1:ncr) = By_yhL*(ecr(1:ncr) - ecr_yL(1:ncr))/(dy) &
               & + Bx_yhL*(ecr_yhL_xR(1:ncr) - ecr_yhL_xL(1:ncr))/(2.0d0*dx)

! ---- On calcule la vitesse de streaming dot u --------------------------------
          rho_yhL = 0.5d0*(uin(l, i, j, k, 1) + uin(l, i, j-1, k, 1)) 
          vsty_yhL(1:ncr) = - By_yhL*Bgradecr_yhL(1:ncr)/(sqrt(rho_yhL)*abs(Bgradecr_yhL(1:ncr)))

! ---- On calcule les flux yh inf --------------------------------------------
          fluxh_L(1:ncr) = gamma_rad(1:ncr)*ecr_yhL(1:ncr)*vsty_yhL(1:ncr)
          myflux(l, i, j, k, 1:ncr) = -fluxh_L(1:ncr)*dt/dy

#endif
#if NDIM==3

#endif

     enddo
  enddo
  enddo
  enddo

end subroutine driftYcrflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine driftZcrflx(Temp,bf,myflux,dx,dy,dz,dt,ngrid,compute,ffdx,Dpara,kperp)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid,compute
  real(dp)::dx,dy,dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp,Dpara,kperp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::Bnorm,fz,oneovertwodx,oneoverfourdx,dx_loc
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fz1,fz2,fz3,fz4
  real(dp)::kpar,oneminuskperp,kparaz1,kparaz2,kparaz3,kparaz4
  real(dp)::kperpz1,kperpz2,kperpz3,kperpz4
  real(dp)::oneminuskperpz1,oneminuskperpz2,oneminuskperpz3,oneminuskperpz4
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_kappa

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t

  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

if(isotrope_cond)then
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        dTdz1  =(Temp(l,i,j,k)-Temp(l,i,j,k-1))/dx
        dx_loc=max(ffdx(l,i,j,k),ffdx(l,i,j,k-1))
        kpar=2d0/(Dpara(l,i,j,k-1)+Dpara(l,i,j,k))
        if(compute.ne.3)then
           fz    =kpar*dTdz1
        else
           fz=kpar/dx
        end if
        fz=fz/dx_loc
        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo

else

  if (.not.slopelim_cond)then
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
#if NDIM==3
!!$          ---------
!!$         / |      /|
!!$        /  |     / |
!!$        --------   |
!!$        |  |   |   |
!!$        | /3   |  /4
!!$        |/     | /
!!$        1------2

        ! Centered symmetric scheme
        if(compute .ne. 3)then
           dTdx1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i  ,j-1,k-1) &
                - Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdx2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i+1,j-1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdx3=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i-1,j+1,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdx4=(Temp(l,i+1,j+1,k  )+Temp(l,i+1,j  ,k  )+Temp(l,i+1,j+1,k-1)+Temp(l,i+1,j  ,k-1) &
                - Temp(l,i  ,j+1,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx

           dTdy1=(Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i-1,j  ,k-1) &
                - Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdy2=(Temp(l,i+1,j  ,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i+1,j-1,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdy3=(Temp(l,i  ,j+1,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i-1,j+1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdy4=(Temp(l,i+1,j+1,k  )+Temp(l,i  ,j+1,k  )+Temp(l,i+1,j+1,k-1)+Temp(l,i  ,j+1,k-1) &
                - Temp(l,i+1,j  ,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i+1,j  ,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx

           dTdz1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i-1,j-1,k  ) &
                - Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdz2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  ) &
                - Temp(l,i+1,j  ,k-1)-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdz3=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i-1,j  ,k  ) &
                - Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdz4=(Temp(l,i+1,j+1,k  )+Temp(l,i+1,j  ,k  )+Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  ) &
                - Temp(l,i+1,j+1,k-1)-Temp(l,i+1,j  ,k-1)-Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx
        end if

        bx1=0.25d0*(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=0.25d0*(bf(l,i+1,j-1,k-1,1)+bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1))
        bx3=0.25d0*(bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j+1,k-1,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1))
        bx4=0.25d0*(bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j+1,k-1,1)+bf(l,i+1,j  ,k  ,1)+bf(l,i+1,j+1,k  ,1))

        by1=0.25d0*(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=0.25d0*(bf(l,i  ,j  ,k-1,2)+bf(l,i+1,j  ,k-1,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2))
        by3=0.25d0*(bf(l,i-1,j+1,k-1,2)+bf(l,i  ,j+1,k-1,2)+bf(l,i-1,j+1,k  ,2)+bf(l,i  ,j+1,k  ,2))
        by4=0.25d0*(bf(l,i  ,j+1,k-1,2)+bf(l,i+1,j+1,k-1,2)+bf(l,i  ,j+1,k  ,2)+bf(l,i+1,j+1,k  ,2))

        bz1=0.25d0*(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=0.25d0*(bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i+1,j-1,k  ,3)+bf(l,i+1,j  ,k  ,3))
        bz3=0.25d0*(bf(l,i-1,j  ,k  ,3)+bf(l,i-1,j+1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3))
        bz4=0.25d0*(bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3)+bf(l,i+1,j  ,k  ,3)+bf(l,i+1,j+1,k  ,3))

        Bnorm=sqrt(bx1*bx1+by1*by1+bz1*bz1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
           bz1=bz1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2+bz2*bz2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
           bz2=bz2/Bnorm
        endif
        Bnorm=sqrt(bx3*bx3+by3*by3+bz3*bz3)
        if(Bnorm.gt.0.0)then
           bx3=bx3/Bnorm
           by3=by3/Bnorm
           bz3=bz3/Bnorm
        endif
        Bnorm=sqrt(bx4*bx4+by4*by4+bz4*bz4)
        if(Bnorm.gt.0.0)then
           bx4=bx4/Bnorm
           by4=by4/Bnorm
           bz4=bz4/Bnorm
        endif

!!$        kpar=8d0/(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
!!$             +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
!!$             +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
!!$             +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
!!$             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
!!$             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
!!$             +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1))
!!$        kpar=8d0/(Dpara(l,i+1,j+1,k  )+Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j+1,k-1) &
!!$             +      Dpara(l,i+1,j  ,k-1)+Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  ) &
!!$             +      Dpara(l,i  ,j+1,k-1)+Dpara(l,i  ,j  ,k-1))
!!$

!!$        if(alfven_diff_coeff)then
           kparaz1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
             +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
             +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
           kparaz2=0.125d0*(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
                +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
           kparaz3=0.125d0*(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
                +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1))
           kparaz4=0.125d0*(Dpara(l,i+1,j+1,k  )+Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j+1,k-1) &
                +      Dpara(l,i+1,j  ,k-1)+Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j+1,k-1)+Dpara(l,i  ,j  ,k-1))

           kperpz1=0.125d0*(kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  )+kperp(l,i  ,j  ,k-1) &
             +      kperp(l,i  ,j-1,k-1)+kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ) &
             +      kperp(l,i-1,j  ,k-1)+kperp(l,i-1,j-1,k-1))
           kperpz2=0.125d0*(kperp(l,i+1,j  ,k  )+kperp(l,i+1,j-1,k  )+kperp(l,i+1,j  ,k-1) &
                +      kperp(l,i+1,j-1,k-1)+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  ) &
             +      kperp(l,i  ,j  ,k-1)+kperp(l,i  ,j-1,k-1))
           kperpz3=0.125d0*(kperp(l,i  ,j+1,k  )+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j+1,k-1) &
                +      kperp(l,i  ,j  ,k-1)+kperp(l,i-1,j+1,k  )+kperp(l,i-1,j  ,k  ) &
                +      kperp(l,i-1,j+1,k-1)+kperp(l,i-1,j  ,k-1))
           kperpz4=0.125d0*(kperp(l,i+1,j+1,k  )+kperp(l,i+1,j  ,k  )+kperp(l,i+1,j+1,k-1) &
                +      kperp(l,i+1,j  ,k-1)+kperp(l,i  ,j+1,k  )+kperp(l,i  ,j  ,k  ) &
                +      kperp(l,i  ,j+1,k-1)+kperp(l,i  ,j  ,k-1))

           oneminuskperpz1 = 1.0d0-kperpz1
           oneminuskperpz2 = 1.0d0-kperpz2
           oneminuskperpz3 = 1.0d0-kperpz3
           oneminuskperpz4 = 1.0d0-kperpz4

!!$        else
!!$           kparaz1 = kpar
!!$           kparaz2 = kpar
!!$           kparaz3 = kpar
!!$           kparaz4 = kpar
!!$           kperpz1 = k_perp
!!$           kperpz2 = k_perp
!!$           kperpz3 = k_perp
!!$           kperpz4 = k_perp
!!$           oneminuskperpz1 = oneminuskperp
!!$           oneminuskperpz2 = oneminuskperp
!!$           oneminuskperpz3 = oneminuskperp
!!$           oneminuskperpz4 = oneminuskperp
!!$        end if

        if(compute .ne. 3)then
           fz1=kparaz1*(bz1*oneminuskperpz1*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)+kperpz1*dTdz1)
           fz2=kparaz2*(bz2*oneminuskperpz2*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)+kperpz2*dTdz2)
           fz3=kparaz3*(bz3*oneminuskperpz3*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)+kperpz3*dTdz3)
           fz4=kparaz4*(bz4*oneminuskperpz4*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)+kperpz4*dTdz4)
        else
           fz1=kparaz1*(bz1*oneminuskperpz1*(bx1+by1+bz1)+kperpz1)
           fz2=kparaz2*(bz2*oneminuskperpz2*(bx2+by2+bz2)+kperpz2)
           fz3=kparaz3*(bz3*oneminuskperpz3*(bx3+by3+bz3)+kperpz3)
           fz4=kparaz4*(bz4*oneminuskperpz4*(bx4+by4+bz4)+kperpz4)
        end if
        fz=0.25d0*(fz1+fz2+fz3+fz4)
#endif

        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo
  endif

  if (slopelim_cond)then
     ! TODO
  endif
endif

end subroutine driftZcrflx



subroutine minmod_cr(a,b,sigma)
  use amr_parameters
  implicit none
  real(dp)::a,b
  real(dp)::sigma
  if (abs(a).gt.abs(b)) sigma=b
  if (abs(b).ge.abs(a)) sigma=a
  if (a*b.le.0.0d0) sigma=0.0d0
end subroutine minmod_cr
!###########################################################
!###########################################################
!###########################################################
!###########################################################
