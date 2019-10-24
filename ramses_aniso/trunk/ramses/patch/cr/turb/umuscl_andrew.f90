! ---------------------------------------------------------------
!  MAG_UNSPLIT Unsplit second order Godunov integrator for
!              polytropic magnetized gas dynamics using 
!              MUSCL-HANCOCK scheme
!              with various Riemann solvers and slope limiters.
!              The sheme follows closely the paper by
!              Londrillo & Del Zanna ApJ 2000, 530, 508, 
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
!  This routine was written by Sebastien Fromang and Patrick Hennebelle
! ----------------------------------------------------------------
! modif nimhd
!subroutine mag_unsplit(uin,gravin,flux,emfx,emfy,emfz,tmp,dx,dy,dz,dt,ngrid)
subroutine mag_unsplit(uin,gravin,flux,emfx,emfy,emfz,tmp,dx,dy,dz,dt,ngrid,ind_grid,jcell)
! fin modif nimhd
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  ! modif nimhd
  integer,dimension(1:nvector) :: ind_grid
  ! fin modif nimhd
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 

  ! Output electromotive force
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfx
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfy
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfz

  ! Output courant vector in the cell
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcell

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin 
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),save::bf  

  ! Cell-centered slopes

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Face-centered slopes
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim),save::dbf

  ! Face-averaged left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp
  
  ! Edge-averaged left-right and top-bottom state arrays
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3),save::qLB

  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   ),save::tx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)       ,save::emf

  ! modif nimhd
  ! WARNING following quantities defined with three components even
  ! if ndim<3 !
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::flxmagx,flxmagy,flxmagz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3),save::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::jcentersquare,jxbsquare
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::jemfx,jemfy,jemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::florentzx,florentzy,florentzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::fluxmd,fluxh,fluxad
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::emfambdiff,fluxambdiff
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::emfohmdiss,fluxohm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::emfhall,fluxhall
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3),save:: fvisco
  integer:: ntest
  !real(dp) :: gammaadbis,densionbis
  !real(dp) :: dttemp, rhotemp,btemp
  !integer :: nbidouille, ncompte
  ! fin modif nimhd
  
  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ! modif nimhd
  bmagij=0.d0
  emfambdiff=0.d0
  fluxambdiff=0.d0
  emfohmdiss=0.d0
  fluxohm=0.d0
  jcentersquare=0.d0
  emfhall=0.d0
  fluxhall=0.d0
  fluxmd=0.d0
  fluxh=0.d0
  fluxad=0.d0
  
  bemfx=0.d0
  bemfy=0.d0
  bemfz=0.d0
  jemfx=0.d0
  jemfy=0.d0
  jemfz=0.d0
  florentzx=0.d0
  florentzy=0.d0
  florentzz=0.d0
  
  jcell=0.0d0
 
  fvisco=0.d0
  ! fin modif nimhd
  
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)


  ! Translate to primative variables, compute sound speeds  
  call ctoprim(uin,qin,bf,gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope(bf,qin,dq,dbf,dx,dt,ngrid)

  ! modif nimhd
  if((nambipolar.eq.1).or.(nmagdiffu.eq.1).or.(nhall.eq.1)) then
  
  ! compute Lorentz Force with magnetic fluxes
  !    call computejb(uin,qin,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad)
  
  !  compute Lorentz Force with current
  call computejb2(uin,qin,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,jcell)
  
  endif
  
  ! AMBIPOLAR DIFFUSION
  
  if(nambipolar.eq.1) then
     call computambip(uin,qin,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare)
  
  endif
  
  ! Hall effect
  
  if(nhall.eq.1)then
     call computehall(uin,qin,ngrid,dx,dy,dz,florentzx,florentzy,florentzz,fluxh,emfhall,fluxhall)
  
  endif
  
  ! OHMIC DISSIPATION
  
  if(nmagdiffu.eq.1) then
     call computdifmag(uin,qin,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare)
  endif
  
  !  END OHMIC DISSIPATION
  ! fin modif nimhd
  
  ! modif cmm
  ! Pseudo viscosity
  
  if(nvisco.eq.1) then
     call computevisco(qin,ngrid,dx,dy,dz,dt,fvisco)
  endif
  ! fin modif cmm

  ! Compute 3D traced-states in all three directions
#if NDIM==1
! modif nimhd
!     call trace1d(qin   ,dq    ,qm,qp                ,dx      ,dt,ngrid)
     call trace1d(qin   ,dq    ,qm,qp                ,dx      ,dt,ngrid,jcentersquare)
! fin modif nimhd
#endif
#if NDIM==2
! modif nimhd
!     call trace2d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy   ,dt,ngrid)
     call trace2d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy   ,dt,ngrid,jcentersquare)
! fin modif nimhd
#endif
#if NDIM==3
! modif nimhd
! call trace3d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid)
     call trace3d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid,jcentersquare)
! fin modif nimhd
#endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,6,7,8,fx,tx,ngrid)
  ! Save flux in output array
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! modif nimhd
  ! Energy flux from ohmic term dB/dt=rot(-eta*J)
  if((nambipolar.eq.1).or.(nmagdiffu.eq.1).or.(nhall.eq.1) .and. (.not.radiative_nimhdheating)) then
  
     ivar=5
     do k=klo,khi
        do j=jlo,jhi
           do i=if1,if2
              do l=1,ngrid
  
                 flux(l,i,j,k,ivar,1)=flux(l,i,j,k,ivar,1)+(fluxambdiff(l,i,j,k,1)+fluxhall(l,i,j,k,1)+fluxohm(l,i,j,k,1) )*dt/dx
  
              end do
           end do
        end do
     end do
     
  endif
  ! fin modif nimhd
  
  ! modif cmm
  ! momentum flux from pseudo-viscosity
  if(nvisco.eq.1) then
  
      do ivar=2,4
         do k=klo,khi
            do j=jlo,jhi
               do i=if1,if2                  
                  do l=1,ngrid
      
                    flux(l,i,j,k,ivar,1)=flux(l,i,j,k,ivar,1)+fvisco(l,i,j,k,ivar-1,1)*dt/dx
      
                  end do
               end do
            end do
         end do
      end do
      
  endif
  ! fin modif cmm

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,7,6,8,fx,tx,ngrid)
  ! Save flux in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do

  ! modif nimhd
  ! Energy flux from ohmic term dB/dt=rot(-eta*J)
  if((nambipolar.eq.1).or.(nmagdiffu.eq.1).or.(nhall.eq.1)) then
  
     ivar=5
     do k=klo,khi
        do j=jf1,jf2
           do i=ilo,ihi
              do l=1,ngrid
  
                 flux(l,i,j,k,ivar,2)=flux(l,i,j,k,ivar,2)+(fluxambdiff(l,i,j,k,2)+fluxhall(l,i,j,k,2)+fluxohm(l,i,j,k,2))*dt/dy
  
              end do
           end do
        end do
     end do
     
  endif
  ! fin modif nimhd
  
  ! modif cmm
  ! momentum flux from pseudo-viscosity
  if(nvisco.eq.1) then
  
      do ivar=2,4
        do k=klo,khi
            do j=jf1,jf2
               do i=ilo,ihi
                  do l=1,ngrid
      
                    flux(l,i,j,k,ivar,2)=flux(l,i,j,k,ivar,2)+ fvisco(l,i,j,k,ivar-1,2)*dt/dy
      
                  end do
               end do
            end do
         end do
      end do
  
  endif
  ! fin modif cmm
#endif

  ! Solve for 1D flux in Z direction
#if NDIM==3
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,8,6,7,fx,tx,ngrid)
  ! Save flux in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
  
  ! modif nimhd
  ! Energy flux from ohmic term dB/dt=rot(-eta*J)
  if((nambipolar.eq.1).or.(nmagdiffu.eq.1).or.(nhall.eq.1)) then
  
     ivar=5
     do k=kf1,kf2
        do j=jlo,jhi
           do i=ilo,ihi
              do l=1,ngrid
  
                 flux(l,i,j,k,ivar,3)=flux(l,i,j,k,ivar,3)+(fluxambdiff(l,i,j,k,3)+fluxhall(l,i,j,k,3)+fluxohm(l,i,j,k,3))*dt/dz
  
              end do
           end do
        end do
     end do
     
  endif
  ! fin modif nimhd
  
  ! modif cmm
  ! momentum flux from pseudo-viscosity
  if(nvisco.eq.1) then
  
      do ivar=2,4
         do k=kf1,kf2
            do j=jlo,jhi
               do i=ilo,ihi
                  do l=1,ngrid
      
                    flux(l,i,j,k,ivar,3)=flux(l,i,j,k,ivar,3)+ fvisco(l,i,j,k,ivar-1,3)*dt/dz
      
                  end do
               end do
            end do
         end do
      end do
      
  endif
  ! fin modif nimhd
#endif

! modif nimhd
! emf rather than fluxes
#if NDIM==1

  do k=kf1,kf2
     do j=jf1,jf2
        do i=ilo,ihi
           do l=1,ngrid
              emfx(l,i,j,k)=( emfambdiff(l,i,j,k,nxx)+emfohmdiss(l,i,j,k,nxx)+emfhall(l,i,j,k,nxx)  )*dt/dx
           end do
        end do
     end do
  end do
  
  do k=kf1,kf2
     do j=jlo,jhi
        do i=if1,if2
           do l=1,ngrid
              emfy(l,i,j,k)=( emfambdiff(l,i,j,k,nyy)+emfohmdiss(l,i,j,k,nyy)+emfhall(l,i,j,k,nyy) )*dt/dx
           end do
        end do
     end do
  end do
  do k=klo,khi
     do j=jf1,jf2
        do i=if1,if2
           do l=1,ngrid
              emfz(l,i,j,k)=( emfambdiff(l,i,j,k,nzz)+emfohmdiss(l,i,j,k,nzz)+emfhall(l,i,j,k,nzz) )*dt/dx
           end do
        end do
     end do
  end do

#endif
! fin modif nimhd

#if NDIM>1
  ! Solve for EMF in Z direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &           qLT,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               if1  ,if2  ,jf1  ,jf2  ,klo  ,khi  , 2,3,4,6,7,8,emf,ngrid)
 
 ! Save vector in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k)*dt/dx
        ! modif nimhd
        emfz(l,i,j,k)=emfz(l,i,j,k)+(emfambdiff(l,i,j,k,nzz)+emfohmdiss(l,i,j,k,nzz)+emfhall(l,i,j,k,nzz) )*dt/dx
        ! fin modif nimhd
     end do
  end do
  end do
  end do
#if NDIM==2
  do k=2,2
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k-1)*dt/dx
     end do
  end do
  end do
  end do
#endif
#endif

#if NDIM>2
  ! Solve for EMF in Y direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1,ju2,ku1+1,ku2+1, &
       &           qLT,iu1  ,iu2  ,ju1,ju2,ku1+1,ku2+1, &
       &           qRB,iu1+1,iu2+1,ju1,ju2,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1,ju2,ku1  ,ku2  , &
       &               if1  ,if2  ,jlo,jhi,kf1  ,kf2  , 4,2,3,8,6,7,emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=if1,if2
     do l=1,ngrid
        emfy(l,i,j,k)=emf(l,i,j,k)*dt/dx
        ! modif nimhd
        emfy(l,i,j,k)=emfy(l,i,j,k)+ ( emfambdiff(l,i,j,k,nyy)+emfohmdiss(l,i,j,k,nyy)+emfhall(l,i,j,k,nyy) )*dt/dx
        ! fin modif nimhd
     end do
  end do
  end do
  end do
  ! Solve for EMF in X direction
  CALL cmp_mag_flx(qRT,iu1,iu2,ju1+1,ju2+1,ku1+1,ku2+1, &
       &           qRB,iu1,iu2,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLT,iu1,iu2,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qLB,iu1,iu2,ju1  ,ju2  ,ku1  ,ku2  , &
       &               ilo,ihi,jf1  ,jf2  ,kf1  ,kf2  , 3,4,2,7,8,6,emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jf1,jf2
  do i=ilo,ihi
     do l=1,ngrid
        emfx(l,i,j,k)=emf(l,i,j,k)*dt/dx
        ! modif nimhd
        emfx(l,i,j,k)=1.d0*emfx(l,i,j,k)+( emfambdiff(l,i,j,k,nxx)+emfohmdiss(l,i,j,k,nxx)+emfhall(l,i,j,k,nxx) )*dt/dx
        ! fin modif nimhd
     end do
  end do
  end do
  end do
#endif

end subroutine mag_unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! modif cmm
subroutine computevisco(q,ngrid,dx,dy,dz,dt,fvisco)

  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3) :: fvisco

  real(dp) :: muvisco

! declare local variables
  INTEGER :: i, j, k, l
  real(dp) :: rhox,rhoy,rhoz


 do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              rhox=0.5d0*(q(l,i,j,k,1)+q(l,i-1,j,k,1))
              rhoy=0.5d0*(q(l,i,j,k,1)+q(l,i,j-1,k,1))
              rhoz=0.5d0*(q(l,i,j,k,1)+q(l,i,j,k-1,1))

! WARNING Flux F defined as dU/dt+dF/dx=0 
              fvisco(l,i,j,k,1,1)=-muvisco(rhox)*(q(l,i,j,k,2)-q(l,i-1,j,k,2))/dx
              fvisco(l,i,j,k,1,2)=-muvisco(rhoy)*(q(l,i,j,k,2)-q(l,i,j-1,k,2))/dy
              fvisco(l,i,j,k,1,3)=-muvisco(rhoz)*(q(l,i,j,k,2)-q(l,i,j,k-1,2))/dz
              fvisco(l,i,j,k,2,1)=-muvisco(rhox)*(q(l,i,j,k,3)-q(l,i-1,j,k,3))/dx
              fvisco(l,i,j,k,2,2)=-muvisco(rhoy)*(q(l,i,j,k,3)-q(l,i,j-1,k,3))/dy
              fvisco(l,i,j,k,2,3)=-muvisco(rhoz)*(q(l,i,j,k,3)-q(l,i,j,k-1,3))/dz
              fvisco(l,i,j,k,3,1)=-muvisco(rhox)*(q(l,i,j,k,4)-q(l,i-1,j,k,4))/dx
              fvisco(l,i,j,k,3,2)=-muvisco(rhoy)*(q(l,i,j,k,4)-q(l,i,j-1,k,4))/dy
              fvisco(l,i,j,k,3,3)=-muvisco(rhoz)*(q(l,i,j,k,4)-q(l,i,j,k-1,4))/dz

           end do
        end do
     end do
  end do

  
end subroutine computevisco
! fin modif cmm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! modif nimhd
subroutine computejb(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad)

  USE amr_parameters
  use hydro_commons
   USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 

  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt


real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxh,fluxad
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij


! declare local variables
  INTEGER ::i, j, k, l, m, n 

real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
 real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::flxmagxx,flxmagxy,flxmagxz,flxmagyx,flxmagyy,flxmagyz,flxmagzx,flxmagzy,flxmagzz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface,fluxbis,fluxter,fluxquat
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bcenter
real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
real(dp)::b12x,b12y,b12z,emag,bsquare
real(dp)::computdivbisx,computdivbisy,computdivbisz
real(dp)::computdxbis,computdybis,computdzbis
real(dp)::crossprodx,crossprody,crossprodz

! magnetic field at center of cells

do k=ku1,ku2
     do j=ju1,ju2
        do i=iu1,iu2
           
           do l=1,ngrid
              bcenter(l,i,j,k,nxx)=q(l,i,j,k,6)
              bcenter(l,i,j,k,nyy)=q(l,i,j,k,7)
              bcenter(l,i,j,k,nzz)=q(l,i,j,k,8)

           end do
        end do
     end do
  end do


!!!!!!!!!!!!!!!!!!
! EMF x
!!!!!!!!!!!!!!!!!!

! magnetic field at location of EMF

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfx(l,i,j,k,1)=0.25d0*( q(l,i,j,k,6)+q(l,i,j-1,k,6)+q(l,i,j,k-1,6)+q(l,i,j-1,k-1,6) )

           end do
        end do
     end do
  end do

 do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfx(l,i,j,k,2)=0.5d0*( u(l,i,j,k,7)+u(l,i,j,k-1,7) )

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid
                            
              bemfx(l,i,j,k,3)=0.5d0*(u(l,i,j,k,8)+u(l,i,j-1,k,8))
              
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
! EMF y
!!!!!!!!!!!!!!!!!!


! magnetic field at location of EMF

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfy(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j,k-1,6))

           end do
        end do
     end do
  end do

 do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bemfy(l,i,j,k,2)=0.25d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7)+q(l,i,j,k-1,7)+q(l,i-1,j,k-1,7))

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
                            
              bemfy(l,i,j,k,3)=0.5d0*(u(l,i-1,j,k,8)+u(l,i,j,k,8))
        
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
! EMF z
!!!!!!!!!!!!!!!!!!

! magnetic field at location of EMF

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfz(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6))

           end do
        end do
     end do
  end do

 do k=ku1,ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bemfz(l,i,j,k,2)=0.5d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7))

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
                            
              bemfz(l,i,j,k,3)=0.25d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8)+q(l,i,j-1,k,8)+q(l,i-1,j-1,k,8))
 
           end do
        end do
     end do
  end do

! bmagij is the value of the magnetic field Bi where Bj 
! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
! and we can write it Bx,y

  do k=ku1,ku2
     do j=ju1,ju2
        do i=iu1,iu2
           do l=1,ngrid
              
              do m=1,3
                
!! m+5 mandatory cf Bx=uin(l,i,j,k,6)
                 bmagij(l,i,j,k,m,m)=u(l,i,j,k,m+5)


              end do
           end do
        end do
     end do
  end do


! case Bx,y

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,1,2)=0.5d0*(q(l,i,j,k,6)+q(l,i,j-1,k,6))

           end do
        end do
     end do
  end do


! case Bx,z

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,1,3)=0.5d0*(q(l,i,j,k,6)+q(l,i,j,k-1,6))

           end do
        end do
     end do
  end do

! case By,x

  do k=ku1,ku2
     do j=ju1,max(1,ju2-1)
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,2,1)=0.5d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7))

           end do
        end do
     end do
  end do

! case By,z

  do k=min(1,ku1+1),ku2
     do j=ju1,max(1,ju2-1)
        do i=iu1,iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,2,3)=0.5d0*(q(l,i,j,k,7)+q(l,i,j,k-1,7))

           end do
        end do
     end do
  end do

! case Bz,x

  do k=ku1,max(1,ku2-1)
     do j=ju1,ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,3,1)=0.5d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8))

           end do
        end do
     end do
  end do

! case Bz,y

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2
        do i=iu1,iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,3,2)=0.5d0*(q(l,i,j,k,8)+q(l,i,j-1,k,8))

           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
!
! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
! Bn at i-1/2,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!

do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2
        do i=iu1,iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,1)=0.25d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6)+u(l,i,j,k-1,6)+u(l,i,j-1,k-1,6))

           end do
        end do
     end do
  end do

! case By for Lorentz force EMF 

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,2)=0.25d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7)+u(l,i,j,k-1,7)+u(l,i-1,j,k-1,7)) 
  
           end do
        end do
     end do
  end do
 
! case Bz for Lorentz force EMF 

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,3)=0.25d0*(u(l,i,j,k,8)+u(l,i-1,j,k,8)+u(l,i,j-1,k,8)+u(l,i-1,j-1,k,8)) 
 
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! computation of the component of j where EMFs are located
! jemfx(l,i,j,k,n) is the component Jn at i,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid

              jemfx(l,i,j,k,1)=(u(l,i,j,k,8)-u(l,i,j-1,k,8))/dy-(u(l,i,j,k,7)-u(l,i,j,k-1,7))/dz 
              jemfx(l,i,j,k,2)=(bmagij(l,i,j,k,1,2)-bmagij(l,i,j,k-1,1,2))/dz- (bmagijbis(l,i+1,j,k,3)-bmagijbis(l,i,j,k,3))/dx
              jemfx(l,i,j,k,3)=(bmagijbis(l,i+1,j,k,2) -bmagijbis(l,i,j,k,2))/dx- (bmagij(l,i,j,k,1,3)-bmagij(l,i,j-1,k,1,3))/dy
              

              jemfy(l,i,j,k,1)=(bmagijbis(l,i,j+1,k,3)-bmagijbis(l,i,j,k,3))/dy-(bmagij(l,i,j,k,2,1) - bmagij(l,i,j,k-1,2,1) )/dz
              jemfy(l,i,j,k,2)=(u(l,i,j,k,6)-u(l,i,j,k-1,6))/dz-(u(l,i,j,k,8)-u(l,i-1,j,k,8))/dx
              jemfy(l,i,j,k,3)=(bmagij(l,i,j,k,2,3)-bmagij(l,i-1,j,k,2,3))/dx-(bmagijbis(l,i,j+1,k,1)-bmagijbis(l,i,j,k,1))/dy


              jemfz(l,i,j,k,1)=(bmagij(l,i,j,k,3,1) -bmagij(l,i,j-1,k,3,1))/dy-(bmagijbis(l,i,j,k+1,2)-bmagijbis(l,i,j,k,2))/dz
              jemfz(l,i,j,k,2)=( bmagijbis(l,i,j,k+1,1)-bmagijbis(l,i,j,k,1))/dz-(bmagij(l,i,j,k,3,2)-bmagij(l,i-1,j,k,3,2))/dx
              jemfz(l,i,j,k,3)=(u(l,i,j,k,7)-u(l,i-1,j,k,7))/dx-(u(l,i,j,k,6)-u(l,i,j-1,k,6))/dy

           end do
        end do
     end do
  end do


if((nambipolar.eq.1).or.(nhall.eq.1)) then

! Fx,x

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              b12x=bmagijbis(l,i,j,k,1)
              b12y=bmagijbis(l,i,j,k,2)
              b12z=bmagijbis(l,i,j,k,3)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)

              flxmagxx(l,i,j,k,1)=b12x*b12x-emag
              
           end do
        end do
     end do
  end do


! Fx,y

  do k=min(1,ku1+1),ku2
     do j=ju1,max(1,ju2-1)
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid
              
              b12x=bmagij(l,i,j,k,1,3)
              b12y=bmagij(l,i,j,k,2,3)

              flxmagxx(l,i,j,k,2)=b12x*b12y

           end do
        end do
     end do
  end do

! Fx,z

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2
        do i=iu1,max(1,iu2-1)

           do l=1,ngrid
              
              b12x=bmagij(l,i,j,k,1,2)
              b12z=bmagij(l,i,j,k,3,2)

              flxmagxx(l,i,j,k,3)=b12x*b12z

           end do
        end do
     end do
  end do

! Fy,x
           
  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              b12x=bmagijbis(l,i,j,k,1)
              b12y=bmagijbis(l,i,j,k,2)
              
              flxmagyx(l,i,j,k,1)=b12y*b12x

           end do
        end do
     end do
  end do

! Fy,y

  do k=min(1,ku1+1),ku2
     do j=ju1,max(1,ju2-1)
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,3)
              b12y=bmagij(l,i,j,k,2,3)
              b12z=bmagij(l,i,j,k,3,3)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                           
              flxmagyx(l,i,j,k,2)=b12y*b12y-emag
              
           end do
        end do
     end do
  end do

! Fy,z

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2
        do i=iu1,max(1,iu2-1)

           do l=1,ngrid
              
              b12y=bmagij(l,i,j,k,2,2)
              b12z=bmagij(l,i,j,k,3,2)
              flxmagyx(l,i,j,k,3)=b12y*b12z
              
           end do
        end do
     end do
  end do

! Fz,x

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              b12x=bmagijbis(l,i,j,k,1)
              b12z=bmagijbis(l,i,j,k,3)
              
              flxmagzx(l,i,j,k,1)=b12z*b12x

           end do
        end do
     end do
  end do

! Fz,y

  do k=min(1,ku1+1),ku2
     do j=ju1,max(1,ju2-1)
        do i=iu1,max(1,iu2-1)

           do l=1,ngrid

              b12y=bmagij(l,i,j,k,2,3)
              b12z=bmagij(l,i,j,k,3,3)
                                         
              flxmagzx(l,i,j,k,2)=b12z*b12y
              
           end do
        end do
     end do
  end do

! Fz,z

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,2) 
              b12y=bmagij(l,i,j,k,2,2)
              b12z=bmagij(l,i,j,k,3,2)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
              flxmagzx(l,i,j,k,3)=b12z*b12z-emag
              
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Fx,x

  do k=min(1,ku1+1),ku2
     do j= ju1,max(1,ju2-1)     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,3)
              b12y=bmagij(l,i,j,k,2,3)
              b12z=bmagij(l,i,j,k,3,3)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)

              flxmagxy(l,i,j,k,1)=b12x*b12x-emag

           end do
        end do
     end do
  end do

! Fx,y

 do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2

           do l=1,ngrid
              
              b12x=bmagijbis(l,i,j,k,1)
              b12y=bmagijbis(l,i,j,k,2)
              flxmagxy(l,i,j,k,2)=b12x*b12y
              
           end do
        end do
     end do
  end do

! Fx,z
  
  do k=ku1,max(1,ku2-1)
    do j=ju1,max(1,ju2-1)
       do i=min(1,iu1+1),iu2

           do l=1,ngrid
              
              b12x=bmagij(l,i,j,k,1,1)
              b12z=bmagij(l,i,j,k,3,1)
              flxmagxy(l,i,j,k,3)=b12x*b12z
              
           end do
        end do
     end do
  end do

! Fy,x


  do k=min(1,ku1+1),ku2
     do j= ju1,max(1,ju2-1)     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,3)
              b12y=bmagij(l,i,j,k,2,3)

              flxmagyy(l,i,j,k,1)=b12y*b12x

           end do
        end do
     end do
  end do

! Fy,y

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
                      
           do l=1,ngrid
              
              b12x=bmagijbis(l,i,j,k,1)
              b12y=bmagijbis(l,i,j,k,2)
              b12z=bmagijbis(l,i,j,k,3)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)

              flxmagyy(l,i,j,k,2)=b12y*b12y-emag
              
           end do
        end do
     end do
  end do

! Fy,z
  
  do k=ku1,max(1,ku2-1)
    do j=ju1,max(1,ju2-1)
       do i=min(1,iu1+1),iu2

           do l=1,ngrid
              
              b12y=bmagij(l,i,j,k,2,1)
              b12z=bmagij(l,i,j,k,3,1)
              flxmagyy(l,i,j,k,3)=b12y*b12z
              
           end do
        end do
     end do
  end do


! Fz,x

  do k=min(1,ku1+1),ku2
     do j= ju1,max(1,ju2-1)     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,3)
              b12z=bmagij(l,i,j,k,3,3)

              flxmagzy(l,i,j,k,1)=b12z*b12x

           end do
        end do
     end do
  end do

! Fz,y

 do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
              
              b12y=bmagijbis(l,i,j,k,2)
              b12z=bmagijbis(l,i,j,k,3)

              flxmagzy(l,i,j,k,2)=b12z*b12y
              
           end do
        end do
     end do
  end do

! Fz,z
  
  do k=ku1,max(1,ku2-1)
    do j=ju1,max(1,ju2-1)
       do i=min(1,iu1+1),iu2

           do l=1,ngrid
              
              b12x=bmagij(l,i,j,k,1,1)
              b12y=bmagij(l,i,j,k,2,1)
              b12z=bmagij(l,i,j,k,3,1)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
              
              flxmagzy(l,i,j,k,3)=b12z*b12z-emag
              
           end do
        end do
     end do
  end do


! Fx,x

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,2)
              b12y=bmagij(l,i,j,k,2,2)
              b12z=bmagij(l,i,j,k,3,2)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)

              flxmagxz(l,i,j,k,1)=b12x*b12x-emag

           end do
        end do
     end do
  end do

! Fx,y

 do k=ku1,max(1,ku2-1)
     do j=ju1,max(1,ju2-1)      
        do i=min(1,iu1+1),iu2 
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,1)
              b12y=bmagij(l,i,j,k,2,1)
             
              flxmagxz(l,i,j,k,2)=b12x*b12y
              
           end do
        end do
     end do
  end do

! Fx,z

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
              
              b12x=bmagijbis(l,i,j,k,1)
              b12z=bmagijbis(l,i,j,k,3)
              flxmagxz(l,i,j,k,3)=b12x*b12z
              
           end do
        end do
     end do
  end do


! Fy,x


  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,2)
              b12y=bmagij(l,i,j,k,2,2)
             
              flxmagyz(l,i,j,k,1)=b12y*b12x

           end do
        end do
     end do
  end do

! Fy,y

 do k=ku1,max(1,ku2-1)
     do j=ju1,max(1,ju2-1)      
        do i=min(1,iu1+1),iu2 
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,1)
              b12y=bmagij(l,i,j,k,2,1)
              b12z=bmagij(l,i,j,k,3,1)

              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
             
              flxmagyz(l,i,j,k,2)=b12y*b12y-emag
              
           end do
        end do
     end do
  end do

! Fy,z

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
              
              b12y=bmagijbis(l,i,j,k,2)
              b12z=bmagijbis(l,i,j,k,3)
              flxmagyz(l,i,j,k,3)=b12y*b12z
              
           end do
        end do
     end do
  end do

! Fz,x


  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2     
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid

              b12x=bmagij(l,i,j,k,1,2)
              b12z=bmagij(l,i,j,k,3,2)
             
              flxmagzz(l,i,j,k,1)=b12z*b12x

           end do
        end do
     end do
  end do

! Fz,y

 do k=ku1,max(1,ku2-1)
     do j=ju1,max(1,ju2-1)      
        do i=min(1,iu1+1),iu2 
           
           do l=1,ngrid

              b12y=bmagij(l,i,j,k,2,1)
              b12z=bmagij(l,i,j,k,3,1)
             
              flxmagzz(l,i,j,k,2)=b12z*b12y
              
           end do
        end do
     end do
  end do

! Fz,z

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
              
              b12x=bmagijbis(l,i,j,k,1)
              b12y=bmagijbis(l,i,j,k,2)
              b12z=bmagijbis(l,i,j,k,3)
              emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)

              flxmagzz(l,i,j,k,3)=b12z*b12z-emag
              
           end do
        end do
     end do
  end do



!!!!!!!!!!!!!!!!!!!!!!!
do k=min(1,ku1+1),max(1,ku2-1)
    do j=min(1,ju1+1),max(1,ju2-1)
       do i=min(1,iu1+1),max(1,iu2-1)

          do l = 1, ngrid

             florentzx(l,i,j,k,1)=computdivbisx(flxmagxx,l,i,j,k,dx,dy,dz)
             florentzx(l,i,j,k,2)=computdivbisx(flxmagyx,l,i,j,k,dx,dy,dz)
             florentzx(l,i,j,k,3)=computdivbisx(flxmagzx,l,i,j,k,dx,dy,dz)

             florentzy(l,i,j,k,1)=computdivbisy(flxmagxy,l,i,j,k,dx,dy,dz)
             florentzy(l,i,j,k,2)=computdivbisy(flxmagyy,l,i,j,k,dx,dy,dz)
             florentzy(l,i,j,k,3)=computdivbisy(flxmagzy,l,i,j,k,dx,dy,dz)

             florentzz(l,i,j,k,1)=computdivbisz(flxmagxz,l,i,j,k,dx,dy,dz)
             florentzz(l,i,j,k,2)=computdivbisz(flxmagyz,l,i,j,k,dx,dy,dz)
             florentzz(l,i,j,k,3)=computdivbisz(flxmagzz,l,i,j,k,dx,dy,dz)

          end do
       end do
    end do
 end do

! end if((nambipolar.eq.1).or.(nhall.eq.1)) then
endif

! computation of current on faces

if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1).or.(nmagdiffu2.eq.1)) then

! face at i-1/2,j,k


do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)           
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,1,1)=computdybis(bemfz,3,l,i,j,k,dy)-computdzbis(bemfy,2,l,i,j,k,dz)

           end do
        end do
     end do
  end do



 do k=min(1,ku1+1),max(1,ku2-1)
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,2,1)=computdzbis(bemfy,1,l,i,j,k,dz)-computdxbis(bcenter,3,l,i-1,j,k,dx)

           end do
        end do
     end do
  end do


  do k=ku1,ku2
     do j=min(1,ju1+1),max(1,ju2-1)       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,3,1)=computdxbis(bcenter,2,l,i-1,j,k,dx)-computdybis(bemfz,1,l,i,j,k,dy)

           end do
        end do
     end do
  end do

! face at i,j-1/2,k

do k=min(1,ku1+1),max(1,ku2-1) 
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid
              
              jface(l,i,j,k,1,2)=computdybis(bcenter,3,l,i,j-1,k,dy)-computdzbis(bemfx,2,l,i,j,k,dz)


           end do
        end do
     end do
  end do

do k=min(1,ku1+1),max(1,ku2-1) 
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),max(1,iu2-1) 
           
           do l=1,ngrid
              
              jface(l,i,j,k,2,2)=computdzbis(bemfx,1,l,i,j,k,dz)-computdxbis(bemfz,3,l,i,j,k,dx)

           end do
        end do
     end do
  end do


do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),max(1,iu2-1) 
           
           do l=1,ngrid

              jface(l,i,j,k,3,2)=computdxbis(bemfz,2,l,i,j,k,dx)-computdybis(bcenter,1,l,i,j-1,k,dy)

           end do
        end do
     end do
  end do

! face at i,j,k-1/2



  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),max(1,ju2-1)        
        do i=iu1,iu2
           
           do l=1,ngrid
                            
              jface(l,i,j,k,1,3)=computdybis(bemfx,3,l,i,j,k,dy)-computdzbis(bcenter,2,l,i,j,k-1,dz)             
              
           end do
        end do
     end do
  end do



do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid
                            
              jface(l,i,j,k,2,3)=computdzbis(bcenter,1,l,i,j,k-1,dz)-computdxbis(bemfy,3,l,i,j,k,dx)             

           end do
        end do
     end do
  end do

do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),max(1,ju2-1)      
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid
                            
              jface(l,i,j,k,3,3)=computdxbis(bemfy,2,l,i,j,k,dx)-computdybis(bemfx,1,l,i,j,k,dx)            

           end do
        end do
     end do
  end do

 do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(jface,bmagij,fluxbis,l,i,j,k)

              fluxmd(l,i,j,k,1)=fluxbis(l,i,j,k,1,1)
              fluxmd(l,i,j,k,2)=fluxbis(l,i,j,k,2,2)
              fluxmd(l,i,j,k,3)=fluxbis(l,i,j,k,3,3)
  
           end do
        end do
     end do
  end do


! end if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1)) then
endif

if((nambipolar.eq.1).or.(nhall.eq.1)) then

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(fluxbis,bmagij,fluxter,l,i,j,k)

              fluxh(l,i,j,k,1)=fluxter(l,i,j,k,1,1)
              fluxh(l,i,j,k,2)=fluxter(l,i,j,k,2,2)
              fluxh(l,i,j,k,3)=fluxter(l,i,j,k,3,3)

   
           end do
        end do
     end do
  end do

! end if((nambipolar.eq.1).or.(nhall.eq.1)) then
endif


if(nambipolar.eq.1)then

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)

              fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
              fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
              fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)

           end do
        end do
     end do
  end do

! end if(nambipolar.eq.1) then
endif


end subroutine computejb
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computejb2(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,jcell)

  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 


  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt


real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxh,fluxad
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcell


! declare local variables
  INTEGER ::i, j, k, l, m, n 

real(dp)::computdx,computdy,computdz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::flxmagxx,flxmagxy,flxmagxz,flxmagyx,flxmagyy,flxmagyz,flxmagzx,flxmagzy,flxmagzz
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bcenter
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::fluxbis,fluxter,fluxquat
real(dp)::b12x,b12y,b12z,emag,bsquare
real(dp)::computdivbisx,computdivbisy,computdivbisz
real(dp)::computdxbis,computdybis,computdzbis

! magnetic field at center of cells
do k=ku1,ku2
     do j=ju1,ju2
        do i=iu1,iu2
           do l=1,ngrid
              bcenter(l,i,j,k,nxx)=q(l,i,j,k,6)
              bcenter(l,i,j,k,nyy)=q(l,i,j,k,7)
              bcenter(l,i,j,k,nzz)=q(l,i,j,k,8)
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
! EMF x
!!!!!!!!!!!!!!!!!!

! magnetic field at location of EMF

  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfx(l,i,j,k,1)=0.25d0*( q(l,i,j,k,6)+q(l,i,j-1,k,6)+q(l,i,j,k-1,6)+q(l,i,j-1,k-1,6) )

           end do
        end do
     end do
  end do

 do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfx(l,i,j,k,2)=0.5d0*( u(l,i,j,k,7)+u(l,i,j,k-1,7) )

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid
                            
              bemfx(l,i,j,k,3)=0.5d0*(u(l,i,j,k,8)+u(l,i,j-1,k,8))
              
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
! EMF y
!!!!!!!!!!!!!!!!!!


! magnetic field at location of EMF

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfy(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j,k-1,6))

           end do
        end do
     end do
  end do

 do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bemfy(l,i,j,k,2)=0.25d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7)+q(l,i,j,k-1,7)+q(l,i-1,j,k-1,7))

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
                            
              bemfy(l,i,j,k,3)=0.5d0*(u(l,i-1,j,k,8)+u(l,i,j,k,8))
        
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
! EMF z
!!!!!!!!!!!!!!!!!!

! magnetic field at location of EMF

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid

              bemfz(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6))

           end do
        end do
     end do
  end do

 do k=ku1,ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bemfz(l,i,j,k,2)=0.5d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7))

           end do
        end do
     end do
  end do

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
                            
              bemfz(l,i,j,k,3)=0.25d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8)+q(l,i,j-1,k,8)+q(l,i-1,j-1,k,8))
 
           end do
        end do
     end do
  end do

! bmagij is the value of the magnetic field Bi where Bj 
! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
! and we can write it Bx,y

  do k=ku1,ku2
     do j=ju1,ju2
        do i=iu1,iu2
           do l=1,ngrid
              
              do m=1,3
                
!! m+5 mandatory cf Bx=uin(l,i,j,k,6)
                 bmagij(l,i,j,k,m,m)=u(l,i,j,k,m+5)

 
              end do
           end do
        end do
     end do
  end do


! case Bx,y

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,1,2)=0.5d0*(q(l,i,j,k,6)+q(l,i,j-1,k,6))

           end do
        end do
     end do
  end do


! case Bx,z

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2
        do i=iu1,max(1,iu2-1)
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,1,3)=0.5d0*(q(l,i,j,k,6)+q(l,i,j,k-1,6))

           end do
        end do
     end do
  end do

! case By,x

  do k=ku1,ku2
     do j=ju1,max(1,ju2-1)
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,2,1)=0.5d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7))

           end do
        end do
     end do
  end do

! case By,z

  do k=min(1,ku1+1),ku2
     do j=ju1,max(1,ju2-1)
        do i=iu1,iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,2,3)=0.5d0*(q(l,i,j,k,7)+q(l,i,j,k-1,7))

           end do
        end do
     end do
  end do

! case Bz,x

  do k=ku1,max(1,ku2-1)
     do j=ju1,ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,3,1)=0.5d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8))

           end do
        end do
     end do
  end do

! case Bz,y

  do k=ku1,max(1,ku2-1)
     do j=min(1,ju1+1),ju2
        do i=iu1,iu2
           
           do l=1,ngrid
               
              bmagij(l,i,j,k,3,2)=0.5d0*(q(l,i,j,k,8)+q(l,i,j-1,k,8))

           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!
!
! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
! Bn at i-1/2,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!

do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2
        do i=iu1,iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,1)=0.25d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6)+u(l,i,j,k-1,6)+u(l,i,j-1,k-1,6))

           end do
        end do
     end do
  end do

! case By for Lorentz force EMF 

  do k=min(1,ku1+1),ku2
     do j=ju1,ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,2)=0.25d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7)+u(l,i,j,k-1,7)+u(l,i-1,j,k-1,7)) 
  
           end do
        end do
     end do
  end do
 
! case Bz for Lorentz force EMF 

  do k=ku1,ku2
     do j=min(1,ju1+1),ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bmagijbis(l,i,j,k,3)=0.25d0*(u(l,i,j,k,8)+u(l,i-1,j,k,8)+u(l,i,j-1,k,8)+u(l,i-1,j-1,k,8)) 
 
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! computation of the component of j where EMFs are located
! jemfx(l,i,j,k,n) is the component Jn at i,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid

              jemfx(l,i,j,k,1)=(u(l,i,j,k,8)-u(l,i,j-1,k,8))/dy-(u(l,i,j,k,7)-u(l,i,j,k-1,7))/dz 
              jemfx(l,i,j,k,2)=(bmagij(l,i,j,k,1,2)-bmagij(l,i,j,k-1,1,2))/dz- (bmagijbis(l,i+1,j,k,3)-bmagijbis(l,i,j,k,3))/dx
              jemfx(l,i,j,k,3)=(bmagijbis(l,i+1,j,k,2) -bmagijbis(l,i,j,k,2))/dx- (bmagij(l,i,j,k,1,3)-bmagij(l,i,j-1,k,1,3))/dy
              


              jemfy(l,i,j,k,1)=(bmagijbis(l,i,j+1,k,3)-bmagijbis(l,i,j,k,3))/dy-(bmagij(l,i,j,k,2,1) - bmagij(l,i,j,k-1,2,1) )/dz
              jemfy(l,i,j,k,2)=(u(l,i,j,k,6)-u(l,i,j,k-1,6))/dz-(u(l,i,j,k,8)-u(l,i-1,j,k,8))/dx
              jemfy(l,i,j,k,3)=(bmagij(l,i,j,k,2,3)-bmagij(l,i-1,j,k,2,3))/dx-(bmagijbis(l,i,j+1,k,1)-bmagijbis(l,i,j,k,1))/dy


              jemfz(l,i,j,k,1)=(bmagij(l,i,j,k,3,1) -bmagij(l,i,j-1,k,3,1))/dy-(bmagijbis(l,i,j,k+1,2)-bmagijbis(l,i,j,k,2))/dz
              jemfz(l,i,j,k,2)=( bmagijbis(l,i,j,k+1,1)-bmagijbis(l,i,j,k,1))/dz-(bmagij(l,i,j,k,3,2)-bmagij(l,i-1,j,k,3,2))/dx
              jemfz(l,i,j,k,3)=(u(l,i,j,k,7)-u(l,i-1,j,k,7))/dx-(u(l,i,j,k,6)-u(l,i,j-1,k,6))/dy

           end do
        end do
     end do
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! computation of the component of j at center of cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    do k=1,1!min(1,ku1+1),max(1,ku2-1)
!       do j=1,1!min(1,ju1+1),max(1,ju2-1)
!          do i=1,1!min(1,iu1+1),max(1,iu2-1)
   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
              jcell(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcell(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcell(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)

            end do
         end do
      end do
   end do


if((nambipolar.eq.1).or.(nhall.eq.1).or.(nambipolar2.eq.1)) then

! EMF x
  
  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid
              
              call crossprod(jemfx,bemfx,florentzx,l,i,j,k)
              call crossprod(jemfy,bemfy,florentzy,l,i,j,k)
              call crossprod(jemfz,bemfz,florentzz,l,i,j,k)
                    
           end do
        end do
     end do
  end do
  

! end if((nambipolar.eq.1).or.(nhall.eq.1)) then
endif


! computation of current on faces

if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1).or.(nambipolar2.eq.1).or.(nmagdiffu2.eq.1)) then

! face at i-1/2,j,k


do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)           
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,1,1)=computdybis(bemfz,3,l,i,j,k,dy)-computdzbis(bemfy,2,l,i,j,k,dz)

           end do
        end do
     end do
  end do



 do k=min(1,ku1+1),max(1,ku2-1)
     do j=ju1,ju2       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,2,1)=computdzbis(bemfy,1,l,i,j,k,dz)-computdxbis(bcenter,3,l,i-1,j,k,dx)

           end do
        end do
     end do
  end do


  do k=ku1,ku2
     do j=min(1,ju1+1),max(1,ju2-1)       
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              jface(l,i,j,k,3,1)=computdxbis(bcenter,2,l,i-1,j,k,dx)-computdybis(bemfz,1,l,i,j,k,dy)

           end do
        end do
     end do
  end do

! face at i,j-1/2,k

do k=min(1,ku1+1),max(1,ku2-1) 
     do j=min(1,ju1+1),ju2       
        do i=iu1,iu2
           
           do l=1,ngrid
              
              jface(l,i,j,k,1,2)=computdybis(bcenter,3,l,i,j-1,k,dy)-computdzbis(bemfx,2,l,i,j,k,dz)


           end do
        end do
     end do
  end do

do k=min(1,ku1+1),max(1,ku2-1) 
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),max(1,iu2-1) 
           
           do l=1,ngrid
              
              jface(l,i,j,k,2,2)=computdzbis(bemfx,1,l,i,j,k,dz)-computdxbis(bemfz,3,l,i,j,k,dx)

           end do
        end do
     end do
  end do


do k=ku1,ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),max(1,iu2-1) 
           
           do l=1,ngrid

              jface(l,i,j,k,3,2)=computdxbis(bemfz,2,l,i,j,k,dx)-computdybis(bcenter,1,l,i,j-1,k,dy)

           end do
        end do
     end do
  end do

! face at i,j,k-1/2



  do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),max(1,ju2-1)        
        do i=iu1,iu2
           
           do l=1,ngrid
                            
              jface(l,i,j,k,1,3)=computdybis(bemfx,3,l,i,j,k,dy)-computdzbis(bcenter,2,l,i,j,k-1,dz)             
              
           end do
        end do
     end do
  end do



do k=min(1,ku1+1),ku2
     do j=ju1,ju2       
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid
                            
              jface(l,i,j,k,2,3)=computdzbis(bcenter,1,l,i,j,k-1,dz)-computdxbis(bemfy,3,l,i,j,k,dx)             

           end do
        end do
     end do
  end do

do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),max(1,ju2-1)      
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid
                            
              jface(l,i,j,k,3,3)=computdxbis(bemfy,2,l,i,j,k,dx)-computdybis(bemfx,1,l,i,j,k,dx)            

           end do
        end do
     end do
  end do


do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(jface,bmagij,fluxbis,l,i,j,k)

              fluxmd(l,i,j,k,1)=fluxbis(l,i,j,k,1,1)
              fluxmd(l,i,j,k,2)=fluxbis(l,i,j,k,2,2)
              fluxmd(l,i,j,k,3)=fluxbis(l,i,j,k,3,3)
  
           end do
        end do
     end do
  end do


! end if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1)) then
endif

if((nambipolar.eq.1).or.(nhall.eq.1).or.(nambipolar2==1)) then

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(fluxbis,bmagij,fluxter,l,i,j,k)

              fluxh(l,i,j,k,1)=fluxter(l,i,j,k,1,1)
              fluxh(l,i,j,k,2)=fluxter(l,i,j,k,2,2)
              fluxh(l,i,j,k,3)=fluxter(l,i,j,k,3,3)

   
           end do
        end do
     end do
  end do

! end if((nambipolar.eq.1).or.(nhall.eq.1)) then
endif


if((nambipolar.eq.1).or.(nambipolar2==1))then

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)

              fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
              fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
              fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)

           end do
        end do
     end do
  end do


! end if(nambipolar.eq.1) then
endif




do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)

              fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
              fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
              fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)

           end do
        end do
     end do
  end do


end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare)

  USE amr_parameters
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  USE const
  use units_commons
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 

  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd

real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3):: emfohmdiss,fluxohm 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare

! declare local variables
  INTEGER ::i, j, k, l, m, n, ht,h


! WARNING following quantities defined with three components even
! if ndim<3 !
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcenter
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemf


real(dp)::computdx,computdy,computdz,tcell,rhocell,bcell
real(dp)::crossprodx,crossprody,crossprodz, Cv

real(dp)::etaohmdiss,etaod2,tcellx,tcelly,tcellz,bsquarex,bsquarey,bsquarez,etaohmdissbricolo,dtlim
real(dp)::pressurex,pressurey,pressurez,rhox,rhoy,rhoz,epsx,epsy,epsz
real(dp)::etaod2x,etaod2y,etaod2z,rhof,pf,bsqf,epsf,tcellf,barotrop1D

integer , dimension(1:3) :: index_i,index_j,index_k


index_i = (/1,0,0/)
index_j = (/0,1,0/)
index_k = (/0,0,1/)

dtlim = dt !neil

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid

              jemf(l,i,j,k,1)=jemfx(l,i,j,k,1)
              jemf(l,i,j,k,2)=jemfy(l,i,j,k,2)
              jemf(l,i,j,k,3)=jemfz(l,i,j,k,3)


              rhox=0.25d0*(u(l,i,j,k,   1)+u(l,i  ,j-1,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i  ,j-1,k-1,   1))
              rhoy=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i-1,j  ,k-1,   1))
              rhoz=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j-1,k  ,   1)+u(l,i-1,j-1,k  ,   1))
              epsx=0.25d0*(u(l,i,j,k,nvar)+u(l,i  ,j-1,k,nvar)+u(l,i,j  ,k-1,nvar)+u(l,i  ,j-1,k-1,nvar))
              epsy=0.25d0*(u(l,i,j,k,nvar)+u(l,i-1,j  ,k,nvar)+u(l,i,j  ,k-1,nvar)+u(l,i-1,j  ,k-1,nvar))
              epsz=0.25d0*(u(l,i,j,k,nvar)+u(l,i-1,j  ,k,nvar)+u(l,i,j-1,k  ,nvar)+u(l,i-1,j-1,k  ,nvar))
               if(nmagdiffu .eq.1)then
                 bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
                 bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
                 bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2
               else if(nmagdiffu2 .eq.1)then
                  bsquarex=u(l,i,j,k,   2)
                  bsquarey=u(l,i,j,k,   2)
                  bsquarez=u(l,i,j,k,   2)
                  epsx=u(l,i,j,k,nvar)
                  epsy=u(l,i,j,k,nvar)
                  epsz=u(l,i,j,k,nvar)
                  rhox=u(l,i,j,k,1)
                  rhoy=u(l,i,j,k,1)
                  rhoz=u(l,i,j,k,1)
                  if(epsx .ne.u(l,i,j,k,3))then
                     ! Attention, on est sur les boundary du domaine, divu et enew ne sont pas connus....
                     bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
                     bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
                     bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2
                  end if
               end if

               if(ntestDADM.eq.1)then
                  tcellx=1.0d0
                  tcelly=1.0d0
                  tcellz=1.0d0
               else
!                  print*,'x',rhox,epsx,u(l,i,j,k,2),bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
!                  if(epsy* scale_d*scale_v**2  .lt. 1.d-16) print*,'y',rhoy,epsy
!                  if(epsz* scale_d*scale_v**2  .lt. 1.d-16) print*,'z',rhoz,epsz
!                  print*,rhox,epsx,rhoy,epsy,rhoz,epsz
                  call temperature_eos(rhox,epsx,tcellx,ht)
                  call temperature_eos(rhoy,epsy,tcelly,ht)
                  call temperature_eos(rhoz,epsz,tcellz,ht)
!!$                  tcelly=10.
!!$                  tcellx=10.
!!$                  tcellz=10.
               endif
               
!               etaod2x=etaohmdiss(rhox,bsquarex,tcellx)
!               etaod2y=etaohmdiss(rhoy,bsquarey,tcelly)
!               etaod2z=etaohmdiss(rhoz,bsquarez,tcellz)
              etaod2x=etaohmdissbricolo(rhox,bsquarex,tcellx,dtlim,dx)
              etaod2y=etaohmdissbricolo(rhoy,bsquarey,tcelly,dtlim,dx)
              etaod2z=etaohmdissbricolo(rhoz,bsquarez,tcellz,dtlim,dx)
              
! WARNING dB/dt=-curl(eta*J)
              emfohmdiss(l,i,j,k,nxx)=-etaod2x*jemf(l,i,j,k,1)
              emfohmdiss(l,i,j,k,nyy)=-etaod2y*jemf(l,i,j,k,2)
              emfohmdiss(l,i,j,k,nzz)=-etaod2z*jemf(l,i,j,k,3)

! !!!!!!!!!!!!!!!!!!!!!!!
! !
! ! compute j at center of cells
! !
! ! mandatory for non isotherm case
! 
! ! bmagij is the value of the magnetic field Bi where Bj 
! ! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
! ! and we can write it Bx,y

              jcenter(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcenter(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcenter(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! ! computation of j at faces
! ! jface is the value of the current where Bj 
! ! is naturally defined; Ex jface(l,i,j,k,1,2) is Jx at i,j-1/2,k
! ! and we can write it Jx,y
              if(nmagdiffu2 .eq.0)then
                 do h = 1,3
                    
                    rhof=0.5d0*(u(l,i,j,k,   1)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),   1))
!                 epsf=u(l,i,j,k,3)
                    epsf=0.5d0*(u(l,i,j,k,nvar)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),nvar))
                    bsqf=bmagij(l,i,j,k,1,h)**2+bmagij(l,i,j,k,2,h)**2+bmagij(l,i,j,k,3,h)**2
                    
                 ! Compute gas temperature in cgs
                    if(eos ) then 
                       call temperature_eos(rhof,epsf,tcellf,ht)
                    elseif(barotrop)then
                       tcellf=barotrop1D(rhof*scale_d)
                    elseif(ntestDADM.eq.1)then
                       tcellf=1.0d0
                    else
                       write(*,*) 'Temperature EOS needs updating!'
                    endif
                    
                    etaod2=etaohmdiss(rhof,bsqf,tcellf)
                    fluxohm(l,i,j,k,h)=etaod2*fluxmd(l,i,j,k,h)
                    
                    !               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
                    !               pf=0.5d0*(q(l,i,j,k,5)+q(l,i-1,j,k,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,1)=etaod2*fluxmd(l,i,j,k,1)
                    
                 enddo
              end if

!            end do
!         end do
!      end do
!   end do
! 
! 
!   do k=min(1,ku1+1),max(1,ku2-1) 
! !     do j=min(1,ju1+1),ju2 
!      do j=min(1,ju1+1),max(1,ju2-1)
!         do i=min(1,iu1+1),max(1,iu2-1) 
!            
!            do l=1,ngrid
! 
!               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
!               pf=0.5d0*(q(l,i,j,k,5)+q(l,i,j-1,k,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,2)=etaod2*fluxmd(l,i,j,k,2)
! 
!            end do
!         end do
!      end do
!   end do
! 
! 
!  do k=min(1,ku1+1),max(1,ku2-1) 
! !     do j=min(1,ju1+1),ju2  
!     do j=min(1,ju1+1),max(1,ju2-1)
!         do i=min(1,iu1+1),max(1,iu2-1) 
!            
!            do l=1,ngrid
! 
!               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
!               pf=0.5d0*(q(l,i,j,k,5)+q(l,i,j,k-1,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,3)=etaod2*fluxmd(l,i,j,k,3)
! 
!            end do
!         end do
!      end do
!   end do


! compute contribution to energy flux +eta*I*B


!            do l=1,ngrid
              if(nmagdiffu2 .eq. 1)then
                 jcentersquare(l,i,j,k)=jcenter(l,i,j,k,1)*jcenter(l,i,j,k,1)+jcenter(l,i,j,k,2)*jcenter(l,i,j,k,2)+jcenter(l,i,j,k,3)*jcenter(l,i,j,k,3)
                 
                 rhocell = u(l,i,j,k,1)
                 bcell   = u(l,i,j,k,2)
                 if(u(l,i,j,k,nvar) .ne.u(l,i,j,k,3))then
                    ! Attention, on est sur les boundary du domaine, divu et enew ne sont pas connus....
                    bcell=(0.5*(u(l,i,j,k,6)+u(l,i,j,k,nvar+1)))**2 + (0.5*(u(l,i,j,k,7)+u(l,i,j,k,nvar+2)))**2 +(0.5*(u(l,i,j,k,8)+u(l,i,j,k,nvar+3)))**2
                 end if

                 if(ntestDADM.eq.1)then
                    tcell=1.0d0
                 else 
                    call temperature_eos(rhocell,u(l,i,j,k,nvar),tcell,ht)
!                    if(nmagdiffu2.eq.1)call temperature_eos(rhocell,u(l,i,j,k,3),tcell,ht)
                    end if
                    
                    !neil jcentersquare(l,i,j,k) = jcentersquare(l,i,j,k)*etaohmdiss(rhocell,bcell,tcell)*dt
                    jcentersquare(l,i,j,k) = jcentersquare(l,i,j,k)*etaohmdissbricolo(rhocell,bcell,tcell,dtlim,dx)*dt
                    
                 end if
              end do
        end do
     end do
  end do
  

end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare)

  use amr_commons
  USE amr_parameters
  use hydro_commons
  USE const
  use units_commons
  IMPLICIT NONE
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  
  INTEGER ::ngrid,ht
  REAL(dp)::dx,dy,dz,dt,dtambdiff2,barotrop1D
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::rhocellmin,bsquaremax
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxad
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::emfambdiff
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxambdiff

! declare local variables
  INTEGER ::i, j, k, l, m, n, ntest,ic,ivar
  real(dp)::computdx,computdy,computdz

  real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
  real(dp)::rhofx,rhofy,rhofz
  real(dp)::bsquarex,bsquarey,bsquarez,bsquare
  real(dp)::bsquarexx,bsquareyy,bsquarezz
  real(dp)::betaad2,betaadbricolo
  real(dp)::rhox,rhoy,rhoz,rhocell,bcell,bcellold,tcell
  real(dp)::dtlim,Cv,eps
  real(dp)::crossprodx,crossprody,crossprodz

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentz

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jxbsquare
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcenter
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jxb


!modif pour voir les lieux du seuil
!  real(dp),dimension(1:3)             :: skip_loc
!  real(dp),dimension(1:twotondim,1:3) :: xc
!  integer                             :: nx_loc
!  real(dp)                            :: scale


!nx_loc = (icoarse_max -icoarse_min+1)
!scale = dble(nx_loc)/boxlen
!print*, dx, 0.5d0**11, 0.5d0**7/scale
!
!do  ind=1,twotondim
!    iz=(ind-1)/4
!    iy=



! do NOT change value below Variation of betaad
! to avoid too small time step allowed
  ntest=0

  dtlim=dt!*coefalfven
!dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

jxb=0.0d0

jxbsquare=0.0d0
jcenter=0.0d0

  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
              rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
              rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))

              rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
              rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
              rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
              
              rhocellmin(l,i,j,k)=min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)

              rhocell = u(l,i,j,k,1)

              ! Compute gas temperature in cgs
              if(ntestDADM.eq.1) then
                 tcell=1.0d0
              else
                 call temperature_eos(u(l,i,j,k,1),u(l,i,j,k,nvar),tcell,ht)
                 if(nambipolar2.eq.1)call temperature_eos(u(l,i,j,k,1),u(l,i,j,k,3),tcell,ht)
              end if
             

              bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
              bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
              bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2


              bsquarexx=bmagij(l,i,j,k,1,1)**2+bmagij(l,i,j,k,2,1)**2+bmagij(l,i,j,k,3,1)**2
              bsquareyy=bmagij(l,i,j,k,1,2)**2+bmagij(l,i,j,k,2,2)**2+bmagij(l,i,j,k,3,2)**2
              bsquarezz=bmagij(l,i,j,k,1,3)**2+bmagij(l,i,j,k,2,3)**2+bmagij(l,i,j,k,3,3)**2

              bsquaremax(l,i,j,k)=max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)
                 
! EMF x
  
              v1x=florentzx(l,i,j,k,1)
              v1y=florentzx(l,i,j,k,2)
              v1z=florentzx(l,i,j,k,3)
              v2x=bemfx(l,i,j,k,1)
              v2y=bemfx(l,i,j,k,2)
              v2z=bemfx(l,i,j,k,3)
             
              emfambdiff(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
              rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
!                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
!!$                 rhox=rhocell
!!$                 rhoy=rhocell
!!$                 rhoz=rhocell
              end if

              if(nambipolar2 .eq. 0)rhocell = rhocellmin(l,i,j,k)
! alfven time alone maybe not correct
!             betaad2=betaadbricolo(rhox,dtlim,bsquare,dx,ntest)
! comparison with hydro+idealMHD
!!$              rhox= u(l,i,j,k,3)
!!$              rhoy= u(l,i,j,k,3)
!!$              rhoz= u(l,i,j,k,3)

              betaad2=betaadbricolo(rhocell,rhox,dtlim,bcell,bcellold,dx,ntest,tcell)
!              betaad2=betaadbricolo(rhocell,rhox,dtlim,bcellold,bcellold,dx,ntest,tcell)

              emfambdiff(l,i,j,k,1)=emfambdiff(l,i,j,k,1)*betaad2 

! EMF y
              v1x=florentzy(l,i,j,k,1)
              v1y=florentzy(l,i,j,k,2)
              v1z=florentzy(l,i,j,k,3)
              v2x=bemfy(l,i,j,k,1)
              v2y=bemfy(l,i,j,k,2)
              v2z=bemfy(l,i,j,k,3)
             
              emfambdiff(l,i,j,k,2)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)

              rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
!                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
              end if

              if(nambipolar2 .eq. 0)rhocell = rhocellmin(l,i,j,k)
! alfven time alone maybe not correct
!             betaad2=betaadbricolo(rhoy,dtlim,bsquare,dx,ntest)
! comparison with hydro+idealMHD 

             betaad2=betaadbricolo(rhocell,rhoy,dtlim,bcell,bcellold,dx,ntest,tcell)
!             betaad2=betaadbricolo(rhocell,rhoy,dtlim,bcellold,bcellold,dx,ntest,tcell)

             emfambdiff(l,i,j,k,2)=emfambdiff(l,i,j,k,2)*betaad2            
                    
! EMF z

             v1x=florentzz(l,i,j,k,1)
             v1y=florentzz(l,i,j,k,2)
             v1z=florentzz(l,i,j,k,3)
             v2x=bemfz(l,i,j,k,1)
             v2y=bemfz(l,i,j,k,2)
             v2z=bemfz(l,i,j,k,3)
             
             emfambdiff(l,i,j,k,3)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)
              rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
 !                bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
              end if
              if(nambipolar2 .eq. 0) rhocell = rhocellmin(l,i,j,k)
             
             betaad2=betaadbricolo(rhocell,rhoz,dtlim,bcell,bcellold,dx,ntest,tcell)
!             betaad2=betaadbricolo(rhocell,rhoz,dtlim,bcellold,bcellold,dx,ntest,tcell)

             emfambdiff(l,i,j,k,3)=emfambdiff(l,i,j,k,3)*betaad2

! energy flux on faces

              v2x=bmagij(l,i,j,k,1,1)
              v2y=bmagij(l,i,j,k,2,1)
              v2z=bmagij(l,i,j,k,3,1)

             bcell = bsquaremax(l,i,j,k)
             if(nambipolar2.eq.1)then
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if


              rhocell = rhocellmin(l,i,j,k)
              rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
              betaad2=betaadbricolo(rhocell,rhofx,dtlim,bcell,bcell,dx,ntest,tcell)
              fluxambdiff(l,i,j,k,1)=-betaad2*fluxad(l,i,j,k,1)

              v2x=bmagij(l,i,j,k,1,2)
              v2y=bmagij(l,i,j,k,2,2)
              v2z=bmagij(l,i,j,k,3,2)

              rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
              bcell = bsquaremax(l,i,j,k)
              if(nambipolar2.eq.1)then
!                 bcell=0.5d0*(u(l,i,j,k,2)+u(l,i,j-1,k,2))
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              betaad2=betaadbricolo(rhocell,rhofy,dtlim,bcell,bcell,dx,ntest,tcell)
              fluxambdiff(l,i,j,k,2)=-betaad2*fluxad(l,i,j,k,2)

              v2x=bmagij(l,i,j,k,1,3)
              v2y=bmagij(l,i,j,k,2,3)
              v2z=bmagij(l,i,j,k,3,3)
!              bsquare=v2x*v2x+v2y*v2y+v2z*v2z
              rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
              bcell = bsquaremax(l,i,j,k)
              if(nambipolar2.eq.1)then
!                 bcell=0.5d0*(u(l,i,j,k,2)+u(l,i,j,k-1,2))
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              betaad2=betaadbricolo(rhocell,rhofz,dtlim,bcell,bcell,dx,ntest,tcell)
              fluxambdiff(l,i,j,k,3)=-betaad2*fluxad(l,i,j,k,3)

              v2x=u(l,i,j,k,6)
              v2y=u(l,i,j,k,7)
              v2z=u(l,i,j,k,8)
             !              bsquare=v2x*v2x+v2y*v2y+v2z*v2z
              bcellold=bcell
             if(nambipolar2.eq.1)then
                 bcellold=u(l,i,j,k,2)
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              jcenter(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcenter(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcenter(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)

              call crossprod(jcenter,u(:,:,:,:,6:8),jxb,l,i,j,k)

              jxbsquare(l,i,j,k)=(jxb(l,i,j,k,1)*jxb(l,i,j,k,1)+jxb(l,i,j,k,2)*jxb(l,i,j,k,2)+jxb(l,i,j,k,3)*jxb(l,i,j,k,3))*&
              & betaadbricolo(u(l,i,j,k,1),u(l,i,j,k,1),dtlim,bcell,bcellold,dx,ntest,tcell)*dtlim
              



           end do
        end do
     end do
  end do

end SUBROUTINE computambip
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computehall(u,q,ngrid,dx,dy,dz,florentzx,florentzy,florentzz,fluxh,emfhall,fluxhall)


  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  
  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz,fluxh
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::emfhall,fluxhall

! declare local variables
  INTEGER ::i, j, k, l
  real(dp)::rhox,rhoy,rhoz,rhofx,rhofy,rhofz,pfx,pfy,pfz
  real(dp)::crossprodx,crossprody,crossprodz,reshall
  real(dp)::computdivbisx,computdivbisy,computdivbisz
  real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

! EMF x
  
  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid            

              emfhall(l,i,j,k,1)=florentzx(l,i,j,k,1)

              rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))

              emfhall(l,i,j,k,1)=-reshall(rhox)*emfhall(l,i,j,k,1)
                    
           end do
        end do
     end do
  end do

! EMF y
  
  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid
              
              emfhall(l,i,j,k,2)=florentzy(l,i,j,k,2)

              rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))

              emfhall(l,i,j,k,2)=-reshall(rhoy)*emfhall(l,i,j,k,2)
                    
           end do
        end do
     end do
  end do

! EMF z
  
  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l = 1, ngrid

              emfhall(l,i,j,k,3)=florentzz(l,i,j,k,3)

              rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))

              emfhall(l,i,j,k,3)=-reshall(rhoz)*emfhall(l,i,j,k,3)
                    
           end do
        end do
     end do
  end do

! energy flux on faces

  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1) 
        do i=min(1,iu1+1),max(1,iu2-1) 
           do l=1,ngrid
              
              rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
!              pf=0.5d0*(q(l,i,j,k,5)+q(l,i-1,j,k,5))
              fluxhall(l,i,j,k,1)=reshall(rhofx)*fluxh(l,i,j,k,1)

              rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
!              pfy=0.5d0*(q(l,i,j,k,5)+q(l,i,j-1,k,5))
              fluxhall(l,i,j,k,2)=reshall(rhofy)*fluxh(l,i,j,k,2)

               rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
!              pfz=0.5d0*(q(l,i,j,k,5)+q(l,i,j,k-1,5))
               fluxhall(l,i,j,k,3)=reshall(rhofz)*fluxh(l,i,j,k,3)

           end do
        end do
     end do
  end do

end SUBROUTINE computehall
! fin modif nimhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==1
! modif nimhd
!SUBROUTINE  trace1d(q,dq,qm,qp,dx,dt,ngrid)
SUBROUTINE  trace1d(q,dq,qm,qp,dx,dt,ngrid,jcentersquare)
! fin modif nimhd
  USE amr_parameters
  USE hydro_parameters
  USE const
  use radiation_parameters,only:small_er,is_radiative_energy
  use units_commons
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx,dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 
  ! modif nimhd
  REAL(dp):: etaohmdiss,bcell,tcell,barotrop1D
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare
! fin modif nimhd

  ! declare local variables
  INTEGER ::i, j, k, l, n, ht
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC
  REAL(dp)::dtdx
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0

#if USE_FLD==1 || USE_M_1==1
  REAL(dp),ALLOCATABLE,DIMENSION(:)::E,dEx,sE0  
  if(ngrp>0)ALLOCATE(E(1:nrad),dEx(1:nrad),sE0(1:nrad))
#endif

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir = 1; iu = 2; iv = 3 ; iw = 4 ; ip = 5; iA = 6; iB = 7; iC = 8

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r = q(l,i,j,k,ir)
              u = q(l,i,j,k,iu)
              v = q(l,i,j,k,iv)
              w = q(l,i,j,k,iw)
              p = q(l,i,j,k,ip)
              A = q(l,i,j,k,iA)
              B = q(l,i,j,k,iB)
              C = q(l,i,j,k,iC)
#if USE_FLD==1 || USE_M_1==1
              do n = 1,nrad
                 E(n) = q(l,i,j,k,iC+n)
              end do
#endif
              ! TVD slopes in X direction
              drx = half*dq(l,i,j,k,ir,1)
              dux = half*dq(l,i,j,k,iu,1)
              dvx = half*dq(l,i,j,k,iv,1)
              dwx = half*dq(l,i,j,k,iw,1)
              dpx = half*dq(l,i,j,k,ip,1)
              dBx = half*dq(l,i,j,k,iB,1)
              dCx = half*dq(l,i,j,k,iC,1)
#if USE_FLD==1 || USE_M_1==1
              do n = 1,nrad
                 dEx(n) = half*dq(l,i,j,k,iC+n,1)
              end do
#endif

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-r*dux
              if(i_scheme .ne. i_induction)then
              su0 = -u*dux-(dpx+B*dBx+C*dCx)/r  
#if USE_FLD==1
              do n = 1,ngrp
                 su0 = su0 - dEx(n)/r/3.0D0
              end do
#endif
              sv0 = -u*dvx+(A*dBx)/r
              sw0 = -u*dwx+(A*dCx)/r
              endif
              sp0 = -u*dpx-gamma*p*dux
              sB0 = -u*dBx+A*dvx-B*dux
              sC0 = -u*dCx+A*dwx-C*dux
#if USE_FLD==1
              do n = 1,ngrp
                 sE0(n) = -u*dEx(n)-4.0d0/3.0d0*E(n)*dux
              end do
#endif
#if USE_M_1==1
              do n = 1,ngrp
                 sE0(n) = -u*dEx(n)
              end do
#endif
              
              ! modif nimhd
              if(magohm.eq.0 )then
                 if(ntestDADM.eq.1) then
                    tcell=1.0d0
                 else
                    call temperature_eos(q(l,i,j,k,1),q(l,i,j,k,nvar),tcell,ht)
                 endif
                 bcell = A*A + B*B + C*C
                 sp0 = sp0 +(gamma-1.d0)*etaohmdiss(q(l,i,j,k,1),bcell,tcell)*jcentersquare(l,i,j,k)
              endif     
              ! fin modif nimhd
              
              ! Cell-centered predicted states
              r = r + sr0*dtdx
              u = u + su0*dtdx
              v = v + sv0*dtdx
              w = w + sw0*dtdx
              p = p + sp0*dtdx
              B = B + sB0*dtdx
              C = C + sC0*dtdx
#if USE_FLD==1 || USE_M_1==1
              do n = 1,ngrp
                 E(n) = E(n) + sE0(n)*dtdx
              end do
#endif

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = A
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
#if USE_FLD==1 || USE_M_1==1
              do n = 1,nrad
                 qp(l,i,j,k,iC+n,1) = E(n) - dEx(n)
                 if(is_radiative_energy(n+1)) qp(l,i,j,k,iC+n,1) = max(small_er, qp(l,i,j,k,iC+n,1))
              end do
#endif

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = A
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
#if USE_FLD==1 || USE_M_1==1
              do n = 1,nrad
                 qm(l,i,j,k,iC+n,1) = E(n) + dEx(n)
                 if(is_radiative_energy(n+1)) qm(l,i,j,k,iC+n,1) = max(small_er, qm(l,i,j,k,iC+n,1))
              end do
#endif
           END DO
        END DO
     END DO
  END DO
  
  ! passive scalars
  DO n = 9+nrad, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 a   = q(l,i,j,k,n )           ! Cell centered values
                 u   = q(l,i,j,k,iu)  
                 dax = half * dq(l,i,j,k,n,1)  ! TVD slope
                 sa0 = -u*dax                  ! Source terms
                 a   = a + sa0*dtdx            ! Predicted state
                 qp(l,i,j,k,n,1) = a - dax     ! Right state
                 qm(l,i,j,k,n,1) = a + dax     ! Left state
              END DO
           END DO
        END DO
     END DO
  END DO


END SUBROUTINE trace1d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==2
! modif nimhd
!SUBROUTINE trace2d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dt,ngrid)
SUBROUTINE trace2d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dt,ngrid,jcentersquare)
! fin modif nimhd
  USE amr_parameters
  USE hydro_parameters
  USE const
  use radiation_parameters,only:small_er
  use units_commons
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx, dy, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB
  ! modif nimhd
  REAL(dp)::etaohmdiss,bcell,tcell,barotrop1D
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare
  ! fin modif nimhd

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ez
  INTEGER ::i, j, k, l, n, ht
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC 
  REAL(dp)::dtdx, dtdy, smallp
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::dry, duy, dvy, dwy, dpy, dAy, dBy, dCy
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0
  REAL(dp)::AL, AR, BL, BR
  REAL(dp)::dALy, dARy, dBLx, dBRx
  REAL(DP)::sAL0, sAR0, sBL0, sBR0

#if USE_FLD==1
  REAL(dp),ALLOCATABLE,DIMENSION(:)::E,dEx,dEy,sE0  
  if(ngrp>0)ALLOCATE(E(1:ngrp),dEx(1:ngrp),dEy(1:ngrp),sE0(1:ngrp))
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  smallp = smallr*smallc**2/gamma

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5; ia=6; ib=7; ic=8 

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           DO l = 1, ngrid
              u = 0.25*(q(l,i-1,j-1,k,iu)+q(l,i-1,j,k,iu)+q(l,i,j-1,k,iu)+q(l,i,j,k,iu))
              v = 0.25*(q(l,i-1,j-1,k,iv)+q(l,i-1,j,k,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k,iv))
              A = 0.5*(bf(l,i,j-1,k,1)+bf(l,i,j,k,1))
              B = 0.5*(bf(l,i-1,j,k,2)+bf(l,i,j,k,2))
              Ez(l,i,j,k)=u*B-v*A
           END DO
        END DO
     END DO
  END DO

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)
              p =    q(l,i,j,k,ip)
              A =    q(l,i,j,k,iA)
              B =    q(l,i,j,k,iB)
              C =    q(l,i,j,k,iC)
#if USE_FLD==1
              do n = 1, ngrp
                 E(n) = q(l,i,j,k,iC+n)
              end do
#endif

              ! Face centered variables
              AL =  bf(l,i  ,j  ,k,1)
              AR =  bf(l,i+1,j  ,k,1)
              BL =  bf(l,i  ,j  ,k,2)
              BR =  bf(l,i  ,j+1,k,2)

              ! Cell centered TVD slopes in X direction
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)
              dBx = half * dq(l,i,j,k,iB,1)
              dCx = half * dq(l,i,j,k,iC,1)
#if USE_FLD==1
              do n = 1, ngrp
                 dEx(n) = half*dq(l,i,j,k,iC+n,1)
              end do
#endif

              ! Cell centered TVD slopes in Y direction
              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)
              dAy = half * dq(l,i,j,k,iA,2)
              dCy = half * dq(l,i,j,k,iC,2)
#if USE_FLD==1
              do n = 1, ngrp
                 dEy(n) = half*dq(l,i,j,k,iC+n,2)
              end do
#endif

              ! Face centered TVD slopes in transverse direction
              dALy = half * dbf(l,i  ,j  ,k,1,1)
              dARy = half * dbf(l,i+1,j  ,k,1,1)
              dBLx = half * dbf(l,i  ,j  ,k,2,1)
              dBRx = half * dbf(l,i  ,j+1,k,2,1)

              ! Edge centered electric field Ez = uB-vA
              ELL = Ez(l,i  ,j  ,k)
              ELR = Ez(l,i  ,j+1,k)
              ERL = Ez(l,i+1,j  ,k)
              ERR = Ez(l,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(ELR-ELL)*dtdy*half
              sAR0 = +(ERR-ERL)*dtdy*half
              sBL0 = -(ERL-ELL)*dtdx*half
              sBR0 = -(ERR-ELR)*dtdx*half

              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              
              ! Source terms (including transverse derivatives)
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy
              if(i_scheme.ne.i_induction)then
              su0 = (-u*dux-(dpx+B*dBx+C*dCx)/r)*dtdx + (-v*duy+B*dAy/r)*dtdy 
#if USE_FLD==1
              do n = 1, ngrp
                 su0 = su0 - (third*dEx(n)/r)*dtdx
              end do
#endif
              sv0 = (-u*dvx+A*dBx/r)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)/r)*dtdy
#if USE_FLD==1
              do n = 1, ngrp
                 sv0 = sv0 - (third*dEy(n)/r)*dtdy
              end do
#endif
              sw0 = (-u*dwx+A*dCx/r)*dtdx + (-v*dwy+B*dCy/r)*dtdy
              endif
              sp0 = (-u*dpx-dux*gamma*p)*dtdx + (-v*dpy-dvy*gamma*p)*dtdy
              sC0 = (-u*dCx-C*dux+A*dwx)*dtdx + (-v*dCy-C*dvy+B*dwy)*dtdy
#if USE_FLD==1
              do n = 1, ngrp
                 sE0(n) = (-u*dEx(n)-4.0d0*third*E(n)*dux)*dtdx &
                      & + (-v*dEy(n)-4.0d0*third*E(n)*dvy)*dtdy
              end do
#endif

              ! modif nimhd
              if(magohm.eq.0 )then
                 if(ntestDADM.eq.1) then
                    tcell=1.0d0
                 else
                    call temperature_eos(q(l,i,j,k,1),q(l,i,j,k,nvar),tcell,ht)
                 endif
                 bcell = A*A + B*B + C*C
                 sp0 = sp0 +(gamma-1.d0)*etaohmdiss(q(l,i,j,k,1),bcell,tcell)*jcentersquare(l,i,j,k)*dt
              endif
              ! fin modif nimhd
              
              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              C = C + sC0
              A = 0.5*(AL+AR)
              B = 0.5*(BL+BR)
#if USE_FLD==1
              do n=1,ngrp
                 E(n)=E(n)+sE0(n)
              end do
#endif
              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = AL     
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = MAX(smallp, qp(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qp(l,i,j,k,iC+n,1) = E(n) - dEx(n)
                 qp(l,i,j,k,iC+n,1) = max(small_er, qp(l,i,j,k,iC+n,1))
              end do
#endif

              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = AR     
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = MAX(smallp, qm(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qm(l,i,j,k,iC+n,1) = E(n) + dEx(n)
                 qm(l,i,j,k,iC+n,1) = max(small_er, qm(l,i,j,k,iC+n,1))
              end do
#endif

              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,iA,2) = A - dAy
              qp(l,i,j,k,iB,2) = BL     
              qp(l,i,j,k,iC,2) = C - dCy
              qp(l,i,j,k,ir,2) = MAX(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = MAX(smallp, qp(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qp(l,i,j,k,iC+n,2) = E(n) - dEy(n)
                 qp(l,i,j,k,iC+n,2) = max(small_er, qp(l,i,j,k,iC+n,2))
              end do
#endif

              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,iA,2) = A + dAy
              qm(l,i,j,k,iB,2) = BR     
              qm(l,i,j,k,iC,2) = C + dCy
              qm(l,i,j,k,ir,2) = MAX(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = MAX(smallp, qm(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qm(l,i,j,k,iC+n,2) = E(n) + dEy(n)
                 qm(l,i,j,k,iC+n,2) = max(small_er, qm(l,i,j,k,iC+n,2))
              end do
#endif

              ! Edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,3) = r + (+drx+dry)
              qRT(l,i,j,k,iu,3) = u + (+dux+duy)
              qRT(l,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(l,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(l,i,j,k,ip,3) = p + (+dpx+dpy)
              qRT(l,i,j,k,iC,3) = C + (+dCx+dCy)
              qRT(l,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(l,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(l,i,j,k,ir,3) = MAX(smallr, qRT(l,i,j,k,ir,3))
              qRT(l,i,j,k,ip,3) = MAX(smallp, qRT(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qRT(l,i,j,k,iC+n,3) = E(n) + (+dEx(n)+dEy(n))
                 qRT(l,i,j,k,iC+n,3) = max(small_er, qRT(l,i,j,k,iC+n,3))
              end do
#endif

              ! Edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,3) = r + (+drx-dry)
              qRB(l,i,j,k,iu,3) = u + (+dux-duy)
              qRB(l,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(l,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(l,i,j,k,ip,3) = p + (+dpx-dpy)
              qRB(l,i,j,k,iC,3) = C + (+dCx-dCy)
              qRB(l,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(l,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(l,i,j,k,ir,3) = MAX(smallr, qRB(l,i,j,k,ir,3))
              qRB(l,i,j,k,ip,3) = MAX(smallp, qRB(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qRB(l,i,j,k,iC+n,3) = E(n) + (+dEx(n)-dEy(n))
                 qRB(l,i,j,k,iC+n,3) = max(small_er, qRB(l,i,j,k,iC+n,3))
              end do
#endif

              ! Edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,3) = r + (-drx+dry)
              qLT(l,i,j,k,iu,3) = u + (-dux+duy)
              qLT(l,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(l,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(l,i,j,k,ip,3) = p + (-dpx+dpy)
              qLT(l,i,j,k,iC,3) = C + (-dCx+dCy)
              qLT(l,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(l,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(l,i,j,k,ir,3) = MAX(smallr, qLT(l,i,j,k,ir,3))
              qLT(l,i,j,k,ip,3) = MAX(smallp, qLT(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qLT(l,i,j,k,iC+n,3) = E(n) + (-dEx(n)+dEy(n))
                 qLT(l,i,j,k,iC+n,3) = max(small_er, qLT(l,i,j,k,iC+n,3))
              end do
#endif

              ! Edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,3) = r + (-drx-dry)
              qLB(l,i,j,k,iu,3) = u + (-dux-duy)
              qLB(l,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(l,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(l,i,j,k,ip,3) = p + (-dpx-dpy)
              qLB(l,i,j,k,iC,3) = C + (-dCx-dCy)
              qLB(l,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(l,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(l,i,j,k,ir,3) = MAX(smallr, qLB(l,i,j,k,ir,3))
              qLB(l,i,j,k,ip,3) = MAX(smallp, qLB(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qLB(l,i,j,k,iC+n,3) = E(n) + (-dEx(n)-dEy(n))
                 qLB(l,i,j,k,iC+n,3) = max(small_er, qLB(l,i,j,k,iC+n,3))
              end do
#endif

           END DO
        END DO
     END DO
  END DO

  ! Passive scalars
  DO n = 9+nrad, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 r   = q(l,i,j,k,n )              ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 drx = half * dq(l,i,j,k,n,1)     ! TVD slopes
                 dry = half * dq(l,i,j,k,n,2)
                 sr0 = -u*drx*dtdx -v*dry*dtdy    ! Source terms
                 r   = r + sr0                    ! Predicted state
                 qp(l,i,j,k,n,1) = r - drx        ! Right state
                 qm(l,i,j,k,n,1) = r + drx        ! Left state
                 qp(l,i,j,k,n,2) = r - dry        ! Top state
                 qm(l,i,j,k,n,2) = r + dry        ! Bottom state
              END DO
           END DO
        END DO
     END DO
  END DO

END SUBROUTINE trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==3
! modif nimhd
!SUBROUTINE trace3d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid)
SUBROUTINE trace3d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid,jcentersquare)
! fin modif nimhd
  USE amr_parameters
  USE hydro_parameters
  USE const
  use radiation_parameters,only:small_er
  use units_commons
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx, dy, dz, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB
  ! modif nimhd
  REAL(dp)::etaohmdiss,bcell,tcell,barotrop1D
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare
  ! fin modif nimhd

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ex
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ey
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ez

  INTEGER ::i, j, k, l, n, ht
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC 
  REAL(dp)::dtdx, dtdy, dtdz, smallp
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::FLL, FLR, FRL, FRR
  REAL(dp)::GLL, GLR, GRL, GRR
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::dry, duy, dvy, dwy, dpy, dAy, dbY, dCy
  REAL(dp)::drz, duz, dvz, dwz, dpz, dAz, dBz, dCz
  REAL(dp)::sr0, su0=0, sv0=0, sw0=0, sp0, sA0, sB0, sC0
  REAL(dp)::AL, AR, BL, BR, CL, CR
  REAL(dp)::dALy, dARy, dALz, dARz
  REAL(dp)::dBLx, dBRx, dBLz, dBRz
  REAL(dp)::dCLx, dCRx, dCLy, dCRy
  REAL(DP)::sAL0, sAR0, sBL0, sBR0, sCL0, sCR0

#if USE_FLD==1
  REAL(dp),ALLOCATABLE,DIMENSION(:)::E,dEx,dEy,dEz,sE0  
  if(ngrp>0)ALLOCATE(E(1:ngrp),dEx(1:ngrp),dEy(1:ngrp),dEz(1:ngrp),sE0(1:ngrp))
#endif
  
  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  smallp = smallr*smallc**2/gamma

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5; ia=6; ib=7; ic=8 

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           DO l = 1, ngrid
              v = 0.25*(q(l,i,j-1,k-1,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k-1,iv)+q(l,i,j,k,iv))
              w = 0.25*(q(l,i,j-1,k-1,iw)+q(l,i,j-1,k,iw)+q(l,i,j,k-1,iw)+q(l,i,j,k,iw))
              B = 0.5*(bf(l,i,j,k-1,2)+bf(l,i,j,k,2))
              C = 0.5*(bf(l,i,j-1,k,3)+bf(l,i,j,k,3))
              Ex(l,i,j,k)=v*C-w*B

              u = 0.25*(q(l,i-1,j,k-1,iu)+q(l,i-1,j,k,iu)+q(l,i,j,k-1,iu)+q(l,i,j,k,iu))
              w = 0.25*(q(l,i-1,j,k-1,iw)+q(l,i-1,j,k,iw)+q(l,i,j,k-1,iw)+q(l,i,j,k,iw))
              A = 0.5*(bf(l,i,j,k-1,1)+bf(l,i,j,k,1))
              C = 0.5*(bf(l,i-1,j,k,3)+bf(l,i,j,k,3))
              Ey(l,i,j,k)=w*A-u*C

              u = 0.25*(q(l,i-1,j-1,k,iu)+q(l,i-1,j,k,iu)+q(l,i,j-1,k,iu)+q(l,i,j,k,iu))
              v = 0.25*(q(l,i-1,j-1,k,iv)+q(l,i-1,j,k,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k,iv))
              A = 0.5*(bf(l,i,j-1,k,1)+bf(l,i,j,k,1))
              B = 0.5*(bf(l,i-1,j,k,2)+bf(l,i,j,k,2))
              Ez(l,i,j,k)=u*B-v*A
           END DO
        END DO
     END DO
  END DO

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)            
              p =    q(l,i,j,k,ip)
              A =    q(l,i,j,k,iA)
              B =    q(l,i,j,k,iB)
              C =    q(l,i,j,k,iC)            
#if USE_FLD==1
              do n = 1, ngrp
                 E(n) = q(l,i,j,k,iC+n)
              end do
#endif

              ! Face centered variables
              AL =  bf(l,i  ,j  ,k  ,1)
              AR =  bf(l,i+1,j  ,k  ,1)
              BL =  bf(l,i  ,j  ,k  ,2)
              BR =  bf(l,i  ,j+1,k  ,2)
              CL =  bf(l,i  ,j  ,k  ,3)
              CR =  bf(l,i  ,j  ,k+1,3)

              ! Cell centered TVD slopes in X, Y and Z directions
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)
              dBx = half * dq(l,i,j,k,iB,1)
              dCx = half * dq(l,i,j,k,iC,1)
#if USE_FLD==1
              do n = 1, ngrp
                 dEx(n) = half*dq(l,i,j,k,iC+n,1)
              end do
#endif

              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)
              dAy = half * dq(l,i,j,k,iA,2)
              dCy = half * dq(l,i,j,k,iC,2)
#if USE_FLD==1
              do n = 1, ngrp
                 dEy(n) = half*dq(l,i,j,k,iC+n,2)
              end do
#endif

              drz = half * dq(l,i,j,k,ir,3)
              duz = half * dq(l,i,j,k,iu,3)
              dvz = half * dq(l,i,j,k,iv,3)
              dwz = half * dq(l,i,j,k,iw,3)
              dpz = half * dq(l,i,j,k,ip,3)
              dAz = half * dq(l,i,j,k,iA,3)
              dBz = half * dq(l,i,j,k,iB,3)
#if USE_FLD==1
              do n = 1, ngrp
                 dEz(n) = half*dq(l,i,j,k,iC+n,3)
              end do
#endif

              ! Face centered TVD slopes in transverse directions
              dALy = half * dbf(l,i  ,j  ,k  ,1,1)
              dARy = half * dbf(l,i+1,j  ,k  ,1,1)
              dALz = half * dbf(l,i  ,j  ,k  ,1,2)
              dARz = half * dbf(l,i+1,j  ,k  ,1,2)

              dBLx = half * dbf(l,i  ,j  ,k  ,2,1)
              dBRx = half * dbf(l,i  ,j+1,k  ,2,1)
              dBLz = half * dbf(l,i  ,j  ,k  ,2,2)
              dBRz = half * dbf(l,i  ,j+1,k  ,2,2)

              dCLx = half * dbf(l,i  ,j  ,k  ,3,1)
              dCRx = half * dbf(l,i  ,j  ,k+1,3,1)
              dCLy = half * dbf(l,i  ,j  ,k  ,3,2)
              dCRy = half * dbf(l,i  ,j  ,k+1,3,2)

              ! Edge centered electric field in X, Y and Z directions
              ELL = Ex(l,i,j  ,k  )
              ELR = Ex(l,i,j  ,k+1)
              ERL = Ex(l,i,j+1,k  )
              ERR = Ex(l,i,j+1,k+1)

              FLL = Ey(l,i  ,j,k  )
              FLR = Ey(l,i  ,j,k+1)
              FRL = Ey(l,i+1,j,k  )
              FRR = Ey(l,i+1,j,k+1)

              GLL = Ez(l,i  ,j  ,k)
              GLR = Ez(l,i  ,j+1,k)
              GRL = Ez(l,i+1,j  ,k)
              GRR = Ez(l,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(GLR-GLL)*dtdy*half -(FLR-FLL)*dtdz*half
              sAR0 = +(GRR-GRL)*dtdy*half -(FRR-FRL)*dtdz*half
              sBL0 = -(GRL-GLL)*dtdx*half +(ELR-ELL)*dtdz*half
              sBR0 = -(GRR-GLR)*dtdx*half +(ERR-ERL)*dtdz*half
              sCL0 = +(FRL-FLL)*dtdx*half -(ERL-ELL)*dtdy*half
              sCR0 = +(FRR-FLR)*dtdx*half -(ERR-ELR)*dtdy*half
              
              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              CL = CL + sCL0
              CR = CR + sCR0
              
              ! Source terms (including transverse derivatives)
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy + (-w*drz-dwz*r)*dtdz 
              if(i_scheme.ne.i_induction)then
              su0 = (-u*dux-(dpx+B*dBx+C*dCx)/r)*dtdx + (-v*duy+B*dAy/r)*dtdy + (-w*duz+C*dAz/r)*dtdz 
#if USE_FLD==1
              do n = 1, ngrp
                 su0 = su0 - (third*dEx(n)/r)*dtdx
!!$!B+M
!!$                 su0 = su0 - (dEx(n)/r)*dtdx
              end do
#endif
              sv0 = (-u*dvx+A*dBx/r)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)/r)*dtdy + (-w*dvz+C*dBz/r)*dtdz
#if USE_FLD==1
              do n = 1, ngrp
                 sv0 = sv0 - (third*dEy(n)/r)*dtdy
!!$                 sv0 = sv0 - (dEy(n)/r)*dtdy
              end do
#endif
              sw0 = (-u*dwx+A*dCx/r)*dtdx + (-v*dwy+B*dCy/r)*dtdy + (-w*dwz-(dpz+A*dAz+B*dBz)/r)*dtdz 
#if USE_FLD==1
              do n = 1, ngrp
                 sw0 = sw0 - (third*dEz(n)/r)*dtdz
!!$                 sw0 = sw0 - (dEz(n)/r)*dtdz
              end do
#endif
              endif
              sp0 = (-u*dpx-dux*gamma*p)*dtdx + (-v*dpy-dvy*gamma*p)*dtdy + (-w*dpz-dwz*gamma*p)*dtdz
#if USE_FLD==1
              do n = 1, ngrp
!!$!B+M
!!$                 sE0(n) = (-1.0d0*E(n)*dux-u*dEx(n))*dtdx &
!!$                      & + (-1.0d0*E(n)*dvy-v*dEy(n))*dtdy &
!!$                      & + (-1.0d0*E(n)*dwz-w*dEz(n))*dtdz
!!$                 sE0(n) = (-2.0d0*E(n)*dux-u*dEx(n))*dtdx &
!!$                      & + (-2.0d0*E(n)*dvy-v*dEy(n))*dtdy &
!!$                      & + (-2.0d0*E(n)*dwz-w*dEz(n))*dtdz
                 sE0(n) = (-4.0d0*third*E(n)*dux-u*dEx(n))*dtdx &
                      & + (-4.0d0*third*E(n)*dvy-v*dEy(n))*dtdy &
                      & + (-4.0d0*third*E(n)*dwz-w*dEz(n))*dtdz
              end do
#endif

              ! modif nimhd
              if(magohm.eq.0 )then
                 if(ntestDADM.eq.1) then
                    tcell=1.0d0
                 else
                    call temperature_eos(q(l,i,j,k,1),q(l,i,j,k,nvar),tcell,ht)
                 endif
                 bcell = A*A + B*B + C*C
                 sp0 = sp0 +(gamma-1.d0)*etaohmdiss(q(l,i,j,k,1),bcell,tcell)*jcentersquare(l,i,j,k)*dt
              endif
              ! fin modif nimhd
              
              ! Cell-centered predicted states
              r = r + sr0
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              A = 0.5*(AL+AR)
              B = 0.5*(BL+BR)
              C = 0.5*(CL+CR)
#if USE_FLD==1
              do n = 1, ngrp
                 E(n) = E(n) + sE0(n)
              end do
#endif

              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = AL
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = MAX(smallp, qp(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qp(l,i,j,k,iC+n,1) = E(n) - dEx(n)
                 qp(l,i,j,k,iC+n,1) = max(small_er, qp(l,i,j,k,iC+n,1))
              end do
#endif
                 
              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = AR
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = MAX(smallp, qm(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qm(l,i,j,k,iC+n,1) = E(n) + dEx(n)
                 qm(l,i,j,k,iC+n,1) = max(small_er, qm(l,i,j,k,iC+n,1))
              end do
#endif
              
              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,iA,2) = A - dAy
              qp(l,i,j,k,iB,2) = BL
              qp(l,i,j,k,iC,2) = C - dCy
              qp(l,i,j,k,ir,2) = MAX(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = MAX(smallp, qp(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qp(l,i,j,k,iC+n,2) = E(n) - dEy(n)
                 qp(l,i,j,k,iC+n,2) = max(small_er, qp(l,i,j,k,iC+n,2))
              end do
#endif
              
              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,iA,2) = A + dAy
              qm(l,i,j,k,iB,2) = BR
              qm(l,i,j,k,iC,2) = C + dCy
              qm(l,i,j,k,ir,2) = MAX(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = MAX(smallp, qm(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qm(l,i,j,k,iC+n,2) = E(n) + dEy(n)
                 qm(l,i,j,k,iC+n,2) = max(small_er, qm(l,i,j,k,iC+n,2))
              end do
#endif
              
              ! Face averaged front state at back interface
              qp(l,i,j,k,ir,3) = r - drz
              qp(l,i,j,k,iu,3) = u - duz
              qp(l,i,j,k,iv,3) = v - dvz
              qp(l,i,j,k,iw,3) = w - dwz
              qp(l,i,j,k,ip,3) = p - dpz
              qp(l,i,j,k,iA,3) = A - dAz
              qp(l,i,j,k,iB,3) = B - dBz
              qp(l,i,j,k,iC,3) = CL
              qp(l,i,j,k,ir,3) = MAX(smallr, qp(l,i,j,k,ir,3))
              qp(l,i,j,k,ip,3) = MAX(smallp, qp(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qp(l,i,j,k,iC+n,3) = E(n) - dEz(n)
                 qp(l,i,j,k,iC+n,3) = max(small_er, qp(l,i,j,k,iC+n,3))
              end do
#endif

              ! Face averaged back state at front interface
              qm(l,i,j,k,ir,3) = r + drz
              qm(l,i,j,k,iu,3) = u + duz
              qm(l,i,j,k,iv,3) = v + dvz
              qm(l,i,j,k,iw,3) = w + dwz
              qm(l,i,j,k,ip,3) = p + dpz
              qm(l,i,j,k,iA,3) = A + dAz
              qm(l,i,j,k,iB,3) = B + dBz
              qm(l,i,j,k,iC,3) = CR
              qm(l,i,j,k,ir,3) = MAX(smallr, qm(l,i,j,k,ir,3))
              qm(l,i,j,k,ip,3) = MAX(smallp, qm(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qm(l,i,j,k,iC+n,3) = E(n) + dEz(n)
                 qm(l,i,j,k,iC+n,3) = max(small_er, qm(l,i,j,k,iC+n,3))
              end do
#endif

              ! X-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,1) = r + (+dry+drz)
              qRT(l,i,j,k,iu,1) = u + (+duy+duz)
              qRT(l,i,j,k,iv,1) = v + (+dvy+dvz)
              qRT(l,i,j,k,iw,1) = w + (+dwy+dwz)
              qRT(l,i,j,k,ip,1) = p + (+dpy+dpz)
              qRT(l,i,j,k,iA,1) = A + (+dAy+dAz)
              qRT(l,i,j,k,iB,1) = BR+ (   +dBRz)
              qRT(l,i,j,k,iC,1) = CR+ (+dCRy   )
              qRT(l,i,j,k,ir,1) = MAX(smallr, qRT(l,i,j,k,ir,1))
              qRT(l,i,j,k,ip,1) = MAX(smallp, qRT(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qRT(l,i,j,k,iC+n,1) = E(n) + (+dEy(n)+dEz(n))
                 qRT(l,i,j,k,iC+n,1) = max(small_er, qRT(l,i,j,k,iC+n,1))
              end do
#endif

              ! X-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,1) = r + (+dry-drz)
              qRB(l,i,j,k,iu,1) = u + (+duy-duz)
              qRB(l,i,j,k,iv,1) = v + (+dvy-dvz)
              qRB(l,i,j,k,iw,1) = w + (+dwy-dwz)
              qRB(l,i,j,k,ip,1) = p + (+dpy-dpz)
              qRB(l,i,j,k,iA,1) = A + (+dAy-dAz)
              qRB(l,i,j,k,iB,1) = BR+ (   -dBRz)
              qRB(l,i,j,k,iC,1) = CL+ (+dCLy   )
              qRB(l,i,j,k,ir,1) = MAX(smallr, qRB(l,i,j,k,ir,1))
              qRB(l,i,j,k,ip,1) = MAX(smallp, qRB(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qRB(l,i,j,k,iC+n,1) = E(n) + (+dEy(n)-dEz(n))
                 qRB(l,i,j,k,iC+n,1) = max(small_er, qRB(l,i,j,k,iC+n,1))
              end do
#endif

              ! X-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,1) = r + (-dry+drz)
              qLT(l,i,j,k,iu,1) = u + (-duy+duz)
              qLT(l,i,j,k,iv,1) = v + (-dvy+dvz)
              qLT(l,i,j,k,iw,1) = w + (-dwy+dwz)
              qLT(l,i,j,k,ip,1) = p + (-dpy+dpz)
              qLT(l,i,j,k,iA,1) = A + (-dAy+dAz)
              qLT(l,i,j,k,iB,1) = BL+ (   +dBLz)
              qLT(l,i,j,k,iC,1) = CR+ (-dCRy   )
              qLT(l,i,j,k,ir,1) = MAX(smallr, qLT(l,i,j,k,ir,1))
              qLT(l,i,j,k,ip,1) = MAX(smallp, qLT(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qLT(l,i,j,k,iC+n,1) = E(n) + (-dEy(n)+dEz(n))
                 qLT(l,i,j,k,iC+n,1) = max(small_er, qLT(l,i,j,k,iC+n,1))
              end do
#endif

              ! X-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,1) = r + (-dry-drz)
              qLB(l,i,j,k,iu,1) = u + (-duy-duz)
              qLB(l,i,j,k,iv,1) = v + (-dvy-dvz)
              qLB(l,i,j,k,iw,1) = w + (-dwy-dwz)
              qLB(l,i,j,k,ip,1) = p + (-dpy-dpz)
              qLB(l,i,j,k,iA,1) = A + (-dAy-dAz)
              qLB(l,i,j,k,iB,1) = BL+ (   -dBLz)
              qLB(l,i,j,k,iC,1) = CL+ (-dCLy   )
              qLB(l,i,j,k,ir,1) = MAX(smallr, qLB(l,i,j,k,ir,1))
              qLB(l,i,j,k,ip,1) = MAX(smallp, qLB(l,i,j,k,ip,1))
#if USE_FLD==1
              do n = 1, ngrp
                 qLB(l,i,j,k,iC+n,1) = E(n) + (-dEy(n)-dEz(n))
                 qLB(l,i,j,k,iC+n,1) = max(small_er, qLB(l,i,j,k,iC+n,1))
              end do
#endif

              ! Y-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,2) = r + (+drx+drz)
              qRT(l,i,j,k,iu,2) = u + (+dux+duz)
              qRT(l,i,j,k,iv,2) = v + (+dvx+dvz)
              qRT(l,i,j,k,iw,2) = w + (+dwx+dwz)
              qRT(l,i,j,k,ip,2) = p + (+dpx+dpz)
              qRT(l,i,j,k,iA,2) = AR+ (   +dARz)
              qRT(l,i,j,k,iB,2) = B + (+dBx+dBz)
              qRT(l,i,j,k,iC,2) = CR+ (+dCRx   )
              qRT(l,i,j,k,ir,2) = MAX(smallr, qRT(l,i,j,k,ir,2))
              qRT(l,i,j,k,ip,2) = MAX(smallp, qRT(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qRT(l,i,j,k,iC+n,2) = E(n) + (+dEx(n)+dEz(n))
                 qRT(l,i,j,k,iC+n,2) = max(small_er, qRT(l,i,j,k,iC+n,2))
              end do
#endif

              ! Y-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,2) = r + (+drx-drz)
              qRB(l,i,j,k,iu,2) = u + (+dux-duz)
              qRB(l,i,j,k,iv,2) = v + (+dvx-dvz)
              qRB(l,i,j,k,iw,2) = w + (+dwx-dwz)
              qRB(l,i,j,k,ip,2) = p + (+dpx-dpz)
              qRB(l,i,j,k,iA,2) = AR+ (   -dARz)
              qRB(l,i,j,k,iB,2) = B + (+dBx-dBz)
              qRB(l,i,j,k,iC,2) = CL+ (+dCLx   )
              qRB(l,i,j,k,ir,2) = MAX(smallr, qRB(l,i,j,k,ir,2))
              qRB(l,i,j,k,ip,2) = MAX(smallp, qRB(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qRB(l,i,j,k,iC+n,2) = E(n) + (+dEx(n)-dEz(n))
                 qRB(l,i,j,k,iC+n,2) = max(small_er, qRB(l,i,j,k,iC+n,2))
              end do
#endif

              ! Y-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,2) = r + (-drx+drz)
              qLT(l,i,j,k,iu,2) = u + (-dux+duz)
              qLT(l,i,j,k,iv,2) = v + (-dvx+dvz)
              qLT(l,i,j,k,iw,2) = w + (-dwx+dwz)
              qLT(l,i,j,k,ip,2) = p + (-dpx+dpz)
              qLT(l,i,j,k,iA,2) = AL+ (   +dALz)
              qLT(l,i,j,k,iB,2) = B + (-dBx+dBz)
              qLT(l,i,j,k,iC,2) = CR+ (-dCRx   )
              qLT(l,i,j,k,ir,2) = MAX(smallr, qLT(l,i,j,k,ir,2))
              qLT(l,i,j,k,ip,2) = MAX(smallp, qLT(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qLT(l,i,j,k,iC+n,2) = E(n) + (-dEx(n)+dEz(n))
                 qLT(l,i,j,k,iC+n,2) = max(small_er, qLT(l,i,j,k,iC+n,2))
              end do
#endif

              ! Y-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,2) = r + (-drx-drz)
              qLB(l,i,j,k,iu,2) = u + (-dux-duz)
              qLB(l,i,j,k,iv,2) = v + (-dvx-dvz)
              qLB(l,i,j,k,iw,2) = w + (-dwx-dwz)
              qLB(l,i,j,k,ip,2) = p + (-dpx-dpz)
              qLB(l,i,j,k,iA,2) = AL+ (   -dALz)
              qLB(l,i,j,k,iB,2) = B + (-dBx-dBz)
              qLB(l,i,j,k,iC,2) = CL+ (-dCLx   )
              qLB(l,i,j,k,ir,2) = MAX(smallr, qLB(l,i,j,k,ir,2))
              qLB(l,i,j,k,ip,2) = MAX(smallp, qLB(l,i,j,k,ip,2))
#if USE_FLD==1
              do n = 1, ngrp
                 qLB(l,i,j,k,iC+n,2) = E(n) + (-dEx(n)-dEz(n))
                 qLB(l,i,j,k,iC+n,2) = max(small_er, qLB(l,i,j,k,iC+n,2))
              end do
#endif

              ! Z-edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,3) = r + (+drx+dry)
              qRT(l,i,j,k,iu,3) = u + (+dux+duy)
              qRT(l,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(l,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(l,i,j,k,ip,3) = p + (+dpx+dpy)
              qRT(l,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(l,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(l,i,j,k,iC,3) = C + (+dCx+dCy)
              qRT(l,i,j,k,ir,3) = MAX(smallr, qRT(l,i,j,k,ir,3))
              qRT(l,i,j,k,ip,3) = MAX(smallp, qRT(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qRT(l,i,j,k,iC+n,3) = E(n) + (+dEx(n)+dEy(n))
                 qRT(l,i,j,k,iC+n,3) = max(small_er, qRT(l,i,j,k,iC+n,3))
              end do
#endif

              ! Z-edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,3) = r + (+drx-dry)
              qRB(l,i,j,k,iu,3) = u + (+dux-duy)
              qRB(l,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(l,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(l,i,j,k,ip,3) = p + (+dpx-dpy)
              qRB(l,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(l,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(l,i,j,k,iC,3) = C + (+dCx-dCy)
              qRB(l,i,j,k,ir,3) = MAX(smallr, qRB(l,i,j,k,ir,3))
              qRB(l,i,j,k,ip,3) = MAX(smallp, qRB(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qRB(l,i,j,k,iC+n,3) = E(n) + (+dEx(n)-dEy(n))
                 qRB(l,i,j,k,iC+n,3) = max(small_er, qRB(l,i,j,k,iC+n,3))
              end do
#endif

              ! Z-edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,3) = r + (-drx+dry)
              qLT(l,i,j,k,iu,3) = u + (-dux+duy)
              qLT(l,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(l,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(l,i,j,k,ip,3) = p + (-dpx+dpy)
              qLT(l,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(l,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(l,i,j,k,iC,3) = C + (-dCx+dCy)
              qLT(l,i,j,k,ir,3) = MAX(smallr, qLT(l,i,j,k,ir,3))
              qLT(l,i,j,k,ip,3) = MAX(smallp, qLT(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qLT(l,i,j,k,iC+n,3) = E(n) + (-dEx(n)+dEy(n))
                 qLT(l,i,j,k,iC+n,3) = max(small_er, qLT(l,i,j,k,iC+n,3))
              end do
#endif

              ! Z-edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,3) = r + (-drx-dry)
              qLB(l,i,j,k,iu,3) = u + (-dux-duy)
              qLB(l,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(l,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(l,i,j,k,ip,3) = p + (-dpx-dpy)
              qLB(l,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(l,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(l,i,j,k,iC,3) = C + (-dCx-dCy)
              qLB(l,i,j,k,ir,3) = MAX(smallr, qLB(l,i,j,k,ir,3))
              qLB(l,i,j,k,ip,3) = MAX(smallp, qLB(l,i,j,k,ip,3))
#if USE_FLD==1
              do n = 1, ngrp
                 qLB(l,i,j,k,iC+n,3) = E(n) + (-dEx(n)-dEy(n))
                 qLB(l,i,j,k,iC+n,3) = max(small_er, qLB(l,i,j,k,iC+n,3))
              end do
#endif

           END DO
        END DO
     END DO
  END DO

  ! Passive scalars
  DO n = 9+nrad, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 r   = q(l,i,j,k,n )            ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 w   = q(l,i,j,k,iw)
                 drx = half * dq(l,i,j,k,n,1)   ! TVD slopes
                 dry = half * dq(l,i,j,k,n,2)
                 drz = half * dq(l,i,j,k,n,3)
                 sr0 = -u*drx*dtdx -v*dry*dtdy -w*drz*dtdz   ! Source terms
                 r   = r + sr0                  ! Predicted state
                 qp(l,i,j,k,n,1) = r - drx      ! Right state
                 qm(l,i,j,k,n,1) = r + drx      ! Left state
                 qp(l,i,j,k,n,2) = r - dry      ! Top state
                 qm(l,i,j,k,n,2) = r + dry      ! Bottom state
                 qp(l,i,j,k,n,3) = r - drz      ! Front state
                 qm(l,i,j,k,n,3) = r + drz      ! Back state
              END DO
           END DO
        END DO
     END DO
  END DO

END SUBROUTINE trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, &
     &                ln ,lt1,lt2,bn ,bt1,bt2, flx,tmp,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2,bn,bt1,bt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp),dimension(1:nvar)::qleft,qright,fgdnv
  REAL(dp)::zero_flux, bn_mean, entho

  xdim=ln-1
  entho=one/(gamma-one)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
           
              ! Enforce continuity for normal magnetic field
              bn_mean = half*(qm(l,i,j,k,bn,xdim)+qp(l,i,j,k,bn,xdim))

              ! Left state
              qleft (1) = qm(l,i,j,k,1  ,xdim)      ! Mass density
              qleft (2) = qm(l,i,j,k,5  ,xdim)      ! Pressure
              qleft (3) = qm(l,i,j,k,ln ,xdim)      ! Normal velocity
              qleft (4) = bn_mean                   ! Normal magnetic field
              qleft (5) = qm(l,i,j,k,lt1,xdim)      ! Tangential velocity 1
              qleft (6) = qm(l,i,j,k,bt1,xdim)      ! Tangential magnetic field 1
              qleft (7) = qm(l,i,j,k,lt2,xdim)      ! Tangential velocity 2
              qleft (8) = qm(l,i,j,k,bt2,xdim)      ! Tangential magnetic field 2

              ! Right state
              qright(1) = qp(l,i,j,k,1  ,xdim)      ! Mass density
              qright(2) = qp(l,i,j,k,5  ,xdim)      ! Pressure
              qright(3) = qp(l,i,j,k,ln ,xdim)      ! Normal velocity
              qright(4) = bn_mean                   ! Normal magnetic field
              qright(5) = qp(l,i,j,k,lt1,xdim)      ! Tangential velocity 1
              qright(6) = qp(l,i,j,k,bt1,xdim)      ! Tangential magnetic field 1
              qright(7) = qp(l,i,j,k,lt2,xdim)      ! Tangential velocity 2
              qright(8) = qp(l,i,j,k,bt2,xdim)      ! Tangential magnetic field 2

              ! Other advected quantities
              do n = 9, nvar
                 qleft (n) = qm(l,i,j,k,n,xdim)
                 qright(n) = qp(l,i,j,k,n,xdim)    
              end do
          
              ! Solve 1D Riemann problem
              zero_flux = one
              IF(i_scheme.NE.i_induction)THEN
              SELECT CASE (i_riemann)
              CASE (i_roe)
                 CALL athena_roe    (qleft,qright,fgdnv,zero_flux)
              CASE (i_llf)
                 CALL lax_friedrich (qleft,qright,fgdnv,zero_flux)
              CASE (i_hll)
                 if( ( (qright(4)**2+qright(6)**2+qright(8)**2)/qright(2) .gt. switch_solv .or. &
                      & (qleft(4)**2+qleft(6)**2+qleft(8)**2)/qleft(2) .gt. switch_solv  ) .or. &
                      & ( qleft(1)/qright(1) .gt. switch_solv) .or. (qleft(1)/qright(1) .lt. 1.0d0/switch_solv) )  then
                    CALL lax_friedrich(qleft,qright,fgdnv,zero_flux) 
                 else
                    CALL hll(qleft,qright,fgdnv)
                 endif
                !CALL hll(qleft,qright,fgdnv)
              CASE (i_hlld)
                 if( ( (qright(4)**2+qright(6)**2+qright(8)**2)/qright(2) .gt. switch_solv .or.&
                      & (qleft(4)**2+qleft(6)**2+qleft(8)**2)/qleft(2) .gt. switch_solv  ) .or.&
                      & ( qleft(1)/qright(1) .gt. switch_solv) .or. (qleft(1)/qright(1) .lt. 1.0d0/switch_solv) )  then
                    CALL lax_friedrich(qleft,qright,fgdnv,zero_flux) 
                 else
                    CALL hlld(qleft,qright,fgdnv)
                 endif
                 !CALL hlld(qleft,qright,fgdnv)
              CASE (i_upwind)
                 CALL lax_friedrich (qleft,qright,fgdnv,zero_flux)
              CASE (i_hydro)
                 CALL hydro_acoustic(qleft,qright,fgdnv)
              CASE DEFAULT
                 write(*,*)'unknown riemann solver'
                 call clean_stop
              END SELECT
              ELSE
                 CALL upwind(qleft,qright,fgdnv,zero_flux)
              ENDIF
           
              ! Output fluxes
              flx(l,i,j,k,1  ) = fgdnv(1)       ! Mass density
              flx(l,i,j,k,5  ) = fgdnv(2)       ! Total energy
              flx(l,i,j,k,ln ) = fgdnv(3)       ! Normal momentum
              flx(l,i,j,k,bn ) = fgdnv(4)       ! Normal magnetic field
              flx(l,i,j,k,lt1) = fgdnv(5)       ! Transverse momentum 1
              flx(l,i,j,k,bt1) = fgdnv(6)       ! Transverse magnetic field 1
              flx(l,i,j,k,lt2) = fgdnv(7)       ! Transverse momentum 2
              flx(l,i,j,k,bt2) = fgdnv(8)       ! Transverse magnetic field 2

              ! Other advected quantities
              do n = 9, nvar
                 flx(l,i,j,k,n) = fgdnv(n)
              end do

              ! Normal velocity estimate
              tmp(l,i,j,k,1) = half*(qleft(3)+qright(3))
              ! Internal energy flux
              if(fgdnv(1)>zero)then
                 tmp(l,i,j,k,2) = qleft (2)/qleft (1)*entho*fgdnv(1)
              else
                 tmp(l,i,j,k,2) = qright(2)/qright(1)*entho*fgdnv(1)
              end if

           end do
        end do
     end do
  end do
  
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cmp_mag_flx(qRT,irt1,irt2,jrt1,jrt2,krt1,krt2, &
       &               qRB,irb1,irb2,jrb1,jrb2,krb1,krb2, &
       &               qLT,ilt1,ilt2,jlt1,jlt2,klt1,klt2, &
       &               qLB,ilb1,ilb2,jlb1,jlb2,klb1,klb2, &
       &                   ilo ,ihi ,jlo ,jhi ,klo ,khi , &
       &                   lp1 ,lp2 ,lor ,bp1 ,bp2 ,bor ,emf,ngrid)
  ! 2D Riemann solver to compute EMF at cell edges
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  ! indices of the 2 planar velocity lp1 and lp2 and the orthogonal one,
  ! lor and idem for the magnetic field
  INTEGER ::lp1,lp2,lor,bp1,bp2,bor
  INTEGER ::irt1,irt2,jrt1,jrt2,krt1,krt2
  INTEGER ::irb1,irb2,jrb1,jrb2,krb1,krb2
  INTEGER ::ilt1,ilt2,jlt1,jlt2,klt1,klt2
  INTEGER ::ilb1,ilb2,jlb1,jlb2,klb1,klb2
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  REAL(dp),DIMENSION(1:nvector,irt1:irt2,jrt1:jrt2,krt1:krt2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,irb1:irb2,jrb1:jrb2,krb1:krb2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,ilt1:ilt2,jlt1:jlt2,klt1:klt2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2,1:nvar,1:3)::qLB

  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2):: emf

  ! local variables
  INTEGER ::i, j, k, n, l, idim, xdim, m
  REAL(dp),DIMENSION(1:nvector,1:nvar)::qLL,qRL,qLR,qRR
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fmean_x,fmean_y,qtmp
  REAL(dp) :: ELL,ERL,ELR,ERR,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  REAL(dp) :: zero_flux,E
  REAL(dp) :: cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  REAL(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  REAL(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  REAL(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR
  REAL(dp) :: ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  REAL(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  REAL(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  REAL(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  REAL(dp) :: AstarT,AstarB,BstarR,BstarL

  xdim = lor - 1

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi

           ! Density
           DO l = 1, ngrid
              qLL (l,1) = qRT(l,i,j,k,1,xdim)
              qRL (l,1) = qLT(l,i,j,k,1,xdim)
              qLR (l,1) = qRB(l,i,j,k,1,xdim)
              qRR (l,1) = qLB(l,i,j,k,1,xdim)
           END DO           

           ! Pressure 
           DO l = 1, ngrid
              qLL (l,2) = qRT(l,i,j,k,5,xdim)
              qRL (l,2) = qLT(l,i,j,k,5,xdim)
              qLR (l,2) = qRB(l,i,j,k,5,xdim)
              qRR (l,2) = qLB(l,i,j,k,5,xdim)
           END DO

           ! First parallel velocity 
           DO l = 1, ngrid
              qLL (l,3) = qRT(l,i,j,k,lp1,xdim)
              qRL (l,3) = qLT(l,i,j,k,lp1,xdim)
              qLR (l,3) = qRB(l,i,j,k,lp1,xdim)
              qRR (l,3) = qLB(l,i,j,k,lp1,xdim)
           END DO

           ! Second parallel velocity 
           DO l = 1, ngrid
              qLL (l,4) = qRT(l,i,j,k,lp2,xdim)
              qRL (l,4) = qLT(l,i,j,k,lp2,xdim)
              qLR (l,4) = qRB(l,i,j,k,lp2,xdim)
              qRR (l,4) = qLB(l,i,j,k,lp2,xdim)
           END DO

           ! First parallel magnetic field (enforce continuity)
           DO l = 1, ngrid
              qLL (l,6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qRL (l,6) = half*(qRT(l,i,j,k,bp1,xdim)+qLT(l,i,j,k,bp1,xdim))
              qLR (l,6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))
              qRR (l,6) = half*(qRB(l,i,j,k,bp1,xdim)+qLB(l,i,j,k,bp1,xdim))
           END DO

           ! Second parallel magnetic field (enforce continuity)
           DO l = 1, ngrid
              qLL (l,7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRL (l,7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))
              qLR (l,7) = half*(qRT(l,i,j,k,bp2,xdim)+qRB(l,i,j,k,bp2,xdim))
              qRR (l,7) = half*(qLT(l,i,j,k,bp2,xdim)+qLB(l,i,j,k,bp2,xdim))
           END DO

           ! Orthogonal velocity 
           DO l = 1, ngrid
              qLL (l,5) = qRT(l,i,j,k,lor,xdim)
              qRL (l,5) = qLT(l,i,j,k,lor,xdim)
              qLR (l,5) = qRB(l,i,j,k,lor,xdim)
              qRR (l,5) = qLB(l,i,j,k,lor,xdim)
           END DO

           ! Orthogonal magnetic Field
           DO l = 1, ngrid
              qLL (l,8) = qRT(l,i,j,k,bor,xdim)
              qRL (l,8) = qLT(l,i,j,k,bor,xdim)
              qLR (l,8) = qRB(l,i,j,k,bor,xdim)
              qRR (l,8) = qLB(l,i,j,k,bor,xdim)
           END DO

#if USE_FLD==1
           ! Radiative energy
           do n = 1, ngrp
              DO l = 1, ngrid
                 qLL (l,8+n) = qRT(l,i,j,k,8+n,xdim)
                 qRL (l,8+n) = qLT(l,i,j,k,8+n,xdim)
                 qLR (l,8+n) = qRB(l,i,j,k,8+n,xdim)
                 qRR (l,8+n) = qLB(l,i,j,k,8+n,xdim)
              END DO
           end do
#endif

           ! Compute final fluxes
            DO l = 1, ngrid 

               ! vx*by - vy*bx at the four edge centers
               ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
               ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
               ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)  
               ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6) 

               if(i_riemann2d == i_hlld)then

                  
                  rLL=qLL(l,1); pLL=qLL(l,2); uLL=qLL(l,3); vLL=qLL(l,4); ALL=qLL(l,6); BLL=qLL(l,7) ; CLL=qLL(l,8) 
                  rLR=qLR(l,1); pLR=qLR(l,2); uLR=qLR(l,3); vLR=qLR(l,4); ALR=qLR(l,6); BLR=qLR(l,7) ; CLR=qLR(l,8) 
                  rRL=qRL(l,1); pRL=qRL(l,2); uRL=qRL(l,3); vRL=qRL(l,4); ARL=qRL(l,6); BRL=qRL(l,7) ; CRL=qRL(l,8) 
                  rRR=qRR(l,1); pRR=qRR(l,2); uRR=qRR(l,3); vRR=qRR(l,4); ARR=qRR(l,6); BRR=qRR(l,7) ; CRR=qRR(l,8) 

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLL(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLR(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRL(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRR(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLL(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLR(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRL(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRR(l,8+n)
                  end do
#endif
                  call find_speed_fast(qtmp,cfastRRy)

                  SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                  SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
                  SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
                  ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)

                  ELL=uLL*BLL-vLL*ALL
                  ELR=uLR*BLR-vLR*ALR
                  ERL=uRL*BRL-vRL*ARL
                  ERR=uRR*BRR-vRR*ARR

                  PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
                  PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
                  PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
                  PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
                  
                  rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
                  rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
                  rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
                  rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
                  
                  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
                  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
                  
                  rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
                  rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
                  rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
                  EstarLLx=ustar*BstarLL-vLL  *ALL
                  EstarLLy=uLL  *BLL    -vstar*AstarLL
                  EstarLL =ustar*BstarLL-vstar*AstarLL

                  rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
                  rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
                  rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
                  EstarLRx=ustar*BstarLR-vLR  *ALR
                  EstarLRy=uLR  *BLR    -vstar*AstarLR
                  EstarLR =ustar*BstarLR-vstar*AstarLR

                  rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
                  rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
                  rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
                  EstarRLx=ustar*BstarRL-vRL  *ARL
                  EstarRLy=uRL  *BRL    -vstar*AstarRL
                  EstarRL =ustar*BstarRL-vstar*AstarRL

                  rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
                  rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
                  rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
                  EstarRRx=ustar*BstarRR-vRR  *ARR
                  EstarRRy=uRR  *BRR    -vstar*AstarRR
                  EstarRR =ustar*BstarRR-vstar*AstarRR

                  calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       &       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
                  calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       &       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
                  calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                       &       abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
                  calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                       &       abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
                  SAL=min(ustar-calfvenL,zero); SAR=max(ustar+calfvenR,zero)
                  SAB=min(vstar-calfvenB,zero); SAT=max(vstar+calfvenT,zero)
                  AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
                  BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

                  if(SB>0d0)then
                     if(SL>0d0)then
                     E=ELL
                     else if(SR<0d0)then
                     E=ERL
                     else
                     E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
                     endif
                  else if (ST<0d0)then
                     if(SL>0d0)then
                     E=ELR
                     else if(SR<0d0)then
                     E=ERR
                     else
                     E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
                     endif
                  else if(SL>0d0)then
                     E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
                  else if(SR<0d0)then
                     E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
                  else
                     E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
                          & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
                  endif
                  
                  emf(l,i,j,k) = E

               else if(i_riemann2d== i_hll)then

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLL(l,8+n)
                  end do
#endif
                  vLLx=qtmp(3); call find_speed_fast(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLR(l,8+n)
                  end do
#endif
                  vLRx=qtmp(3); call find_speed_fast(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRL(l,8+n)
                  end do
#endif
                  vRLx=qtmp(3); call find_speed_fast(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRR(l,8+n)
                  end do
#endif
                  vRRx=qtmp(3); call find_speed_fast(qtmp,cRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLl(l,8+n)
                  end do
#endif
                  vLLy=qtmp(3); call find_speed_fast(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLR(l,8+n)
                  end do
#endif
                  vLRy=qtmp(3); call find_speed_fast(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRL(l,8+n)
                  end do
#endif
                  vRLy=qtmp(3); call find_speed_fast(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRR(l,8+n)
                  end do
#endif
                  vRRy=qtmp(3); call find_speed_fast(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else if (i_riemann2d== i_hlla)then

                  ! Compute 4 Alfven velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLL(l,8+n)
                  end do
#endif
                  vLLx=qtmp(3); call find_speed_alfven(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qLR(l,8+n)
                  end do
#endif
                  vLRx=qtmp(3); call find_speed_alfven(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRL(l,8+n)
                  end do
#endif
                  vRLx=qtmp(3); call find_speed_alfven(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qtmp(8+n)=qRR(l,8+n)
                  end do
#endif
                  vRRx=qtmp(3); call find_speed_alfven(qtmp,cRRx)

                  ! Compute 4 Alfven relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  vLLy=qtmp(3); call find_speed_alfven(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  vLRy=qtmp(3); call find_speed_alfven(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  vRLy=qtmp(3); call find_speed_alfven(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  vRRy=qtmp(3); call find_speed_alfven(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else

                  ! find the average value of E
                  E = forth*(ELL+ERL+ELR+ERR)
                  
                  ! call the first solver in the x direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qLR(l,1))
                  qright(1) = half*(qRR(l,1)+qRL(l,1))

                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qLR(l,2))
                  qright(2) = half*(qRR(l,2)+qRL(l,2))
                  
                  ! vt1 becomes normal velocity
                  qleft (3) = half*(qLL(l,3)+qLR(l,3))
                  qright(3) = half*(qRR(l,3)+qRL(l,3))
                  
                  ! bt1 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,6)+qLR(l,6))
                  qright(4) = half*(qRR(l,6)+qRL(l,6))
                  
                  ! vt2 becomes transverse velocity field
                  qleft (5) = half*(qLL(l,4)+qLR(l,4))
                  qright(5) = half*(qRR(l,4)+qRL(l,4))
                  
                  ! bt2 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,7)+qLR(l,7))
                  qright(6) = half*(qRR(l,7)+qRL(l,7))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qLR(l,5))
                  qright(7) = half*(qRR(l,5)+qRL(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qLR(l,8))
                  qright(8) = half*(qRR(l,8)+qRL(l,8))

#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qleft (8+n) = half*(qLL(l,8+n)+qLR(l,8+n))
                     qright(8+n) = half*(qRR(l,8+n)+qRL(l,8+n))
                  end do
                  ! Internal energy (dummy argument)
                  qleft (nvar) = 0.0d0!half*(qLL(l,8+n)+qRL(l,8+n))
                  qright(nvar) = 0.0d0!half*(qRR(l,8+n)+qLR(l,8+n))
#endif
                  
                  zero_flux = 0.0
                  SELECT CASE (i_riemann2d)
                  CASE (i_roe)
                     CALL athena_roe   (qleft,qright,fmean_x,zero_flux)
                  CASE (i_llf)
                     CALL lax_friedrich(qleft,qright,fmean_x,zero_flux)
                  CASE (i_upwind)
                     CALL upwind       (qleft,qright,fmean_x,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     call clean_stop
                  END SELECT

                  ! call the second solver in the y direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qRL(l,1))
                  qright(1) = half*(qRR(l,1)+qLR(l,1))
                  
                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qRL(l,2))
                  qright(2) = half*(qRR(l,2)+qLR(l,2))
                  
                  ! vt2 becomes normal velocity
                  qleft (3) = half*(qLL(l,4)+qRL(l,4))
                  qright(3) = half*(qRR(l,4)+qLR(l,4))
                  
                  ! bt2 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,7)+qRL(l,7))
                  qright(4) = half*(qRR(l,7)+qLR(l,7))
                  
                  ! vt1 becomes transverse velocity field 
                  qleft (5) = half*(qLL(l,3)+qRL(l,3))
                  qright(5) = half*(qRR(l,3)+qLR(l,3))
                  
                  ! bt1 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,6)+qRL(l,6))
                  qright(6) = half*(qRR(l,6)+qLR(l,6))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qRL(l,5))
                  qright(7) = half*(qRR(l,5)+qLR(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qRL(l,8))
                  qright(8) = half*(qRR(l,8)+qLR(l,8))
                  
#if USE_FLD==1
                  ! Radiative energy
                  do n = 1, ngrp
                     qleft (8+n) = half*(qLL(l,8+n)+qRL(l,8+n))
                     qright(8+n) = half*(qRR(l,8+n)+qLR(l,8+n))
                  end do
                  ! Internal energy (dummy argument)
                  qleft (nvar) = 0.0d0!half*(qLL(l,8+n)+qRL(l,8+n))
                  qright(nvar) = 0.0d0!half*(qRR(l,8+n)+qLR(l,8+n))
#endif

                  zero_flux = 0.
                  SELECT CASE (i_riemann2d)
                  CASE (i_roe)
                     CALL athena_roe   (qleft,qright,fmean_y,zero_flux)
                  CASE (i_llf)
                     CALL lax_friedrich(qleft,qright,fmean_y,zero_flux)
                  CASE (i_upwind)
                     CALL upwind       (qleft,qright,fmean_y,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     call clean_stop
                  END SELECT
                  
                  ! compute the final value of E including the 2D diffusive
                  ! terms that ensure stability 
                  emf(l,i,j,k) = E + (fmean_x(6) - fmean_y(6))
                  
               endif

            END DO
        END DO
     END DO
  END DO


END SUBROUTINE cmp_mag_flx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,bf,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use radiation_parameters,only:small_er,is_radiative_energy
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  

  integer ::i, j, k, l, n, idim
  real(dp)::eint, smalle, smallp, etot
  real(dp),dimension(1:nvector),save::eken,emag
#if USE_FLD==1 || USE_M_1==1
  real(dp),dimension(1:nvector),save::erad_loc
#endif

  ! EOS
  real(dp)  :: pp_eos

  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma

  ! Store face centered magnetic field
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           DO l = 1, ngrid
              if(i<=iu2)then
                 bf(l,i,j,k,1) = uin(l,i,j,k,6)
              else
                 bf(l,i,j,k,1) = uin(l,i-1,j,k,nvar+1)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           DO l = 1, ngrid
              if(j<=ju2)then
                 bf(l,i,j,k,2) = uin(l,i,j,k,7)
              else
                 bf(l,i,j,k,2) = uin(l,i,j-1,k,nvar+2)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           DO l = 1, ngrid
              if(k<=ku2)then
                 bf(l,i,j,k,3) = uin(l,i,j,k,8)
              else
                 bf(l,i,j,k,3) = uin(l,i,j,k-1,nvar+3)
              endif
           END DO
        end do
     end do
  end do

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do

           ! Debug
           if(debug)then
              do l = 1, ngrid
                 if(uin(l,i,j,k,1).le.smallr)then
                    write(*,*)'negative density'
                    write(*,*)uin(l,i,j,k,1)
                    stop
                 end if
              end do
           end if

           ! Compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/uin(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/uin(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/uin(l,i,j,k,1)
           end do

           ! Compute cell centered magnetic field
           DO l = 1, ngrid
              q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
              q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
              q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
           END DO

#if USE_FLD==1 || USE_M_1==1
           ! Compute radiative energy and flux
           do l = 1, ngrid
              erad_loc(l)=0.0d0
              do n = 1, nrad
                 q(l,i,j,k,8+n) = uin(l,i,j,k,8+n)
                 if(is_radiative_energy(n+1)) then
                    q(l,i,j,k,8+n) = max(q(l,i,j,k,8+n),small_er)
                    erad_loc(l) = erad_loc(l) + q(l,i,j,k,8+n)
                 endif
              end do
           end do
#endif
           
           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute thermal pressure through EOS
           do l = 1, ngrid
              etot = uin(l,i,j,k,5) - emag(l)
#if USE_FLD==1 || USE_M_1==1
              etot = etot - erad_loc(l)
#endif
              eint = etot/uin(l,i,j,k,1)-eken(l)
              eint = eint*q(l,i,j,k,1)   ! volumic
              if(energy_fix)eint=uin(l,i,j,k,nvar)
              call pressure_eos(q(l,i,j,k,1),eint,pp_eos)
              q(l,i,j,k,5)=MAX(pp_eos,smallp)
!!$              if(eos)then
!!$              eint = eint*q(l,i,j,k,1)   ! volumic
!!$              call pressure_eos(q(l,i,j,k,1),eint,pp_eos)
!!$                 q(l,i,j,k,5)=MAX(pp_eos,smallp)
!!$              else
!!$                 q(l,i,j,k,5)=MAX((gamma-one)*q(l,i,j,k,1)*eint,smallp)
!!$              endif
           end do

           ! Gravity predictor step
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) + gravin(l,i,j,k,idim)*dt*half
              end do
           end do

        end do
     end do
  end do

  ! Passive scalar
  do n = 9+nrad, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/uin(l,i,j,k,1)
              end do
           end do
        end do
     end do
  end do
 
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(bf,q,dq,dbf,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,xslope_type
  real(dp)::dpos_lim, dpos_rlim, dpos_ulim, dpos_plim
  real(dp)::beta_x, beta_y, beta_z, beta_abs, dtdx, dxdt, tgp1
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

#if NDIM==1
  if(slope_type==0)then
    dq=zero
  else if(slope_type==1.or.slope_type==2)then  ! minmod or average
    do n = 1, nvar
       do k = klo, khi
          do j = jlo, jhi
             do i = ilo, ihi
                do l = 1, ngrid
                   dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                   drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                   dcen = half*(dlft+drgt)/slope_type
                   dsgn = sign(one, dcen)
                   slop = min(abs(dlft),abs(drgt))
                   dlim = slop
                   if((dlft*drgt)<=zero)dlim=zero
                   dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                end do
             end do
          end do
       end do     
    end do
  else
     write(*,*)'Unknown slope type'
     stop
  end if
  
  if (positivity_type==0) then
     continue ! No positivity limiting
  else if (positivity_type==1) then
     dpos_ulim = min(dx / ((one + gamma)*dt), dx / (2.0d0 * dt))
     dpos_plim = one / (one + gamma)
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
              ! Extra limiting to ensure positivity preservation for
              ! density and pressure during MUSCL-Hancock predictor
              ! step; following K Waagan 2009, JCompPhys
              ! May require CFL number less than 1/ndim?
              ! Density - |drho| must be less than rho/2
              ! for positive pressure/density
              ! PS based on dimensionally split code and I missed the last
              ! (complicated) condition
              dlim = q(l,i,j,k,1)/2.d0
              dsgn = sign(one, dq(l,i,j,k,1,1))
              dq(l,i,j,k,1,1) = dsgn * min(abs(dq(l,i,j,k,1,1)), dlim)
              ! Velocity - |du| must be less than:
              !   0.5 dx/dt for positive density
              !   dx/((1+gamma) dt) for positive pressure
              dlim = dpos_ulim
              dsgn = sign(one, dq(l,i,j,k,2,1))
              dq(l,i,j,k,2,1) = dsgn * min(abs(dq(l,i,j,k,2,1)), dlim)
              ! Pressure - |dp| must be less than p / (1+gamma)
              dlim = dpos_plim * q(l,i,j,k,5)
              dsgn = sign(one, dq(l,i,j,k,5,1))
              dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dlim)
              end do
           end do
        end do
     end do
  else if (positivity_type==2) then
     ! Note this is as restrictive as the 3D version (intentionally)
     dtdx = dt/dx
     dxdt = dx/dt
     tgp1 = (gamma * 3.d0) + 1.d0
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 ! More restrictive limiting following the logic of
                 ! Waagan 2009. Full limiting for unsplit scheme (sort of)
                 ! on density, positivity on edge reconstruction only for
                 ! pressure. Considers the actual speeds rather than just
                 ! assuming the courant condition (less limiting)
                 beta_abs = abs(q(l,i,j,k,2) * dtdx)
                 dpos_rlim = q(l,i,j,k,1) / (beta_abs + 3.d0)
                 dpos_lim = 1.0d0 / (beta_abs + tgp1)
                 dpos_ulim = dpos_lim * dxdt
                 dpos_plim = 2.d0 * q(l,i,j,k,5) * dpos_lim
                 
                 ! Limit the density
                 dsgn = sign(one, dq(l,i,j,k,1,1))
                 dq(l,i,j,k,1,1) = dsgn*min(abs(dq(l,i,j,k,1,1)), dpos_rlim)
                 
                 ! Limit the velocity
                 dsgn = sign(one, dq(l,i,j,k,2,1))
                 dq(l,i,j,k,2,1) = dsgn*min(abs(dq(l,i,j,k,2,1)), dpos_ulim)
                 
                 ! Limit the pressure (guarantees edge reconstruction only)
                 dsgn = sign(one, dq(l,i,j,k,5,1))
                 dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dpos_plim)
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown positivity slope type'
     stop
  end if
#endif

#if NDIM==2   
  if(slope_type==0)then
    dq=zero           
  else if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)
                    
                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    
                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
  ! 1D transverse TVD slopes for face-centered magnetic fields
  ! Bx along direction Y
  if (slope_mag_type==0) then
    dbf=zero
  else if (slope_mag_type==1 .or. slope_mag_type==2) then
    do k = klo, khi
       do j = jlo, jhi
          do i = ilo, ihi+1 ! WARNING HERE
             do l = 1, ngrid
                dlft = slope_mag_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
                drgt = slope_mag_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
                dcen = half*(dlft+drgt)/slope_mag_type
                dsgn = sign(one, dcen)
                slop = min(abs(dlft),abs(drgt))
                dlim = slop
                if((dlft*drgt)<=zero)dlim=zero
                dbf(l,i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
             end do
          enddo
       end do
    end do
    ! By along direction X
    do k = klo, khi
       do j = jlo, jhi+1 ! WARNING HERE
          do i = ilo, ihi
             do l = 1, ngrid
                dlft = slope_mag_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
                drgt = slope_mag_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
                dcen = half*(dlft+drgt)/slope_mag_type
                dsgn = sign(one, dcen)
                slop = min(abs(dlft),abs(drgt))
                dlim = slop
                if((dlft*drgt)<=zero)dlim=zero
                dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
             end do
          enddo
       end do
    end do
  else
     write(*,*)'Unknown mag. slope type'
     stop
  endif
  
  if (positivity_type==0) then
     continue ! No positivity limiting
  else if (positivity_type==1) then
     dpos_ulim = min(dx / ((one + gamma)*dt), dx / (2.0d0 * dt))
     dpos_plim = one / (one + gamma)
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
              ! Extra limiting to ensure positivity preservation for
              ! density and pressure during MUSCL-Hancock predictor
              ! step; following K Waagan 2009, JCompPhys
              ! May require CFL number less than 1/ndim?
              ! Density - |drho| must be less than rho/2
              ! for positive pressure/density
              ! PS based on dimensionally split code and I missed the last
              ! (complicated) condition
              dlim = q(l,i,j,k,1)/2.d0
              dsgn = sign(one, dq(l,i,j,k,1,1))
              dq(l,i,j,k,1,1) = dsgn * min(abs(dq(l,i,j,k,1,1)), dlim)
              dsgn = sign(one, dq(l,i,j,k,1,2))
              dq(l,i,j,k,1,2) = dsgn * min(abs(dq(l,i,j,k,1,2)), dlim)
              ! Velocity - |du| must be less than:
              !   0.5 dx/dt for positive density
              !   dx/((1+gamma) dt) for positive pressure
              dlim = dpos_ulim
              do n=2,3
                 dsgn = sign(one, dq(l,i,j,k,n,1))
                 dq(l,i,j,k,n,1) = dsgn * min(abs(dq(l,i,j,k,n,1)), dlim)
                 dsgn = sign(one, dq(l,i,j,k,n,2))
                 dq(l,i,j,k,n,2) = dsgn * min(abs(dq(l,i,j,k,n,2)), dlim)
              end do
              ! Pressure - |dp| must be less than p / (1+gamma)
              dlim = dpos_plim * q(l,i,j,k,5)
              dsgn = sign(one, dq(l,i,j,k,5,1))
              dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dlim)
              dsgn = sign(one, dq(l,i,j,k,5,2))
              dq(l,i,j,k,5,2) = dsgn * min(abs(dq(l,i,j,k,5,2)), dlim)
              end do
           end do
        end do
     end do
  else if (positivity_type==2) then
     ! Note this is as restrictive as the 3D version (intentionally)
     dtdx = dt/dx
     dxdt = dx/dt
     tgp1 = (gamma * 3.d0) + 1.d0
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 ! More restrictive limiting following the logic of
                 ! Waagan 2009. Full limiting for unsplit scheme (sort of)
                 ! on density, positivity on edge reconstruction only for
                 ! pressure. Considers the actual speeds rather than just
                 ! assuming the courant condition (less limiting)
                 beta_x = q(l,i,j,k,2) * dtdx
                 beta_y = q(l,i,j,k,3) * dtdx
                 beta_abs = abs(beta_x) + abs(beta_y)
                 dpos_rlim = q(l,i,j,k,1) / (beta_abs + 3.d0)
                 dpos_lim = 1.0d0 / (beta_abs + tgp1)
                 dpos_ulim = dpos_lim * dxdt
                 dpos_plim = 2.d0 * q(l,i,j,k,5) * dpos_lim
                 
                 ! Limit the density
                 dsgn = sign(one, dq(l,i,j,k,1,1))
                 dq(l,i,j,k,1,1) = dsgn*min(abs(dq(l,i,j,k,1,1)), dpos_rlim)
                 dsgn = sign(one, dq(l,i,j,k,1,2))
                 dq(l,i,j,k,1,2) = dsgn*min(abs(dq(l,i,j,k,1,2)), dpos_rlim)
                 
                 ! Limit the velocity
                 do n=2,3
                    dsgn = sign(one, dq(l,i,j,k,n,1))
                    dq(l,i,j,k,n,1) = dsgn*min(abs(dq(l,i,j,k,n,1)), dpos_ulim)
                    dsgn = sign(one, dq(l,i,j,k,n,2))
                    dq(l,i,j,k,n,2) = dsgn*min(abs(dq(l,i,j,k,n,2)), dpos_ulim)
                 end do
                 
                 ! Limit the pressure (guarantees edge reconstruction only)
                 dsgn = sign(one, dq(l,i,j,k,5,1))
                 dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dpos_plim)
                 dsgn = sign(one, dq(l,i,j,k,5,2))
                 dq(l,i,j,k,5,2) = dsgn * min(abs(dq(l,i,j,k,5,2)), dpos_plim)
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown positivity slope type'
     stop
  end if
#endif

#if NDIM==3
  if(slope_type==0)then
    dq=zero
  else if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
  
  if (positivity_type==0) then
     continue ! No positivity limiting
  else if (positivity_type==1) then
     dpos_ulim = min(dx / ((one + gamma)*dt), dx / (2.0d0 * dt))
     dpos_plim = one / (one + gamma)
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
              ! Extra limiting to ensure positivity preservation for
              ! density and pressure during MUSCL-Hancock predictor
              ! step; following K Waagan 2009, JCompPhys
              ! May require CFL number less than 1/ndim?
              ! Density - |drho| must be less than rho/2
              ! for positive pressure/density
              ! PS based on dimensionally split code and I missed the last
              ! (complicated) condition
              dlim = q(l,i,j,k,1)/2.d0
              dsgn = sign(one, dq(l,i,j,k,1,1))
              dq(l,i,j,k,1,1) = dsgn * min(abs(dq(l,i,j,k,1,1)), dlim)
              dsgn = sign(one, dq(l,i,j,k,1,2))
              dq(l,i,j,k,1,2) = dsgn * min(abs(dq(l,i,j,k,1,2)), dlim)
              dsgn = sign(one, dq(l,i,j,k,1,3))
              dq(l,i,j,k,1,3) = dsgn * min(abs(dq(l,i,j,k,1,3)), dlim)
              ! Velocity - |du| must be less than:
              !   0.5 dx/dt for positive density
              !   dx/((1+gamma) dt) for positive pressure
              dlim = dpos_ulim
              do n=2,4
                 dsgn = sign(one, dq(l,i,j,k,n,1))
                 dq(l,i,j,k,n,1) = dsgn * min(abs(dq(l,i,j,k,n,1)), dlim)
                 dsgn = sign(one, dq(l,i,j,k,n,2))
                 dq(l,i,j,k,n,2) = dsgn * min(abs(dq(l,i,j,k,n,2)), dlim)
                 dsgn = sign(one, dq(l,i,j,k,n,3))
                 dq(l,i,j,k,n,3) = dsgn * min(abs(dq(l,i,j,k,n,3)), dlim)
              end do
              ! Pressure - |dp| must be less than p / (1+gamma)
              dlim = dpos_plim * q(l,i,j,k,5)
              dsgn = sign(one, dq(l,i,j,k,5,1))
              dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dlim)
              dsgn = sign(one, dq(l,i,j,k,5,2))
              dq(l,i,j,k,5,2) = dsgn * min(abs(dq(l,i,j,k,5,2)), dlim)
              dsgn = sign(one, dq(l,i,j,k,5,3))
              dq(l,i,j,k,5,3) = dsgn * min(abs(dq(l,i,j,k,5,3)), dlim)
              end do
           end do
        end do
     end do
  else if (positivity_type==2) then
     dtdx = dt/dx
     dxdt = dx/dt
     tgp1 = (gamma * 3.d0) + 1.d0
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 ! More restrictive limiting following the logic of
                 ! Waagan 2009. Full limiting for unsplit scheme (sort of)
                 ! on density, positivity on edge reconstruction only for
                 ! pressure. Considers the actual speeds rather than just
                 ! assuming the courant condition (less limiting)
                 beta_x = q(l,i,j,k,2) * dtdx
                 beta_y = q(l,i,j,k,3) * dtdx
                 beta_z = q(l,i,j,k,4) * dtdx
                 beta_abs = abs(beta_x) + abs(beta_y) + abs(beta_z)
                 dpos_rlim = q(l,i,j,k,1) / (beta_abs + 3.d0)
                 dpos_lim = 1.0d0 / (beta_abs + tgp1)
                 dpos_ulim = dpos_lim * dxdt
                 dpos_plim = 2.d0 * q(l,i,j,k,5) * dpos_lim
                 
                 ! Limit the density
                 dsgn = sign(one, dq(l,i,j,k,1,1))
                 dq(l,i,j,k,1,1) = dsgn*min(abs(dq(l,i,j,k,1,1)), dpos_rlim)
                 dsgn = sign(one, dq(l,i,j,k,1,2))
                 dq(l,i,j,k,1,2) = dsgn*min(abs(dq(l,i,j,k,1,2)), dpos_rlim)
                 dsgn = sign(one, dq(l,i,j,k,1,3))
                 dq(l,i,j,k,1,3) = dsgn*min(abs(dq(l,i,j,k,1,3)), dpos_rlim)
                 
                 ! Limit the velocity
                 do n=2,4
                    dsgn = sign(one, dq(l,i,j,k,n,1))
                    dq(l,i,j,k,n,1) = dsgn*min(abs(dq(l,i,j,k,n,1)), dpos_ulim)
                    dsgn = sign(one, dq(l,i,j,k,n,2))
                    dq(l,i,j,k,n,2) = dsgn*min(abs(dq(l,i,j,k,n,2)), dpos_ulim)
                    dsgn = sign(one, dq(l,i,j,k,n,3))
                    dq(l,i,j,k,n,3) = dsgn*min(abs(dq(l,i,j,k,n,3)), dpos_ulim)
                 end do
                 
                 ! Limit the pressure (guarantees edge reconstruction only)
                 dsgn = sign(one, dq(l,i,j,k,5,1))
                 dq(l,i,j,k,5,1) = dsgn * min(abs(dq(l,i,j,k,5,1)), dpos_plim)
                 dsgn = sign(one, dq(l,i,j,k,5,2))
                 dq(l,i,j,k,5,2) = dsgn * min(abs(dq(l,i,j,k,5,2)), dpos_plim)
                 dsgn = sign(one, dq(l,i,j,k,5,3))
                 dq(l,i,j,k,5,3) = dsgn * min(abs(dq(l,i,j,k,5,3)), dpos_plim)
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown positivity slope type'
     stop
  end if

  ! 2D transverse TVD slopes for face-centered magnetic fields
  xslope_type=MIN(slope_mag_type,2)

  if(xslope_type==0)then
    dbf=zero
  else if(xslope_type==1.or.xslope_type==2)then  ! minmod or average
     ! Bx along direction Y and Z
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1 ! WARNING HERE
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
                 drgt = xslope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j,k  ,1) - bf(l,i,j,k-1,1))
                 drgt = xslope_type*(bf(l,i,j,k+1,1) - bf(l,i,j,k  ,1))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,1,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do

     ! By along direction X and Z
     do k = klo, khi
        do j = jlo, jhi+1 ! WARNING HERE
           do i = ilo, ihi
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
                 drgt = xslope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j,k  ,2) - bf(l,i,j,k-1,2))
                 drgt = xslope_type*(bf(l,i,j,k+1,2) - bf(l,i,j,k  ,2))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,2,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do

     ! Bz along direction X and Y
     do k = klo, khi+1 ! WARNING HERE
        do j = jlo, jhi
           do i = ilo, ihi
              ! slopes in first coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i  ,j,k,3) - bf(l,i-1,j,k,3))
                 drgt = xslope_type*(bf(l,i+1,j,k,3) - bf(l,i  ,j,k,3))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one, dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,3,1) = dsgn*min(dlim,abs(dcen))
              end do
              ! slopes in second coordinate direction
              do l = 1, ngrid
                 dlft = xslope_type*(bf(l,i,j  ,k,3) - bf(l,i,j-1,k,3))
                 drgt = xslope_type*(bf(l,i,j+1,k,3) - bf(l,i,j  ,k,3))
                 dcen = half*(dlft+drgt)/xslope_type
                 dsgn = sign(one,dcen)
                 slop = min(abs(dlft),abs(drgt))
                 dlim = slop
                 if((dlft*drgt)<=zero)dlim=zero
                 dbf(l,i,j,k,3,2) = dsgn*min(dlim,abs(dcen))
              end do
           end do
        end do
     end do
  endif
#endif
  
end subroutine uslope

! modif nimhd
! fonctions de produits vectoriels et coef nimhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdx(vec,n2,n3,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
real(dp)::dx
integer::n2,n3,l,i,j,k

computdx=(vec(l,i+1,j,k,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end function computdx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdy(vec,n2,n3,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
real(dp)::dx
integer::n2,n3,l,i,j,k

computdy=(vec(l,i,j+1,k,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end function computdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdz(vec,n2,n3,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
real(dp)::dx
integer::n2,n3,l,i,j,k

computdz=(vec(l,i,j,k+1,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdxbis(vec,n2,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
real(dp)::dx
integer::n2,l,i,j,k

computdxbis=(vec(l,i+1,j,k,n2)-vec(l,i,j,k,n2))/dx

end function computdxbis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdybis(vec,n2,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
real(dp)::dx
integer::n2,l,i,j,k

computdybis=(vec(l,i,j+1,k,n2)-vec(l,i,j,k,n2))/dx

end function computdybis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision  function computdzbis(vec,n2,l,i,j,k,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
real(dp)::dx
integer::n2,l,i,j,k

computdzbis=(vec(l,i,j,k+1,n2)-vec(l,i,j,k,n2))/dx

end function computdzbis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function gammaadbis(rhon,BBcell,BBcellold,temper)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters
use radiation_parameters,only:mu_gas
use cooling_module,ONLY:mH
use units_commons

implicit none  
real(dp)::rhon,rhoH,n_H_max,BBcell,temper,BBcellold
real(dp)::eta_AD_chimie

! function which computes the coefficient gamma which
! appears in ambipolar diffusion dB/dt=1/(gamma*rhoi*rhon)curl*(j*B)*B)+...
! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (6)
! WARNING no mu_0 needed here

n_H_max = 2.5d+17

! C shock Duffin et Pudritz
! gammaadbis in CGS
!gammaadbis=gammaAD

!rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

if(rhoH < n_H_max)then
   gammaadbis=eta_AD_chimie(rhoH,BBcell,BBcellold,temper)
else
   gammaadbis=eta_AD_chimie(n_H_max,BBcell,BBcellold,temper)
endif

gammaadbis=gammaadbis*scale_t*scale_d ! in code units

! test
!gammaadbis=gammaAD

end function gammaadbis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sig_x2d(ll,ii,j,k,lb,ib,sigO,sigH,sigP,bsquare)
use amr_parameters,    only : dp
use hydro_commons ,    only : resistivite_chimie_x
use variables_X
use amr_commons, only : myid
implicit none

integer, intent(in)             :: j,k,ib
real(dp)                        :: B,nH,temper,sigav
real(dp)                        :: j_dp,k_dp,b_dp
real(dp), dimension(nvarchimie) :: x
real(dp), intent(in)            :: ll,ii,lb,bsquare
real(dp), intent(out)           :: sigO,sigH,sigP
integer                         :: i,kk

j_dp = real(j,dp)
kk=min(k,tchimie-1)
k_dp = real(kk,dp)
b_dp = real(ib,dp)


x(1:3)=(1.d0-(ll-j_dp))*(1.d0-(ii-k_dp))*(1.d0-(lb-b_dp))*(resistivite_chimie(1:3,j,kk,ib))+&
           &((ll-j_dp))*(1.d0-(ii-k_dp))*(1.d0-(lb-b_dp))*(resistivite_chimie(1:3,j+1,kk,ib))+&
           &(1.d0-(ll-j_dp))*((ii-k_dp))*(1.d0-(lb-b_dp))*(resistivite_chimie(1:3,j,kk+1,ib))+&
                &((ll-j_dp))*((ii-k_dp))*(1.d0-(lb-b_dp))*(resistivite_chimie(1:3,j+1,kk+1,ib))+&
              (1.d0-(ll-j_dp))*(1.d0-(ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j,kk,ib+1))+&
                  &((ll-j_dp))*(1.d0-(ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j+1,kk,ib+1))+&
                  &(1.d0-(ll-j_dp))*((ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j,kk+1,ib+1))+&
                       &((ll-j_dp))*((ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j+1,kk+1,ib+1))
              
sigP= 10.0d0**x(1)
sigO= 10.0d0**x(2)

! modification since x(3) can be negative we simply use the sign of the leftmost
! point. If there is a sign inversion, we set it to zero.
! If you are using Hall resisitvities, this could be improved by using a linear
! interpolation instead of log.
sigH=(10.0d0**x(3))*resistivite_chimie(0,j,kk,ib)
sigav = sum(resistivite_chimie(0,j:j+1,kk:kk+1,ib:ib+1)) / 8.0d0
if(sigav .ne. resistivite_chimie(0,j,kk,ib))then
   sigH = 0.0_dp
   !if(myid==1) write(*,*)'Sign inversion in Hall resistivity'
endif

return

end subroutine sig_x2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function eta_AD_chimie(rhon,BBcell,BBcellold,temper)
use hydro_commons
use units_commons
use cooling_module ,only : clight
use variables_x,ONLY:dtchimie,dnchimie,nminchimie,tminchimie,&
                    &dbchimie,bminchimie,pi
implicit none

real(dp)     :: sigO,sigH,sigP,densionbis,BBcgs, bbcell,BBcellold
real(dp)::inp,ll,rhon,ii,temper,lb,j_dp
integer :: j,i,ib

!inp=rhon
!ll=(1.d0+(log10(inp)-log10(300.d0))/(15.d0/50.d0))
!ll=(1.d0+(log10(inp)-log10(nminchimie))/dnchimie)
!j=dble(floor(ll))

if(use_res==1)then
   inp=rhon
   ll=(1.d0+(log10(inp)-log10(nminchimie))/dnchimie)
   j=floor(ll)
   j_dp=real(j,dp)
!   ll=(1.d0+(log10(inp)-log10(300.d0))/(17.d0/35.d0))
!   j=dble(floor(ll))
   eta_AD_chimie=(ll-j_dp)*log10(resistivite_chimie_res(6,j+1))+(1.d0-(ll-j_dp))*log10(resistivite_chimie_res(6,j))
   eta_AD_chimie=10**eta_AD_chimie
!   print*, rhon,temper,eta_AD_chimie
!   print*, eta_AD_chimie, inp,(1.43d-7*sqrt(inp))**2
!   stop
else if(use_x2d==1)then
   inp=rhon
   ll=(1.d0+(log10(inp)-log10(nminchimie))/dnchimie)
   j=floor(ll)
   inp=temper
   ii=(1.d0+(log10(inp)-log10(tminchimie))/dtchimie)
   ii=max(ii,1.0d0)
!   ii=(1.d0+(log10(inp)-log10(5.d0))/(3.d0/50.d0))
   i=floor(ii)
   BBcgs=sqrt(BBcellold*(4.d0*pi*scale_d*(scale_v)**2))
!!$   bbcgs=1.43d-7*sqrt(rhon/2.d0/H2_fraction)

!!$   print*, bbcgs, sqrt(BBcellold*(4.d0*pi*scale_d*(scale_v)**2)),rhon
   inp=BBcgs
   lb=(1.d0+(log10(inp)-log10(bminchimie))/dbchimie)
   ib=floor(lb)

   call sig_x2d(ll,ii,j,i,lb,ib,sigO,sigH,sigP,BBcgs) 
!   inp=rhon/xmolaire/H2_fraction     ! inp is neutrals.cc, to fit densionbis
   inp=rhon/2.d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
   eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1.d0/sigP)   ! resistivity in s
!   print*,   eta_AD_chimie,inp*xmolaire*H2_fraction,BBcgs,inp,densionbis(inp),scale_d
!   print*, eta_AD_chimie,inp*xmolaire*H2_fraction,BBcgs

   BBcgs=sqrt(BBcell*(4.d0*pi*scale_d*(scale_v)**2))
!!$   (eta_AD_chimie = max(eta_AD_chimie * (1.0d0-tanh(ii/(dble(tchimie))), 1.d-36)

   eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*clight**2)  ! need B in G, output is gammaad in cgs

!print*,eta_ad_chimie
!   print*, rhon,temper,eta_AD_chimie
!   print*,   eta_AD_chimie,inp*xmolaire*H2_fraction,BBcgs,inp,densionbis(inp),scale_d
!   stop
!   print *,  eta_AD_chimie,inp,ll,j,sigO,sigH,sigP,BBcell,scale_d,densionbis(inp),clight
!   print *, 'biiiiii'
! print*,eta_AD_chimie, ll,j,ii,i,lb,ib,rhon,temper,bbcgs

endif
! stop

!!print*, inp, ll, j,resistivite_chimie(1,1),resistivite_chimie(1,2),resistivite_chimie(6,1),resistivite_chimie(1,35)
!eta_AD_chimie=(ll-j)*log10(resistivite_chimie(6,j+1))+(1.d0-(ll-j))*log10(resistivite_chimie(6,j))
!!print*, eta_AD_chimie
!eta_AD_chimie=10**eta_AD_chimie

! Ad-hoc modification to ensure that the ambipolar resistivity falls to zero when the density exceeds 5.0e1eta_AD_chimie * (1.0d0-tanh(rhon/5.0d13))

end function eta_AD_chimie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function eta_ohm_chimie(rhon,BBcell,temper)
use hydro_commons
use units_commons
use cooling_module, only : clight
use variables_x,ONLY:dtchimie,dnchimie,nminchimie,tminchimie,&
                    &dbchimie,bminchimie,pi
implicit none

real(dp) :: inp,ll,ii,lb,rhon,BBcell
real(dp) :: temper,sigO,sigH,sigP,BBcgs
real(dp) :: j_dp
integer  :: j,i,ib

if(use_res==1)then
   inp=rhon
   ll=(1.d0+(log10(inp)-log10(nminchimie))/dnchimie)
   j=floor(ll)
   j_dp=real(j,dp)
   eta_ohm_chimie=(ll-j_dp)*log10(resistivite_chimie_res(7,j+1))+(1.d0-(ll-j_dp))*log10(resistivite_chimie_res(7,j))
   eta_ohm_chimie=10.0d0**eta_ohm_chimie
else if(use_x2d==1)then
   inp=rhon
   ll=(1.d0+(log10(inp)-log10(nminchimie))/dnchimie)
   j=floor(ll)
   inp=temper
   ii=(1.d0+(log10(inp)-log10(tminchimie))/dtchimie)
   ii=max(ii,1.0d0)
   i=floor(ii)
   BBcgs=sqrt(BBcell*(4.d0*pi*scale_d*(scale_v)**2))
   inp=BBcgs
   lb=(1.d0+(log10(inp)-log10(bminchimie))/dbchimie)
   ib=floor(lb)
   call sig_x2d(ll,ii,j,i,lb,ib,sigO,sigH,sigP,BBcgs)
   eta_ohm_chimie = (1.d0 / sigP) * clight * clight / (4.0_dp*pi)
endif

! Ad-hoc modification to ensure that the ohmic resistivity falls to zero when the density exceeds 1.0e15
! when alkali metals are ionized.
!eta_ohm_chimie = eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15))
eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15)), 1.d-36)

end function eta_ohm_chimie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function densionbis(rhon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters
use units_commons

implicit none 
real(dp)::rhon
real(dp)::xn, rhoncgs

! density of neutral in g/cm3  
rhoncgs=rhon*scale_d

! function which computes the density in g/cm3 of ions 
! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (14)

! density of neutral in number per cm3
!xn=rhoncgs/xmneutre

! density of ions in g/cm3 
!densionbis=densionbis*xmion


! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
! WARNING 3.d-16 si in cgs
densionbis=coefionis*sqrt(rhoncgs)

!!!!!!!!!!!!!!!!!!!!!!!!!
! densionbis in USER UNITS !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

! transformation coefionis in user units
!densionbis=densionbis/sqrt(scale_d)

! back in user units
densionbis=densionbis/scale_d

! test C shock Duffin et Pudritz
if(ntestDADM.eq.1) then
   densionbis=rhoi0
endif

end function densionbis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdxvx(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz

computdxvx=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx 


end function computdxvx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdyvy(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz

computdyvy=(vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy 

end function computdyvy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdzvz(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz

computdzvz=(vec(l,i,j,k+1,nzz)-vec(l,i,j,k,nzz))/dz

end function computdzvz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdiv(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz

if(ndim.eq.1) then
computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx 
endif
if(ndim.eq.2) then
computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx + (vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy 
endif
if(ndim.eq.3) then
computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx + (vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy + (vec(l,i,j,k+1,nzz)-vec(l,i,j,k,nzz))/dz
endif


end function computdiv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdivbisx(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz
real(dp)::computdxvx,computdyvy,computdzvz

if(ndim.eq.1) then
computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)
endif
if(ndim.eq.2) then
computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)
endif
if(ndim.eq.3) then
computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)+computdzvz(vec,l,i,j,k-1,dx,dy,dz)
endif


end function computdivbisx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdivbisy(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz
real(dp)::computdxvx,computdyvy,computdzvz

if(ndim.eq.1) then
computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)
endif
if(ndim.eq.2) then
computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j,k,dx,dy,dz)
endif
if(ndim.eq.3) then
computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j,k,dx,dy,dz)+computdzvz(vec,l,i,j,k-1,dx,dy,dz)
endif

end function computdivbisy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function computdivbisz(vec,l,i,j,k,dx,dy,dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
integer l,i,j,k
real(dp)::dx,dy,dz
real(dp)::computdxvx,computdyvy,computdzvz

if(ndim.eq.1) then
computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)
endif
if(ndim.eq.2) then
computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)
endif
if(ndim.eq.3) then
computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)+computdzvz(vec,l,i,j,k,dx,dy,dz)
endif

end function computdivbisz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine crossprodbis(vec1,vec2,v1crossv2,l,i,j,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec1,vec2,v1crossv2
integer ::l,i,j,k 

real(dp)::v1x,v1y,v1z,v2x,v2y,v2z,crossprodx,crossprody,crossprodz

integer::n

do n=1,3

   v1x=vec1(l,i,j,k,1,n)
   v1y=vec1(l,i,j,k,2,n)
   v1z=vec1(l,i,j,k,3,n)
   
   v2x=vec2(l,i,j,k,1,n)
   v2y=vec2(l,i,j,k,2,n)
   v2z=vec2(l,i,j,k,3,n)
   
   v1crossv2(l,i,j,k,1,n)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
   v1crossv2(l,i,j,k,2,n)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)
   v1crossv2(l,i,j,k,3,n)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

end do

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine crossprod(vec1,vec2,v1crossv2,l,i,j,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec1,vec2,v1crossv2
integer ::l,i,j,k 

real(dp)::v1x,v1y,v1z,v2x,v2y,v2z,crossprodx,crossprody,crossprodz

v1x=vec1(l,i,j,k,1)
v1y=vec1(l,i,j,k,2)
v1z=vec1(l,i,j,k,3)

v2x=vec2(l,i,j,k,1)
v2y=vec2(l,i,j,k,2)
v2z=vec2(l,i,j,k,3)

v1crossv2(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
v1crossv2(l,i,j,k,2)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)
v1crossv2(l,i,j,k,3)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

end subroutine crossprod


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function  crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function which gives the x component of a cross product of two
! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z

use hydro_parameters

implicit none

real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

crossprodx=v1y*v2z-v1z*v2y

end function crossprodx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function crossprody(v1x,v1y,v1z,v2x,v2y,v2z)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function which gives the y component of a cross product of two
! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z

use hydro_parameters
implicit none

real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

crossprody=v1z*v2x-v1x*v2z

end function crossprody

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! function which gives the z component of a cross product of two
! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z

use hydro_parameters
implicit none

real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

crossprodz=v1x*v2y-v2x*v1y

end function crossprodz



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function reshall(rhon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use hydro_parameters
  use units_commons

  implicit none 
  real(dp) ::rhon
  real(dp)::rhocgs,ni
  real(dp)::densionbis

! function which computes the coefficient Rh which
! appears in ohmic dissipation dB/dt=-curl(-1/Rh*J*B)+...
! Rh=1/(Zen_i) n_i in cm-3

! ions density in g/cm3
 ni=densionbis(rhon)
! convert to CGS
ni=ni*scale_d

! convert to cm-3
ni=ni/xmion

! electric elementary charge in cgs : e=4.803d-10
! ions with one elementary charge Rh=1/(Z*e*ni)
resHall=1.d0/(4.803d-10*ni)

! convert to user units : Rh in cm/sqrt(g/cm3)
resHall=resHall*sqrt(scale_d)/scale_l

  if(ntestDADM.eq.1) then

     resHall=rHall
     
  endif

 
end function reshall


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function muvisco(rhon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 

real(dp) ::rhon

muvisco=visco

if(ntestDADM.eq.1) then
   muvisco=visco
endif

end function muvisco

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function etaohmdiss(rhon,BBcell,temper)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_commons
use radiation_parameters,only:mu_gas
use cooling_module,ONLY:mH
use units_commons

implicit none 
real(dp) ::rhon,xpressure,rhoH,rhotemp,BBcell
real(dp)::gammaadbis,densionbis
real(dp)::xionisation,temper,scale_p,xpcgs,rhocgs,xnbcgs,n_H_max
real(dp)::eta_ohm_chimie

if(ntestDADM.eq.0) then

   ! function which computes the coefficient eta which
   ! appears in ohmic dissipation dB/dt=-curl(eta*curl(B))+...
   ! see Machida, Inutsuka, Matsumoto, ApJ, 670,1198-1213, 2007

   n_H_max = 2.5d+17

   ! convert to CGS

   ! scale_p = scale_d*(scale_v**2.)
   ! xpcgs=xpressure*scale_p
   ! rhocgs=rhon*scale_d
   ! ! nb per cm3
   ! xnbcgs=rhocgs/xmneutre
   ! ! temperature in cgs
   ! temper=xpcgs*xmolaire/(rhocgs*rperfectgaz)
   ! !write(*,*)'temper',temper
   ! 
   ! ! degree of ionisation
   ! ! Machida et al 2007 
   ! xionisation=5.7d-4/(xnbcgs)
   ! ! Shu 1987 27.7 p 363 and p 361 m_n=2.33 m_i=29
   ! !xionisation=2.33d0/29.d0*densionbis(rhon)/rhon
   ! 
   ! ! Machida et al 2007 : etaMD=740
   ! !etaohmdiss=etaMD*sqrt(temper/10.d0)*(1.d0-tanh(xnbcgs/1.d15))/xionisation
   ! ! Dapp & Basu 2010
   ! ! etaohmdiss=etaMD*1.3d18*(xnbcgs/1.d12)*sqrt(temper/10.d0)*(1.d0-tanh(xnbcgs/1.d15))
   ! ! back to user units
   ! !print*, etaohmdiss
   ! !stop

   !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

   rhotemp = MAX(rhoH,rho_threshold)

   if(rhotemp < n_H_max)then
      etaohmdiss=eta_ohm_chimie(rhotemp,BBcell,temper)
   else
      etaohmdiss=eta_ohm_chimie(n_H_max,BBcell,temper)
   endif
   
   etaohmdiss=etaohmdiss*scale_t/(scale_l)**2

elseif(ntestDADM.eq.1) then

   ! test Alfven Lessaffre
   !etaohmdiss=2.d-2

   ! test heat diffusion
   !etaohmdiss=1.d0

   ! test oblique shock
   !etaohmdiss=0.15d0

   etaohmdiss=etaMD

endif

 
end function etaohmdiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function etaohmdissbricolo(rhon,BBcell,temper,dtlim,dx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_commons
use radiation_parameters,only:mu_gas
use cooling_module,ONLY:mH
use units_commons

implicit none 
real(dp) ::rhon,xpressure,rhoH,rhotemp,BBcell
real(dp)::gammaadbis,densionbis
real(dp)::xionisation,temper,scale_p,xpcgs,rhocgs,xnbcgs,n_H_max
real(dp)::eta_ohm_chimie,dx,dtlim,xx,dtt

if(ntestDADM.eq.0) then

   ! function which computes the coefficient eta which
   ! appears in ohmic dissipation dB/dt=-curl(eta*curl(B))+...
   ! see Machida, Inutsuka, Matsumoto, ApJ, 670,1198-1213, 2007

   n_H_max = 2.5d+17

   ! convert to CGS

   ! scale_p = scale_d*(scale_v**2.)
   ! xpcgs=xpressure*scale_p
   ! rhocgs=rhon*scale_d
   ! ! nb per cm3
   ! xnbcgs=rhocgs/xmneutre
   ! ! temperature in cgs
   ! temper=xpcgs*xmolaire/(rhocgs*rperfectgaz)
   ! !write(*,*)'temper',temper
   ! 
   ! ! degree of ionisation
   ! ! Machida et al 2007 
   ! xionisation=5.7d-4/(xnbcgs)
   ! ! Shu 1987 27.7 p 363 and p 361 m_n=2.33 m_i=29
   ! !xionisation=2.33d0/29.d0*densionbis(rhon)/rhon
   ! 
   ! ! Machida et al 2007 : etaMD=740
   ! !etaohmdiss=etaMD*sqrt(temper/10.d0)*(1.d0-tanh(xnbcgs/1.d15))/xionisation
   ! ! Dapp & Basu 2010
   ! ! etaohmdiss=etaMD*1.3d18*(xnbcgs/1.d12)*sqrt(temper/10.d0)*(1.d0-tanh(xnbcgs/1.d15))
   ! ! back to user units
   ! !print*, etaohmdiss
   ! !stop

   !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

   rhotemp = MAX(rhoH,rho_threshold)

   if(rhotemp < n_H_max)then
      etaohmdissbricolo=eta_ohm_chimie(rhotemp,BBcell,temper)
   else
      etaohmdissbricolo=eta_ohm_chimie(n_H_max,BBcell,temper)
   endif

   etaohmdissbricolo=etaohmdissbricolo*scale_t/(scale_l)**2

   ! robbery to avoid too small time step
   if(nminitimestep.eq.1 .and. nmagdiffu2.eq.0) then

      if(dtlim.ne.0.d0) then
         xx=etaohmdissbricolo
         if(xx.ne.0.d0) then
          dtt=coefohm*dx*dx/xx   !dtohm pour la cellule
      !    if (myid ==1) print*, dtt,bsquare,betaadbricolo,betaadbricolotemp
         else
            dtt=1.d39
         endif
         if (dtt.le.dtlim) then
            etaohmdissbricolo=coefohm*dx*dx/(dtlim)
         endif
      endif

   endif
   
!   etaohmdissbricolo=etaohmdissbricolo*scale_t/(scale_l)**2

elseif(ntestDADM.eq.1) then

   ! test Alfven Lessaffre
   !etaohmdiss=2.d-2

   ! test heat diffusion
   !etaohmdiss=1.d0

   ! test oblique shock
   !etaohmdiss=0.15d0

   etaohmdissbricolo=etaMD

endif

 
end function etaohmdissbricolo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function betaad(rhon,bsquare,temper)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters

implicit none 
real(dp) ::rhon,rhotemp,bsquare,temper
real(dp)::gammaadbis,densionbis

real(dp)::xx

if(ntestDADM.eq.0) then

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   ! Warning gammaadbis and densionbis already in user units
   ! but NOT rhon/xmneutre

   !betaad=1.4d0/(gammaadbis(rhon)*densionbis(rhon)*rhon/xmneutre )
   ! no xmneutre for Duffin and Pudritz

   rhotemp = MAX(rhon,rho_threshold)

   !xx=gammaadbis(rhotemp,bsquare,temper)*densionbis(rhon)*rhon
   xx=gammaadbis(rhotemp,bsquare,bsquare,temper)*densionbis(rhotemp)*rhotemp

   !xx=gammaadbis(rhotemp)*densionbis(rhotemp)*rhotemp
   !write(*,*)'gammaadbis',gammaadbis(rhon),densionbis(rhon),rhon
   if(xx.ne.0.d0) then
      betaad=1.d0/xx 
   else
      betaad=1.d39
      if(rhotemp < 1.0d+14)then
         write(*,*)'WARNING gammaadbis(rhotemp,bsquare,temper)*densionbis(rhon)*rhon equal 0',gammaadbis(rhotemp,bsquare,bsquare,temper),densionbis(rhotemp),rhotemp
      endif
   endif

   ! Barenblatt
   !betaad=1.d0

elseif(ntestDADM.eq.1) then

   ! test Barenblatt
   !betaadbricolo=1.d0
   ! test C shock
      betaad=1.d0/(gammaAD*rhoi0*rhon)
   !betaadbricolo=0.d0
   
endif

!rhon, gammaadbis(rhon) and densionbis(rhon) already in user units
!!betaad=betaad/scale_d

end function betaad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function betaadbricolo(rhocelln,rhon,dtlim,bsquare,bsquareold,dx,ntest,temper)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use hydro_parameters
use amr_commons
use cooling_module
use variables_X,only:pi
use units_commons

implicit none 

integer :: ntest
real(dp) ::rhocelln,rhon,betaadbricolotemp,dtlim,bsquare,bsquareold,dx,temper
real(dp)::gammaadbis,densionbis,rhotemp,rhotemp_cell

real(dp)::xx,dtt,bbcgs

if(ntestDADM.eq.0) then

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   ! Warning gammaadbis and densionbis already in user units
   ! but NOT rhon/xmneutre

   !betaad=1.4d0/(gammaadbis(rhon)*densionbis(rhon)*rhon/xmneutre )
   ! no xmneutre for Duffin and Pudritz

   rhotemp = MAX(rhon,rho_threshold)
   rhotemp_cell = MAX(rhocelln,rho_threshold)
   !if (myid ==1) then
   !   print*,densionbis(rhocelln),rhocelln
   !   print*,densionbis(rhotemp_cell),rhotemp_cell
   !   print*, 'rho thres', rho_threshold
   !end if


!    xx=gammaadbis(rhotemp_cell,bsquare,temper)*densionbis(rhocelln)*rhocelln  ! dans la cellule

   xx=gammaadbis(rhotemp_cell,bsquare,bsquareold,temper)*densionbis(rhotemp_cell)*rhotemp_cell  ! dans la cellule

   if(xx.ne.0.d0) then
      betaadbricolo=1.d0/xx 
   else
      betaadbricolo=1.d39
      if(rhotemp < 1.0d+14)then
         write(*,*)'WARNING gammaadbis(rhocelln,bsquare,bsquareold,temper)*densionbis(rhocelln)*rhocelln in the cell equals 0',gammaadbis(rhotemp_cell,bsquare,bsquareold,temper),densionbis(rhocelln),rhocelln,bsquare,bsquareold,temper
      endif
   endif

   !xx=gammaadbis(rhotemp,bsquare,bsquareold,temper)*densionbis(rhon)*rhon   ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille

   xx=gammaadbis(rhotemp,bsquare,bsquareold,temper)*densionbis(rhotemp)*rhotemp  

 ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille

   if(xx.ne.0.d0) then
      betaadbricolotemp=1.d0/xx 
   else
      betaadbricolotemp=1.d39
      if(rhotemp < 1.0d+14)then
         write(*,*)'WARNING gammaadbis(rhon,bsquare,bsquareold,temper)*densionbis(rhon)*rhon at the interface equals 0',gammaadbis(rhotemp,bsquare,bsquareold,temper),densionbis(rhon),rhon
      endif
   endif


   ! robbery to avoid too small time step
   if(nminitimestep.eq.1 .and. nambipolar2.eq.0) then

      if(dtlim.ne.0.d0) then
         xx=bsquare*betaadbricolo
         if(xx.ne.0.d0) then
          dtt=coefad*dx*dx/xx   !dtAD pour la cellule
      !    if (myid ==1) print*, dtt,bsquare,betaadbricolo,betaadbricolotemp
         else
            dtt=1.d39
         endif
         if (dtt.le.dtlim) then   ! on compare bien dtAD calcule pour la cellule (rhocelln) avec le temps de la simu
            betaadbricolo=coefad*dx*dx/(dtlim*bsquare)
      !      write(*,*) 'la où ça seuille rho et B valent : ', rhocelln, bsquare
            !ici dtlim est le dt le plus petit : normal, ou seuillé si besoin est.
         else
            betaadbricolo=betaadbricolotemp  ! le betaad normal calcule avec rho a l'interface
         endif
      endif

   endif

elseif(ntestDADM.eq.1) then
   ! test Barenblatt
   !betaadbricolo=1.d0
   ! test C shock
      betaadbricolo=1.d0/(gammaAD*rhoi0*rhon)
   !betaadbricolo=0.d0
endif

end function betaadbricolo
! fin modif nimhd

