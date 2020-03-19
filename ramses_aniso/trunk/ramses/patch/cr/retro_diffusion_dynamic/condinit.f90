!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,8+ivar)=q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+ivar)
  enddo
#endif
#if NVAR>8+NENER
  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega
     
     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine PcrLoann(x, Pcr, r_center, Pcr_0, Rsnr, dx, t, ncell)!(x,v,dx,t,ncell)
     use amr_parameters
     use hydro_parameters  
     implicit none
     integer ::ncell                         ! Size of input arrays
     real(dp)::dx                            ! Cell size
     real(dp)::t                             ! Current time

     !real(dp),dimension(1:nvector,1:3)::v    ! Velocity field

     real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
     real(dp),dimension(1:nvector)::Pcr      ! CR Pressure field 
     real(dp),dimension(1:3)::r_center       ! Center position of the SNR
     real(dp)::Pcr_0                         ! CR Pressure amplitude of the step function 
     real(dp)::Rsnr                          ! Radius of the CR Cloud 
     real(dp)::xx,yy,zz

     integer::i

     do i=1,ncell
          xx=x(i,1)
#if NDIM > 1
               yy=x(i,2)
#endif
#if NDIM > 2
               zz=x(i,3)
#endif
          
! Pour le moment j'injecte un carré de RCs
#if NDIM == 1
               if (xx >= r_center(1) - Rsnr .and. xx <= r_center(1) + Rsnr) then
                    Pcr(i) = Pcr_0
               else
                    Pcr(i) = 0
               endif 
#endif 

! Injection carrée 
!#if NDIM == 2
!               if ((xx >= r_center(1) - Rsnr .and. xx <= r_center(1) + Rsnr) .and. (yy >= r_center(2) - Rsnr .and. yy <= r_center(2) + Rsnr)) then
!                    Pcr(i) = Pcr_0
!               else
!                    Pcr(i) = 0
!               endif 
!#endif 

! Injection disque 
#if NDIM == 2
               if (sqrt((r_center(1) - xx)**2 + (r_center(2) - yy)**2) <= Rsnr) then
                    Pcr(i) = Pcr_0
               else
                    Pcr(i) = 0
               endif 
#endif 

#if NDIM == 3
               !write(*,*) "xx = ",xx, " yy = ",yy, " zz = ",zz
               !write(*,*)"r_center(1) = ",r_center(1),", r_center(2) = ",r_center(2),", r_center(3) = ",r_center(3)
               !write(*,*)"sqrt = ",sqrt((r_center(1) - xx)**2 + (r_center(2) - yy)**2 + (r_center(3) - zz)**2),", Rsnr = ",Rsnr
               if (sqrt((r_center(1) - xx)**2 + (r_center(2) - yy)**2 + (r_center(3) - zz)**2) <= Rsnr) then
                    Pcr(i) = Pcr_0
               else
                    Pcr(i) = 0
               endif 

               if (Pcr(i) > 0.) write(*,*) "Pcr _in = ",Pcr(i)

               !write(*,*) "xx = ",xx, "yy = ",yy, "zz = ",zz, "rx_center = ",r_center(1), "ry_center = ",r_center(2), "rz_center = ",r_center(3), "Rsnr = ",Rsnr, "Pcr(i) = ",Pcr(i)
#endif 

     end do 
     !================================================================
     ! This routine computes the user defined velocity fields.
     ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
     ! v(i,1:3) is the imposed 3-velocity in user units.
     !================================================================
   end subroutine PcrLoann