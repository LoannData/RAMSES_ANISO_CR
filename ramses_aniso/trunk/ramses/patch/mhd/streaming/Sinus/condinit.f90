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
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xl,yl,xr,yr,rl,rr,xx,yy,zz,ttmin,ttmax,pi,xcenter,Al,Ar,R0,A0,theta

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  pi=acos(-1d0)
  ttmin=-pi/8d0
  ttmax= pi/8d0
  R0=0.3
  A0=1d0/sqrt(2d0)
  do i=1,nn
     q(i,1  )=1.0d0
     q(i,2:4)=0.0d0
     q(i,5  )=10.0d0

#if NDIM==1
     q(i,6  )=1.0d0
     q(i,7:8)=0.0d0
     q(i,nvar+1)=1.0d0
     q(i,nvar+2:nvar+3)=0.0d0
     q(i,nvar)=(gamma_rad(1)-1d0)*(1d0+0.5d0*sin(2d0*acos(-1d0)*x(i,1)*boxlen))
#endif

#if NDIM==2
!!$     q(i,6  )=1.0d0
!!$     q(i,7:8)=0.0d0
!!$     q(i,nvar+1)=1.0d0
!!$     q(i,nvar+2:nvar+3)=0.0d0
!!$     q(i,6:7)=sqrt(2d0)
!!$     q(i,8  )=0.0d0
!!$     q(i,nvar+1:nvar+2)=sqrt(2d0)
!!$     q(i,nvar+3)=0.0d0
!!$     q(i,nvar)=(gamma_rad(1)-1d0)*(1d0+0.5d0*sin(2d0*acos(-1d0)*x(i,2)*boxlen))

     xx=x(i,1)-0.5d0*boxlen
     yy=x(i,2)-0.5d0*boxlen
     rr=sqrt(xx**2+yy**2)
     if(xx .ge. 0.0d0 .and. yy .ge. 0.0d0)theta=atan(yy/xx)
     if(xx .lt. 0.0d0 .and. yy .ge. 0.0d0)theta=pi/2d0    +atan(abs(xx/yy))
     if(xx .lt. 0.0d0 .and. yy .lt. 0.0d0)theta=pi        +atan(abs(yy/xx))
     if(xx .ge. 0.0d0 .and. yy .lt. 0.0d0)theta=pi*3d0/2d0+atan(abs(xx/yy))

!!$     if(rr.gt.0.2d0*boxlen .and. rr.lt.0.3d0*boxlen .and. theta.gt.ttmin .and. theta.lt.ttmax .and. xx.lt.0d0)then
     if(rr.gt.0.15d0*boxlen .and. rr.lt.0.35d0*boxlen)then
        q(i,nvar)=(gamma_rad(1)-1d0)*(1d0+0.5d0*cos(2d0*theta))
!!$        q(i,nvar)=(gamma_rad(1)-1d0)*1d0
     else
        q(i,nvar)=1d-2
     endif
     xl=x(i,1)-0.5*dx-boxlen/2.0
     xr=x(i,1)+0.5*dx-boxlen/2.0
     yl=x(i,2)-0.5*dx-boxlen/2.0
     yr=x(i,2)+0.5*dx-boxlen/2.0

     Ar = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,6)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     q(i,nvar+1)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,7)=(Al-Ar)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     q(i,nvar+2)=(Al-Ar)/dx
     q(i,8)=0.0
     q(i,nvar+3)=0.0

#endif

#if NDIM==3
     xx=x(i,1)-0.5d0*boxlen
     yy=x(i,2)-0.5d0*boxlen
     zz=x(i,3)-0.5d0*boxlen
     rr=sqrt(xx**2+yy**2)
     if(xx .ge. 0.0d0 .and. yy .ge. 0.0d0)theta=atan(yy/xx)
     if(xx .lt. 0.0d0 .and. yy .ge. 0.0d0)theta=pi/2d0    +atan(abs(xx/yy))
     if(xx .lt. 0.0d0 .and. yy .lt. 0.0d0)theta=pi        +atan(abs(yy/xx))
     if(xx .ge. 0.0d0 .and. yy .lt. 0.0d0)theta=pi*3d0/2d0+atan(abs(xx/yy))

     if(rr.gt.0.15d0*boxlen .and. rr.lt.0.35d0*boxlen .and. abs(zz).lt.0.1)then
        q(i,nvar)=(gamma_rad(1)-1d0)*(1d0+0.5d0*cos(2d0*theta))
     else
        q(i,nvar)=1d-2
     endif
     xl=x(i,1)-0.5*dx-boxlen/2.0
     xr=x(i,1)+0.5*dx-boxlen/2.0
     yl=x(i,2)-0.5*dx-boxlen/2.0
     yr=x(i,2)+0.5*dx-boxlen/2.0

     Ar = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,6)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     q(i,nvar+1)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,7)=(Al-Ar)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     q(i,nvar+2)=(Al-Ar)/dx
     q(i,8)=0.0
     q(i,nvar+3)=0.0
#endif

  end do

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
