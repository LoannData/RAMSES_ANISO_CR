!################################################################
! TAKEN FROM PATRICK HENNEBELLE'S sink_particle.f90 AND PUT HERE
! BY SAM GEEN, JANUARY 2015 
!################################################################
!################################################################
!################################################################
!################################################################
subroutine star_jet_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine checks jet events in cells where a
  ! star particle has been spawned.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,isink,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar,temp_blast
  real(dp)::T2_AGN,T2_min,T2_max,delta_mass_max
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering star_jet_feedback'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
!  nx_loc=(icoarse_max-icoarse_min+1)
!  scale=boxlen/dble(nx_loc)
!  dx_min=(0.5D0**nlevelmax)*scale
!  vol_min=dx_min**ndim


  ! Compute the direction of the magnetic field in the sink
!  call compute_magnetic_field


  ! Compute the direction of the angular momentum around the sink
!  call compute_angular_momentum


  ! Compute the grid discretization effects
  call average_jet

!#ifndef WITHOUTMPI
!  call MPI_ALLREDUCE(ok_blast_agn,ok_blast_agn_all,nsink,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
!  ok_blast_agn=ok_blast_agn_all
!#endif

  ! Modify hydro quantities to account for the AGN blast
  call put_jet


  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar + 3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine star_jet_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_jet
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  !created by PH 02/2014
  ! This routine average the hydro quantities inside the jet location
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,j,isink,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,drr,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m_ms
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector),save::ok

  real(dp):: bx,by,bz,dxtot,cos_b_r,btot,lx,ly,lz,ltot
  logical::first_write

  integer,dimension(1:nsinkmax)::ntot_vol,ntot_vol_all


  

  first_write=.true.

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m_ms = scale_d*scale_l**3d0 / 2.d33


  if(nsink==0)return
  if(verbose)write(*,*)'Entering average_jet'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax


  ! Initialize the averaged variables
  vol_gas_jet=0.0
  ntot_vol=0
  ntot_vol_all=0  




   dx=0.5D0**nlevelmax 
   dx_min=dx*scale
   ! Maximum radius of the ejecta equal to twice the accretion radius
   rmax=4.0*dble(ir_cloud)*dx_min
   rmax2=rmax*rmax
   vol_loc=dx_min**ndim

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 



     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do isink=1,nsink

                  if( msink(isink)*scale_m_ms .lt. mass_jet_sink  .or. delta_mass(isink)/msink(isink) .lt. 1.d-3 ) cycle
                   ! Check if the cell lies within the sink radius
                   dxx=x-xsink(isink,1)
                   if(dxx> 0.5*scale)then
                       dxx=dxx-scale
                   endif
                   if(dxx<-0.5*scale)then
                       dxx=dxx+scale
                   endif
                   dyy=y-xsink(isink,2)
                   if(dyy> 0.5*scale)then
                       dyy=dyy-scale
                   endif
                   if(dyy<-0.5*scale)then
                       dyy=dyy+scale
                   endif
                   dzz=z-xsink(isink,3)
                   if(dzz> 0.5*scale)then
                       dzz=dzz-scale
                   endif
                   if(dzz<-0.5*scale)then
                       dzz=dzz+scale
                   endif
                   drr=dxx*dxx+dyy*dyy+dzz*dzz


                   !magnetic field/angular momentum around the sink
!                   bx = bsink(isink,1)
!                   by = bsink(isink,2)
!                   bz = bsink(isink,3)


                   lx = lsink(isink,1)
                   ly = lsink(isink,2)
                   lz = lsink(isink,3)

!                   bx = ang_mom_sink(isink,1)
!                   by = ang_mom_sink(isink,2)
!                   bz = ang_mom_sink(isink,3)

                   ltot = sqrt(lx**2+ly**2+lz**2)


                   if(ltot .ne. 0) then
                      lx=lx/ltot
                      ly=ly/ltot
                      lz=lz/ltot
                   else
                      if(myid .eq. 1. .and. first_write) then
                        write(*,*) 'problem angular momentum is 0, choose arbitrary direction'
                        write(*,*) isink,msink(isink),scale_m_ms,mass_jet_sink
                        first_write=.false.
                      endif
                      lx=1.
                      ly=0.
                      lz=0.
                   endif


!                      lx=1.
!                      ly=0.
!                      lz=0.

                   dxtot = sqrt(dxx**2+dyy**2+dzz**2)
  
                   cos_b_r = (lx*dxx + ly*dyy + lz*dzz) / dxtot


!!!! ******************

!                   write(*,*) drr,rmax2,cos_b_r,dxtot,btot


                    !we select the cells which are  outside the sink particle and aligned with the magnetic field in the sink
                    if(drr.lt.rmax2 .and. drr .gt. rmax2/4. .and. abs(cos_b_r) .gt. sqrt(3.)/2.)then
                       vol_gas_jet(isink)=vol_gas_jet(isink)+vol_loc
                       ntot_vol(isink) = ntot_vol(isink) + 1
!                       mass_gas_agn(isink)=mass_gas_agn(isink)+vol_loc*uold(ind_cell(i),1)
                    endif



!                    if(dr_cell.le.dx_loc/2.0)then
!                       ind_blast_agn(isink)=ind_cell(i)
!                       vol_blast_agn(isink)=vol_loc
!                       mass_blast_agn(isink)=vol_loc*uold(ind_cell(i),1)
!                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas_jet,vol_gas_jet_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!  call MPI_ALLREDUCE(mass_gas_agn,mass_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_jet=vol_gas_jet_all
!  mass_gas_agn=mass_gas_agn_all

  call MPI_ALLREDUCE(ntot_vol,ntot_vol_all,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
!  call MPI_ALLREDUCE(mass_gas_agn,mass_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_vol = ntot_vol_all

#endif


!  do isink=1,nsink
!           if (myid .eq. 1) write(*,*) 'isink, vol_gas_jet, ntot_vol ', isink,vol_gas_jet(isink),ntot_vol(isink)
!  end do


  if(verbose)write(*,*)'Exiting average_jet'

end subroutine average_jet
!################################################################
!################################################################
!################################################################
!################################################################
subroutine put_jet
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine adds the jet velocity around the sink particle 
  ! and in the direction of the magnetic field within the sink
  ! created by PH 02/2014
  !------------------------------------------------------------------------
  integer::ilevel,j,isink,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,drr,d,u,v,w,ek,u_r,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,T2_AGN,T2_max
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector),save::ok

  real(dp):: bx,by,bz,dxtot,btot,lx,ly,lz,ltot

  real(dp)::dm_ej,frac_acc_ej,e,emag
  real(dp)::d_jet,u_jet,v_jet,w_jet,Utot
  real(dp)::cos_b_r,scale_m_ms,coef
  real(dp)::bx1,bx2,by1,by2,bz1,bz2

  integer::ivar

  real(dp),dimension(1:nsinkmax,1:ndim)::mom_tot,mom_tot_all
  logical::first_loc

  logical::first_write

  integer,dimension(1:nsinkmax)::ntot_vol,ntot_vol_all

  first_write=.true.


  ntot_vol = 0
  ntot_vol_all = 0

!  write(*,*) 'enter put jet'

  mom_tot=0. ; mom_tot_all=0.

  first_loc = .true.

  frac_acc_ej = 1./3.
 

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Utot = Ucoef*1.6d6 ! 16 km/s from Wang et al. 2010
  Utot = Utot / scale_v !in code units

  scale_m_ms = scale_d*scale_l**3d0 / 2.d33


  if(nsink==0)return
  if(verbose)write(*,*)'Entering put_jet'

!write(*,*)'Entering put_jet'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
 
   dx=0.5D0**nlevelmax 
   dx_min=dx*scale
   ! Maximum radius of the ejecta equal to twice the accretion radius
   rmax=4.0*dble(ir_cloud)*dx_min
   rmax2=rmax*rmax
   vol_loc=dx_min**ndim

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 



     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do isink=1,nsink

                  if( msink(isink)*scale_m_ms .lt. mass_jet_sink  .or. delta_mass(isink)/msink(isink) .lt. 1.d-3 .or. vol_gas_jet(isink) .eq. 0.) cycle
                   ! Check if the cell lies within the sink radius
                   dxx=x-xsink(isink,1)
                   if(dxx> 0.5*scale)then
                       dxx=dxx-scale
                   endif
                   if(dxx<-0.5*scale)then
                       dxx=dxx+scale
                   endif
                   dyy=y-xsink(isink,2)
                   if(dyy> 0.5*scale)then
                       dyy=dyy-scale
                   endif
                   if(dyy<-0.5*scale)then
                       dyy=dyy+scale
                   endif
                   dzz=z-xsink(isink,3)
                   if(dzz> 0.5*scale)then
                       dzz=dzz-scale
                   endif
                   if(dzz<-0.5*scale)then
                       dzz=dzz+scale
                   endif
                   drr=dxx*dxx+dyy*dyy+dzz*dzz


                   !magnetic field/angular momentum around the sink
!                   bx = bsink(isink,1)
!                   by = bsink(isink,2)
!                   bz = bsink(isink,3)


                   lx = lsink(isink,1)
                   ly = lsink(isink,2)
                   lz = lsink(isink,3)

!                   bx = ang_mom_sink(isink,1)
!                   by = ang_mom_sink(isink,2)
!                   bz = ang_mom_sink(isink,3)

                   ltot = sqrt(lx**2+ly**2+lz**2)


                   if(ltot .ne. 0) then
                      lx=lx/ltot
                      ly=ly/ltot
                      lz=lz/ltot
                   else
                      if(myid .eq. 1. .and. first_write) then
                        write(*,*) 'problem angular momentum is 0, choose arbitrary direction'
                        write(*,*) isink,msink(isink),scale_m_ms,mass_jet_sink
                        first_write=.false.
                      endif
                      lx=1.
                      ly=0.
                      lz=0.
                   endif


!                      lx=1.
!                      ly=0.
!                      lz=0.

                   dxtot = sqrt(dxx**2+dyy**2+dzz**2)
  
                   cos_b_r = (lx*dxx + ly*dyy + lz*dzz) / dxtot

                   !put the jet velocity in the radial direction
                   if(rad_jet) then 
                     lx = dxx / dxtot
                     ly = dyy / dxtot
                     lz = dzz / dxtot
                   endif


!!!! ******************

!                   write(*,*) drr,rmax2,cos_b_r,dxtot,btot


                    !we select the cells which are  outside the sink particle and aligned with the magnetic field in the sink
                    if(drr.lt.rmax2 .and. drr .gt. rmax2/4. .and. abs(cos_b_r) .gt. sqrt(3.)/2.)then


                       ntot_vol(isink) = ntot_vol(isink) + 1


                     !we eject a fraction of the mass since last accretion
                     !equal to delta_mass 
                       dm_ej =  frac_acc_ej * max(delta_mass(isink),0.) / vol_gas_jet(isink)

                       d = uold(ind_cell(i),1)

                       !new density in the cell
                       d_jet = d + dm_ej

                       !add a velocity along the sink magnetic field and in the frame of the sink 
                       !the sqrt(msink) is taken from Wang et al. 2010, APJ 709
                       !however since we have big sink (due to coarse resolution) which are 
                       !subfragmented, we assume that the sink reprensent sqrt(M) of mass sqrt(M)
                       !thus the 1/4 instead of 1/2.
                       coef = Utot * (msink(isink)*scale_m_ms)**(0.25) 
                       if( .not. rad_jet) coef=coef * cos_b_r / abs(cos_b_r)

                       u = uold(ind_cell(i),2) / d                 
                       v = uold(ind_cell(i),3) / d                
                       w = uold(ind_cell(i),4) / d                

                       u_jet = d/d_jet * u + (1.-d/d_jet) * (coef * lx + vsink(isink,1) )
                       v_jet = d/d_jet * v + (1.-d/d_jet) * (coef * ly + vsink(isink,2) ) 
                       w_jet = d/d_jet * w + (1.-d/d_jet) * (coef * lz + vsink(isink,3) ) 

                       mom_tot(isink,1) = mom_tot(isink,1) + dm_ej * coef * lx
                       mom_tot(isink,2) = mom_tot(isink,2) + dm_ej * coef * ly
                       mom_tot(isink,3) = mom_tot(isink,3) + dm_ej * coef * lz

!                       if(first_loc) then 
!                         write(*,*) 'dir angular mom ',lx,ly,lz
!                         write(*,*) 'vol_gas_jet ',vol_gas_jet(isink)
!                         write(*,*) 'coef ',coef,'d ',d, 'd_jet ', d_jet, 'msink ',msink(isink)*scale_m_ms
!                         write(*,*) 'delta_mass', delta_mass(isink),'u ',u, 'u_jet ',u_jet,'v ',v, 'v_jet ',v_jet,'w ',w, 'w_jet ',w_jet
!                         first_loc = .false.
!                       endif



!                       u_jet = ( coef * bx + vsink(isink,1) )
!                       v_jet = ( coef * by + vsink(isink,2) ) 
!                       w_jet = ( coef * bz + vsink(isink,3) ) 

!          vsink_new(isink,1)=vsink_new(isink,1)+acc_mass*u
!          vsink_new(isink,2)=vsink_new(isink,2)+acc_mass*v
!          vsink_new(isink,3)=vsink_new(isink,3)+acc_mass*w

!                       write(*,*) 'myid ',myid, 'isink ',isink, 'idp, j', ind_part(j),j
!                       write(*,*) 'd, d_jet', d, d_jet
!                       write(*,*) 'u_jet, v_jet, w_jet', u_jet, v_jet, w_jet, cos_b_r




                       bx1 = uold(ind_cell(i),6)
                       by1 = uold(ind_cell(i),7)
                       bz1 = uold(ind_cell(i),8)

                       bx2 = uold(ind_cell(i),nvar+1)
                       by2 = uold(ind_cell(i),nvar+2)
                       bz2 = uold(ind_cell(i),nvar+3)

                       emag = 0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d_jet


                       e = uold(ind_cell(i),5) - 0.5d0*(u**2+v**2+w**2)*d - emag
                       e = e / d * d_jet
                       e = e + 0.5d0*(u_jet**2+v_jet**2+w_jet**2)*d_jet + emag

                       uold(ind_cell(i),1) = d_jet
                       uold(ind_cell(i),2) = d_jet*u_jet
                       uold(ind_cell(i),3) = d_jet*v_jet
                       uold(ind_cell(i),4) = d_jet*w_jet
                       uold(ind_cell(i),5) = e
!                       do ivar=imetal,nvar
!                          uold(ind_cell(i),ivar) = d_jet*z(ivar)
!                       end do

                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels


  ! Reset accreted mass
  do isink=1,nsink
!           if (myid .eq. 1) write(*,*) 'msink ', msink(isink), 'delta_mass ', delta_mass(isink)

           !withdraw from the sink mass what has been put in the cells
           msink(isink) = msink(isink) - delta_mass(isink) * frac_acc_ej
           !set delta_mass to zero
           delta_mass(isink)= 0. 
  end do


  call MPI_ALLREDUCE(mom_tot,mom_tot_all,nsink*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!  call MPI_ALLREDUCE(mass_gas_agn,mass_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  mom_tot = mom_tot_all


  call MPI_ALLREDUCE(ntot_vol,ntot_vol_all,nsink,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
!  call MPI_ALLREDUCE(mass_gas_agn,mass_gas_agn_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_vol = ntot_vol_all

!  do isink=1,nsink
!        if (myid .eq. 1) write(*,*) 'isink, ntot ',isink, ntot_vol(isink)
!        write(*,*) 'mom_tot ',mom_tot(isink,1),mom_tot(isink,2),mom_tot(isink,3)
!  end do


  if(verbose)write(*,*)'Exiting put_jet'

end subroutine put_jet
!###########################################################
!###########################################################
!###########################################################
!###########################################################
