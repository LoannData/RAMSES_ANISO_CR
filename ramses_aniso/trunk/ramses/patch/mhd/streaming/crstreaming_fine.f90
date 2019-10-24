!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_dflux_crstreaming(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets the cr streaming fluxes to 0 before calling
  ! the flux scheme.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,icr
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Set flux_cr = 0 for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do icr=1,ncr
        do i=1,active(ilevel)%ngrid
           dflux_cr(active(ilevel)%igrid(i)+iskip,icr) = 0.0d0
        end do
     end do
  end do

  ! Set fcr to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do icr=1,ncr
        do i=1,reception(icpu,ilevel)%ngrid
          dflux_cr(reception(icpu,ilevel)%igrid(i)+iskip,icr)= 0.0d0
        end do
     end do
  end do
  end do

111 format('   Entering set_dflux_crstreaming for level ',i2)

end subroutine set_dflux_crstreaming
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold_cr(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold =unew for the quantities that have been 
  ! updated because of CR streaming
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,icr

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,5) = unew(active(ilevel)%igrid(i)+iskip,5)
           do icr=1,ncr
              uold(active(ilevel)%igrid(i)+iskip,8+icr) = unew(active(ilevel)%igrid(i)+iskip,8+icr)
           end do
        end do
  end do

111 format('   Entering set_uold_dust for level ',i2)

end subroutine set_uold_cr
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crstreaming_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the thermal conduction scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,icr
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweepsuold
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call crstreamingfine1(ind_grid,ngrid,ilevel)
  end do

  do icr=1,ncr
      call make_virtual_reverse_dp(dflux_cr(1,icr),ilevel)
  end do

111 format('   Entering crstreaming_fine for level ',i2)
end subroutine crstreaming_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crstreamingfine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
!  use poisson_commons
!  use cooling_module
  implicit none
  integer::ilevel,ncache,icr
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! cr streaming solver that compute cr energy flux. This flux is zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array unew(:), both at the current level and at the
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uuloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uuuloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  !real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::xloc=0.0d0
!!$  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),save::residu_loc=0.0d0

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
!!$  real(dp),dimension(1:nvector),save:: residu

  integer::neul=5
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dflux_x,dflux_y,dflux_z
  real(dp)::dx,scale,erad,etot
  real(dp)::dflux,weight,oneontwotondim

!!$  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_omega,scale_tau,scale_kappa
!!$  real(dp)::vx,vy,vz,dens,bnorm,bx,by,bz,vnorm,va,Ma,Dpara,kperp,kpar

  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  flux = 0.0d0
  uloc = 0.0d0
  uuloc= 0.0d0
  uuuloc= 0.0d0
  facdx= 0.0d0
  ok   = .false.

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  oneontwotondim=1.0d0/dble(twotondim)


  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     end if
     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        !Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i), i3, j3, k3, ivar) = uold(ind_cell(i), ivar)
              facdx(ind_exist(i),i3,j3,k3)=1.0d0
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do
        do i=1,nbuffer
           if(interpol_type.gt.0)then
              facdx(ind_nexist(i),i3,j3,k3)=1.0d0
           else
              facdx(ind_nexist(i),i3,j3,k3)=1.5d0
           endif
        end do
        !Computes and stores (in uuuloc) the density
        do i=1,nexist
           erad=0.0
#if NCR>0
           do irad=1,ncr
              erad=erad+uold(ind_cell(i),8+irad)
           end do
#endif
           etot=uold(ind_cell(i),5)
           uuuloc(ind_exist(i),i3,j3,k3,1)=etot-erad
        end do

     end do
     end do
     end do
     ! End loop over cells
  end do
  end do
  end do
  ! End loop over neighboring grids

  !-----------------------------------------------
  ! Compute energy flux due to thermal conduction
  !-----------------------------------------------
  call crstreaming_split(uloc,flux,dx,dx,dx,dtnew(ilevel),ncache,facdx)

 !Reset fluxes at refined interfaces
  do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do icr=1,ncr
            do i=1,ncache
               if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                  flux(i,i3,j3,k3,icr,idim)=0.0d0
               end if
            end do
         end do
      end do
      end do
      end do
   end do
! -------- A traiter a la fin ---------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------
  !   !Conservative  update at level ilevel for the cr streaming
  !-------------------------------------------------------------
  do idim=1,ndim
  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           i3=1+i2
           j3=1+j2
           k3=1+k2
           ! Update cr energy
           do icr = 1, ncr
           do i=1,ncache
           if(son(ind_cell(i))==0)then
!!$              write(*,'(A,4e11.3)')'before uuloc,uold,fl,fr ',uuloc(i,i3,j3,k3,icr),uold(ind_cell(i),9)&
!!$                   &,flux(i,i3,j3,k3,icr,idim),flux(i,i3+i0,j3+j0,k3+k0,icr,idim)
              uuloc(i,i3,j3,k3,icr)= uuloc(i,i3,j3,k3,icr)&
                   &+(flux(i,i3,j3,k3,icr,idim)&
                   &-flux(i,i3+i0,j3+j0,k3+k0,icr,idim))
!!$              write(*,'(A,e11.3)')'after uuloc              ',uuloc(i,i3,j3,k3,icr)
           end if
          end do
       end do
    end do
 end do
end do
end do
!Update conservative variables new state vector
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     ind_son=1+i2+2*j2+4*k2
     iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              unew(ind_cell(i),5) =  uuuloc(i,i3,j3,k3,1)
              do icr=1,ncr
                 ! Update cr energy
                 unew(ind_cell(i),8+icr)=uuloc(i,i3,j3,k3,icr)+uold(ind_cell(i),8+icr) &
                      &+dflux_cr(ind_cell(i),icr)
                 unew(ind_cell(i),5) = unew(ind_cell(i),5)+unew(ind_cell(i),8+icr)
              enddo
           end if
        end do
     end do
  end do
end do


  if(ilevel>levelmin)then
  !-----------------------------------------------------
  ! update at level ilevel-1
  !-----------------------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     !----------------------
     ! Left flux at boundary
     !----------------------
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux
     do icr=1,ncr
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              dflux_cr(ind_buffer(i),icr)=dflux_cr(ind_buffer(i),icr) &
                   &-flux(ind_cell(i),i3,j3,k3,icr,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do

     !-----------------------
     ! Right flux at boundary
     !-----------------------
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux
     do icr=1,ncr
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              dflux_cr(ind_buffer(i),icr)=dflux_cr(ind_buffer(i),icr) &
                   &+flux(ind_cell(i),i3+i0,j3+j0,k3+k0,icr,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
  end do
  ! End loop over dimensions
end if

end subroutine crstreamingfine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
