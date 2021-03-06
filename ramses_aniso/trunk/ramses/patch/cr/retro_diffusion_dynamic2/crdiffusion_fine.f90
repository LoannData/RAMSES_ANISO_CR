!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdiff_fine(ilevel,compute)
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
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,compute,igroup
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call crdifffine1(ind_grid,ngrid,ilevel,compute,igroup)
  end do

111 format('   Entering conduction_fine for level ',i2)

end subroutine crdiff_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdifffine1(ind_grid,ncache,ilevel,compute,igroup)
  use amr_commons
  use hydro_commons
  use poisson_commons
!  use radiation_parameters,ONLY:dt_imp
  use cooling_module
  implicit none
  integer::ilevel,ncache,compute,igroup
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! thermal conduction solver that compute energy flux. This flux is zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::indpass ! YD:test
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+7),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+7),save::upass
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp

  !real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::xloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),save::residu_loc=0.0d0

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  real(dp),dimension(1:nvector),save:: residu

  integer::neul=5
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim,length,twodx,twodx2,bnorm2,DCRmax_code
  real(dp)::eth,ekin,emag,erad,bb
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_omega,scale_tau,scale_kappa
  real(dp)::vx,vy,vz,dens,bnorm,bx,by,bz,vnorm,va,Ma,Dpara,kperp,kpar,deCRL,deCRR,deCRC,deCR
  real(dp)::deCRx,deCRy,deCRz,bdeCR,bdeCRx_l,bdeCRx_r,bdeCRy_l,bdeCRy_r
  integer::ind_cr1,l
  real(dp)::kparasub,kperpsub,Temp,gradPcr, gradPcr_x, gradPcr_y, gradPcr_z

!!$  scale_tau=1d0/(scale_t/scale_l**3) ! Time in code units

  ind_cr1=9
  if(twotemp)ind_cr1=10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t
  kpar=Dcr/scale_kappa
  kpar=1.d10*flinj*boxlen*scale_l/scale_kappa ! Linj*c/3
  DCRmax_code=Dmax/scale_kappa

  oneontwotondim = 1.d0/dble(twotondim)

  residu     = 0.0d0
  residu_loc = 0.0d0

  uloc =0.0d0
  upass=0.0d0

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  twodx =2.0d0*dx
  twodx2=twodx**2
 
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
        call interpol_hydro_cond(u1,ind1,u2,nbuffer)
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
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              indpass(ind_exist(i),i3,j3,k3)=ind_cell(i) !! YD: test
              uloc (ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              upass(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              facdx(ind_exist(i),i3,j3,k3)=1.0d0

              if(compute==2 .and. ivar==1)then 
                 upass(ind_exist(i),i3,j3,k3,ivar)=unew(ind_cell(i),3) ! y of PBiCGSTAB
                 ! neighbor cell at ilevel+1, put value to zero in vector
                 ! because it is considered as a Dirichlet BC, i.e. in the RHS
                 if(son(ind_cell(i))>0) upass(ind_exist(i),i3,j3,k3,ivar)=0.0d0
              endif
              if(compute==4 .and. ivar==2)then 
                 upass(ind_exist(i),i3,j3,k3,ivar)=unew(ind_cell(i),8)  ! z of PBiCGSTAB
                 ! neighbor cell at ilevel+1, put value to zero in vector
                 ! because it is considered as a Dirichlet BC, i.e. in the RHS
                 if(son(ind_cell(i))>0) upass(ind_exist(i),i3,j3,k3,ivar)=0.0d0
              endif

           end do
           do i=1,nbuffer
              indpass(ind_nexist(i),i3,j3,k3)=ind_son !! YD: test
              uloc (ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              upass(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              ! neighbor cell at ilevel+1, put value to zero in vector
              ! because it is considered as a Dirichlet BC, i.e. in the RHS
              if(compute==2 .and. ivar==1)upass(ind_nexist(i),i3,j3,k3,ivar)=0.0
              if(compute==4 .and. ivar==2)upass(ind_nexist(i),i3,j3,k3,ivar)=0.0
           end do
        end do
        
        do i=1,nbuffer
           if(interpol_type_cond.gt.0)then
              facdx(ind_nexist(i),i3,j3,k3)=1.0d0
           else
              facdx(ind_nexist(i),i3,j3,k3)=1.5d0
           endif
        end do

        !if(streaming_diffusion)then
        !   do i=1,nexist
        !      uloc(ind_exist(i),i3,j3,k3,5)=divu(ind_cell(i))              
        !  enddo
        !   do i=1,nbuffer
        !     uloc(ind_nexist(i),i3,j3,k3,5)=uloc(ind_nexist(i),i3,j3,k3,igroup)
        !   enddo
        !endif

        if(alfven_diff_coeff)then
           ! Compute Alfvenic Mach number = V/V_A and store it into uloc(:,3)
           do i=1,nexist
              dens   = uold(ind_cell(i),1)
              vx     = uold(ind_cell(i),2)/dens
              vy     = uold(ind_cell(i),3)/dens
              vz     = uold(ind_cell(i),4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx     = 0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              by     = 0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              bz     = 0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              bnorm  = (bx**2+by**2+bz**2)**0.5
              
              va = bnorm/(dens**0.5)
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/(Ma**3)
                 kperp = 1.0d0
              else
                 Dpara = kpar*(Ma**2)
                 kperp=Ma**4
              end if
              upass(ind_exist(i),i3,j3,k3,3)=Dpara
              upass(ind_exist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
           end do
           do i=1,nbuffer
              dens= u2(i,ind_son,1)
              vx = u2(i,ind_son,2)/dens
              vy = u2(i,ind_son,3)/dens
              vz = u2(i,ind_son,4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx = 0.5*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              by = 0.5*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              bz = 0.5*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              bnorm = (bx**2+by**2+bz**2)**0.5
              va = bnorm/dens**0.5
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/Ma**3
                 kperp = 1.0d0
              else
                 Dpara = kpar*Ma**2
                 kperp=Ma**4
              end if
              
              upass(ind_nexist(i),i3,j3,k3,3)=Dpara
              upass(ind_nexist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
           end do
        endif

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids
  if(variable_diff_coeff)then
   ! Compute diffusion coefficient due to streaming instability and store it into uloc(:,5,6)

#if NDIM==1
   ! Yo: It does not work with NCR>1
   j1=1;k1=1
   do l=1,ncache
   do i1=0,3
       dens   = uloc(l,i1,j1,k1,1)
       vx     = uloc(l,i1,j1,k1,2)/dens
       vy     = uloc(l,i1,j1,k1,3)/dens
       vz     = uloc(l,i1,j1,k1,4)/dens
       bx     = 0.5*(uloc(l,i1,j1,k1,6) + uloc(l,i1,j1,k1,nvar+1))
       by     = 0.5*(uloc(l,i1,j1,k1,7) + uloc(l,i1,j1,k1,nvar+2))
       bz     = 0.5*(uloc(l,i1,j1,k1,8) + uloc(l,i1,j1,k1,nvar+3))

       bnorm  = sqrt(bx**2+by**2+bz**2)

       va = bnorm/(dens**0.5)
       
       ! compute temperature
       ekin=0.5*dens*(vx**2+vy**2+vz**2)
       emag=0.0d0
       do idim=1,3
           bb=uloc(l,i1,j1,k1,5+idim)+uloc(l,i1,j1,k1,nvar+idim)
           emag=emag+bb*bb
       end do
       emag=emag*0.125d0
       erad=0.0d0
#if NENER>0
       do irad=1,nener
           erad=erad+uloc(l,i1,j1,k1,8+irad)
       end do
#endif
       
       eth=uloc(l,i1,j1,k1,5)-ekin-emag-erad
       
       temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!
       
       gradPcr_x=(uloc(l,i1+1,j1,k1,igroup)-uloc(l,i1-1,j1,k1,igroup))/twodx
       gradPcr = (gradPcr_x*bx)*(gamma_rad(irad)-1.0d0)/bnorm

       call subgridcr_diffusion(dens, Temp, bnorm, gradPcr, kparasub, kperpsub)

       upass(l,i1,j1,k1,nvar+5)=kparasub
       upass(l,i1,j1,k1,nvar+6)=kperpsub

       !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
   enddo
   end do
#endif  
#if NDIM==2
   ! Yo: It does not work with NCR>1
   k1=1
   do l=1,ncache
   do j1=0,3
       do i1=0,3
           dens   = uloc(l,i1,j1,k1,1)
           vx     = uloc(l,i1,j1,k1,2)/dens
           vy     = uloc(l,i1,j1,k1,3)/dens
           vz     = uloc(l,i1,j1,k1,4)/dens
           bx     = 0.5*(uloc(l,i1,j1,k1,6) + uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5*(uloc(l,i1,j1,k1,7) + uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5*(uloc(l,i1,j1,k1,8) + uloc(l,i1,j1,k1,nvar+3))
           
           bnorm  = sqrt(bx**2+by**2+bz**2)
           
           va = bnorm/(dens**0.5)
           
           ! compute temperature
           ekin=0.5*dens*(vx**2+vy**2+vz**2)
           emag=0.0d0
           do idim=1,3
               bb=uloc(l,i1,j1,k1,5+idim)+uloc(l,i1,j1,k1,nvar+idim)
               emag=emag+bb*bb
           end do
           emag=emag*0.125d0
           erad=0.0d0
#if NENER>0
           do irad=1,nener
               erad=erad+uloc(l,i1,j1,k1,8+irad)
           end do
#endif
           
           eth=uloc(l,i1,j1,k1,5)-ekin-emag-erad
           
           temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!
           
           gradPcr_x=(uloc(l,i1+1,j1  ,k1  ,igroup)&
               &    -uloc(l,i1-1,j1  ,k1  ,igroup))/twodx
           gradPcr_y=(uloc(l,i1  ,j1+1,k1  ,igroup)&
               &    -uloc(l,i1  ,j1-1,k1  ,igroup))/twodx
           gradPcr = (gradPcr_x*bx+gradPcr_y*by)*(gamma_rad(irad)-1.0d0)/bnorm
           
           !write(*,*) "Bnorm_outside_before_2 = ",bnorm
           call subgridcr_diffusion(dens, Temp, bnorm, gradPcr, kparasub, kperpsub)
           !write(*,*) "Bnorm_outside_after_2 = ",bnorm

           upass(l,i1,j1,k1,nvar+5)=kparasub
           upass(l,i1,j1,k1,nvar+6)=kperpsub

           !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
           
       end do
   end do
   end do
#endif  
#if NDIM==3
   ! Yo: It does not work with NCR>1
   do l=1,ncache
   do k1=0,3
   do j1=0,3
       do i1=0,3
           dens   = uloc(l,i1,j1,k1,1)
           vx     = uloc(l,i1,j1,k1,2)/dens
           vy     = uloc(l,i1,j1,k1,3)/dens
           vz     = uloc(l,i1,j1,k1,4)/dens
           bx     = 0.5*(uloc(l,i1,j1,k1,6) + uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5*(uloc(l,i1,j1,k1,7) + uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5*(uloc(l,i1,j1,k1,8) + uloc(l,i1,j1,k1,nvar+3))
           
           bnorm  = sqrt(bx**2+by**2+bz**2)

           !write(*,*) bnorm
           
           va = bnorm/(dens**0.5)
           
           ! compute temperature
           ekin=0.5*dens*(vx**2+vy**2+vz**2)
           emag=0.0d0
           do idim=1,3
               bb=uloc(l,i1,j1,k1,5+idim)+uloc(l,i1,j1,k1,nvar+idim)
               emag=emag+bb*bb
           end do
           emag=emag*0.125d0
           erad=0.0d0
#if NENER>0
           do irad=1,nener
               erad=erad+uloc(l,i1,j1,k1,8+irad)
           end do
#endif
           
           eth=uloc(l,i1,j1,k1,5)-ekin-emag-erad
           
           temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!

           !write(*,*) "temp : ",temp, eth, scale_t2, mu_gas
           
           !write (*,*) "igroup = ",igroup," nvar - igroup = ",nvar-igroup
           gradPcr_x=(uloc(l,i1+1,j1  ,k1  ,igroup)&
               &    -uloc(l,i1-1,j1  ,k1  ,igroup))/twodx
           gradPcr_y=(uloc(l,i1  ,j1+1,k1  ,igroup)&
               &    -uloc(l,i1  ,j1-1,k1  ,igroup))/twodx
           gradPcr_z=(uloc(l,i1  ,j1  ,k1+1,igroup)&
               &    -uloc(l,i1  ,j1  ,k1-1,igroup))/twodx
           gradPcr = (gradPcr_x*bx+gradPcr_y*by+gradPcr_z*bz)*(gamma_rad(irad)-1.0d0)/bnorm


           !write(*,*) "gradPcr_x = ",gradPcr_x," gradPcr_y = ",gradPcr_y," gradPcr_z = ",gradPcr_z," gradPcr =",gradPcr

           !write(*,*) "out : ",temp
           !write(*,*) gradPcr
           !write(*,*) "out : ",dens, bnorm
           !write(*,*) "Bnorm_outside_before_3 = ",bnorm
           !write(*,*) gradPcr_x, gradPcr_y, gradPcr_z, uloc(l,i1+1,j1  ,k1  ,igroup), uloc(l,i1-1,j1  ,k1  ,igroup)
           
           if (gradPcr .ne. 0.) then

           !write(*,*)"dens = ",dens,", temp = ",temp,", bnorm = ",bnorm,"gradPcr = ",gradPcr
           !write(*,*)"kparasub = ",kparasub,", kperpsub = ",kperpsub 
           call subgridcr_diffusion(dens, temp, bnorm, gradPcr, kparasub, kperpsub)
           !write(*,*) "Bnorm_outside_after_3 = ",bnorm
           !write(*,*) kparasub
           upass(l,i1,j1,k1,nvar+5)=kparasub
           upass(l,i1,j1,k1,nvar+6)=kperpsub

           !write(*,*) "kparasub = ",kparasub*scale_kappa," DCR = ",DCR
           !write(*,*) "gradPcr = ", gradPcr !* scale_d*scale_v**2/scale_l
           
           else 
           ! In the case where the CRs gradient is null, the self-generated diffusion coefficient is 
           ! infinite (high) in order to make the background diffusion coefficient dominate 
               upass(l,i1,j1,k1,nvar+5)=Dcr*1d4/scale_kappa
               upass(l,i1,j1,k1,nvar+6)=Dcr*1d4/scale_kappa
           endif 

           !if (kparasub < 0. .or. kperpsub < 0.) then 
           !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
           !endif
           
           
           
       end do
   end do
   end do
end do
   
#endif  
end if



  if(streaming_diffusion)then

#if NDIM==1
     ! Yo: It does not work with NCR>1
     j1=1;k1=1
     do l=1,ncache
        do i1=0,3
           deCR=uloc(l,i1+1,j1,k1,ind_cr1)-uloc(l,i1-1,j1,k1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCR) / twodx / sqrt(bnorm2)

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
     enddo
#endif     
#if NDIM==2
     ! Yo: It does not work with NCR>1
     k1=1
     do l=1,ncache
        do j1=0,3
        do i1=0,3
!!$           deCRx=uloc(l,i1+1,j1,k1,5)-uloc(l,i1-1,j1,k1,5)
!!$           deCRy=uloc(l,i1,j1+1,k1,5)-uloc(l,i1,j1-1,k1,5)
           deCRx=uloc(l,i1+1,j1,k1,ind_cr1)-uloc(l,i1-1,j1,k1,ind_cr1)
           deCRy=uloc(l,i1,j1+1,k1,ind_cr1)-uloc(l,i1,j1-1,k1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCRx + by*deCRy) / twodx / sqrt(bnorm2)
!!$           bdeCR= sqrt(deCRx*deCRx + deCRy*deCRy) / twodx

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
!!$           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,5)/DCRmax_code)
!!$           length=MIN(uloc(l,i1,j1,k1,5)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
        enddo
     enddo
#endif     
#if NDIM==3
     ! Yo: It does not work with NCR>1
     do l=1,ncache
        do k1=0,3
        do j1=0,3
        do i1=0,3
           deCRx=uloc(l,i1+1,j1  ,k1  ,ind_cr1)-uloc(l,i1-1,j1  ,k1  ,ind_cr1)
           deCRy=uloc(l,i1  ,j1+1,k1  ,ind_cr1)-uloc(l,i1  ,j1-1,k1  ,ind_cr1)
           deCRz=uloc(l,i1  ,j1  ,k1+1,ind_cr1)-uloc(l,i1  ,j1  ,k1-1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCRx + by*deCRy + bz*deCRz) / twodx / sqrt(bnorm2)
!!$           bdeCR= sqrt(deCRx**2 + deCRy**2 + deCRz**2) / twodx

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
        enddo
        enddo
     enddo
#endif     

  endif

  !-----------------------------------------------
  ! Compute energy flux due to thermal conduction
  !-----------------------------------------------
  call crdiff_split(upass,flux,dx,dx,dx,dt_imp,ncache,compute,facdx,igroup)
!!$  call crdiff_split(upass,flux,dx,dx,dx,dt_imp,ncache,compute,facdx,igroup,indpass)

  !-----------------------------------------------------
  ! update at level ilevel
  !-----------------------------------------------------
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
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then        
                 if(compute==1)then
                    residu_loc(i,i3   ,j3   ,k3   )=0.0d0
                 else if(compute==2)then
                    residu_loc(i,i3   ,j3   ,k3   )=unew(ind_cell(i),3)
                 else if(compute==3)then
                    residu_loc(i,i3   ,j3   ,k3   )=1.0d0
                 else if(compute==4)then
                    residu_loc(i,i3   ,j3   ,k3   )=unew(ind_cell(i),8)
                 end if
              
              endif
           end do
        end do
     end do
  end do

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
        ! Update conservative variables new state vector
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              residu_loc(i,i3   ,j3   ,k3   )=residu_loc(i,i3   ,j3   ,k3   )+ &
                   & (flux(i,i3   ,j3   ,k3   ,5,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,5,idim))
           endif
        end do
     end do
     end do
     end do
  end do

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
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 if(compute==1)then
                    unew(ind_cell(i),1     ) = -residu_loc(i,i3   ,j3   ,k3  ) ! r_0 of PBiCGSTAB
                    unew(ind_cell(i),nvar+2) = -residu_loc(i,i3   ,j3   ,k3  ) ! rzero of PBiCGSTAB
                    if(variable_diff_coeff)then
                     if(alfven_diff_coeff)then

                        !upass(ind_exist(i),i3,j3,k3,3)=Dpara
                        !upass(ind_exist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara

                        Dcrdiff(ind_cell(i),1)=1d0/(1d0/(upass(ind_exist(i),i3,j3,k3,3))+1d0/upass(i,i3,j3,k3,nvar+5))*scale_kappa
                        Dcrdiff(ind_cell(i),2)=MAX(upass(ind_exist(i),i3,j3,k3,4)*upass(ind_exist(i),i3,j3,k3,3), upass(i,i3,j3,k3,nvar+6))*scale_kappa
                        !Dcrdiff(ind_cell(i),2)=1d0/(1d0/(upass(ind_exist(i),i3,j3,k3,4)*upass(ind_exist(i),i3,j3,k3,3))+1d0/upass(i,i3,j3,k3,nvar+6))*scale_kappa

                         !Dcrdiff(ind_cell(i),1)=0.
                         !Dcrdiff(ind_cell(i),2)=0. !Do somethin !!!!!
                     else
                         Dcrdiff(ind_cell(i),1)=1d0/(1d0/(DCR/scale_kappa)+1d0/upass(i,i3,j3,k3,nvar+5))*scale_kappa
                         Dcrdiff(ind_cell(i),2)=1d0/(1d0/(k_perp*DCR/scale_kappa)+1d0/upass(i,i3,j3,k3,nvar+6))*scale_kappa
                         ! Tester avec le cas MAX(Dperp, Dself_perp)
                         
                         !if (DCR > uloc(i,i3,j3,k3,nvar+5)*scale_kappa) then
                         !write(*,*) "DCR = ",DCR," Dself,par = ",uloc(i,i3,j3,k3,nvar+5)*scale_kappa
                         !endif

                         !Dcrdiff(ind_cell(i),1)=DCR
                         !Dcrdiff(ind_cell(i),2)=k_perp*DCR
                     endif

                     !write(*,*) "Dcrdiff(ind_cell(i),1) = ",uloc(i,i3,j3,k3,nvar+5)*scale_kappa
                     !write(*,*) "Dcrdiff(ind_cell(i),2) = ",uloc(i,i3,j3,k3,nvar+6)*scale_kappa
                 end if
                  else if(compute==2)then
                    unew(ind_cell(i),6     ) = residu_loc(i,i3   ,j3   ,k3   ) ! vi of PBiCGSTAB
                 else if(compute==3)then ! Compute 1/A
                    unew(ind_cell(i),4     ) = 1.0d0!/residu(i)
                 else if(compute==4)then
                    unew(ind_cell(i),nvar+1) = residu_loc(i,i3   ,j3   ,k3   ) ! t of PBiCGSTAB
                 end if
              endif
           end do
        end do
     end do
  end do


end subroutine crdifffine1

! -------------- Subroutines for the calculation of the sub-resolution diffusion coefficient ----------------
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine subgridcr_diffusion(rho_sim, T_sim, B0_sim, gradPcr, kpara, kperp)
   use cooling_module
   use hydro_parameters
   implicit none
   ! - Some constants ---------------------------------------
   real::mp, e, GeV, pi 
   
   ! - Variables for NN interpolation -----------------------
   integer::k, n, kvalues  
   real(dp)::X_test, Y_test, rho_test, T_test
   real(dp)::calc_value1, calc_value2, calc_value3, betacr 
   real(dp),dimension(1:7)::mi,ni,mn,nn,Xion,T,rho
   real(dp),dimension(1:7)::X_train, Y_train, distance  
   real(dp),dimension(1:7):: f_value1, f_value2, f_value3
   real(dp),dimension(1:7)::dist_k, x_k, y_k, val1, val2, val3
   
   ! - Variables from RAMSES-ISM ----------------------------
   real(dp)::rho_sim, T_sim, B0_sim, gradPcr, scale_m, scale_n
   real(dp)::ECR_sim, VT, Linj, theta, Dfg, Dnll
   ! - Variables for RAMSES-ISM -----------------------------
   real(dp)::Gamma_in,rho_i,rho_n,rmi,rmn,chi,rXion,xi_n,rnn,rni,rg
   real(dp)::nuin,Va,Vai,Va_temp, Va_final,Turb
   !parameter(k_reel = selected_real_kind(6,90))
   real(dp),intent(out)::kpara,kperp
   real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_omega,scale_tau,scale_kappa
   
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   scale_kappa = scale_l**2/scale_t
   scale_m = scale_d*scale_l**3
   scale_n = scale_l**3
   
   ! - Some constants 
   mp = 1.6726d-24 ! Proton mass [g]
   e = 4.80326d-10 ! Electron charge [statcoulomb]
   GeV = 0.00160218 ! GeV in Erg
   pi = acos(-1.0d0)
   
   rho_sim = rho_sim * scale_d
   gradPcr = gradPcr * scale_d*scale_v**2/scale_l
   
   !if (gradPcr .ne. 0.) then 
   !    write(*,*) "gradPcr = ",gradPcr," erg/cm^3"
   !endif
   
   
   !write(*,*) "Bnorm_inside_before = ",B0_sim
   B0_sim  = B0_sim  * sqrt(4.*pi*scale_d*scale_v**2)
   !write(*,*) "pi = ",pi," scale_d = ",scale_d," scale_v = ",scale_v
   !write(*,*) "Bnorm_inside_after = ",B0_sim
   
   !write(*,*) "In : ",T_sim
   
   
   
   ! - kNN interpolation algorithm parameters ----------------
   k = 7  ! Number of nearest neighbors (Only working for k = 1 for instance !!!)
   n = 20 ! Distance ponderation parameter 
   
   ! - Variables from RAMSES-ISM ----------------------------
   !!$rho_sim = 1.d-21 ! g.cm^-3 
   !!$T_sim   = 50.    ! Kelvin
   !!$B0_sim = 10.d-6   ! Gauss 
   ECR_sim = 1.*GeV ! Erg     
   !!$gradPcr = 1d-29  ! Erg/cm^4
   
   ! - Observationnal and simulation values ------------------
   !T    = [  8000.,     50.,     50.,     30.,     20.]
   !mi   = [ 1.0*mp, 12.0*mp, 12.0*mp, 29.0*mp, 29.0*mp]
   !mn   = [1.21*mp, 1.21*mp, 1.67*mp, 2.12*mp, 2.12*mp]
   !ni   = [  0.007,  2.4e-2,    0.03,    3e-2,    0.01]
   !Xion = [   0.02,    8e-4,    1e-4,    1e-5,    1e-6]
   !-----------------------------------------------------------
   !T    = [1000000.,  8000.,     50.,     50.,     30.,     20.]
   !mi   = [ 1.0*mp, 1.0*mp, 12.0*mp, 12.0*mp, 29.0*mp, 29.0*mp]
   !mn   = [1.21*mp,1.21*mp, 1.21*mp, 1.67*mp, 2.12*mp, 2.12*mp]
   !ni   = [   0.01,  0.007,  2.4e-2,    0.03,    3e-2,    0.01]
   !Xion = [    0.99,  0.02,    8e-4,    1e-4,    1e-5,    1e-6]
   T    = [1000000.,  10000., 8000.,     50.,     50.,     30.,     20.]
   mi   = [ 1.0*mp, 1.0*mp, 1.0*mp, 12.0*mp, 12.0*mp, 29.0*mp, 29.0*mp]
   mn   = [1.21*mp, 1.21*mp, 1.21*mp, 1.21*mp, 1.67*mp, 2.12*mp, 2.12*mp]
   ni   = [   0.01,  99., 0.007,  2.4e-2,    0.03,    3e-2,    0.01]
   Xion = [    0.99, 0.99, 0.02,    8e-4,    1e-4,    1e-5,    1e-6]
   
   do kvalues=1,7
       rho(kvalues) = (mi(kvalues)*ni(kvalues) + mn(kvalues)*ni(kvalues)*(Xion(kvalues)**(-1) - 1.))
       !rho(kvalues) = ni(kvalues)/Xion(kvalues)*mi(kvalues)
       X_train(kvalues) = log10(T(kvalues))
       Y_train(kvalues) = log10(rho(kvalues))
       !write(*,*) "rho = ",rho(kvalues),", log10(T) = ",X_train(kvalues)," log10(rho) = ",Y_train(kvalues)
   end do 
   
   X_test = log10(T_sim)
   Y_test = log10(rho_sim)

   
   
   ! - Algorithm ---------------------------------------------
   do kvalues=1,7
       f_value1(kvalues) = log10(Xion(kvalues))
       f_value2(kvalues) = mn(kvalues)
       f_value3(kvalues) = mi(kvalues)
   end do 
   
   !write(*,*) X_test
   call getValue_linear(X_train, Y_train, distance, f_value1, f_value2, f_value3, X_test, Y_test, k , dist_k, x_k, y_k, val1, val2, val3, n, calc_value1, calc_value2, calc_value3)
   !write(*,*) X_test
   rXion = 10.**(calc_value1) ! From log axis to lin axis 
   rmn = calc_value2
   rmi = calc_value3

   !write(*,*) "Tsim = ",T_sim,", rho_sim = ",rho_sim/mp,", Xion = ",rXion
   !write(*,*) "==================================================================================="
   
   !write(*,*) rXion, rmn/mp, rmi/mp


   ! Case where we are in young SNRs -> Dself = Dbohm 
   rg    = ECR_sim/(e*B0_sim)
   if (T_sim > 1.0D6) then 

      !kpara = 4*clight/(3*3.14720)*rg / scale_kappa 
      kpara = Dcr / 10. / scale_kappa 
      !write(*,*) "Kappa = ",kpara
      kperp = kpara*1.d-2 ! In order to avoid scheme divergences
   
   ! If we are outside young SNRs 
   else 
      ! In fully ionized phases 
      if (rXion > 0.9) then 
         if (abs(gradPcr) > 1e-50) then
            Va  = B0_sim/sqrt(4.*pi*rho_sim)
            VT = 10*1e5 ! [10 km/s in cm/s]
            Linj = 50*3.086e18 ! [pc]
            theta = 9.2e-14*T_sim
            Dfg  = 0.9*(VT/Va)**(3./2)*(clight*rg**0.5/(abs(gradPcr)*Linj**0.5))
            Dnll = 0.3*theta**0.25*(clight/Va)**0.5*clight*rg**0.5/((abs(gradPcr))**0.5)
            !write(*,*) "VT = ",VT,", Linj = ",Linj,", theta = ",theta,", Va_final = ",Va
            !write(*,*) "Dfg = ",Dfg,", Dnll = ",Dnll,", Xion = ",rXion,", rg = ",rg
            kpara = (Dfg**(-1) + Dnll**(-1))**(-1) / scale_kappa 
            kperp = kpara*1.d-2 ! In order to avoid scheme divergences
         else 
            kpara = 1e30 / scale_kappa
            kperp = kpara*1.d-2 ! In order to avoid scheme divergences
         end if 
      ! In partially ionized phases
      else 
         ! - Basic calculations 
         chi   = (rmn/rmi)*(rXion**(-1) - 1.) 
         rho_n = rho_sim/(1. + chi**(-1))
         rho_i = chi**(-1)/(1. + chi**(-1))*rho_sim
         xi_n  = 1./(1 + chi**(-1))
         rnn   = rho_n/rmn 
         rni   = rho_i/rmi 
         

         ! - nu_in and Gamma_in calculation 
         Gamma_in = 0.
         nuin  = 2*rnn*8.4d-9*(T_sim/1.d4)**(0.4) 
         Vai = B0_sim/sqrt(4.*pi*rho_i)
         Gamma_in = -nuin/2.
         Va_final = Vai
         
         ! - Energy density of waves (Turbulence level)
         Turb = abs(Va_final*abs(gradPcr)/(2.*Gamma_in)/(B0_sim**2/(8.*pi)))
         betacr = sqrt(1 - (1/(1 + ECR_sim/(mp*clight)**2))**2)

         ! - Self-generated diffusion coefficient 
         kpara = 4.*pi*rg*betacr*clight/(3.*Turb) / scale_kappa
         kperp = kpara*1.d-2 ! In order to avoid scheme divergences 

         !write(*,*) "kpara = ",kpara*scale_kappa
      end if 
   end if 



   






   
   ! - Basic calculations 
   !chi   = (rmn/rmi)*(rXion**(-1) - 1.) 
   !rho_n = rho_sim/(1. + chi**(-1))
   !rho_i = chi**(-1)/(1. + chi**(-1))*rho_sim
   !xi_n  = 1./(1 + chi**(-1))
   !rnn   = rho_n/rmn 
   !rni   = rho_i/rmi 
   !rg    = ECR_sim/(e*B0_sim)
   
   ! - nu_in and Gamma_in calculation 
   !Gamma_in = 0.
   !nuin  = 2*rnn*8.4d-9*(T_sim/1.d4)**(0.4) 
   !Va  = B0_sim/sqrt(4.*pi*rho_sim)
   !Vai = B0_sim/sqrt(4.*pi*rho_i)
   !Va_temp = B0_sim/sqrt(4.*pi*rho_i)
   !Va_temp = B0_sim/(4.*pi*rho_sim)

   !Gamma_in = -nuin/2.
   !Va_final = Vai
   
   ! - Energy density of waves (Turbulence level)
   !Turb = abs(Va_final*abs(gradPcr)/(2.*Gamma_in)/(B0_sim**2/(8.*pi)))
   


   !! - Diffusion coefficients 
   !!betacr = sqrt(1 - (1/(1 + ECR_sim/(mp*clight)**2))**2)
   !!kpara = 4.*pi*rg*betacr*clight/(3.*Turb) / scale_kappa

   ! - Diffusion coefficients 
   !if (rXion < 0.9) then 
   !   !write(*,*)  4.*pi*rg*clight/(3.*Turb)
   !   betacr = sqrt(1 - (1/(1 + ECR_sim/(mp*clight)**2))**2)
   !   kpara = 4.*pi*rg*betacr*clight/(3.*Turb) / scale_kappa
   !   !write(*,*) "Hello !!"
   !else 
   !   if (abs(gradPcr) > 1e-50) then 
   !      VT = 10*1e5 ! [10 km/s in cm/s]
   !      Linj = 50*3.086e18 ! [pc]
   !      theta = 9.2e-14*T_sim
   !      Dfg  = 0.9*(VT/Va)**(3./2)*(clight*rg**0.5/(abs(gradPcr)*Linj**0.5))
  !       Dnll = 0.3*theta**0.25*(clight/Va)**0.5*clight*rg**0.5/((abs(gradPcr))**0.5)
   !      !write(*,*) "VT = ",VT,", Linj = ",Linj,", theta = ",theta,", Va_final = ",Va
   !      !write(*,*) "Dfg = ",Dfg,", Dnll = ",Dnll,", Xion = ",rXion,", rg = ",rg
   !      kpara = (Dfg**(-1) + Dnll**(-1))**(-1) / scale_kappa 
   !      !write(*,*) "kpara = ",kpara*scale_kappa 
   !   else 
   !      kpara = 1e30 / scale_kappa
   !   end if 
   !end if 


   !write(*,*) "kpara = ",kpara*scale_kappa

      
   
   !write(*,*) "betacr = ",betacr
   !kperp = kpara*Turb**2
   
   
   !write(*,*) "-----------------------------------------------------------------"
   !write(*,*) "Tested values : rho = ",rho_sim," g.cm^-3, T = ",T_sim," Kelvin" 
   !write(*,*) "Result : Xion = ",rXion," mi = ",rmi/mp," mp, mn = ",rmn/mp," mp" 
   !write(*,*) "chi = ",chi,", rho_n = ",rho_n," g.cm^-3, rho_i = ",rho_i," g.cm^-3" 
   !write(*,*) "Xi_n = ",xi_n,", nn = ",rnn," cm^-3, ni = ",rni," cm^-3" 
   !write(*,*) "rg = ",rg," cm, nu_in = ",nuin," s^-1, Va_final",Va_final," cm.s^-1" 
   !write(*,*) "Gamma_in = ",Gamma_in," s^-1, I(k) = ",Turb
   !write(*,*) "kpara = ",kpara*scale_kappa," cm^2.s^-1, kperp = ",kperp*scale_kappa," cm^2.s^-1" 
   !write(*,*) "-----------------------------------------------------------------"
   
   end subroutine subgridcr_diffusion
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine getValue_linear(X_train, Y_train, distance, value1, value2, value3, X_test, Y_test, k , dist_kk, x_k, y_k, val1, val2, val3, n, calc_value1, calc_value2, calc_value3)
       use cooling_module
   implicit none
   ! - Input values -----------------------------------------
   integer::k, n 
   real(dp),intent(in)::X_test, Y_test
   real(dp),dimension(1:7)::X_train, Y_train, value1, value2, value3, distance 
   real(dp),dimension(1:k)::x_k, y_k, val1, val2, val3, dist_kk
   
   ! - Output values ----------------------------------------
   real(dp)::calc_value1, calc_value2, calc_value3 
   
   ! - InCode values ----------------------------------------
   real(dp)::distinv
   real(dp)::temp_value1
   real(dp)::temp_value2
   real(dp)::temp_value3
   integer::ivalues
   
   distinv = 0.0D0 
   temp_value1 = 0.0D0
   temp_value2 = 0.0D0
   temp_value3 = 0.0D0 
   
   
   !write(*,*) "Valeur X : ",X_test, "Valeur Y :", log10(10**(Y_test)/1.6726d-24),", k = ",k
   calc_value1 = 0.
   calc_value2 = 0.
   calc_value3 = 0.
   !write(*,*) "Before : ",calc_value1, calc_value2, calc_value3
   
   call selectKNN(X_train, Y_train, X_test, Y_test, k, value1, value2, value3, dist_kk, x_k, y_k, val1, val2, val3)

   !do ivalues=1,k 
   !   write(*,*) "dist_k = ",dist_kk(ivalues)
   !end do 
   
   
   
   if (k.EQ.1) then 
       calc_value1 = val1(1)
       !write(*,*) "val1 = ",val1(1)
       calc_value2 = val2(1)
       calc_value3 = val3(1)
   end if 
   if (k > 1) then 
       do ivalues=1,k 
           distinv = distinv + 1./(dist_kk(ivalues)**n) 
           temp_value1 = temp_value1 + val1(ivalues)/(dist_kk(ivalues)**n)
           !write(*,*) "temp_value1 = ",temp_value1, ", other = ",val1(ivalues)/(dist_k(ivalues)**n)
           !write(*,*) "dist_k = ",dist_kk(ivalues)
           temp_value2 = temp_value2 + val2(ivalues)/(dist_kk(ivalues)**n)
           temp_value3 = temp_value3 + val3(ivalues)/(dist_kk(ivalues)**n)
       end do 
       calc_value1 = temp_value1/distinv 
       calc_value2 = temp_value2/distinv 
       calc_value3 = temp_value3/distinv
   end if 
   
   
 
   
   !write(*,*) "After : ",calc_value1, calc_value2, calc_value3
   
   
   
   end subroutine getValue_linear 
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine selectKNN(X_train, Y_train, X_test, Y_test, k, value1, value2, value3, dist_k, x_k, y_k, val1, val2, val3) 
       use cooling_module
   implicit none
   ! - Input values ------------------------------------------
   integer::k, iivalues 
   real(dp),dimension(1:7)::X_train, Y_train, value1, value2, value3, distance  
   real(dp)::X_test, Y_test
   
   ! - Output values ----------------------------------------- 
   real(dp),dimension(1:k)::dist_k, x_k, y_k, val1, val2, val3 
   
   
   call getDistance(X_train, Y_train, X_test, Y_test, value1, value2, value3, distance)
   
   !write(*,*) "========================================================================"
   !write(*,*) "k = ",k
   do iivalues=1,k 
       dist_k(iivalues) = distance(iivalues)
       x_k(iivalues)    = X_train(iivalues)
       y_k(iivalues)    = Y_train(iivalues)
       val1(iivalues)   = value1(iivalues)
       val2(iivalues)   = value2(iivalues)
       val3(iivalues)   = value3(iivalues)

       !write(*,*) "sub dist_k = ",distance(ivalues)
   
       !write(*,*) "dist_k = ",dist_k(ivalues)," X_train = ", X_train(ivalues),", Y_train = ",Y_train(ivalues), "Xval = ",val1(ivalues)
       !write(*,*) "vals = ",val1(ivalues), val2(ivalues), val3(ivalues)
   end do 
   
   end subroutine selectKNN
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine getDistance(X_train, Y_train, X_test, Y_test, value1, value2, value3, distance) 
       use cooling_module
   implicit none 
   ! - Input values ------------------------------------------
   integer::ivalues,jvalues
   real(dp),dimension(1:7)::X_train, Y_train, value1, value2, value3  
   real(dp)::X_test,Y_test
   
   ! - Output values -----------------------------------------
   real(dp),dimension(1:7)::distance
   integer,dimension(1:7)::order 
   real(dp)::buffer_value
   
   do ivalues = 1,7 
       distance(ivalues) = sqrt((X_train(ivalues)-X_test)**2. + (Y_train(ivalues)-Y_test)**2.)
       !write(*,*) ivalues," : ",distance(ivalues)
       !write(*,*) "Tsim = ",X_test,"  Ttrain = ",X_train(ivalues)
       !write(*,*) "rhosim = ",Y_test," rhotrain = ",Y_train(ivalues)
       !write(*,*) distance(ivalues)
       order(ivalues) = ivalues
   end do
   
   do ivalues = 1,7
       do jvalues = 1,7
           if (ivalues.NE.jvalues) then 
               if ((distance(ivalues) <= distance(jvalues)).AND.(ivalues > jvalues)) then
   
                   buffer_value = distance(ivalues)
                   distance(ivalues) = distance(jvalues)
                   distance(jvalues) = buffer_value
   
                   buffer_value = X_train(ivalues)
                   X_train(ivalues) = X_train(jvalues)
                   X_train(jvalues) = buffer_value
   
                   buffer_value = Y_train(ivalues)
                   Y_train(ivalues) = Y_train(jvalues)
                   Y_train(jvalues) = buffer_value
   
                   buffer_value = order(ivalues)
                   order(ivalues) = order(jvalues)
                   order(jvalues) = buffer_value

                   buffer_value = value1(ivalues)
                   value1(ivalues) = value1(jvalues)
                   value1(jvalues) = buffer_value

                   buffer_value = value2(ivalues)
                   value2(ivalues) = value2(jvalues)
                   value2(jvalues) = buffer_value

                   buffer_value = value3(ivalues)
                   value3(ivalues) = value3(jvalues)
                   value3(jvalues) = buffer_value
               end if 
           endif
       end do
   end do
   end subroutine getDistance

