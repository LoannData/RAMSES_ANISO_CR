!*************************************************************************
SUBROUTINE sink_wind_feedback(ilevel, dt)

! This routine adds photons from radiating sinks to appropriate cells in 
! the  hydro grid.
! ilegel =>  grid level in which to perform the feedback
! ti     =>  initial time for the timestep (code units)
! dt     =>  real timestep length in code units
!-------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use rt_parameters
  implicit none
  integer:: ilevel
  real(dp):: dt
  integer:: igrid, jgrid, ipart, jpart, next_part
  integer:: i, ig, ip, npart1, npart2, icpu
  integer,dimension(1:nvector),save:: ind_grid, ind_part, ind_grid_part
!-------------------------------------------------------------------------
  if(nsink .le. 0 ) return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)
        npart2=0
        if(npart1 > 0)then
          ipart = headp(igrid) 
           ! Loop over particles
           do jpart = 1, npart1
              next_part = nextp(ipart)
              ! only sink cloud particles
              if(idp(ipart) .lt. 0) then 
                 npart2 = npart2+1     
              endif
              ipart = next_part
           end do
        endif

        ! Gather sink and cloud particles within the grid
        if(npart2 > 0)then        
           ig = ig+1
           ind_grid(ig) = igrid
           ipart = headp(igrid)
           ! Loop over particles
           do jpart = 1, npart1
              next_part = nextp(ipart)
              if(idp(ipart) .lt. 0) then
                 if(ig==0)then
                    ig=1      
                    ind_grid(ig)=igrid
                 end if
                 ip = ip+1
                 ind_part(ip) = ipart
                 ind_grid_part(ip) = ig
              endif
              if(ip == nvector)then
                 call sink_RT_vsweep( &
                         ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
                 ip = 0
                 ig = 0
              end if
              ipart = next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid = next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip > 0) then
        call sink_RT_vsweep( &
             ind_grid,ind_part,ind_grid_part,ig,ip,dt,ilevel)
     endif
  end do 
  ! End loop over cpus

111 format('   Entering sink_rt_feedback DEBUG for level ',I2)

END SUBROUTINE sink_wind_feedback

!*************************************************************************
SUBROUTINE sink_wind_vsweep(ind_grid,ind_part,ind_grid_part,ng,np,dt,ilevel)

! This routine is called by subroutine sink_rt_feedback.
! Each sink and cloud  particle dumps a number of photons into the nearest grid cell
! using array rtunew.
! Radiation is injected into cells at level ilevel, but it is important
! to know that ilevel-1 cells may also get some radiation. This is due
! to sink and cloud particles that have just crossed to a coarser level. 
!
! ind_grid      =>  grid indexes in amr_commons (1 to ng)
! ind_part      =>  sink indexes in pm_commons(1 to np)
! ind_grid_part =>  points from star to grid (ind_grid) it resides in
! ng            =>  number of grids
! np            =>  number of sink and cloud particles
! dt            =>  timestep length in code units
! ilevel        =>  amr level at which we're adding radiation
!-------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use rt_hydro_commons
  use rt_parameters
  !use rt_cooling_module, only:iIR
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  real(dp)::dt, sourcemass
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ip,isink
  real(dp)::dx,dx_loc,scale,vol_loc,vol_cgs
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,ngroups),save::part_NpInp
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,3),save::id=0,igd=0,icd=0
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::Ep2Np
  ! units and temporary quantities
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, scale_Np  & 
            , scale_Fp, age, z, scale_inp, scale_Nphot, dt_Gyr           &
            , dt_loc_Gyr, scale_msun
  real(dp),parameter::vol_factor=2**ndim   ! Vol factor for ilevel-1 cells
  ! changed:
  real(dp),dimension(1:nvector),save::dn
  real(dp)::hnu,scale_energy,ener
  ! Emission law based on Vacca 1996
  real(dp),parameter::qh_power=3.184
  real(dp),parameter::qh_const=10d0**43.97
  ! Sink cluster totals
  real(dp)::scluster,mcluster,wcluster,ssink
  ! Particle position
  real(dp),dimension(1:3)::rvec 

!-------------------------------------------------------------------------
! if(.not. metal) z = log10(max(z_ave*0.02, 10.d-5))![log(m_metals/m_tot)]
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_Np, scale_Fp)
  dt_Gyr = dt*scale_t*sec2Gyr
  ! Mesh spacing in ilevel
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1) = dble(icoarse_min)
  if(ndim>1)skip_loc(2) = dble(jcoarse_min)
  if(ndim>2)skip_loc(3) = dble(kcoarse_min)
  scale = boxlen/dble(nx_loc) 
  dx_loc = dx*scale
  vol_loc = dx_loc**ndim
  scale_inp = rt_esc_frac * scale_d / scale_np / vol_loc / m_sun    
  scale_nPhot = vol_loc * scale_np * scale_l**ndim / 1.d50
  scale_msun = scale_d * scale_l**ndim / m_sun    
  vol_cgs = (dx_loc*scale_l)**ndim



  ! Ep2Np=(scale_d * scale_v**2)/( scale_Np * group_egy(iIR) * ev_to_erg)


  ! Lower left corners of 3x3x3 grid-cubes (with given grid in center)
  do idim = 1, ndim
     do i = 1, ng
        x0(i,idim) = xg(ind_grid(i),idim) - 3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i) = father(ind_grid(i))
  end do
  call get3cubefather(&
          ind_cell, nbors_father_cells, nbors_father_grids, ng, ilevel)

  ! Rescale position of stars to positions within 3x3x3 cell supercube
  do idim = 1, ndim
     do j = 1, np
        x(j,idim) = xp(ind_part(j),idim)/scale + skip_loc(idim)
        x(j,idim) = x(j,idim) - x0(ind_grid_part(j),idim)
        x(j,idim) = x(j,idim)/dx 
     end do
  end do

   ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim) = x(j,idim)
     end do
  end do

  ! Compute parent grids
  do idim = 1, ndim
     do j = 1, np
        igd(j,idim) = id(j,idim)/2
     end do
  end do
  do j = 1, np
     kg(j) = 1 + igd(j,1) + 3*igd(j,2) + 9*igd(j,3) ! 1 to 27
  end do
  do j = 1, np
     igrid(j) = son(nbors_father_cells(ind_grid_part(j),kg(j))) 
  end do

  ! Check if particles are entirely in level ilevel.
  ok(1:np) = .true.
  do j = 1, np
     ok(j) = ok(j) .and. igrid(j) > 0
  end do

  ! Compute parent cell position within it's grid
  do idim = 1, ndim
     do j = 1, np
        if( ok(j) ) then
           icd(j,idim) = id(j,idim) - 2*igd(j,idim)
        end if
     end do
  end do
  do j = 1, np
     if( ok(j) ) then
        icell(j) = 1 + icd(j,1) + 2*icd(j,2) + 4*icd(j,3)
     end if
  end do

  ! Compute parent cell adress
  do j = 1, np
     if( ok(j) )then
        indp(j) = ncoarse + (icell(j)-1)*ngridmax + igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
     end if
  end do
  
  ! Calculate photon emission totals for model
  scluster = 0d0
  mcluster = 0d0
  wcluster = 0d0
  do isink=1,nsink
     ! Handwavy fit assuming dominant star in sink is msink/3
     call vaccafit(0.3*msink(isink)*scale_msun,ssink)
     mcluster = mcluster + msink(isink)
     wcluster = wcluster + ssink
  end do
  ! Analytic fit to emission from IMF using Vacca 1996
  scluster = 4.6d46 * (mcluster*scale_msun)**1.0746

  ! Run through each particle and deposit stuff onto the grid
  do j=1,np
     if( ok(j) ) then                                      !   ilevel cell
        ! NEED:
        ! MASS LOSS RATE IN g / s
        ! LUMINOSITY IN erg / s
        ! WIND TEMPERATURE -> ENERGY
        ! NOTE THAT EWIND HAS TO INCLUDE PWIND
        ! Get sink index
        isink=-idp(ind_part(j))
        ! Use weighted fit
        ! sink emission = scluster * vacca(msink/3) / sum(vacca(msink/3))
        call vaccafit(0.3*msink(isink)*scale_msun,ssink)
        dn(j) = scluster * ssink / wcluster
        ! Get normal vector from sink to cloud particle
        rvec = xp(ind_part(j),:) - xsink(isink)
        rvec /= sqrt(sum(rvec**2))
        ! Deposit stuff onto the grid
        unew(indp(j),1)=unew(indp(j),1)+mwind/vol_loc
        unew(indp(j),2:4)=unew(indp(j),2:4)+pwind/vol_loc*rvec
        unew(indp(j),5)=unew(indp(j),5)+ewind/vol_loc
                        
        ! Scale that photon emission rate right up
        !dn(j) = dn(j) * dt*scale_t / dble(ncloud_sink) / vol_cgs / scale_Np
        !rtunew(indp(j),1)=rtunew(indp(j),1)+dn(j)
     endif
  end do

END SUBROUTINE sink_wind_vsweep

SUBROUTINE vaccafit(M,S)
  ! M - stellar mass / solar masses
  ! S - photon emission rate in / s

  use amr_parameters

  real(dp),intent(in)::M
  real(dp),intent(out)::S
  real(dp),parameter::stf_K=9.634642584812752d48 ! s**(-1)
  real(dp),parameter::stf_m0=2.728098824280431d1 ! Msun
  real(dp),parameter::stf_a=6.840015602892084d0
  real(dp),parameter::stf_b=4.353614230584390d0
  real(dp),parameter::stf_c=1.142166657042991d0 
  
  S = stf_K * (M / stf_m0)**stf_a / (1 + (M / stf_m0)**stf_b)**stf_c

END SUBROUTINE
