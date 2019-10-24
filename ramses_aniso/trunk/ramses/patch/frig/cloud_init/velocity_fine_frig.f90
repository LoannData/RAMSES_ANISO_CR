!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
!  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,j
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,Cwnm
  real(dp)::dx_min, fact, Emag,Emag0

! STG HACK - ignore if not MHD
! TODO: Take boundary cleaner and use for non-MHD solver
#ifndef SOLVERmhd
  return
#endif 

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Cwnm = sqrt(8000./scale_T2)

  pi=ACOS(-1.0d0)

  time = t * Cwnm / boxlen



  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  dx_min = (0.5D0**levelmin)/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        

       do i=1,ngrid



        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),1) .lt. 0d0) then
           write(*,*) "DENSITY < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TEMPERATURE < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN BEFORE VELOCITY_FINE, OH NO", ind_cell(i),j,uold(ind_cell(i),1)
              call clean_stop
           end if
        end do


        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .lt. 2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),2))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+1) = uold(ind_cell_vois,6)
 

              uold(ind_cell(i),6)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+1) =  uold(ind_cell_vois,6) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),6) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           endif




           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif



        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .gt. boxlen-2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),1))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
             ind_cell_vois = ind_cell_vois + ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
              uold(ind_cell(i),6) = uold(ind_cell_vois,nvar+1)
 
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),9) of the preceeding case 
              uold(ind_cell(i),6) =  uold(ind_cell(i),7) + uold(ind_cell(i),8)  + uold(ind_cell_vois,nvar+1) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif





        !impose vanishing gradient conditions at the y  faces
        if(  xx(i,2) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),4))


           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (1) is 4 (2)
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+2) = uold(ind_cell_vois,7)
 
              uold(ind_cell(i),7)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),nvar+1 ) + uold(ind_cell_vois,7) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),7) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        if(  xx(i,2) .gt. boxlen-2.*dx_min ) then 
             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),3))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (4) is 1 (2)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
             ind_cell_vois = ind_cell_vois + 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
              uold(ind_cell(i),7) = uold(ind_cell_vois,nvar+2)
 
              uold(ind_cell(i),nvar+2)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),10) of the preceeding case 
              uold(ind_cell(i),7) =  uold(ind_cell(i),6 ) + uold(ind_cell_vois,nvar+2) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif



        if(  xx(i,3) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),6))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 5 (6) is 1 (2)
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+3) = uold(ind_cell_vois,8)
 
              uold(ind_cell(i),8)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),7) 
           else
              !should be equal to uold(ind_cell(i),8) of the preceeding case 
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i), nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell_vois,8) - uold(ind_cell(i),6)  - uold(ind_cell(i),7) 

              !ensure div B
              uold(ind_cell(i),8) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),7) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))


           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif


        if(  xx(i,3) .gt. boxlen-2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),5))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 1 (2) is 5 (6)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
             ind_cell_vois = ind_cell_vois + 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
              uold(ind_cell(i),8) = uold(ind_cell_vois,nvar+3)
 
              uold(ind_cell(i),nvar+3)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 
           else
              !should be equal to uold(ind_cell(i),nvar+3) of the preceeding case 
              uold(ind_cell(i),8) =  uold(ind_cell(i), 6) + uold(ind_cell(i),7) + uold(ind_cell_vois,nvar+3) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+2) 

              !ensure div B
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TEMPERATURE < 0 AFTER VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN AFTER VELOCITY_FINE, OH NO", ind_cell(i),j
              call clean_stop
           end if
        end do


       enddo



       
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
