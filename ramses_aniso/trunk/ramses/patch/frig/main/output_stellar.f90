subroutine backup_stellar(filename)
    use amr_commons
    use pm_commons

    use feedback_module

    implicit none

#ifndef WITHOUTMPI
    include 'mpif.h'
#endif 

    character(len=80):: filename

    integer:: ilun
    character(len=80):: fileloc
    character(len=5):: nchar
    real(dp), allocatable, dimension(:):: xdp
    integer, allocatable, dimension(:):: xin
    integer, parameter:: tag = 1135
    integer:: dummy_io, info2

    integer:: nstellar_var, idim

    if(.not. stellar) return

    if(verbose) write(*,*) 'Entering backup_stellar'

    nstellar_var = ndim + 3 ! positions, mass, birth and life times

    ilun = 4*ncpu + myid + 11

    call title(myid, nchar)
    fileloc = TRIM(filename) // TRIM(nchar)

    ! Wait for the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE > 0) then
        if(mod(myid-1, IOGROUPSIZE) /= 0) then
            call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag, &
                & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
        end if
    end if
#endif

    open(unit=ilun, file=TRIM(fileloc), form='unformatted')
    rewind(ilun)

    write(ilun) nstellar_var
    write(ilun) nstellar

!added by PH 
!    write(ilun) nstellar_tot

    if(nstellar > 0) then
        allocate(xdp(1:nstellar))
        allocate(xin(1:nstellar))

        ! Write stellar object position
        do idim = 1, ndim
            xdp = xstellar(1:nstellar, idim)
            write(ilun) xdp
        end do

        ! Write stellar object mass
        xdp = mstellar(1:nstellar)
        write(ilun) xdp

        ! Write stellar object birth time
        xdp = tstellar(1:nstellar)
        write(ilun) xdp

        ! Write stellar object life time
        xdp = ltstellar(1:nstellar)
        write(ilun) xdp

        ! Write stellar object parent sink id
        xin = id_stellar(1:nstellar)
        write(ilun) xin
    end if

    close(ilun)

    ! Send the token
#ifndef WITHOUTMPI
    if(IOGROUPSIZE > 0) then
        if(mod(myid, IOGROUPSIZE) /= 0 .and. (myid < ncpu)) then
            dummy_io = 1
            call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid-1+1, tag, &
                & MPI_COMM_WORLD, info2)
        end if
    end if
#endif

end subroutine backup_stellar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_stellar_csv(filename)
  use amr_commons
  use pm_commons

  use feedback_module

  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,icpu,istellar

  if(verbose)write(*,*)'Entering output_stellar_csv'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write stellar properties
  !======================
  do istellar=1,nstellar

     write(ilun,'(I10,6(A1,ES20.10))')id_stellar(istellar),',',mstellar(istellar),&
          ',',xstellar(istellar,1),',',xstellar(istellar,2),',',xstellar(istellar,3),&
          ',',tstellar(istellar),',',ltstellar(istellar)
  end do

  close(ilun)

end subroutine output_stellar_csv
