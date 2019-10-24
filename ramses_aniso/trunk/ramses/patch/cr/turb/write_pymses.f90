SUBROUTINE write_pymses(output_dir)
   use amr_commons
   use hydro_parameters
   implicit none
   
   integer              :: ilun           ! I/O unit number
   integer              :: id             ! Data id counter
   integer              :: i              ! Loop counter
   integer              :: npassive       ! Number of passive quantities
   integer              :: num_list(1:20) ! Number list
   character(len=*), intent(in) :: output_dir ! Output directory
   character(len=50000) :: file_buffer    ! Buffer for header file
   character(len=1000)  :: item_buffer    ! Buffer for data item details
   character(len=3)     :: str_num        ! Number in string format
   character(len=255)   :: str_num_list   ! Number list in string format
   
   ilun = myid+10
   
   open(ilun,file=trim(output_dir)//'data_info.txt', status="unknown",&
        form="formatted")
  
   file_buffer = '{"hydro" : [ '
   
   id = 0
   write(str_num, '(I0)') id
   item_buffer = '("Scalar", "rho", '//trim(str_num)//')'
   file_buffer = trim(file_buffer) // trim(item_buffer)
   id = id + 1
   
   do i=1,ndim
      num_list(i) = id
      id = id + 1
   end do
   call list_of_numbers(ndim, num_list, str_num_list)
   item_buffer = '("Vector", "vel", '//trim(str_num_list)//')'
   file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   
   do i=1,ndim
      num_list(i) = id
      id = id + 1
   end do
   call list_of_numbers(ndim, num_list, str_num_list)
   item_buffer = '("Vector", "B left", '//trim(str_num_list)//')'
   file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   
   do i=1,ndim
      num_list(i) = id
      id = id + 1
   end do
   call list_of_numbers(ndim, num_list, str_num_list)
   item_buffer = '("Vector", "B right", '//trim(str_num_list)//')'
   file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   
   write(str_num, '(I0)') id
   item_buffer = '("Scalar", "P", '//trim(str_num)//')'
   file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   id = id + 1
   
!!$   if (nrad == 1) then
!!$      write(str_num, '(I0)') id
!!$      item_buffer = '("Scalar", "E_{rad}", '//trim(str_num)//')'
!!$      file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
!!$      id = id + 1
!!$   else if (nrad > 1) then
!!$      do i=1,nrad
!!$         num_list(i) = id
!!$         id = id + 1
!!$      end do
!!$      call list_of_numbers(nrad, num_list, str_num_list)
!!$      item_buffer = '("Vector", "E_{rad}", '//trim(str_num_list)//')'
!!$      file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
!!$   end if
   
   write(str_num, '(I0)') id
   item_buffer = '("Scalar", "E_{int}", '//trim(str_num)//')'
   file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   id = id + 1
   
   npassive = nvar - (9)!+nrad)
   if (npassive == 1) then
      write(str_num, '(I0)') id
      item_buffer = '("Scalar", "passive", '//trim(str_num)//')'
      file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
      id = id + 1
   else if (npassive > 1) then
      do i=1,npassive
         num_list(i) = id
         id = id + 1
      end do
      call list_of_numbers(npassive, num_list, str_num_list)
      item_buffer = '("Vector", "passive", '//trim(str_num_list)//')'
      file_buffer = trim(file_buffer) // ', ' // trim(item_buffer)
   end if
   
   file_buffer = trim(file_buffer) // ']'
   
! 0 - density
! 1, 2, 3 - velocity
! 4, 5, 6 - left B field
! 7, 8, 9 - right B field
! 10 - pressure
! 11... (9, 8+nrad) - radiative energy
! ... (9+nrad,nvar) - passive scalars
   
!    {
! ... "hydro" : [ Scalar("rho", 0), Vector("vel", [1, 2, 3]), Scalar("P", 4), Scalar("Z", 5), Scalar("HCO", 6) ],
! ... "grav"  : [ Vector("g", [0, 1, 2]) ]
! ... }
   
   if (poisson) then
      call list_of_numbers(3, (/(i,i=0,ndim-1)/), str_num_list)
      file_buffer = trim(file_buffer)//&
      & ', "grav" : [ ("Vector", "g", '//trim(str_num_list)//') ]'
   end if
   file_buffer = trim(file_buffer) // '}'
   
   write(ilun,'(A)') trim(file_buffer)
   close(ilun)
   
   return

   contains
   
   subroutine list_of_numbers(N, numbers, string_out)
      integer, intent(in)           :: N           ! Number of numbers
      integer, intent(in)           :: numbers(:)  ! Array of integers
      character(LEN=*), intent(out) :: string_out  ! String containing list
      
      integer                       :: i           ! Loop counter
      character(LEN=20)             :: str_num     ! Number as string
      
      string_out = '['
      
      if (N==0) then
         continue
      else if (N==1) then
         write(str_num, '(I0)') numbers(1)
         string_out = trim(string_out) // str_num
      else
         do i=1,N-1
            write(str_num, '(I0)') numbers(i)
            string_out = trim(string_out) // trim(str_num) // ', '
         end do
         write(str_num, '(I0)') numbers(N)
         string_out = trim(string_out) // trim(str_num)
      end if
      
      string_out = trim(string_out) // ']'
      return
   
   end subroutine list_of_numbers
  
END SUBROUTINE write_pymses
