module cloud_module
  ! EMPTY MODULE TO ALLOW COMPATABILITY WITH NON-CLOUD VERSIONS
  ! SEE cloudinit/ FOR FULL CLOUD INITIAL CONDITIONS MODULE

  real(dp):: switch_solv=1.d20

end module cloud_module

subroutine read_cloud_params(nml_ok)

  use cloud_module
  implicit none
  logical::nml_ok

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cloud_params/switch_solv

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist


end subroutine read_cloud_params

