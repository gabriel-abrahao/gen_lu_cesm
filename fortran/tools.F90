module tools
contains

! Get the dimension sizes of the 4D reference variable
subroutine get_ref_dimzises(fname,varname,nlat,nlon,ntim,npft)
  use netcdf
  implicit none
  character(*), INTENT(IN) :: fname,varname
  integer ncid
  integer status, varid
  integer nlat,nlon,ntim,npft
  integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids

  write(*,*) "PASSEI 1"
  write(*,*) "varname=",trim(ADJUSTL(varname))

  status = nf90_open(trim(ADJUSTL(fname)), NF90_NOWRITE, ncid)
  status = nf90_inq_varid(ncid,trim(ADJUSTL(varname)),varid)
  status = nf90_inquire_variable(ncid,varid,dimids = dimids)
  status = nf90_inquire_dimension(ncid,dimids(1),len = nlon)
  status = nf90_inquire_dimension(ncid,dimids(2),len = nlat)
  status = nf90_inquire_dimension(ncid,dimids(3),len = npft)
  status = nf90_inquire_dimension(ncid,dimids(4),len = ntim)
  status = nf90_close(ncid)
  write(*,*) "LALALA", fname,varname

  return
end subroutine get_ref_dimzises

subroutine get_ref_grid(fname,nlat,nlon,lats,latn,lonw,lone)
  use netcdf
  implicit none
  character(*), INTENT(IN) :: fname
  integer ncid
  integer status, varid
  integer nlat,nlon,ntim,npft
  integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
  real*8, ALLOCATABLE, DIMENSION(:,:) :: lats,latn,lonw,lone

  ALLOCATE(lats(nlon,nlat),latn(nlon,nlat),lonw(nlon,nlat),lone(nlon,nlat))

  ! status = nf90_open(fname, NF90_NOWRITE, ncid)
  ! status = nf90_inq_varid(ncid,"PCT_PFT",varid)
  ! status = nf90_inquire_variable(ncid,varid,dimids = dimids)
  ! status = nf90_inquire_dimension(ncid,dimids(1),len = nlon)
  ! status = nf90_inquire_dimension(ncid,dimids(2),len = nlat)
  ! status = nf90_inquire_dimension(ncid,dimids(3),len = npft)
  ! status = nf90_inquire_dimension(ncid,dimids(4),len = ntim)
  ! status = nf90_close(ncid)

  return
end subroutine get_ref_grid

end module
