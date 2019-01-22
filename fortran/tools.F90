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

    status = nf90_open(trim(ADJUSTL(fname)), NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,trim(ADJUSTL(varname)),varid)

    status = nf90_inquire_variable(ncid,varid,dimids = dimids)

    status = nf90_inquire_dimension(ncid,dimids(1),len = nlon)
    status = nf90_inquire_dimension(ncid,dimids(2),len = nlat)
    status = nf90_inquire_dimension(ncid,dimids(3),len = npft)
    status = nf90_inquire_dimension(ncid,dimids(4),len = ntim)

    status = nf90_close(ncid)

    return
  end subroutine get_ref_dimzises

  ! Get the 2D lat lon bounds of the 4D reference variable
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

    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"LATS",varid)
    status = nf90_get_var(ncid,varid,lats)

    status = nf90_inq_varid(ncid,"LATN",varid)
    status = nf90_get_var(ncid,varid,latn)

    status = nf90_inq_varid(ncid,"LONW",varid)
    status = nf90_get_var(ncid,varid,lonw)

    status = nf90_inq_varid(ncid,"LONE",varid)
    status = nf90_get_var(ncid,varid,lone)

    status = nf90_close(ncid)

    return
  end subroutine get_ref_grid
! Get the dimension sizes of a 2D file
  subroutine get_2d_dimsizes(fname,varname,nlat,nlon)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname,varname
    integer ncid
    integer status, varid
    integer nlat,nlon
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids

    status = nf90_open(trim(ADJUSTL(fname)), NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,trim(ADJUSTL(varname)),varid)

    status = nf90_inquire_variable(ncid,varid,dimids = dimids)

    status = nf90_inquire_dimension(ncid,dimids(1),len = nlon)
    status = nf90_inquire_dimension(ncid,dimids(2),len = nlat)

    status = nf90_close(ncid)

    return
  end subroutine get_2d_dimsizes

! Get the 2D lat lon bounds of a 2D file, that has 2D centerpoint variables, assuming a regular grid
  subroutine get_veg_grid(fname,nlat,nlon,lats,latn,lonw,lone)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,ntim,npft
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    real*8, ALLOCATABLE, DIMENSION(:,:) :: latixy,longxy,lats,latn,lonw,lone

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ALLOCATE(lats(nlon,nlat),latn(nlon,nlat),lonw(nlon,nlat),lone(nlon,nlat))
    ALLOCATE(latixy(nlon,nlat),longxy(nlon,nlat))

    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"LATIXY",varid)
    status = nf90_get_var(ncid,varid,latixy)

    status = nf90_inq_varid(ncid,"LONGXY",varid)
    status = nf90_get_var(ncid,varid,longxy)

    status = nf90_close(ncid)

    ! Assuming a regular grid
    reslon = longxy(2,1) - longxy(1,1)
    reslat = latixy(1,2) - latixy(1,1)

    lonw = longxy - (reslon/2.0d0)
    lone = longxy + (reslon/2.0d0)
    lats = latixy - (reslat/2.0d0)
    latn = latixy + (reslat/2.0d0)

    return
  end subroutine get_veg_grid

end module
