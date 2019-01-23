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

  ! Read the potential vegetation data (lon,lat,pft)
  subroutine read_veg_data(fname,nlat,nlon,npft,data)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    real*8, ALLOCATABLE, DIMENSION(:,:,:) :: data

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ALLOCATE(data(nlon,nlat,npft))

    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"PCT_PFT",varid)
    status = nf90_get_var(ncid,varid,data)

    status = nf90_close(ncid)

    return
  end subroutine read_veg_data

  subroutine dum_write_3d(fname,data,nlat,nlon,npft)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft
    integer i
    integer, dimension(3) :: dimids ! To store the dimension ids
    real*8, DIMENSION(:,:,:), INTENT(IN) :: data

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ncid = 987

    do i = 1,size(dimids)
      dimids(i) = i
    end do
    varid = i+1

    status = nf90_create(fname, NF90_CLOBBER, ncid)

    status = nf90_def_dim(ncid,"lon",nlon,dimids(1))
    status = nf90_def_dim(ncid,"lat",nlat,dimids(2))
    status = nf90_def_dim(ncid,"pft",npft,dimids(3))
    !write(*,*) dimids
    status = nf90_def_var(ncid,"dummy",nf90_double,dimids(1:3),varid)

    status = nf90_enddef(ncid)

    status = nf90_put_var(ncid,varid,data,(/1,1,1/),(/nlon,nlat,npft/))
    write(*,*) status

    ! status = nf90_inq_varid(ncid,"PCT_PFT",varid)
    ! status = nf90_get_var(ncid,varid,data)
    !
    status = nf90_close(ncid)

    return
  end subroutine dum_write_3d

  subroutine dum_write_2d(fname,data,nlat,nlon)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft
    integer i
    integer, dimension(2) :: dimids ! To store the dimension ids
    real*8, DIMENSION(:,:), INTENT(IN) :: data

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ncid = 987

    do i = 1,size(dimids)
      dimids(i) = i
    end do
    varid = i+1

    status = nf90_create(fname, NF90_CLOBBER, ncid)

    status = nf90_def_dim(ncid,"lon",nlon,dimids(1))
    status = nf90_def_dim(ncid,"lat",nlat,dimids(2))
    !status = nf90_def_dim(ncid,"pft",npft,dimids(3))
    !write(*,*) dimids
    status = nf90_def_var(ncid,"dummy",nf90_double,dimids,varid)

    status = nf90_enddef(ncid)

    status = nf90_put_var(ncid,varid,data,(/1,1/),(/nlon,nlat/))
    write(*,*) status

    ! status = nf90_inq_varid(ncid,"PCT_PFT",varid)
    ! status = nf90_get_var(ncid,varid,data)
    !
    status = nf90_close(ncid)

    return
  end subroutine dum_write_2d

  ! Flip the longitude variables and data (global 3D)
  subroutine flip_lon_global_3d(data,nlat,nlon,npft,lats,latn,lonw,lone)
    use netcdf
    implicit none
    integer nlat,nlon,npft
    integer midl,midr ! Middle (left and right) longitude indices
    real*8, DIMENSION(:,:) :: lats,latn,lonw,lone
    real*8, DIMENSION(:,:,:) :: data
    real*8, ALLOCATABLE, DIMENSION(:,:,:) :: dum3d
    real*8, ALLOCATABLE, DIMENSION(:,:) :: dum2d

    ! We'll allocate in turns to save memory
    ALLOCATE(dum2d(nlon,nlat))

    ! We are assuming a regular longitude grid here, with a nonodd number of points
    midl = nlon/2
    midr = midl+1

    dum2d(1:midl,:) = lonw(midr:nlon,:)
    dum2d(midr:nlon,:) = lonw(1:midl,:) + 360.0d0
    lonw = dum2d

    dum2d(1:midl,:) = lone(midr:nlon,:)
    dum2d(midr:nlon,:) = lone(1:midl,:) + 360.0d0
    lone = dum2d

    ! Deallocating 2D to save memory
    deallocate(dum2d)
    ALLOCATE(dum3d(nlon,nlat,npft))

    dum3d(1:midl,:,:) = data(midr:nlon,:,:)
    dum3d(midr:nlon,:,:) = data(1:midl,:,:)
    data = dum3d

    return
  end subroutine flip_lon_global_3d

! Checks if a little interval is at least partially inside a big interval
  function is_inside_vec(biglo,bighi,litlo,lithi)
    real*8 biglo,bighi,litlo,lithi
    logical is_inside_vec

    if ( (litlo .ge. biglo .and. litlo .le. bighi) .or. (lithi .ge. biglo .and. lithi .le. bighi)) then
      is_inside_vec = .true.
    else
      is_inside_vec = .false.
    end if

    return
  end function is_inside_vec

! Checks if a little interval is completely inside a big interval
  function is_contained_vec(biglo,bighi,litlo,lithi)
    real*8 biglo,bighi,litlo,lithi
    logical is_contained_vec

    if (.not.is_inside_vec(biglo,bighi,litlo,lithi)) then
      is_contained_vec = .false.
      return
    else
      if ( litlo .ge. biglo .and. lithi .le. bighi ) then
        is_contained_vec = .true.
      else
        is_contained_vec = .false.
      end if
    end if

    return
  end function is_contained_vec


    function find_bound_inds_vec(biglo,bighi,litlovec,lithivec)
      implicit none
      real*8 biglo,bighi
      real*8, dimension(:) :: litlovec,lithivec

      integer, dimension(2) :: find_bound_inds_vec

      integer n
      integer i

      n = size(litlovec)

      i = n/2
      write(*,*) i
      ! do while (.not.(is_inside_vec(biglo,bighi,litlovec(i),lithivec(i))))
      ! end do

      return
    end function find_bound_inds_vec

end module
