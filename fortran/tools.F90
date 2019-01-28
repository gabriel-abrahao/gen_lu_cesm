module tools
contains

  ! Count the lines in a text file
  function count_lines(fname)
    implicit none
    character(*), INTENT(IN) ::fname
    integer count_lines
    integer status
    open(999,file = trim(ADJUSTL(fname)))
    count_lines = -1
    do while (status.ge.0)
      count_lines = count_lines+1
      read(999,*,iostat = status)
    end do
    if (status.gt.0) then
      write(*,*) "Error in count_lines(",fname,")"
      write(*,*) "status: ",status
    end if
    close(999)
    return
  end function count_lines

  subroutine read_codes(fname,codes,classes,ncod)
    use iso_fortran_env
    implicit none
    character(*), INTENT(IN) ::fname
    integer ncod
    INTEGER, ALLOCATABLE, DIMENSION(:) :: codes
    CHARACTER*200, ALLOCATABLE, DIMENSION(:) :: classes
    integer i

    allocate(codes(ncod),classes(ncod))

    open(999,file = fname)

    do i = 1,ncod
      read(999,*) codes(i),classes(i)
    end do

    close(999)
    return
  end subroutine read_codes

  function ispan(start,end,stride)
    implicit none
    integer start,end,stride,n,i,count
    INTEGER, ALLOCATABLE,DIMENSION(:) :: ispan

    n = CEILING(real(end-start+1)/real(stride))

    allocate(ispan(n))

    count = 0
    do i = start,end,stride
      count = count + 1
      ispan(count) = i
    end do

    return
  end function ispan

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

  ! Get the dimension sizes of a 3D file
  subroutine get_3d_dimsizes(fname,varname,nlat,nlon,ntim)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname,varname
    integer ncid
    integer status, varid
    integer nlat,nlon,ntim
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids

    status = nf90_open(trim(ADJUSTL(fname)), NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,trim(ADJUSTL(varname)),varid)

    status = nf90_inquire_variable(ncid,varid,dimids = dimids)

    status = nf90_inquire_dimension(ncid,dimids(1),len = nlon)
    status = nf90_inquire_dimension(ncid,dimids(2),len = nlat)
    status = nf90_inquire_dimension(ncid,dimids(3),len = ntim)

    status = nf90_close(ncid)

    return
  end subroutine get_3d_dimsizes

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

  ! Get the 2D lat lon bounds of a 2D file, that has 2D centerpoint variables, assuming a regular grid
  subroutine get_landuse_grid(fname,nlat,nlon,lats,latn,lonw,lone)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,ntim,npft
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    real*8, ALLOCATABLE, DIMENSION(:) :: latc,lonc,lats,latn,lonw,lone

    INTEGER i
    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ALLOCATE(lats(nlat),latn(nlat),lonw(nlon),lone(nlon))
    ALLOCATE(latc(nlat),lonc(nlon))

    ncid = 999
    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"lat",varid)
    status = nf90_get_var(ncid,varid,latc)

    status = nf90_inq_varid(ncid,"lon",varid)
    status = nf90_get_var(ncid,varid,lonc)

    status = nf90_close(ncid)

    ! Assuming a regular grid
    reslon = lonc(2) - lonc(1)
    reslat = latc(2) - latc(1)

    lonw = lonc - (reslon/2.0d0)
    lone = lonc + (reslon/2.0d0)
    lats = latc - (reslat/2.0d0)
    latn = latc + (reslat/2.0d0)

    return
  end subroutine get_landuse_grid

  subroutine get_landuse_years(fname,years,ntim)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,ntim,npft
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    real*8, ALLOCATABLE, DIMENSION(:) :: years

    INTEGER i

    allocate(years(ntim))

    ncid = 999
    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"time",varid)
    status = nf90_get_var(ncid,varid,years)

    write(*,*) years

    status = nf90_close(ncid)
    return
  end subroutine get_landuse_years

  ! Read the potential vegetation data (lon,lat,pft)
  subroutine read_landuse_data_timestep(fname,data,nlat,nlon,timestep)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft,ntim,timestep
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    ! real*8, ALLOCATABLE, DIMENSION(:,:) :: data
    integer, ALLOCATABLE, DIMENSION(:,:) :: data

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ALLOCATE(data(nlon,nlat))

    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"landuse",varid)
    status = nf90_get_var(ncid,varid,data,(/1,1,timestep/),(/nlon,nlat,1/))

    status = nf90_close(ncid)

    return
  end subroutine read_landuse_data_timestep

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

  subroutine read_veg_mask(fname,nlat,nlon,data)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft
    integer, dimension(nf90_max_var_dims) :: dimids ! To store the dimension ids
    real*8, ALLOCATABLE, DIMENSION(:,:) :: data

    ALLOCATE(data(nlon,nlat))

    status = nf90_open(fname, NF90_NOWRITE, ncid)

    status = nf90_inq_varid(ncid,"LANDMASK",varid)
    status = nf90_get_var(ncid,varid,data)

    status = nf90_close(ncid)

    return
  end subroutine read_veg_mask

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

    ! write(*,*) lats(1,1),latn(1,1)
    ! write(*,*) lats(nlon,1),latn(nlon,1)
    ! write(*,*) lats(1,nlat),latn(1,nlat)
    ! write(*,*) lats(nlon,nlat),latn(nlon,nlat)
    !
    ! write(*,*) lonw(1,1),lone(1,1)
    ! write(*,*) lonw(nlon,1),lone(nlon,1)
    ! write(*,*) lonw(1,nlat),lone(1,nlat)
    ! write(*,*) lonw(nlon,nlat),lone(nlon,nlat)



    return
  end subroutine flip_lon_global_3d

  ! Flip the longitude variable (vecto longitude, NO CROSSING GREENWICH!!)
  subroutine flip_lon_nocross_vec(nlon,lonw,lone)
    use netcdf
    implicit none
    integer nlat,nlon,npft
    integer midl,midr ! Middle (left and right) longitude indices
    real*8, DIMENSION(:) :: lonw,lone

    lonw = lonw + 360.0d0
    lone = lone + 360.0d0

    return
  end subroutine flip_lon_nocross_vec

  ! Flip the longitude variables and data (global 2D)
  subroutine flip_lon_global_2d_nometa(data,nlat,nlon)
    use netcdf
    implicit none
    integer nlat,nlon
    integer midl,midr ! Middle (left and right) longitude indices
    real*8, DIMENSION(:,:) :: data
    real*8, ALLOCATABLE, DIMENSION(:,:) :: dum2d

    ! We'll allocate in turns to save memory
    ALLOCATE(dum2d(nlon,nlat))

    ! We are assuming a regular longitude grid here, with a nonodd number of points
    midl = nlon/2
    midr = midl+1

    dum2d(1:midl,:) = data(midr:nlon,:)
    dum2d(midr:nlon,:) = data(1:midl,:)
    data = dum2d

    ! Deallocating 2D to save memory
    deallocate(dum2d)

    return
  end subroutine flip_lon_global_2d_nometa

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

    integer, dimension(2) :: bnds,find_bound_inds_vec

    integer n
    integer i, maxi, mini, fnd
    integer iter,maxiter

    n = size(litlovec)

    ! The out-of-range case
    if (biglo.ge.lithivec(n) .or. bighi.le.litlovec(1)) then
      write(*,*) "find_bound_inds_vec: OUT OF RANGE, you should never have called me"
      write(*,*) "biglo,bighi,litlovec(1),lithivec(n)"
      write(*,*) biglo,bighi,litlovec(1),lithivec(n)
      stop
    end if

    maxiter = int(log(real(n))) * 2 + 3

    ! ! write(*,*) "maxiter = ",maxiter
    ! write(*,*) litlovec(n-5:n)
    ! write(*,*) "biglo,bighi = ",biglo,bighi

    ! write(*,*) litlovec(179:181)
    ! stop
    i = n/2
    maxi = n
    mini = 1
    iter = 1
    ! write(*,*) i
    do while (.not.(is_inside_vec(biglo,bighi,litlovec(i),lithivec(i))))
      iter = iter + 1
      !maxi = max(n,i,maxi)
      !mini = min(1,i,mini)
      ! write(*,*) i,mini,maxi, iter
      ! write(*,*) litlovec(i),biglo
      if (litlovec(i).ge.biglo) then
        maxi = i
        i = i - ((i-mini)/2)
        if ((maxi-mini).eq.1) then ! We are by one
          i = mini
        end if
        ! write(*,*) "PASSEI ESQ"
      else
        mini = i
        i = i + ((maxi-i)/2)
        if ((maxi-mini).eq.1) then ! We are by one
          i = maxi
        end if
        ! write(*,*) "PASSEI DIR"
      end if
      ! write(*,*) i
      if (iter.ge.maxiter) then
        write(*,*) "ERROR: find_bound_inds_vec did not find any matches"
        stop
      end if
      !exit
    end do

    fnd = i

    ! write(*,*) "FOUND ",fnd
    ! write(*,*) litlovec(i),lithivec(i)

    bnds(:) = fnd

    i = bnds(1)
    do while (is_inside_vec(biglo,bighi,litlovec(i),lithivec(i)))
      bnds(1) = i
      i = i-1
      if (i.lt.1) exit
    end do
    i = bnds(2)
    do while (is_inside_vec(biglo,bighi,litlovec(i),lithivec(i)) .or. i.lt.1 .or. i.gt.n)
      bnds(2) = i
      i = i+1
      if (i.gt.n) exit
    end do

    ! write(*,*) "biglo,bighi = ",biglo,bighi
    ! write(*,*) "bnds = ",bnds
    ! write(*,*) "Pixels inside:"
    ! do i = bnds(1),bnds(2)
    !   write(*,*) litlovec(i),lithivec(i)
    ! end do

    find_bound_inds_vec = bnds
    return
  end function find_bound_inds_vec

subroutine write_pft_data(fname,data,nlat,nlon,npft,lats,latn,lonw,lone)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,npft
    real*8, ALLOCATABLE, DIMENSION(:) :: lat,lon
    real*8, DIMENSION(:,:) :: lats,latn,lonw,lone
    ! real*8, DIMENSION(:,:) :: mask
    real*8, DIMENSION(:,:,:) :: data
    integer i
    integer, dimension(3) :: dimids,coordids ! To store the dimension ids
    integer latsid,latnid,lonwid,loneid

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ncid = 9876

    ! Define ids explicitly
    do i = 1,size(dimids)
      dimids(i) = i
    end do
    do i = 1,size(coordids)
      coordids(i) = i + dimids(size(dimids))
    end do
    varid = coordids(size(coordids))+1
    latsid = varid + 1
    latnid = latsid + 1
    lonwid = latnid + 1
    loneid = lonwid + 1

    ! Get 1D versions of coordvars for ease of plotting later
    allocate(lat(nlat),lon(nlon))
    lat = (lats(1,:)+latn(1,:))/2.0d0
    lon = (lonw(:,1)+lone(:,1))/2.0d0

    !NF90_CLOBBER allows overwriting
    status = nf90_create(fname, NF90_CLOBBER, ncid)

    status = nf90_def_dim(ncid,"lon",nlon,dimids(1))
    status = nf90_def_dim(ncid,"lat",nlat,dimids(2))
    status = nf90_def_dim(ncid,"pft",npft,dimids(3))

    status = nf90_def_var(ncid,"lon",nf90_double,dimids(1),coordids(1))
    status = nf90_def_var(ncid,"lat",nf90_double,dimids(2),coordids(2))
    status = nf90_put_att(ncid, coordids(1), "long_name", "lon")
    status = nf90_put_att(ncid, coordids(1), "units", "degrees east")
    status = nf90_put_att(ncid, coordids(2), "long_name", "lat")
    status = nf90_put_att(ncid, coordids(2), "units", "degrees north")

    !write(*,*) dimids
    status = nf90_def_var(ncid,"PCT_PFT",nf90_double,dimids(1:3),varid)
    status = nf90_put_att(ncid, varid, "long_name", "percent_pft")
    status = nf90_put_att(ncid, varid, "units", "unitless")

    status = nf90_def_var(ncid,"LATS",nf90_double,dimids(1:2),latsid)
    status = nf90_def_var(ncid,"LATN",nf90_double,dimids(1:2),latnid)
    status = nf90_def_var(ncid,"LONW",nf90_double,dimids(1:2),lonwid)
    status = nf90_def_var(ncid,"LONE",nf90_double,dimids(1:2),loneid)
    status = nf90_put_att(ncid, latsid, "long_name", "latitude of southern edge")
    status = nf90_put_att(ncid, latnid, "long_name", "latitude of northern edge")
    status = nf90_put_att(ncid, lonwid, "long_name", "longitude of western edge")
    status = nf90_put_att(ncid, loneid, "long_name", "longitude of eastern edge")
    status = nf90_put_att(ncid, latsid, "units", "degrees north")
    status = nf90_put_att(ncid, latnid, "units", "degrees north")
    status = nf90_put_att(ncid, lonwid, "units", "degrees east")
    status = nf90_put_att(ncid, loneid, "units", "degrees east")

    ! Exit define mode
    status = nf90_enddef(ncid)

    ! Put the variables
    status = nf90_put_var(ncid,varid,data,(/1,1,1/),(/nlon,nlat,npft/))

    status = nf90_put_var(ncid,latsid,lats,(/1,1/),(/nlon,nlat/))

    ! Put coordinate variables
    status = nf90_put_var(ncid,coordids(1),lon)
    status = nf90_put_var(ncid,coordids(2),lat)



    ! Close the file
    status = nf90_close(ncid)

    return
end subroutine write_pft_data

subroutine write_rochedo_data(fname,data,nlat,nlon,ncod,ntim,lats,latn,lonw,lone)
    use netcdf
    implicit none
    character(*), INTENT(IN) :: fname
    integer ncid
    integer status, varid
    integer nlat,nlon,ncod,ntim
    real*8, ALLOCATABLE, DIMENSION(:) :: lat,lon
    real*8, DIMENSION(:,:) :: lats,latn,lonw,lone
    ! real*8, DIMENSION(:,:) :: mask
    real*8, DIMENSION(:,:,:,:) :: data
    integer i
    integer, dimension(4) :: dimids,coordids ! To store the dimension ids
    integer latsid,latnid,lonwid,loneid

    real*8 reslat,reslon ! The lat and lon resolutions of a regular grid

    ncid = 9876

    ! Define ids explicitly
    do i = 1,size(dimids)
      dimids(i) = i
    end do
    do i = 1,size(coordids)
      coordids(i) = i + dimids(size(dimids))
    end do
    varid = coordids(size(coordids))+1
    latsid = varid + 1
    latnid = latsid + 1
    lonwid = latnid + 1
    loneid = lonwid + 1

    ! Get 1D versions of coordvars for ease of plotting later
    allocate(lat(nlat),lon(nlon))
    lat = (lats(1,:)+latn(1,:))/2.0d0
    lon = (lonw(:,1)+lone(:,1))/2.0d0

    !NF90_CLOBBER allows overwriting
    status = nf90_create(fname, NF90_CLOBBER, ncid)

    status = nf90_def_dim(ncid,"lon",nlon,dimids(1))
    status = nf90_def_dim(ncid,"lat",nlat,dimids(2))
    status = nf90_def_dim(ncid,"pft",ncod,dimids(3))
    status = nf90_def_dim(ncid,"time",ntim,dimids(4))

    status = nf90_def_var(ncid,"lon",nf90_double,dimids(1),coordids(1))
    status = nf90_def_var(ncid,"lat",nf90_double,dimids(2),coordids(2))
    status = nf90_put_att(ncid, coordids(1), "long_name", "lon")
    status = nf90_put_att(ncid, coordids(1), "units", "degrees east")
    status = nf90_put_att(ncid, coordids(2), "long_name", "lat")
    status = nf90_put_att(ncid, coordids(2), "units", "degrees north")

    !write(*,*) dimids
    status = nf90_def_var(ncid,"PCT_PFT",nf90_double,dimids(1:4),varid)
    status = nf90_put_att(ncid, varid, "long_name", "percent_pft")
    status = nf90_put_att(ncid, varid, "units", "unitless")

    status = nf90_def_var(ncid,"LATS",nf90_double,dimids(1:2),latsid)
    status = nf90_def_var(ncid,"LATN",nf90_double,dimids(1:2),latnid)
    status = nf90_def_var(ncid,"LONW",nf90_double,dimids(1:2),lonwid)
    status = nf90_def_var(ncid,"LONE",nf90_double,dimids(1:2),loneid)
    status = nf90_put_att(ncid, latsid, "long_name", "latitude of southern edge")
    status = nf90_put_att(ncid, latnid, "long_name", "latitude of northern edge")
    status = nf90_put_att(ncid, lonwid, "long_name", "longitude of western edge")
    status = nf90_put_att(ncid, loneid, "long_name", "longitude of eastern edge")
    status = nf90_put_att(ncid, latsid, "units", "degrees north")
    status = nf90_put_att(ncid, latnid, "units", "degrees north")
    status = nf90_put_att(ncid, lonwid, "units", "degrees east")
    status = nf90_put_att(ncid, loneid, "units", "degrees east")

    ! Exit define mode
    status = nf90_enddef(ncid)

    ! Put the variables
    status = nf90_put_var(ncid,varid,data,(/1,1,1,1/),(/nlon,nlat,ncod,ntim/))
    write(*,*) status
    write(*,*) shape(data)
    write(*,*) (/nlon,nlat,ncod,ntim/)

    status = nf90_put_var(ncid,latsid,lats,(/1,1/),(/nlon,nlat/))
    status = nf90_put_var(ncid,latnid,latn,(/1,1/),(/nlon,nlat/))
    status = nf90_put_var(ncid,lonwid,lonw,(/1,1/),(/nlon,nlat/))
    status = nf90_put_var(ncid,loneid,lone,(/1,1/),(/nlon,nlat/))

    ! Put coordinate variables
    status = nf90_put_var(ncid,coordids(1),lon)
    status = nf90_put_var(ncid,coordids(2),lat)



    ! Close the file
    status = nf90_close(ncid)

    return
end subroutine write_rochedo_data

end module
