program remap_rochedo
  ! use fpl !FPL module
  use netcdf
  use tools
  implicit none
  character*200 rootfolder,inputfolder,outputfolder,reffname,inpfname,codfname,outfname,outfpref
  integer ncid
  integer status, varid
  integer, TARGET :: inpnlat,inpnlon,refnlat,refnlon,ntim,npft,ncod
  integer,ALLOCATABLE, DIMENSION(:) :: codes
  character*200,ALLOCATABLE, DIMENSION(:) :: classes
  real*8, TARGET, ALLOCATABLE, DIMENSION(:,:) :: reflats,reflatn,reflonw,reflone,veglats,veglatn,veglonw,veglone

  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: vegdata
  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: outdata

  ! real*8, ALLOCATABLE, DIMENSION(:,:) :: inpdata
  integer, ALLOCATABLE, DIMENSION(:,:) :: inpdata

  real*8, pointer :: inplats(:,:),inplatn(:,:),inplonw(:,:),inplone(:,:)

  integer :: outnlon,outnlat

  real*8, TARGET, ALLOCATABLE :: outlats(:,:),outlatn(:,:),outlonw(:,:),outlone(:,:)
  real*8, POINTER :: vecoutlats(:),vecoutlatn(:),vecoutlonw(:),vecoutlone(:)

  real*8, ALLOCATABLE,DIMENSION(:) :: vecinplats(:),vecinplatn(:),vecinplonw(:),vecinplone(:)
  real*8, ALLOCATABLE,DIMENSION(:) :: years

  integer outi,outj,i,j,k,ii,jj, lasti,lastj, timestep

  integer latbnds(2),lonbnds(2)

  real*8 latsize,lonsize,latwgt,lonwgt,wgt,totwgt
  real*8, ALLOCATABLE, dimension(:) ::  outval ! The output value of a pixel in every PFT

  rootfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/"

  inputfolder = trim(ADJUSTL(rootfolder))//"input/"
  outputfolder = trim(ADJUSTL(rootfolder))//"out_rochedo/"

  ! reffname = trim(ADJUSTL(inputfolder))//"surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
  reffname = trim(ADJUSTL(inputfolder))//"min_surfdata.nc"
  inpfname = trim(ADJUSTL(inputfolder))//"weg_comp_2013_2015.nc"
  codfname = trim(ADJUSTL(rootfolder))//"codes_rochedo.csv"

  outfpref = trim(ADJUSTL(outputfolder))//"weg_remap_"

  ! Read the codes from the file
  ncod = count_lines(codfname)
  call read_codes(codfname,codes,classes,ncod)


  ! Get the sizes of the 4D reference file
  call get_ref_dimzises(reffname,"PCT_PFT",outnlat,outnlon,ntim,npft)
  write(*,*) outnlon,outnlat,npft,ntim

  ! Allocate the output variable
  allocate(outdata(outnlon,outnlat,ncod))

  ! Get the lat lon 2D bounds from the reference file
  call get_ref_grid(reffname,outnlat,outnlon,outlats,outlatn,outlonw,outlone)

! Get the sizes of the input landuse flie
  call get_3d_dimsizes(inpfname,"landuse",inpnlat,inpnlon,ntim)
  write(*,*) inpnlon,inpnlat,ntim
! Get the years in the time dimension
  call get_landuse_years(inpfname,years,ntim)

!
! Get the lat lon 2D bounds, assuming a regular grid
  call get_landuse_grid(inpfname,inpnlat,inpnlon,vecinplats,vecinplatn,vecinplonw,vecinplone)
! Flip the longitude variable, assuming its not crossing Greenwich
  call flip_lon_nocross_vec(inpnlon,vecinplonw,vecinplone)


! TODO: PUT TIME LOOP HERE!!!!
  timestep = 1
  ! Read a single year of data
  write(*,*) ">>>>>>>>>>>>>>> Reading timestep ",timestep
  call read_landuse_data_timestep(inpfname,inpdata,inpnlat,inpnlon,timestep)
  write(*,*) ">>>>>>>>>>>>>>> Finished reading timestep ",timestep

! call dum_write_2d("dummy.nc",outlats,outnlat,outnlon)

!   ! Use pointers to more generic names here (TODO: Refactor the code)
!   inpnlon => vegnlon
!   inpnlat => vegnlat
!
!   inplats => veglats
!   inplatn => veglatn
!   inplonw => veglonw
!   inplone => veglone
!
!   outnlon => refnlon
!   outnlat => refnlat
!
!   outlats => reflats
!   outlatn => reflatn
!   outlonw => reflonw
!   outlone => reflone
!
!   ! FIXME: For now, we'll simply assume the grid is regular
!   vecinplats => inplats(1,:)
!   vecinplatn => inplatn(1,:)
!   vecinplonw => inplonw(:,1)
!   vecinplone => inplone(:,1)
!
  vecoutlats => outlats(1,:)
  vecoutlatn => outlatn(1,:)
  vecoutlonw => outlonw(:,1)
  vecoutlone => outlone(:,1)

  ! Preallocate outdata with zeros
  outdata(:,:,:) = 0.0d0

  ! Allocate outval and fill it with zeros
  allocate(outval(ncod))
  outval(:) = 0.0d0



  ! Put output loop here
  ! outi = 67
  ! outj = 232
    do outi = 80,80
      do outj = 244,244
      ! do outi = 1,1
      !   do outj = 1,1

!   do outi = 1,outnlat
    ! write(*,*) "Running latitude ",outi," of ",outnlat
!     do outj = 1,outnlon
!
      write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< outi,outj >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) outi,outj
      write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      ! Out of bounds, continue the loop
      if (vecoutlats(outi).ge.vecinplatn(inpnlat) .or. vecoutlatn(outi).le.vecinplats(1) .or. vecoutlonw(outj).ge.vecinplone(inpnlon) .or. vecoutlone(outj).le.vecinplonw(1)) then
        cycle
      end if


      latbnds = find_bound_inds_vec(vecoutlats(outi),vecoutlatn(outi),vecinplats,vecinplatn)
      lonbnds = find_bound_inds_vec(vecoutlonw(outj),vecoutlone(outj),vecinplonw,vecinplone)
      write(*,*) "latbnds = ",latbnds
      write(*,*) "lonbnds = ",lonbnds

      write(*,*) "======================== vecout =========================="
      write(*,*) vecoutlats(outi),vecoutlatn(outi),vecoutlonw(outj),vecoutlone(outj)
      write(*,*) "=========================================================="

      outval(:) = 0.0d0
      totwgt = 0.0d0
      do i = latbnds(1),latbnds(2)
        do j = lonbnds(1),lonbnds(2)

          ! Calculate lat weight for that pixel
          if (is_contained_vec(vecoutlats(outi),vecoutlatn(outi),vecinplats(i),vecinplatn(i))) then
            latwgt = 1.0d0
          else

            latsize = vecinplatn(i)-vecinplats(i)
            if (vecinplatn(i) .ge. vecoutlats(outi) .and. vecinplatn(i) .le. vecoutlatn(outi)) then
              latwgt = (vecinplatn(i)-vecoutlats(outi))/latsize
            else
              latwgt = (vecoutlatn(outi) - vecinplats(i))/latsize
            end if
          end if

          ! Calculate lon weight for that pixel
          if (is_contained_vec(vecoutlonw(outj),vecoutlone(outj),vecinplonw(j),vecinplone(j))) then

            lonwgt = 1.0d0
          else
            lonsize = vecinplone(j)-vecinplonw(j)
            if (vecinplone(j) .ge. vecoutlonw(outj) .and. vecinplone(j) .le. vecoutlone(outj)) then
              lonwgt = (vecinplone(j)-vecoutlonw(outj))/lonsize
            else
              lonwgt = (vecoutlone(outj) - vecinplonw(j))/lonsize
            end if
          end if

          wgt = lonwgt*latwgt

          ! ! Check if it's inside the mask
          ! ! if (vegmask(j,i).ne.1) wgt = 0.0d0
          !
          totwgt = totwgt + wgt

          ! Accumulate val*wgt in outval
          do k = 1,ncod
            if (inpdata(j,i).eq.codes(k)) then
              outval(k) = outval(k) + wgt
            end if
          end do !k, cods

          write(*,*) i,j,wgt,inpdata(j,i)

        end do ! j, lonbnds
      end do !i, latbnds

      do k = 1,ncod
          write(*,*) codes(k),outval(k),100.0d0*outval(k)/totwgt
      end do

      ! Now divide by the weights
      outval = outval/totwgt

      ! Write to the variable
      outdata(outj,outi,:) = outval(:)

    end do !outj, outnlat
  end do !outi, outnlon
!

  ! Write dataset. The division by zero (wgt=0 when mask=0) leads to missing values in the output
  ! call write_pft_data(outfname,outdata,outnlat,outnlon,npft,outlats,outlatn,outlonw,outlone) ! Checking


  ! write(*,*) reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2)
  ! write(*,*) is_inside_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))
  ! write(*,*) is_contained_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))


end program remap_rochedo
