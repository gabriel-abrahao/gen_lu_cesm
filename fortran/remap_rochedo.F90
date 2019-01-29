program remap_rochedo
  ! use fpl !FPL module
  use netcdf
  use tools
  implicit none
  character*200 rootfolder,inputfolder,outputfolder,reffname,inpfname,codfname,outfname
  character*1000 dumchar
  integer ncid
  integer status, varid
  integer, TARGET :: inpnlat,inpnlon,refnlat,refnlon,ntim,npft,ncod,refntim
  integer,ALLOCATABLE, DIMENSION(:) :: codes
  character*200,ALLOCATABLE, DIMENSION(:) :: classes

  real*8, TARGET, ALLOCATABLE, DIMENSION(:,:) :: reflats,reflatn,reflonw,reflone,veglats,veglatn,veglonw,veglone

  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: vegdata
  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: outdata
  real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: fulloutdata

  ! real*8, ALLOCATABLE, DIMENSION(:,:) :: inpdata
  integer, ALLOCATABLE, DIMENSION(:,:) :: inpdata

  real*8, pointer :: inplats(:,:),inplatn(:,:),inplonw(:,:),inplone(:,:)

  integer :: outnlon,outnlat

  real*8, TARGET, ALLOCATABLE :: outlats(:,:),outlatn(:,:),outlonw(:,:),outlone(:,:)
  real*8, POINTER :: vecoutlats(:),vecoutlatn(:),vecoutlonw(:),vecoutlone(:)

  real*8, ALLOCATABLE,DIMENSION(:) :: vecinplats(:),vecinplatn(:),vecinplonw(:),vecinplone(:)
  real*8, ALLOCATABLE,DIMENSION(:) :: years
  integer, ALLOCATABLE,DIMENSION(:) :: intyears
  integer syear
  integer misscod,missind

  integer outi,outj,i,j,k,ii,jj, lasti,lastj, timestep

  integer latbnds(2),lonbnds(2)

  real*8 latsize,lonsize,latwgt,lonwgt,wgt,totwgt
  real*8, ALLOCATABLE, dimension(:) ::  outval ! The output value of a pixel in every PFT

  rootfolder = "../"

  ! Folders relative to rootfolder
  inputfolder = "input/"
  outputfolder = "out_rochedo/"

  ! Input file names (relative to input folder)
  reffname = "min_surfdata.nc"
  inpfname = "weg_comp_2013_2015.nc"

  codfname = trim(ADJUSTL(rootfolder))//"codes_rochedo.csv"

  misscod = 0 ! Code that refers to missing cells in Rochedo
  syear = 2013 ! Used to recreate the time coordvar

  outfname = "fraction_"//trim(ADJUSTL(inpfname))
! ################################ END OF INPUTS ###############################

  ! Resolving full paths
  inputfolder = trim(ADJUSTL(rootfolder))//inputfolder
  outputfolder = trim(ADJUSTL(rootfolder))//outputfolder
  reffname = trim(ADJUSTL(inputfolder))//reffname
  inpfname = trim(ADJUSTL(inputfolder))//inpfname
  outfname = trim(ADJUSTL(outputfolder))//outfname

  write(*,*) "rootfolder  : ", trim(ADJUSTL(rootfolder))
  write(*,*) "inputfolder  : ", trim(ADJUSTL(inputfolder))
  write(*,*) "outputfolder  : ", trim(ADJUSTL(outputfolder))
  write(*,*) "reffname  : ", trim(ADJUSTL(reffname))
  write(*,*) "inpfname  : ", trim(ADJUSTL(inpfname))
  write(*,*) "outfname  : ", trim(ADJUSTL(outfname))


  ! Read the codes from the file
  ncod = count_lines(codfname)
  write(*,*) ncod
  call read_codes(codfname,codes,classes,ncod)
  missind = minloc(abs(codes-misscod),misscod) ! COD index of the missing value

  dumchar = format_classes(codes,classes)
  ! write(*,*) dumchar
  stop

  ! Get the sizes of the 4D reference file
  call get_ref_dimzises(reffname,"PCT_PFT",outnlat,outnlon,refntim,npft)
  ! write(*,*) outnlon,outnlat,npft,ntim


  ! Get the lat lon 2D bounds from the reference file
  call get_ref_grid(reffname,outnlat,outnlon,outlats,outlatn,outlonw,outlone)

  ! Get the sizes of the input landuse flie
  call get_3d_dimsizes(inpfname,"landuse",inpnlat,inpnlon,ntim)

  write(*,*) inpnlon,inpnlat,ntim

  ! Get the years in the time dimension, but use syear to define integer years
  call get_landuse_years(inpfname,years,ntim)
  intyears = ispan(syear,syear+ntim-1,1)

  write(*,*) intyears


  write(*,*) ncod
  ! Allocate the output variables
  allocate(outdata(outnlon,outnlat,ncod))
  allocate(fulloutdata(outnlon,outnlat,ncod,ntim))
  write(*,*) ncod
  ! write(*,*) shape(outdata)

  write(*,*) shape(fulloutdata)
  !
  ! Get the lat lon 2D bounds, assuming a regular grid
  call get_landuse_grid(inpfname,inpnlat,inpnlon,vecinplats,vecinplatn,vecinplonw,vecinplone)
  ! Flip the longitude variable, assuming its not crossing Greenwich
  call flip_lon_nocross_vec(inpnlon,vecinplonw,vecinplone)

  write(*,*) shape(fulloutdata)


  do timestep = 1,ntim
    ! Read a single year of data
    write(*,*) ">>>>>>>>>>>>>>> Reading timestep ",timestep
    call read_landuse_data_timestep(inpfname,inpdata,inpnlat,inpnlon,timestep)
    write(*,*) ">>>>>>>>>>>>>>> Finished reading timestep ",timestep

    vecoutlats => outlats(1,:)
    vecoutlatn => outlatn(1,:)
    vecoutlonw => outlonw(:,1)
    vecoutlone => outlone(:,1)

    ! Preallocate outdata with zeros
    outdata(:,:,:) = 0.0d0

    ! Allocate outval and fill it with zeros
    if (.not.ALLOCATED(outval)) then
      allocate(outval(ncod))
    end if
    outval(:) = 0.0d0

    do outi = 1,outnlat
      ! write(*,*) "Running latitude ",outi," of ",outnlat
      write(*,'(a,a,1i10,a,1i10)',advance="no") char(13),"Running latitude ",outi," of ",outnlat
      do outj = 1,outnlon
exit
        ! write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< outi,outj >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ! write(*,*) outi,outj
        ! write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ! Out of bounds, continue the loop and give 100% to the missing value cod
        if (vecoutlats(outi).ge.vecinplatn(inpnlat) .or. vecoutlatn(outi).le.vecinplats(1) .or. vecoutlonw(outj).ge.vecinplone(inpnlon) .or. vecoutlone(outj).le.vecinplonw(1)) then
          outdata(outj,outi,:) = 0.0d0
          outdata(outj,outi,missind) = 100.0d0
          cycle
        end if


        latbnds = find_bound_inds_vec(vecoutlats(outi),vecoutlatn(outi),vecinplats,vecinplatn)
        lonbnds = find_bound_inds_vec(vecoutlonw(outj),vecoutlone(outj),vecinplonw,vecinplone)
        ! write(*,*) "latbnds = ",latbnds
        ! write(*,*) "lonbnds = ",lonbnds
        !
        ! write(*,*) "======================== vecout =========================="
        ! write(*,*) vecoutlats(outi),vecoutlatn(outi),vecoutlonw(outj),vecoutlone(outj)
        ! write(*,*) "=========================================================="

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

            ! write(*,*) i,j,wgt,inpdata(j,i)

          end do ! j, lonbnds
        end do !i, latbnds

        ! do k = 1,ncod
        !     write(*,*) codes(k),outval(k),100.0d0*outval(k)/totwgt
        ! end do

        ! Now divide by the weights
        outval = 100.0d0*outval/totwgt

        ! Write to the variable
        outdata(outj,outi,:) = outval(:)

      end do !outj, outnlat
    end do !outi, outnlon

    fulloutdata(:,:,:,timestep) = outdata
    !
  end do !timestep, ntim

  call write_rochedo_data(outfname,fulloutdata,outnlat,outnlon,ncod,ntim,outlats,outlatn,outlonw,outlone,intyears,codes,classes)



end program remap_rochedo
