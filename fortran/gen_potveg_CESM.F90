program gen_potveg_CESM
  ! use fpl !FPL module
  use netcdf
  use tools
  implicit none
  character*200 inputfolder,outputfolder,reffname,vegfname,outfname
  integer ncid
  integer status, varid
  integer, TARGET :: vegnlat,vegnlon,refnlat,refnlon,ntim,npft
  real*8, TARGET, ALLOCATABLE, DIMENSION(:,:) :: reflats,reflatn,reflonw,reflone,veglats,veglatn,veglonw,veglone
  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: vegdata
  real*8, ALLOCATABLE, DIMENSION(:,:,:) :: outdata

  real*8, pointer :: inplats(:,:),inplatn(:,:),inplonw(:,:),inplone(:,:)
  real*8, pointer :: outlats(:,:),outlatn(:,:),outlonw(:,:),outlone(:,:)
  real*8, pointer :: vecinplats(:),vecinplatn(:),vecinplonw(:),vecinplone(:)
  real*8, pointer :: vecoutlats(:),vecoutlatn(:),vecoutlonw(:),vecoutlone(:)
  integer, pointer :: inpnlon,inpnlat
  integer, pointer :: outnlon,outnlat

  integer outi,outj,i,j,k,ii,jj, lasti,lastj

  integer latbnds(2),lonbnds(2)

  real*8 latsize,lonsize,latwgt,lonwgt,wgt,totwgt
  real*8, ALLOCATABLE, dimension(:) ::  outval ! The output value of a pixel in every PFT

  inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
  outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

  ! reffname = trim(ADJUSTL(inputfolder))//"surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
  reffname = trim(ADJUSTL(inputfolder))//"min_surfdata.nc"
  vegfname = trim(ADJUSTL(inputfolder))//"mksrf_pft_potv_CN.nc"

  outfname = trim(ADJUSTL(outputfolder))//"mksrf_pft_potv_CN_0.9x1.25.nc"

! Get the sizes of the 4D reference file
  call get_ref_dimzises(reffname,"PCT_PFT",refnlat,refnlon,ntim,npft)
  write(*,*) refnlon,refnlat,npft,ntim

  ! Allocate the output variable
  allocate(outdata(refnlon,refnlat,npft))

! Get the lat lon 2D bounds from the reference file
  call get_ref_grid(reffname,refnlat,refnlon,reflats,reflatn,reflonw,reflone)

! Get the sizes of the potveg flie
  call get_2d_dimsizes(vegfname,"PCT_PFT",vegnlat,vegnlon)
  write(*,*) vegnlon,vegnlat

! Get the lat lon 2D bounds, assuming a regular grid
  call get_veg_grid(vegfname,vegnlat,vegnlon,veglats,veglatn,veglonw,veglone)

! Read the 3d (npft) potential vegetation file
  call read_veg_data(vegfname,vegnlat,vegnlon,npft,vegdata)

! Flip the longitude variables (TODO: This assumes a lot about the dataset as is, make it more generic)
  call flip_lon_global_3d(vegdata,vegnlat,vegnlon,npft,veglats,veglatn,veglonw,veglone)

  ! call dum_write_2d("dummy.nc",veglone,vegnlat,vegnlon) ! Checking
  ! call dum_write_3d("dummy.nc",vegdata,vegnlat,vegnlon,npft) ! Checking

  ! Use pointers to more generic names here (TODO: Refactor the code)
  inpnlon => vegnlon
  inpnlat => vegnlat

  inplats => veglats
  inplatn => veglatn
  inplonw => veglonw
  inplone => veglone

  outnlon => refnlon
  outnlat => refnlat

  outlats => reflats
  outlatn => reflatn
  outlonw => reflonw
  outlone => reflone

  ! FIXME: For now, we'll simply assume the grid is regular
  vecinplats => inplats(1,:)
  vecinplatn => inplatn(1,:)
  vecinplonw => inplonw(:,1)
  vecinplone => inplone(:,1)

  vecoutlats => outlats(1,:)
  vecoutlatn => outlatn(1,:)
  vecoutlonw => outlonw(:,1)
  vecoutlone => outlone(:,1)

  ! Preallocate outdata with zeros
  outdata(:,:,:) = 0.0d0

  ! Allocate outval and fill it with zeros
  allocate(outval(npft))
  outval(:) = 0.0d0



  ! Put output loop here
  ! outi = 67
  ! outj = 232
  do outi = 1,outnlat
    write(*,*) "Running latitude ",outi," of ",outnlat
    do outj = 1,outnlon

      ! write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< outi,outj >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      ! write(*,*) outi,outj
      ! write(*,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      ! Out of bounds, continue the loop
      if (vecoutlats(outi).ge.vecinplatn(inpnlat) .or. vecoutlatn(outi).le.vecinplats(1) .or. vecoutlonw(outj).ge.vecinplone(inpnlon) .or. vecoutlone(outj).le.vecinplonw(1)) then
        continue
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
          totwgt = totwgt + wgt


          ! Accumulate val*wgt in outval
          do k = 1,npft
            outval(k) = outval(k) + vegdata(j,i,k)*wgt
          end do !k, pft

        end do ! j, lonbnds
      end do !i, latbnds

      ! Now divide by the weights
      outval = outval/totwgt

      ! Write to the variable
      outdata(outj,outi,:) = outval(:)

    end do !outj, outnlat
  end do !outi, outnlon

  !TODO: CORRECT FRACTIONS TO 100% USING LANDFRAC_PFT and PFTDATA_MASK

  call dum_write_3d("dummy.nc",outdata,outnlat,outnlon,npft) ! Checking


  ! write(*,*) reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2)
  ! write(*,*) is_inside_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))
  ! write(*,*) is_contained_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))


end program gen_potveg_CESM
