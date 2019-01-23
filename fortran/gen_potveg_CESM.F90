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
  integer, pointer :: inpnlon,inpnlat
  integer, pointer :: outnlon,outnlat

  logical poi

  inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
  outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

  ! reffname = trim(ADJUSTL(inputfolder))//"surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
  reffname = trim(ADJUSTL(inputfolder))//"min_surfdata.nc"
  vegfname = trim(ADJUSTL(inputfolder))//"mksrf_pft_potv_CN.nc"

  outfname = trim(ADJUSTL(outputfolder))//"mksrf_pft_potv_CN_0.9x1.25.nc"

! Get the sizes of the 4D reference file
  call get_ref_dimzises(reffname,"PCT_PFT",refnlat,refnlon,ntim,npft)
  write(*,*) refnlon,refnlat,npft,ntim

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

  write(*,*) reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2)
  write(*,*) is_inside_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))
  write(*,*) is_contained_vec(reflats(1,2),reflatn(1,2),veglats(1,2),veglatn(1,2))


end program gen_potveg_CESM
