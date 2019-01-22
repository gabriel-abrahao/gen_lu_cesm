program gen_potveg_CESM
  ! use fpl !FPL module
  use netcdf
  use tools
  implicit none
  character*200 inputfolder,outputfolder,reffname,vegfname,outfname
  integer ncid
  integer status, varid
  integer vegnlat,vegnlon,refnlat,refnlon,ntim,npft
  real*8, ALLOCATABLE, DIMENSION(:,:) :: reflats,reflatn,reflonw,reflone,veglats,veglatn,veglonw,veglone
  real*8, ALLOCATABLE, DIMENSION(:,:) :: vegdata
  !real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: vardata

  inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
  outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

  reffname = trim(ADJUSTL(inputfolder))//"surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
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

  call read_veg_data(vegfname,vegnlat,vegnlon,vegdata)


end program gen_potveg_CESM
