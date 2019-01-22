program gen_potveg_CESM
  ! use fpl !FPL module
  use netcdf
  use tools
  implicit none
  character*200 inputfolder,outputfolder,reffname,vegfname,outfname
  integer ncid
  integer status, varid
  integer nlat,nlon,ntim,npft
  real*8, ALLOCATABLE, DIMENSION(:,:) :: lats,latn,lonw,lone
  real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: vardata

  inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
  outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

  reffname = trim(ADJUSTL(inputfolder))//"surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
  vegfname = trim(ADJUSTL(inputfolder))//"mksrf_pft_potv_CN.nc"

  outfname = trim(ADJUSTL(outputfolder))//"mksrf_pft_potv_CN_0.9x1.25.nc"

  call get_ref_dimzises(reffname,"PCT_PFT",nlat,nlon,ntim,npft)

  write(*,*) nlon,nlat,npft,ntim

  call get_ref_grid(reffname,nlat,nlon,lats,latn,lonw,lone)


end program gen_potveg_CESM
