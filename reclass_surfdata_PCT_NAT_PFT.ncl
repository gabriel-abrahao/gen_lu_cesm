; Reclassifies a fractional Rochedo et al. (2018) classes file into a preexisting CESM surfdata file
; Apparently the RCP surfdata files have a PCT_NAT_PFT variable that sums to 100 in all PFTs, including crop (15)
; Therefore, the approach here ignores the water and urban fractions from rochedo

load "~/lib_abrahao.ncl"

begin

rootfolder      = "./"

;Folders relative to rootfolder
inputfolder     = "input/"
rochedofolder   = "out_rochedo/"
potvegfolder    = "out_potveg/"
outputfolder    = "out_surfdata/"

; weg or seg
scenario        = "weg"

inpfname        = "surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_c130524.nc"
outfname        = "surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_rochedo_" + scenario + "_c130524.nc"

rocfname        = "fraction_" + scenario + "_comp_2012-2050.nc"

potfname        = "mksrf_pft_potv_CN_0.9x1.25.nc"

syear           = 2012
eyear           = 2050

forcenewcopy    = False ; Force recopying the reference file, even if the output exists

inpfpath = inputfolder + inpfname
outfpath = outputfolder + outfname
rocfpath = rochedofolder + rocfname
potfpath = potvegfolder + potfname

; Make a copy of inpfname to edit on. Just copy if it doesnt exist
if (.not.forcenewcopy) then
  print("The output file already exists? (erase it if its wrong or set forcenewcopy) ")
  system("[ -e " + outfpath + " ] && echo YES")
  system("! [ -e " + outfpath + " ] && cp " + inpfpath + " " + outfpath + " && 'Making an exact copy of inputfile'")
else
    print("Forcing the copy of a new output file")
    system("cp " + inpfpath + " " + outfpath)
end if

; Open the input file
arqinp = addfile(inpfpath,"r")

; Open the output file, which should have been sometime initialized with the input one
arqout = addfile(outfpath,"w")

; Open the rochedo file from remap_rochedo
arqroc = addfile(rocfpath,"r")

; Read the potential vegetation variable, and make sure it has a coordinate variable
arqpot = addfile(potfpath,"r")
potveg = arqpot->PCT_PFT
dimspotveg = dimsizes(potveg)
npftpot = dimspotveg(0)
potveg&pft = ispan(0,npftpot-1,1)

; TIME LOOP HERE
; year = 2015
do year = syear,eyear
  print("Running year " + year)

  roc_PCT_PFT = arqroc->PCT_PFT({year},:,:,:)
  dimsroc = dimsizes(roc_PCT_PFT)
  npftroc = dimsroc(0)

  ; The missing fraction of each point to WEIGHT WITH THE ORIGINAL surfdata LATER
  roc_missfrac = roc_PCT_PFT(0,:,:)

  ; Set up a mask w/o almost 100% missing, water or urban values
  roc_mask = where(roc_PCT_PFT(0,:,:)+roc_PCT_PFT(1,:,:)+roc_PCT_PFT(2,:,:) .ge. 99.9,False,True)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 1: Renormalize all values, removing the missing and water fractions.
  ; To avoid divisions by zero, we first apply the mask
  masked_roc_PCT_PFT = mask(roc_PCT_PFT,roc_mask,True)
  copy_VarMeta(roc_PCT_PFT,masked_roc_PCT_PFT)
  roc_norm = masked_roc_PCT_PFT
  roc_norm({0:2},:,:) = 0.0
  do i = 3, npftroc-1 ; Start after urban (2)
    roc_norm(i,:,:) = (100.0*masked_roc_PCT_PFT(i,:,:)) / (100-(masked_roc_PCT_PFT({0},:,:)+masked_roc_PCT_PFT({1},:,:)+masked_roc_PCT_PFT({2},:,:)))
  end do

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 2: Separate in macroclasses 0: Primary veg. 1: Planted forest 2: Crop 3: Pasture
  ;Since we don't have information on secondary vegetation, all but the planted forests will be considered primary
  nmacro = 4
  roc_macro = new((/nmacro, dimsroc(1), dimsroc(2)/), typeof(roc_norm))
  copy_VarCoords_exl(roc_norm,roc_macro)
  roc_macro!0 = "cod"
  roc_macro&cod = ispan(0,nmacro-1,1)

  ; Primary are only savanas, savanas_em_AP, florestas, florestas_em_AP. And class_9, that is all zero anyway
  roc_macro(0,:,:) = (/ dim_sum_n(roc_norm({5:9},:,:),0) /)
  ; Planted forest is only planted forests (30)
  roc_macro(1,:,:) = (/ roc_norm({30},:,:) /)
  ; Crops has to skip the planted fores (30) code
  roc_macro(2,:,:) = (/ dim_sum_n(roc_norm({11:26},:,:),0) + dim_sum_n(roc_norm({35:40},:,:),0)/)
  ; Pastures are only pastagem and pastagem_em_AP
  roc_macro(3,:,:) = (/ dim_sum_n(roc_norm({3:4},:,:),0) /)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 3: Calculate the final values as a sum of PFT fraction matrices weighted by the macroclasses
  ; Macroclass 1, planted forest, represented here as all broadleaf_evergreen_temperate trees (5)
  all_forest = potveg
  all_forest(:,:,:) = 0.0
  all_forest({5},:,:) = 100.0
  ; Macroclass 2, all crops, represented as the corn PFT (15), the only working crop PFT in the CLM setup for the RCPs
  all_crop = potveg
  all_crop(:,:,:) = 0.0
  all_crop({15},:,:) = 100.0
  ; Macroclass 3, all pasture, represented as C4 grasses (14) for all of Brazil
  all_past = potveg
  all_past(:,:,:) = 0.0
  all_past({14},:,:) = 100.0
  ; Turn roc_macro to weights (0-1) with higher dimensionality
  wgt_roc_macro = new((/nmacro, npftpot, dimspotveg(1), dimspotveg(2)/),typeof(roc_macro))
  do i = 0, npftpot-1
    wgt_roc_macro(:,i,:,:) = roc_macro/100.0
  end do
  ; Apply the weights
  roc_out_pft = wgt_roc_macro(0,:,:,:)*potveg + wgt_roc_macro(1,:,:,:)*all_forest + wgt_roc_macro(2,:,:,:)*all_crop + wgt_roc_macro(3,:,:,:)*all_past
  copy_VarMeta(potveg,roc_out_pft)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 4: To account for the Rochedo border fractions, we get the values from CESMs original file and weight it with the results
  inp_PCT_NAT_PFT = arqinp->PCT_NAT_PFT({year},:,:,:)
  ; Turn roc_missfrac to weights
  wgt_roc_missfrac = new((/npftpot, dimspotveg(1), dimspotveg(2)/),typeof(roc_missfrac))
  do i = 0, npftpot-1
    wgt_roc_missfrac(i,:,:) = roc_missfrac/100.0
  end do
  out_PCT_NAT_PFT = (wgt_roc_missfrac/100.0)*inp_PCT_NAT_PFT + (1-wgt_roc_missfrac/100.0)*roc_out_pft
  out_PCT_NAT_PFT = where(ismissing(out_PCT_NAT_PFT),inp_PCT_NAT_PFT,out_PCT_NAT_PFT)
  out_PCT_NAT_PFT = where(isnan_ieee(out_PCT_NAT_PFT),inp_PCT_NAT_PFT,out_PCT_NAT_PFT) ; Not sure where those NaNs come from, appear to be a border issue
  copy_VarMeta(inp_PCT_NAT_PFT,out_PCT_NAT_PFT)

  ; Write to output "frame"
  arqout->PCT_NAT_PFT({year},:,:,:) = (/out_PCT_NAT_PFT/)

end do ; years

; dum_write(out_PCT_NAT_PFT,"out_PCT_NAT_PFT")
; dum_write(dim_sum_n(out_PCT_NAT_PFT,0),"sum_out_PCT_NAT_PFT")

; dum_write(roc_out_pft,"roc_out_pft")
; dum_write(dim_sum_n(roc_out_pft,0),"sum_roc_out_pft")

; PCT_NAT_PFT = arqout->PCT_NAT_PFT({year},:,:,:)
; PCT_NAT_PFT_sum = dim_sum_n_Wrap(PCT_NAT_PFT,0)
;
; dum_write(PCT_NAT_PFT_sum,"PCT_NAT_PFT_sum")
; dum_write(dim_sum_n_Wrap(roc_PCT_PFT,0),"rochedo")

end ;Script
