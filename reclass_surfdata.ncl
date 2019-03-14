; Reclassifies a fractional Rochedo et al. (2018) classes file into a preexisting CESM surfdata file
; Apparently the RCP surfdata files have a PCT_NAT_PFT variable that sums to 100 in all PFTs, including crop (15),

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

syear           = 2012
eyear           = 2050

inpfpath = inputfolder + inpfname
outfpath = outputfolder + outfname
rocfpath = rochedofolder + rocfname

; Make a copy of inpfname to edit on
system("! [ -e " + outfpath + " ] && cp " + inpfpath + " " + outfpath + " && 'Making an exact copy of inputfile'")

; Open the output file, which is also serves as the input one
arqout = addfile(outfpath,"w")

arqroc = addfile(rocfpath,"r")

; TIME LOOP HERE
year = 2015

roc_PCT_PFT = arqroc->PCT_PFT({year},:,:,:)
dimsroc = dimsizes(roc_PCT_PFT)
npftroc = dimsroc(0)
print(npftroc)
printVarSummary(roc_PCT_PFT)

; The missing fraction of each point to WEIGHT WITH THE ORIGINAL surfdata LATER
roc_missfrac = roc_PCT_PFT(0,:,:)

; Set up a mask w/o almost 100% missing or water values
roc_mask = where(roc_PCT_PFT(0,:,:)+roc_PCT_PFT(1,:,:) .ge. 99.9,False,True)

; 1: Renormalize all values, removing the missing and water fractions.
; To avoid divisions by zero, we first apply the mask
masked_roc_PCT_PFT = mask(roc_PCT_PFT,roc_mask,True)
copy_VarMeta(roc_PCT_PFT,masked_roc_PCT_PFT)
roc_norm = masked_roc_PCT_PFT
roc_norm({0:1},:,:) = 0.0
do i = 2, npftroc-1 ; Start after the water
  roc_norm(i,:,:) = (100.0*masked_roc_PCT_PFT(i,:,:)) /(100-(masked_roc_PCT_PFT({0},:,:)+masked_roc_PCT_PFT({1},:,:)))
end do

dum_write(roc_norm,"roc_norm")
dum_write(dim_sum_n_Wrap(roc_norm,0),"sum_roc_norm")

; PCT_NAT_PFT = arqout->PCT_NAT_PFT({year},:,:,:)
; PCT_NAT_PFT_sum = dim_sum_n_Wrap(PCT_NAT_PFT,0)
;
; dum_write(PCT_NAT_PFT_sum,"PCT_NAT_PFT_sum")
; dum_write(dim_sum_n_Wrap(roc_PCT_PFT,0),"rochedo")

end ;Script
