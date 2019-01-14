load "~/lib_abrahao.ncl"
; Generate a potential vegetation file in the reference grid using a high resolution map
begin

inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

reffname = inputfolder + "surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
vegfname = inputfolder + "vegtype_5min.nc"

; Flip the longitude of the veg map?
flipveg = True 

; PFT code range
minpft = 1
maxpft = 15
npft = maxpft - minpft + 1

reffile = addfile(reffname,"r")
vegfile = addfile(vegfname,"r")

; Flip the longitude of the veg file if needed
if (flipveg) then
	dumvegvar = vegfile->vegtype(0,0,:,:)
	vegvar = lonFlip(dumvegvar)
else
	vegvar = vegfile->vegtype(0,0,:,:)
end if

refvar = reffile->PCT_PFT(0,0,:,:)

lats = reffile->LATS
latn = reffile->LATN
lone = reffile->LONE
lonw = reffile->LONW

printVarSummary(vegvar)
printVarSummary(refvar)

dims = dimsizes(refvar)
nlatref = dims(0)
nlonref = dims(1)

;print(min(lats))
;print(min(latn))
;
i = 91
j = 239

print(lats(i,j))
print(latn(i,j))
print(lonw(i,j))
print(lone(i,j))


vegpix = vegvar({lats(i,j):latn(i,j)},{lonw(i,j):lone(i,j)})
print(vegpix)

end
