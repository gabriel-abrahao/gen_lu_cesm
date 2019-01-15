load "~/lib_abrahao.ncl"
; Generate a potential vegetation file in the reference grid using a high resolution map
begin

inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

reffname = inputfolder + "surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
vegfname = inputfolder + "vegtype_5min.nc"

outfname = outputfolder + "vegtype_0.9x1.25.nc"

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

; Reference variable, for the grid
refvar = reffile->PCT_PFT(0,0,:,:)
refvar = refvar@_FillValue

; 2D vars with the coordinates of the boundaries
lats = reffile->LATS
latn = reffile->LATN
lone = reffile->LONE
lonw = reffile->LONW

printVarSummary(vegvar)
printVarSummary(refvar)

dims = dimsizes(refvar)
nlatref = dims(0)
nlonref = dims(1)

; Create the output (pft, lat, lon) variable with PFT fractions
outvar = new((/npft, nlatref, nlonref/), typeof(refvar))
outvar!0 = "pft"
outvar&pft = ispan(minpft,maxpft,1)
outvar@comment = "Fraction refers only to nonmissing pixels in the original map, and as such should add up to 100 where there is some land and 0 where theres none"

;print(min(lats))
;print(min(latn))
;
do i = 0, nlatref-1
	print(i + " of " + (nlatref-1))
	do j = 0, nlonref-1
		;i = 91
		;j = 239
		;
		;print("lats(i,j) = " + lats(i,j))
		;print("latn(i,j) = " + latn(i,j))
		;print("lonw(i,j) = " + lonw(i,j))
		;print("lone(i,j) = " + lone(i,j))


		vegpix = vegvar({lats(i,j):latn(i,j)},{lonw(i,j):lone(i,j)})
		;print(dimsizes(vegpix))
		;print(vegpix)
		if (all(ismissing(vegpix))) then
			outvar(:,i,j) = 0.0
		else
			do k = minpft,maxpft	
				;k = 1
				;print(num(vegpix.eq.k))
				;print(num(.not.ismissing(vegpix)))
				;print( 100.0*(todouble(num(vegpix.eq.k)) / todouble( num(.not.ismissing(vegpix))) ))
				outvar({k},i,j) = (100.0*(todouble(num(vegpix.eq.k)) / todouble( num(.not.ismissing(vegpix))) ))
			end do ;k, pfts
		end if ; all(ismissing(vegpix))
		delete(vegpix) ; Pixel size can vary
	end do ;j, lons
end do ;i, lats

; Write output
system("rm " + outfname)
outfile = addfile(outfname,"c")
outfile->vegtype = outvar

;; This is to check that, in the CESM file, PFT fractions dont add up to 100
;poi = reffile->PCT_PFT((/0,155/),:,:,:)
;printVarSummary(poi)
;soma = dim_sum_n_Wrap(poi,1)
;printVarSummary(soma)
;dum_write(soma,"soma")

end
