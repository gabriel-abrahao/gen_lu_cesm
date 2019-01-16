load "~/lib_abrahao.ncl"
; Generate a potential vegetation file in the reference grid using the high resolution map derived using the CESM tools in:
; https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/rawdata/dynlnduse_code.050415/
begin

inputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/input/"
outputfolder = "/home/gabriel/transicao/doutorado/gen_lu_cesm/out_potveg/"

reffname = inputfolder + "surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
vegfname = inputfolder + "mksrf_pft_potv_CN.nc"

outfname = outputfolder + "mksrf_pft_potv_CN_0.9x1.25.nc"

; Flip the longitude of the veg map?
flipveg = True 

; PFT code range
minpft = 0
maxpft = 16
npft = maxpft - minpft + 1

reffile = addfile(reffname,"r")
vegfile = addfile(vegfname,"r")

; Put proper coordinate variables
vegvar = vegfile->PCT_PFT
vegvar&lat = vegfile->LAT
vegvar&lon = vegfile->LON

; Flip the longitude of the veg file if needed
if (flipveg) then
	dumvegvar = vegvar
	delete(vegvar)
	vegvar = lonFlip(dumvegvar)
	delete(dumvegvar)
end if

; Get coordinate variables with coordinate indexing capabilities
veglat = vegvar&lat
veglon = vegvar&lon
veglat!0 = "lat"
veglon!0 = "lon"
veglat&lat = veglat
veglon&lon = veglon

; Reference variable, for the grid
refvar = reffile->PCT_PFT(0,0,:,:)
refvar = refvar@_FillValue

; 2D vars with the coordinates of the boundaries
lats = reffile->LATS
latn = reffile->LATN
lone = reffile->LONE
lonw = reffile->LONW

lonvar = (lone(0,:) + lonw(0,:))/2.0d0 
latvar = (lats(:,0) + latn(:,0))/2.0d0 

printVarSummary(vegvar)
printVarSummary(refvar)

dims = dimsizes(refvar)
nlatref = dims(0)
nlonref = dims(1)

; Assume all pixels in the subpixel have the same width/height
veglonres = (veglon(dimsizes(veglon)-1) - veglon(0))/(dimsizes(veglon)-1)
veglonw = veglon - (veglonres/2.0)
veglone = veglon + (veglonres/2.0)
veglonw!0 = "lon"
veglone!0 = "lon"
veglonw&lon = veglonw
veglone&lon = veglone
veglatres = (veglat(dimsizes(veglat)-1) - veglat(0))/(dimsizes(veglat)-1)
veglats = veglat - (veglatres/2.0)
veglatn = veglat + (veglatres/2.0)
veglats!0 = "lat"
veglatn!0 = "lat"
veglats&lat = veglats
veglatn&lat = veglatn

; Create the output (pft, lat, lon) variable with PFT fractions
outvar = new((/npft, nlatref, nlonref/), typeof(refvar))
copy_VarCoords_1(vegvar,outvar) ;Copy PFT definition
copy_VarAtts(vegvar,outvar) ;Copy attributes
outvar@comment = "Fraction refers only to nonmissing pixels in the original map, and as such should add up to 100 where there is some land and 0 where theres none"
outvar!1 = "lat"
outvar!2 = "lon"
outvar&lat = latvar
outvar&lon = lonvar

printVarSummary(outvar)

;print(min(lats))
;print(min(latn))
;
;do i = 0, nlatref-1
;	print(i + " of " + (nlatref-1))
;	do j = 0, nlonref-1
i = 91
	j = 239
		
;		latrange = latn(i,j)-lats(i,j)
;		lonrange = lone(i,j)-lonw(i,j)
;		lonepix = veglone({lonw(i,j):lone(i,j)})
;		lonwpix = veglonw({lonw(i,j):lone(i,j)})
;		latnpix = veglatn({lats(i,j):latn(i,j)})
;		latspix = veglats({lats(i,j):latn(i,j)})
;		print(latwpix)
;		print(latepix)
;		exit
print("lats(i,j) = " + lats(i,j))
print("latn(i,j) = " + latn(i,j))
print("lonw(i,j) = " + lonw(i,j))
print("lone(i,j) = " + lone(i,j))
print("latrange = " + latrange)
print("lonrange = " + lonrange)
print(lonpix)
exit

		vegpix = vegvar(:,{lats(i,j):latn(i,j)},{lonw(i,j):lone(i,j)})
		print(dimsizes(vegpix))
		print(vegpix)
		print(dim_sum_n_Wrap(vegpix,0))

		k = 0



;		if (all(ismissing(vegpix))) then
;			outvar(:,i,j) = 0.0
;		else
;			do k = minpft,maxpft	
;				;k = 1
;				;print(num(vegpix.eq.k))
;				;print(num(.not.ismissing(vegpix)))
;				;print( 100.0*(todouble(num(vegpix.eq.k)) / todouble( num(.not.ismissing(vegpix))) ))
;				outvar({k},i,j) = (100.0*(todouble(num(vegpix.eq.k)) / todouble( num(.not.ismissing(vegpix))) ))
;			end do ;k, pfts
;		end if ; all(ismissing(vegpix))

exit

		delete(vegpix) ; Pixel size can vary
;	end do ;j, lons
;end do ;i, lats

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
