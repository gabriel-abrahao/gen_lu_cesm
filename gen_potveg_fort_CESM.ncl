external EX01 "./get_pixel_weights.so"
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

outalldims = dimsizes(refvar)
outnlat = outalldims(0)
outnlon = outalldims(1)

lonvar = (lone(0,:) + lonw(0,:))/2.0d0
latvar = (lats(:,0) + latn(:,0))/2.0d0


; printVarSummary(vegvar)
; printVarSummary(refvar)

dims = dimsizes(refvar)
nlatref = dims(0)
nlonref = dims(1)

; Assume all pixels in the subpixel have the same width/height
inpalldims = dimsizes(vegvar)
inpnlat = inpalldims(1)
inpnlon = inpalldims(2)
veglonres = (veglon(dimsizes(veglon)-1) - veglon(0))/(dimsizes(veglon)-1)
vecveglonw = veglon - (veglonres/2.0)
vecveglone = veglon + (veglonres/2.0)
veglatres = (veglat(dimsizes(veglat)-1) - veglat(0))/(dimsizes(veglat)-1)
vecveglats = veglat - (veglatres/2.0)
vecveglatn = veglat + (veglatres/2.0)
veglats = new((/inpnlat,inpnlon/),double)
veglatn = new((/inpnlat,inpnlon/),double)
veglonw = new((/inpnlat,inpnlon/),double)
veglone = new((/inpnlat,inpnlon/),double)

veglats = conform(veglats,vecveglats,0)
veglatn = conform(veglatn,vecveglatn,0)
veglonw = conform(veglonw,vecveglonw,1)
veglone = conform(veglone,vecveglone,1)


; Create the output (pft, lat, lon) variable with PFT fractions
outvar = new((/npft, nlatref, nlonref/), typeof(refvar))
copy_VarCoords_1(vegvar,outvar) ;Copy PFT definition
copy_VarAtts(vegvar,outvar) ;Copy attributes
outvar@comment = "Fraction refers only to nonmissing pixels in the original map, and as such should add up to 100 where there is some land and 0 where theres none"
outvar!1 = "lat"
outvar!2 = "lon"
outvar&lat = latvar
outvar&lon = lonvar

; Weights of a single pixel
wgt = new((/inpnlat,inpnlon/),double)

; printVarSummary(outvar)

;print(min(lats))
;print(min(latn))
;
do i = 0, nlatref-1
	print(i + " of " + (nlatref-1))
	do j = 0, nlonref-1
;i = 91
;	j = 239
; i = 2
; 	j = 2

		; Zero all weights first
		wgt = 0.0
		; print("================== Start Call...")
		EX01::get_pixel_weights(wgt,(i+1),(j+1),veglats,veglatn,veglonw,veglone,inpnlat,inpnlon,lats,latn,lonw,lone,outnlat,outnlon)
		; print("================== End Call...")



		do k = 0,npft-1
			outvar(k,i,j) = sum(wgt*vegvar(k,i,j))/sum(wgt)
		end do

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
