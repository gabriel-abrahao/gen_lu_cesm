; Burns full-potential vegetation values into a preexisting CESM surfdata file

; Optionally, smooth the transition by interpolating PFT data in a smooth period

load "~/lib_abrahao.ncl"

begin

rootfolder      = "./"

;Folders relative to rootfolder
inputfolder     = "input/"
; rochedofolder   = "out_rochedo/"
potvegfolder    = "out_potveg/"
outputfolder    = "out_surfdata/"

scenario        = "veg"

inpfname        = "surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc"
outfname        = "surfdata.pftdyn_0.9x1.25_simyr1850-2005_" + scenario + "_c091008.nc"

; rocfname        = "fraction_" + scenario + "_comp_2012-2050.nc"

;potfname        = "mksrf_pft_potv_CN_0.9x1.25.nc"
potfname        = "mksrf_pft_potveg.c081009_0.9x1.25.nc"

; Reclassification period
syear           = 1950
eyear           = 2005

; Smoothing period, must contain some year in the reclassification period to make sense
smsyear		= 2005
smeyear		= 2012

forcenewcopy    = False ; Force recopying the reference file, even if the output exists
checksums	= True	; Check differences between the PFT sums in input and output files

reclass		= True ; Do the reclassification step
smooth		= False ; Do the smoothing step

inpfpath = inputfolder + inpfname
outfpath = outputfolder + outfname
; rocfpath = rochedofolder + rocfname
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
; arqroc = addfile(rocfpath,"r")

if (reclass) then
	; Read the potential vegetation variable, and make sure it has a coordinate variable
	arqpot = addfile(potfpath,"r")
	potveg = arqpot->PCT_PFT
	dimspotveg = dimsizes(potveg)
	dimsveg = dimspotveg
	npftpot = dimspotveg(0)
	potveg&pft = ispan(0,npftpot-1,1)

	do year = syear,eyear
		print("Reclassifying year " + year)
		;
		; roc_PCT_PFT = arqroc->PCT_PFT({year},:,:,:)
		; dimsroc = dimsizes(roc_PCT_PFT)
		; npftroc = dimsroc(0)
		;
		; ; The missing fraction of each point to WEIGHT WITH THE ORIGINAL surfdata LATER
		; roc_missfrac = roc_PCT_PFT(0,:,:)
		;
		; ; Set up a mask w/o almost 100% missing, water or urban values
		; roc_mask = where(roc_PCT_PFT(0,:,:)+roc_PCT_PFT(1,:,:)+roc_PCT_PFT(2,:,:) .ge. 99.9,False,True)
		;
		; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; ; 1: Renormalize all values, removing the missing and water fractions.
		; ; To avoid divisions by zero, we first apply the mask
		; masked_roc_PCT_PFT = mask(roc_PCT_PFT,roc_mask,True)
		; copy_VarMeta(roc_PCT_PFT,masked_roc_PCT_PFT)
		; roc_norm = masked_roc_PCT_PFT
		; roc_norm({0:2},:,:) = 0.0
		; do i = 3, npftroc-1 ; Start after urban (2)
		; 	roc_norm(i,:,:) = (100.0*masked_roc_PCT_PFT(i,:,:)) / (100-(masked_roc_PCT_PFT({0},:,:)+masked_roc_PCT_PFT({1},:,:)+masked_roc_PCT_PFT({2},:,:)))
		; end do
		veg_missfrac = new((/dimsveg(1), dimsveg(2)/), typeof(potveg))
		veg_missfrac = 0.0
		veg_mask = new((/dimsveg(1), dimsveg(2)/), logical)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; 2: Separate in macroclasses 0: Primary veg. 1: Planted forest 2: Crop 3: Pasture
		;Since we don't have information on secondary vegetation, all but the planted forests will be considered primary
		nmacro = 4
		veg_macro = new((/nmacro, dimsveg(1), dimsveg(2)/), typeof(potveg))
		copy_VarCoords_exl(potveg,veg_macro)
		veg_macro!0 = "cod"
		veg_macro&cod = ispan(0,nmacro-1,1)

		; Primary are only savanas, savanas_em_AP, florestas, florestas_em_AP. And class_9, that is all zero anyway
		veg_macro(0,:,:) = 100.0
		; Planted forest is only planted forests (30)
		veg_macro(1,:,:) = 0.0
		; Crops has to skip the planted fores (30) code
		veg_macro(2,:,:) = 0.0
		; Pastures are only pastagem and pastagem_em_AP
		veg_macro(3,:,:) = 0.0


		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; 3: Calculate the final values as a sum of PFT fraction matrices weighted by the macveglasses
		; Macveglass 1, planted forest, represented here as all broadleaf_evergreen_temperate trees (5)
		all_forest = potveg
		all_forest(:,:,:) = 0.0
		all_forest({5},:,:) = 100.0
		; Macveglass 2, all crops, represented as the corn PFT (15), the only working crop PFT in the CLM setup for the RCPs
		all_crop = potveg
		all_crop(:,:,:) = 0.0
		all_crop({15},:,:) = 100.0
		; Macveglass 3, all pasture, represented as C4 grasses (14) for all of Brazil
		all_past = potveg
		all_past(:,:,:) = 0.0
		all_past({14},:,:) = 100.0
		; Turn veg_macro to weights (0-1) with higher dimensionality
		wgt_veg_macro = new((/nmacro, npftpot, dimspotveg(1), dimspotveg(2)/),typeof(veg_macro))
		do i = 0, npftpot-1
			wgt_veg_macro(:,i,:,:) = veg_macro/100.0
		end do
		; Apply the weights
		veg_out_pft = wgt_veg_macro(0,:,:,:)*potveg + wgt_veg_macro(1,:,:,:)*all_forest + wgt_veg_macro(2,:,:,:)*all_crop + wgt_veg_macro(3,:,:,:)*all_past
		copy_VarMeta(potveg,veg_out_pft)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; 4: To account for the Rochedo border fractions, we get the values from CESMs original file and weight it with the results
		inp_PCT_PFT = arqinp->PCT_PFT({year},:,:,:)
		; Turn veg_missfrac to weights
		wgt_veg_missfrac = new((/npftpot, dimspotveg(1), dimspotveg(2)/),typeof(veg_missfrac))
		do i = 0, npftpot-1
			wgt_veg_missfrac(i,:,:) = veg_missfrac/100.0
		end do
		out_PCT_PFT = (wgt_veg_missfrac/100.0)*inp_PCT_PFT + (1-wgt_veg_missfrac/100.0)*veg_out_pft
		out_PCT_PFT = where(ismissing(out_PCT_PFT),inp_PCT_PFT,out_PCT_PFT)
		out_PCT_PFT = where(isnan_ieee(out_PCT_PFT),inp_PCT_PFT,out_PCT_PFT) ; Not sure where those NaNs come from, appear to be a border issue
		copy_VarMeta(inp_PCT_PFT,out_PCT_PFT)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		; 5: Re-normalize to the actual sum of PCT_PFT, that does not add up to 100 because of urban, wetland and other fractions that are not accounted for

		sum_inp_PCT_PFT = dim_sum_n(inp_PCT_PFT,0)
		sum_out_PCT_PFT = dim_sum_n(out_PCT_PFT,0)
		sum_out_PCT_PFT = where(sum_out_PCT_PFT.eq.0,sum_out_PCT_PFT@_FillValue,sum_out_PCT_PFT)
		do i = 0, npftpot-1
			;out_PCT_PFT(i,:,:) = out_PCT_PFT(i,:,:)*sum_inp_PCT_PFT/100.0
			out_PCT_PFT(i,:,:) = where(.not.ismissing(sum_out_PCT_PFT),(out_PCT_PFT(i,:,:)/sum_out_PCT_PFT)*sum_inp_PCT_PFT,inp_PCT_PFT(i,:,:))
		end do

		if (checksums) then
			sum_inp_PCT_PFT = dim_sum_n(inp_PCT_PFT,0)
			sum_out_PCT_PFT = dim_sum_n(out_PCT_PFT,0)
		       	diffsum = sum_out_PCT_PFT - sum_inp_PCT_PFT
			print("Checksum error: "+max((diffsum)))
			dum_write(diffsum,"diffsum")
		end if


		; Write to output "frame"
		print("Writing output frame " + year)
		arqout->PCT_PFT({year},:,:,:) = (/out_PCT_PFT/)

	end do ; years, reclass
end if ; reclass

if (smooth) then
	nysmooth = smeyear - smsyear
	sta_PCT_PFT = arqout->PCT_PFT({smsyear},:,:,:)
	end_PCT_PFT = arqout->PCT_PFT({smeyear},:,:,:)
	dimsin = dimsizes(sta_PCT_PFT)
	npft = dimsin(0)
	printVarSummary(sta_PCT_PFT)
	printVarSummary(end_PCT_PFT)
	delta = (end_PCT_PFT - sta_PCT_PFT)/nysmooth

	; START TIME LOOP
	c = 0
	do year = smsyear+1,smeyear-1
		c = c + 1
		print("Smoothing year " + year + " (from " + smsyear + " to " + smeyear + ")")
		inp_PCT_PFT = arqout->PCT_PFT({year},:,:,:)
		; Add the delta to the starting value
		out_PCT_PFT = sta_PCT_PFT + delta*c
		; Make sure its within range
		out_PCT_PFT = where(out_PCT_PFT.gt.100.0,100.0,where(out_PCT_PFT.lt.0.0,0.0,out_PCT_PFT))
		; Make sure it adds up to the same amount as the input
		sum_inp_PCT_PFT = dim_sum_n(inp_PCT_PFT,0)
		sum_out_PCT_PFT = dim_sum_n(out_PCT_PFT,0)
		sum_out_PCT_PFT = where(sum_out_PCT_PFT.eq.0,sum_out_PCT_PFT@_FillValue,sum_out_PCT_PFT)
		do i = 0, npft-1
			;out_PCT_PFT(i,:,:) = out_PCT_PFT(i,:,:)*sum_inp_PCT_PFT/100.0
			out_PCT_PFT(i,:,:) = where(.not.ismissing(sum_out_PCT_PFT),(out_PCT_PFT(i,:,:)/sum_out_PCT_PFT)*sum_inp_PCT_PFT,inp_PCT_PFT(i,:,:))
		end do

		if (checksums) then
			sum_inp_PCT_PFT = dim_sum_n(inp_PCT_PFT,0)
			sum_out_PCT_PFT = dim_sum_n(out_PCT_PFT,0)
		       	diffsum = sum_out_PCT_PFT - sum_inp_PCT_PFT
			print("Checksum error: "+max((diffsum)))
			dum_write(diffsum,"diffsum")
		end if

		arqout->PCT_PFT({year},:,:,:) = (/out_PCT_PFT/)
	end do
end if

; Add some comments as a global attribute
globalatts = True
globalatts@reclass_veghedo_comments = arqout@reclass_veghedo_comments + "\n" +systemfunc("date") + ": Ran " + get_script_name() + " with reclass=" + reclass + "(" + syear + "-" + eyear + ") " +" and smooth=" + smooth + "(" + smsyear + "-" + smeyear + ")"

fileattdef(arqout,globalatts)

; dum_write(inp_PCT_PFT,"inp_PCT_PFT")
; dum_write(dim_sum_n(inp_PCT_PFT,0),"sum_inp_PCT_PFT")
; dum_write(out_PCT_PFT,"out_PCT_PFT")
; dum_write(dim_sum_n(out_PCT_PFT,0),"sum_out_PCT_PFT")

; dum_write(out_PCT_NAT_PFT,"out_PCT_NAT_PFT")
; dum_write(dim_sum_n(out_PCT_NAT_PFT,0),"sum_out_PCT_NAT_PFT")

; dum_write(veg_out_pft,"roc_out_pft")
; dum_write(dim_sum_n(roc_out_pft,0),"sum_roc_out_pft")

; PCT_NAT_PFT = arqout->PCT_NAT_PFT({year},:,:,:)
; PCT_NAT_PFT_sum = dim_sum_n_Wrap(PCT_NAT_PFT,0)
;
; dum_write(PCT_NAT_PFT_sum,"PCT_NAT_PFT_sum")
; dum_write(dim_sum_n_Wrap(roc_PCT_PFT,0),"rochedo")

end ;Script
