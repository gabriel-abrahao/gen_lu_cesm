load "~/lib_abrahao.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;**********************************************************************
; Function to retrieve the location of a plot so we can calculate
; a location for a title on the left side of the plot.
;**********************************************************************
function get_title_position(plot)
local vpy, vpx, vph
begin
getvalues plot
"vpYF"       : vpy
"vpXF"       : vpx
"vpHeightF"  : vph
end getvalues
;pos = (/vpx*0.75,vpy-(vph/2.)/)
pos = (/vpx,vpy/)
return(pos)
end
begin

inpfile       =	"../out_surfdata/surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_rochedo_weg_c100319.nc"
varname       =	"PCT_PFT"

;syear		=	2005
;eyear		=	2014
years		=	ispan(2000,2050,5)

rcpstr		=	systemfunc("echo " + inpfile + "| grep -o 'rcp[0-9]\.[0-9]'")
scenariostr	=	systemfunc("echo " + inpfile + "| grep -o 'rochedo_[a-zA-Z]\{3\}'")

;foutname      =	"transition_frac_" + scenariostr + "_" + rcpstr + "_" + syear + "_" + eyear + ".nc"
foutname      =	"transition_frac_" + scenariostr + "_" + rcpstr + "_" + min(years) + "_" + max(years) + ".nc"

minval		=	0.
maxval		=	100.
spacing		=	10.

; cmap		=	"BlueRed"
cmap		=	"matlab_jet"

fillmode	=	"RasterFill"
;fillmode	=	"AreaFill"
;fillmode	=	"CellFill"

ftablepft	=	"../codes_cesm.csv"
usepfts		=	(/4,14,15/)


tablepft = asciiread(ftablepft,(/17/),"string")
names = str_get_field(tablepft,2,"	")



arqinp = addfile(inpfile,"r")
varinp = arqinp->$varname$({years},:,:,:)
;varinp = arqinp->$varname$({syear:eyear},:,:,:)

varinp&lsmlat = arqinp->LATIXY(:,0)
varinp&lsmlon = arqinp->LONGXY(0,:)

varinp&lsmlat@units = "degrees_north"
varinp&lsmlon@units = "degrees_east"


dims = dimsizes(varinp)
npft = dims(1)
nusepft = dimsizes(usepfts)
;nyears = eyear-syear+1
nyears = dimsizes(years)
printVarSummary(varinp)

;dum_write(varinp,"varinp")

wks_type = "pdf"
wks_type@wkOrientation = "Landscape"
; wks_type@wkPaperWidthF  =  8.5  ; in inches
; wks_type@wkPaperHeightF = 100.0  ; in inches
wks = gsn_open_wks(wks_type,foutname)

gsn_define_colormap(wks,cmap)


plots = new((/nusepft,nyears/),graphic)
onedplots = new((/nusepft*nyears/),graphic)
;plotsinp = new(npft,graphic)
;plotsref = new(npft,graphic)

res			     = True

res@gsnDraw = False
res@gsnFrame = False
; res@gsnMaximize = True

res@cnFillOn               = True               ; color fill
res@cnLinesOn              = False              ; no contour lines
res@cnLineLabelsOn         = False              ; no contour labels
res@cnInfoLabelOn          = False             ; no contour info label
res@vpHeightF              = 0.6
;res@lbBoxCount = 8

;  res@tmXBLabelFontHeightF  = 0.025
;  res@tmYLLabelFontHeightF  = 0.025
;  res@lbLabelFontHeightF   = 0.02

res@lbLabelBarOn = False

; res@tmYLLabelFont = 25
; res@tmXBLabelFont = 25
;  res@lbLabelFont   = 25

res@tmXBOn = False
res@tmXTOn = False
res@tmYLOn = False
res@tmYROn = False

 res@mpOutlineBoundarySets = "allBoundaries"    ; divisao dos paises
 ;res@cnFillMode       = "RasterFill"   	 ; nao interpolar..
 res@cnFillMode       =  fillmode	   	 ; nao interpolar..

res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed
; res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks


;res@gsnAddCyclic           = False              ; regional data
res@gsnAddCyclic           = True              ; regional data
res@mpDataBaseVersion     = "MediumRes"
;*****************define scala da legenda***************************

 ;res@cnLevelSelectionMode = "ExplicitLevels"
 ;res@cnLevels             = (/1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0/)

;****************define scala lat/lon******************************

res@mpMaxLatF                   = 7.5
res@mpMinLatF                   = -33.5
res@mpMaxLonF                   = -35.5
res@mpMinLonF                   = -80.5

;************** cor no fundo************************************

;  getvalues wks
;    "wkColorMap" : cmap
;  end getvalues
;  cmap(0,:) = (/1.,1.,1./)
;  cmap(1,:) = (/0.,0.,0./)
;  gsn_define_colormap(wks, cmap)

;*****************titulo e unidade***********************************

; res@tiMainString      = title
res@gsnLeftString      = ""
res@gsnRightString      = ""
;PRE1@long_name = " "
;PRE1@units     = " "

res@mpDataSetName = "Earth..4"
res@mpOutlineSpecifiers = (/"Brazil:states"/)
res@mpDataBaseVersion = "MediumRes"
res@mpOutlineBoundarySets = "National"
res@mpFillOn = False

;****************************legenda*********************************

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = minval
res@cnMaxLevelValF = maxval
res@cnLevelSpacingF = spacing

res@gsnRightStringFontHeightF = 0.04
res@gsnLeftStringFontHeightF = res@gsnRightStringFontHeightF
res@gsnCenterStringFontHeightF = res@gsnRightStringFontHeightF
 
;t = -1
;do year = syear,eyear
;	t = t+1
do t = 0, nyears-1 
	do i = 0,nusepft-1
		pft = usepfts(i)
		year = years(t)
		;res@gsnLeftString = names(usepfts(i))
		res@gsnCenterString = year + "~C~" + str_get_field(names(usepfts(i)),1,"_")
		;if (i.eq.0) then
		;	res@gsnCenterString = year + "~C~" + names(usepfts(i))
		;else
		;	res@gsnCenterString = ""

		;end if
		plots(i,t) = gsn_csm_contour_map(wks,varinp({year},pft,:,:),res)
	end do
end do
;c = 0
;t = -1
;do year = syear,eyear
;	t = t+1
;	do i = 0,nusepft-1
;		pft = usepfts(i)
;		;res@gsnLeftString = names(usepfts(i))
;		res@gsnCenterString = year + "~C~" + names(usepfts(i))
;		;if (i.eq.0) then
;		;	res@gsnCenterString = year + "~C~" + names(usepfts(i))
;		;else
;		;	res@gsnCenterString = ""
;
;		;end if
;		onedplots(c) = gsn_csm_contour_map(wks,varinp({year},pft,:,:),res)
;		print("|year:"+year+"|pft:"+pft+"|"+avg(varinp({year},pft,{-11.80},{360-55.5})))
;		c = c+1
;	end do
;end do

 
 ;do i = 0,npft-1
;  res@gsnRightString = i
;  res@gsnLeftString = "INP"
;  plotsinp(i) = gsn_csm_contour_map(wks,varinp(i,:,:),res)
;  res@gsnLeftString = "REF"
;  plotsref(i) = gsn_csm_contour_map(wks,varref(i,:,:),res)
;end do
;
;plotsl1 = new(16,graphic)
;plotsl2 = new(16,graphic)
;do i = 0,7
;  plotsl1(i) = plotsref(i)
;  plotsl1(i+8) = plotsinp(i)
;  plotsl2(i) = plotsref(i+8)
;  plotsl2(i+8) = plotsinp(i+8)
;end do
;
;respref = True
;; respref@gsnDraw = False
;respref@gsnFrame = False
;; respref@gsnPanelLeft = 0.5
;respref@gsnPanelBottom = 0.5
;respref@gsnPanelLabelBar = True
;respref@lbTitleString = year
;respref@lbTitlePosition = "Bottom"
;
;respinp = respref
;delete(respinp@gsnPanelLeft)
;delete(respinp@gsnPanelBottom)
;; respinp@gsnPanelRight = 0.5
;respinp@gsnPanelTop = 0.5
;
; plotsperpage = 4
; gsn_panel(wks,plotsinp,(/plotsperpage,1/),respinp)
; gsn_panel(wks,plotsref,(/plotsperpage,1/),respref)
;gsn_panel(wks,plotsl1,(/2,8/),respref)
;gsn_panel(wks,plotsl2,(/2,8/),respinp)

resp = True
; resp@gsnDraw = False
resp@gsnFrame = False

resp@gsnPanelMainString = rcpstr + " " + scenariostr

;resp@gsnPanelRight = 0.9
;resp@gsnPanelBottom = 0.1
resp@gsnPanelLabelBar = True
;resp@lbTitleString = year
resp@lbTitlePosition = "Bottom"
;resp@lbOrientation = "Vertical"
resp@lbJustification = "TopCenter"
resp@lbTitleString = "PCT_PFT"

resp@gsnMaximize = True
gsn_panel(wks,ndtooned(plots),(/nusepft,nyears/),resp)
;gsn_panel(wks,onedplots,(/nyears,nusepft/),resp)
;onedplots = ndtooned(plots)
;gsn_panel(wks,onedplots,(/nusepft,nyears/),resp)
;
;pos = new((/nusepft,2/),float)
;
;do i = 0,nusepft-1
;	pos(i,:) = get_title_position(plots(i,0))
;end do
;print(pos)
;print(get_title_position(onedplots(3)))
;txres               = True
;txres@txFontHeightF = 0.012
;txres@txAngleF      = 90.
;
;do i = 0,nusepft-1
;	gsn_text_ndc(wks,names(usepfts(i)),   pos(i,0), pos(i,1), txres)
;end do


;drawNDCGrid(wks)

; maximize_output(wks,False)
frame(wks)

end
