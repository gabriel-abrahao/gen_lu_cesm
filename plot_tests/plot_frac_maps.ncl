load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

inpfile       =	"../out_surfdata/surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_rochedo_weg_c130524.nc"
reffile       = "../input/surfdata.pftdyn_0.9x1.25_rcp8.5_simyr1850-2100_c130524.nc"
varname       =	"PCT_NAT_PFT"
year          = 2050

foutname      =	"compare_frac_" + year

minval		=	0.
maxval		=	100.
spacing		=	10.

; cmap		=	"BlueRed"
cmap		=	"matlab_jet"

fillmode	=	"RasterFill"
;fillmode	=	"AreaFill"

arqinp = addfile(inpfile,"r")
arqref = addfile(reffile,"r")
varinp = arqinp->$varname$({year},:,:,:)
varref = arqref->$varname$({year},:,:,:)

varref&lsmlat = arqref->LATIXY(:,0)
varref&lsmlon = arqref->LONGXY(0,:)

varref&lsmlat@units = "degrees_north"
varref&lsmlon@units = "degrees_east"

copy_VarCoords(varref,varinp)

dims = dimsizes(varinp)
npft = dims(0)
printVarSummary(varinp)
printVarSummary(varref)
; exit

wks_type = "pdf"
wks_type@wkOrientation = "Landscape"
; wks_type@wkPaperWidthF  =  8.5  ; in inches
; wks_type@wkPaperHeightF = 100.0  ; in inches
wks = gsn_open_wks(wks_type,foutname)

gsn_define_colormap(wks,cmap)

plotsinp = new(npft,graphic)
plotsref = new(npft,graphic)

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


res@gsnAddCyclic           = False              ; regional data
res@mpDataBaseVersion     = "MediumRes"
;*****************define scala da legenda***************************

 ;res@cnLevelSelectionMode = "ExplicitLevels"
 ;res@cnLevels             = (/1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0/)

;****************define scala lat/lon******************************

res@mpMaxLatF                   = 13.5
res@mpMinLatF                   = -54.5
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

res@gsnRightStringFontHeightF = 0.03
res@gsnLeftStringFontHeightF = res@gsnRightStringFontHeightF
do i = 0,npft-1
  res@gsnRightString = i
  res@gsnLeftString = "INP"
  plotsinp(i) = gsn_csm_contour_map(wks,varinp(i,:,:),res)
  res@gsnLeftString = "REF"
  plotsref(i) = gsn_csm_contour_map(wks,varref(i,:,:),res)
end do
plotsl1 = new(16,graphic)
plotsl2 = new(16,graphic)
do i = 0,7
  plotsl1(i) = plotsref(i)
  plotsl1(i+8) = plotsinp(i)
  plotsl2(i) = plotsref(i+8)
  plotsl2(i+8) = plotsinp(i+8)
end do

respref = True
; respref@gsnDraw = False
respref@gsnFrame = False
; respref@gsnPanelLeft = 0.5
respref@gsnPanelBottom = 0.5
respref@gsnPanelLabelBar = True
respref@lbTitleString = year
respref@lbTitlePosition = "Bottom"

respinp = respref
delete(respinp@gsnPanelLeft)
delete(respinp@gsnPanelBottom)
; respinp@gsnPanelRight = 0.5
respinp@gsnPanelTop = 0.5

; plotsperpage = 4
; gsn_panel(wks,plotsinp,(/plotsperpage,1/),respinp)
; gsn_panel(wks,plotsref,(/plotsperpage,1/),respref)
gsn_panel(wks,plotsl1,(/2,8/),respref)
gsn_panel(wks,plotsl2,(/2,8/),respinp)



; maximize_output(wks,False)
frame(wks)

end
