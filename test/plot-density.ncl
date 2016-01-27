;***********************************************
; plot-density.ncl
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

;************************************************
; read in netCDF file
;************************************************
  density_file = addfile("pgbhnl.gdas.all_density.nc","r")

  density = density_file->density

;************************************************
; plot resources [options]
;************************************************
  plot = new(1,graphic)                ; define a panel plot

; Resources:
  res                     = True

  res@gsnDraw             = False
  res@gsnFrame            = False

  res@vpKeepAspect        = True
  res@tiMainFontHeightF   = 0.02
  res@vpWidthF            = 0.8
  res@vpHeightF           = 0.4
;  res@gsnMaximize         = True
;  res@gsnPaperOrientation = "portrait"

  wks = gsn_open_wks("eps","plot-density")   

; choose color table:
;  gsn_define_colormap(wks,"cosam12")             ; not enough colors
  gsn_define_colormap(wks,"WhBlGrYeRe")
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; favorite right now
;  gsn_define_colormap(wks,"BlGrYeOrReVi200")
;  gsn_define_colormap(wks,"BlueRed")              ; mid-light blue to red
;  gsn_define_colormap(wks,"BlRe")                ; too dark
;  gsn_define_colormap(wks,"wh-bl-gr-ye-re")
;  gsn_define_colormap(wks,"WhBlGrYeRe")
  gsn_merge_colormaps(wks,"gui_default","BlueDarkRed18")


  res@cnFillOn          = True
  res@gsnSpreadColors   = True

  res@cnLinesOn         = False
  res@cnLineLabelsOn    = False
  res@lbLabelAutoStride = True
  res@cnInfoLabelOn     = False
  res@tmXBLabelsOn      = True
  res@tmXBOn            = True

  res@tmXTLabelsOn      = False
  res@tmXTOn            = False

  res@gsnContourZeroLineThicknessF = 0
  res@gsnContourNegLineDashPattern = 6

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 10.0
  res@cnMinLevelValF = 0.0
  res@cnMaxLevelValF = 160.0

  res@lbLabelBarOn = True
  res@lbLabelFontHeightF = 0.024

  res@tiMainFontHeightF = 0.024

  res@gsnSpreadColorStart = 2
  res@gsnSpreadColorEnd = 23

  res@tiMainString = "Warm-Core Cyclone Density from CFSR"
  plot(0) = gsn_csm_contour_map(wks, density, res)

  pres = True
  gsn_panel(wks, plot, (/1,1/), pres)

end
