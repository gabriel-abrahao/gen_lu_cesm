function get_pixel_weights(inplats,inplatn,inplonw,inplone,inpnlat,inpnlon,outlats,outlatn,outlonw,outlone,outnlat,outnlon)
  integer inpnlat,inpnlon,outnlat,outnlon
  !real*8 inplats(inpnlat,inpnlon),inplatn(inpnlat,inpnlon),inplonw(inpnlat,inpnlon),inplone(inpnlat,inpnlon)
  !real*8 outlats(outnlat,outnlon),outlatn(outnlat,outnlon),outlonw(outnlat,outnlon),outlone(outnlat,outnlon)
  real*8 inplats(inpnlon,inpnlat),inplatn(inpnlon,inpnlat),inplonw(inpnlon,inpnlat),inplone(inpnlon,inpnlat)
  real*8 outlats(outnlon,outnlat),outlatn(outnlon,outnlat),outlonw(outnlon,outnlat),outlone(outnlon,outnlat)
!  integer poi(2,2)
  !poi(1,1) = 9876
  !poi(1,2) = 9876
  !poi(2,1) = 9876
  !poi(2,2) = 9876
  !get_pixel_weights = poi
  return
end function
