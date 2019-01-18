subroutine get_pixel_weights(wgt,outi,outj,inplats,inplatn,inplonw,inplone, &
  inpnlat,inpnlon,outlats,outlatn,outlonw,outlone,outnlat,outnlon)
  implicit none
  integer inpnlat,inpnlon,outnlat,outnlon
  integer outi,outj
  !real*8 inplats(inpnlat,inpnlon),inplatn(inpnlat,inpnlon),inplonw(inpnlat,inpnlon),inplone(inpnlat,inpnlon)
  !real*8 outlats(outnlat,outnlon),outlatn(outnlat,outnlon),outlonw(outnlat,outnlon),outlone(outnlat,outnlon)
  real*8 wgt(inpnlon,inpnlat)

  real*8 inplats(inpnlon,inpnlat),inplatn(inpnlon,inpnlat),inplonw(inpnlon,inpnlat),inplone(inpnlon,inpnlat)
  real*8 outlats(outnlon,outnlat),outlatn(outnlon,outnlat),outlonw(outnlon,outnlat),outlone(outnlon,outnlat)

  ! The point boundaries
  real*8 tlats,tlatn,tlonw,tlone
  ! Dimensional weights
  real*8 latwgt, lonwgt

  ! Cell size
  real*8 latsize, lonsize


  ! Counters
  integer i,j,ii,jj

  tlats = outlats(outj,outi)
  tlatn = outlatn(outj,outi)
  tlone = outlone(outj,outi)
  tlonw = outlonw(outj,outi)
  ! write(*,*) tlats, tlatn, tlonw, tlone

  ! i = 5
  ! j = 720
  !
  ! write(*,*) "inplats(j,i) = ",inplats(j,i)
  ! write(*,*) "inplatn(j,i) = ",inplatn(j,i)
  ! write(*,*) "inplonw(j,i) = ",inplonw(j,i)
  ! write(*,*) "inplone(j,i) = ",inplone(j,i)
  ! stop
  wgt(:,:) = 0.0d0

  do j = 1,inpnlon
    do i = 1,inpnlat
      ! write(*,*) i,j
      ! Latitudes are inside the box
      if ( (inplats(j,i) .ge. tlats .and. inplats(j,i) .le. tlatn) .or. &
      (inplatn(j,i) .ge. tlats .and. inplatn(j,i) .le. tlatn)) then
        ! Longitudes are inside the box
        if ( (inplonw(j,i) .ge. tlonw .and. inplonw(j,i) .le. tlone) .or. &
        (inplone(j,i) .ge. tlonw .and. inplone(j,i) .le. tlone) ) then

          ! Check if it's fully inside the box (lat)
          if ( inplats(j,i) .ge. tlats .and. inplatn(j,i) .le. tlatn ) then
            latwgt = 1.0d0
          else ! It's ONLY partially inside (LAT)
            latsize = inplatn(j,i)-inplats(j,i)
            if (inplatn(j,i) .ge. tlats .and. inplatn(j,i) .le. tlatn) then
              latwgt = (inplatn(j,i)-tlats)/latsize
            else
              latwgt = (tlatn - inplats(j,i))/latsize
            end if
          end if

          ! Check if it's fully inside the box (LON)
          if ( inplonw(j,i) .ge. tlonw .and. inplone(j,i) .le. tlone ) then
            lonwgt = 1.0d0
          else ! It's ONLY partially inside (LON)
            lonsize = inplone(j,i)-inplonw(j,i)
            if (inplone(j,i) .ge. tlonw .and. inplone(j,i) .le. tlone) then
              lonwgt = (inplone(j,i)-tlonw)/lonsize
            else
              lonwgt = (tlone - inplonw(j,i))/lonsize
            end if
          end if

          ! write(*,*) i,j
          ! write(*,*) "inplats(j,i) = ",inplats(j,i)
          ! write(*,*) "inplatn(j,i) = ",inplatn(j,i)
          ! write(*,*) "inplonw(j,i) = ",inplonw(j,i)
          ! write(*,*) "inplone(j,i) = ",inplone(j,i)
          ! write(*,*) "latwgt = ",latwgt
          ! write(*,*) "lonwgt = ",lonwgt
          ! write(*,*) "latwgt*lonwgt = ",latwgt*lonwgt


          wgt(j,i) = latwgt*lonwgt

    end if ! lon part of box
  end if ! lat part of box
end do !i, lats
end do !j, lons

! write(*,*) maxval(wgt)

return

end subroutine
