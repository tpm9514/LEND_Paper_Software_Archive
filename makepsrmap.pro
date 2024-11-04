; FUNCTION SEGMENTREGIONS
;   This code segments permanently shadowed regions from binary PSR map
;   Image has 0 values where not PSR, 1 = PSR. Method is a region growing
;   traverse.  Assigns a unique numeric identifier to all pixel clusters, e.g. in each
;   PSR, minarea is the minimum area PSR to report in pixels
;   Method is developed in Sonka et al., 1999 
;          
;   input = binary PSR map from 0% illumination
;   output = PSR map with spatially independent pixel aggregates.  each pixel aggregate has a unique id
;  
Function SegmentRegions,im, minarea

  dimx = n_elements(im(*,0))-1
  dimy = n_elements(im(0,*))-1
  imstore = float(im)
  imstore(*,*) = 0.

  xchk = [-1, 0, 1, 0,-1, 1,-1, 1]
  ychk = [ 0, 1, 0,-1,-1, 1, 1,-1]

  id = long(1)

  for a = 1, dimx-2 do begin
    for b = 1, dimy-2 do begin
      if imstore(a,b) eq 0 and im(a,b) eq 1 then begin
        xlist = [long(a)]
        ylist = [long(b)]
        imstore(a,b) = id
        strt = long(0)
        xcnt = long(0)
        lend = long(1)
        while strt lt lend do begin
          xloc = xlist[strt]
          yloc = ylist[strt]
          for k = 0,7 do begin
            i = xloc+xchk(k)
            j = yloc+ychk(k)
            if i ge 0 and i le dimx and j ge 0 and j le dimy then begin
              if imstore(i,j) eq 0 and im(i,j) eq 1 then begin
                xlist = [xlist, i]
                ylist = [ylist, j]
                imstore(i,j) = id
                lend = lend + 1
              endif
            endif
          endfor
          strt = strt + 1
        endwhile

        if n_elements(xlist) ge MinArea then begin
          id = id + 1
        endif else begin
          imstore(xlist,ylist) = -1
        endelse
      endif

    endfor
    if a mod 400 eq 0 then print, a
  endfor

  print, 'Segmented Spots: ',max(imstore)  
  return, imstore
end



Pro MakeMaps

@~/merrimac/LENDproc/LEND_Paper/Software/dirsets

common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

restore,SaveSets + 'OceanHaline.sav'

device,retain=2
device,decomposed=0
COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr
loadct, 33
r_curr(255) = 255
g_curr(255) = 255
b_curr(255) = 255
r_curr(0) = 0
tvlct, r_curr,g_curr,b_curr
; Read in LRO maps
;
illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
SlopeAzi = read_tiff(LRO_Maps+'Slope_Azi_Angle.tiff')
Slope = read_tiff(LRO_Maps+'TopoSlope.tiff')

; Diviner Max Temp
;
restore,SaveSets+'SPDiviner_HiRes.sav'
tmax = tempmax(50:7249,50:7249)
tmax = congrid(tmax, 1440,1440,/center)
dt = where(latitgrid le -80)

adim = 1440

illum = congrid(illum,1440,1440,/center)
latitgrid = congrid(latitgrid,1440,1440,/center)
longigrid = congrid(longigrid, 1440,1440, /center)

dem = congrid(dem,1440,1440,/center)
slopeazi = congrid(slopeazi, 1440,1440,/center)
slope = congrid(slope, 1440,1440,/center)
tempmax = congrid(tmax, 1440, 1440,/center)
illum = congrid(illum, 1440,1440, /center)

dt = where(illum eq 0 and latitgrid lt -80.)

ill1 = illum
ill1(*) = 0
ill1(dt) = 1

; The line below calls the SegmentRegions function to generate the 
; independent PSRs based on pixels adjacency.  The line below
; displays the segmented PSRs and their ID's
; 
; Input = binary PSR map, minimum PSR area in pixels, 1 is single pixel resolution
illum1 = segmentregions(ill1,1)
tvscl, illum1

; Example code to list PSR areas and diameters
; Total numbers of PSRs, PSR Count begins at 1, PSR 0 = non-PSR
print, 'Total Numbers of PSRs > 2km pixels: Max Illum1', max(illum1)
psrarr = fltarr(max(illum1)+1)
psrdiam = psrarr
print, 'PSR ID, PSR Mean Internal Latitude (deg), PSR Mean Internal Longitude (deg), PSR Area (km2), PSR diameter (km)'
print, '--------------------------------------------------------------------------------------------------------------'
for i = 1,max(illum1) do begin
  pid = where(illum1 eq i,area)
  psrarr(i) = n_elements(pid)*4.             ; multiply by 4, each pixel = 4 km2
  psrdiam(i) = sqrt(psrarr(i)/!pi)*2.
; print out the PSR id, Mean Latit, Mean Longi, area, diameter

  print, i, mean(latitgrid(pid)), mean(longigrid(pid)), psrarr(i), psrdiam(i)

endfor

hist_PSR = histogram(psrdiam, binsize=2.5, min=0.,max=40.)
psrd = plot(findgen(17)*2.5 + 1.25, hist_psr,title='PSR Diameter Distribution, > 80 S',xtitle='PSR Diameter',yrange=[0.9,250], ytitle='Count')
psrd.save,Figures+'PSR_Diameter_Dist.png'

print, '-------------------'
print, 'PSR Diam (Mean, StdDev)', mean(psrdiam(1:468)),stddev(psrdiam(1:468))

read,xxx

psrs = where(illum1 gt 0.)
lowlat = where(latitgrid ge -80)
ill = illum1
ill(*) = 255
ill(psrs) = 20
lb = where(latitgrid gt -80)
latit = latitgrid
latit(*) = 0
k = fltarr(5,5)
k(*) = 1
k(4,0) = 0
k(0,4) = 0
k(4,4) = 0
k(0,0) = 0
latit(lb) = 1
latit1 = erode(latit,k)
latit = latit-latit1
dt = where(latit eq 1)
ill(dt) = 0
window,xsize=1000,ysize=1000

keep = where(latitgrid le -80.)
kill = where(latitgrid gt -80.)

tmax(kill) = max(tmax(keep))
dem(kill) = max(dem(keep))
slopeazi(kill) = max(slopeazi(keep))

; The line below displays the binary illumination map
; PSRs are blue, non-PSR = white
print, 'PSR = 0% Illumination Map > 80 S, units = 0 to 100'
tvscl, ill
read,xxx

print, 'Diviner Max Temperature Map > 80 S, units = K'
tvscl, tmax
read,xxx

print, 'LOLA topography > 80 S, units = altitude as km deviation from 1737.4 km average lunar radius'
print, 'Altitude, mean, stddev:', mean(dem(keep)), stddev(dem(keep))
tvscl, dem
read,xxx

print, 'SLope Azimuth Angle > 80 S, units = degrees, Poleward-facing = 0, equator-facing = 180, East = West = 90'
tvscl, slopeazi
read,xxx


return
end
