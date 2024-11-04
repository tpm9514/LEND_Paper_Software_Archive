; Program:   GenSlopeAziAngleMap.pro
; Comment:   The program generates the south polar slope azimuth angle map used in the study
;
; Input:  South Polar Topography in polar stereographic format
; Output:  Generates Slope Azimuth Angle map
;
; Written by:   Timothy P. McClanahan, NASA GSFC              March 14, 2024
; Contact:   timothy.p.mcclanahan@nasa.gov
;
; The software is provided as is.   It is part of a data and software repository supporting a
; Planetary Science Journal publication:
; “Evidence for Widespread Hydrogen Sequestration within the Moon’s South Polar Cold Traps”

; The software is designated by NASA as a Creative Commons Zero repository.   NASA’s Software
; Release #:  GSC-19191-1.

; Repository  DOI:  10.5281/zenodo.10027812

; The software written in the Interactive Data Language 8.8.3, which requires a
; commercial license to run available at:   https://www.nv5geospatialsoftware.com/Products/IDL



; FUNCTION DOT
;   Derive dot product from two points
Function Dot, scx,scy, snx,sny
  v1 = [scx, scy]
  v2 = [snx, sny]
  theta = acos(total((v1 / sqrt(total(v1^2)))*(v2 / sqrt(total(v2^2)))))
  return,theta/!dtor
end


; FUNCTION MAG
;   Derive the magnitude of a vector
function mag,v1,v2
  return, sqrt(v1^2.+v2^2.)
end


Pro GenSlopeAziAngleMap

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets


; Read in 7200 x 7200 pixel maps

  latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
  longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
  dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
  SlopeAziAngle1 = read_tiff(LRO_Maps+'Slope_Azi_Angle.tiff')

px = 3600.5
py = px
center = [px,py]
k = fltarr(3)
k = 1./[sqrt(2.), 1., sqrt(2.)]


SlopeAziAngle = dem
SlopeAziAngle(*) = 0.
dim = 7200.
EW_Angle = SlopeAziAngle
EWVec = fltarr(2)
Svec = fltarr(2)

X_grad = dem
X_grad(*) = 0.
Y_grad = X_grad

print, ' ' 
print, "Gen:  X and Y gradient maps'
for i = long(1),dim-2 do begin
  for j = long(1),dim-2 do begin
    X_Grad(i,j) = mean(dem(i-1,j-1:j+1)*k) - mean(dem(i+1,j-1:j+1)*k)
    Y_Grad(i,j) = mean(dem(i-1:i+1,j-1)*k) - mean(dem(i-1:i+1,j+1)*k)
  endfor
  if i mod 400 eq 0 then print,"Working pixel row: ",i
endfor
Grad = sqrt(x_grad^2.+y_Grad^2.)

print, ' '
print, 'Generating Slope Azi maps from x,y, gradient maps'

for i = 2,dim-3 do begin
  for j = 2,dim-3 do begin
    ii = float(i)
    jj = float(j)
    pvec = center - [ii,jj]
    pvec = pvec / sqrt(total(pvec^2.))
; Following 3 lines makes east west map
    ewvec(0) = pvec(1)*(-1.)
    ewvec(1) = pvec(0)
;    svec(*) = grad(i,j,1:2)
    svec(0) = x_grad(i,j)
    svec(1) = y_grad(i,j)
    svec = svec / sqrt(total(svec^2.))
    SlopeAziAngle(i,j) = acos(total(pvec*svec))*(180./!pi)
    EW_Angle(i,j) = acos(total(ewvec*svec))*(180./!pi)
; test
;    chk = acos(total(pvec*ewvec))*(180./!pi)
;    if chk ne 90. then read,xxx
  endfor
 if i mod 400 eq 0 then print, 'Working row, i = ',i
endfor

read,xxx
return
end