; PRO WEH_n_MaxTemp_vs_PSRPredict
;
; Written by:   Timothy P. McClanahan, NASA GSFC              October 31, 2024
; Contact:   timothy.p.mcclanahan@nasa.gov

;
; Input:  Collimated WEH, Topography, Binary PSR map, Maximum Temperature
;
; The software is provided as is. It is part of a data and software repository supporting a
;
; Planetary Science Journal publication:
; “Evidence for Widespread Hydrogen Sequestration within the Moon’s South Polar Cold Traps”

; The software is designated by NASA as a Creative Commons Zero repository.   NASA’s Software
; Release #:  GSC-19191-1.

; Repository  DOI:  10.5281/zenodo.10027812
;

; The software written in the Interactive Data Language 8.8.3, which requires a
; commercial license to run available at:   https://www.nv5geospatialsoftware.com/Products/IDL


; KS test:  Kolmogorov Smirnoff test of two populations CDFs
; Test assumes the populations are not normally distributed.
;
pro prob_ks, D, N_eff, probks
  ;+
  ; NAME:
  ;       PROB_KS
  ;
  ; Source:   NASA Goddard Space Flight Center Astronomy Users Library
  ;    https://asd.gsfc.nasa.gov/archive/idlastro/ftp/pro/math/prob_ks.pro
  ;
  ; PURPOSE:
  ;       Return the significance of the Kolmogoroff-Smirnov statistic
  ; EXPLANATION:
  ;       Returns the significance level of an observed value of the
  ;       Kolmogorov-Smirnov statistic D for an effective number of data points
  ;       N_eff.   Called by KSONE and KSTWO
  ;
  ; CALLING SEQUENCE:
  ;       prob_ks, D, N_eff, probks
  ;
  ; INPUT PARAMETERS:
  ;       D -  Kolmogorov statistic, floating scalar, always non-negative
  ;       N_eff - Effective number of data points, scalar.   For a 2 sided test
  ;               this is given by (N1*N2)/(N1+N2) where N1 and N2 are the number
  ;               of points in each data set.
  ;
  ; OUTPUT PARAMETERS:
  ;       probks - floating scalar between 0 and 1 giving the significance level of
  ;               the K-S statistic.   Small values of PROB suggest that the
  ;               distribution being tested are not the same
  ;
  ; REVISION HISTORY:
  ;       Written     W. Landsman                August, 1992
  ;       Corrected typo (termbv for termbf)    H. Ebeling/W.Landsman  March 1996
  ;       Probably did not affect numeric result, but iteration went longer
  ;       than necessary
  ;       Converted to IDL V5.0   W. Landsman   September 1997
  ;-
  On_error,2

  if N_params() LT 3 then begin
    print,'Syntax - prob_ks, D, N_eff, prob'
    print,'  D - Komolgorov-Smirnov statistic, input'
    print,'  N_eff - effective number of data points, input'
    print,'  prob - Significance level of D, output'
    return
  endif

  eps1 = 0.001    ;Stop if current term less than EPS1 times previous term
  eps2 = 1.e-8    ;Stop if current term changes output by factor less than EPS2

  en = sqrt( N_eff )
  lambda = (en + 0.12 + 0.11/en)*D

  a2 = -2.*lambda^2
  probks = 0.
  termbf = 0.
  sign = 1.

  for j = 1,100 do begin

    term = sign*2*exp(a2*j^2)
    probks = probks + term

    if ( abs(term) LE eps1*termbf ) or $
      ( abs(term) LE eps2*probks ) then return

    sign = -sign                  ;Series alternates in sign
    termbf = abs(term)
  endfor

  probks = 1.          ;Sum did not converge after 100 iterations
  return

end


pro kstwo, data1, data2, D, prob
  ;+
  ; NAME:
  ;       KSTWO
  ; 
  ; Source:   NASA Goddard Space Flight Center Astronomy Users Library
  ;    https://asd.gsfc.nasa.gov/archive/idlastro/ftp/pro/math/kstwo.pro
  ; 
  ; PURPOSE:
  ;       Return the two-sided Kolmogorov-Smirnov statistic
  ; EXPLANATION:
  ;       Returns the Kolmogorov-Smirnov statistic and associated probability
  ;       that two arrays of data values are drawn from the same distribution
  ;       Algorithm taken from procedure of the same name in "Numerical
  ;       Recipes" by Press et al., 2nd edition (1992), Chapter 14
  ;
  ; CALLING SEQUENCE:
  ;       kstwo, data1, data2, D, prob
  ;
  ; INPUT PARAMETERS:
  ;       data1 -  vector of data values, at least 4 data values must be included
  ;               for the K-S statistic to be meaningful
  ;       data2 -  second set of data values, does not need to have the same
  ;               number of elements as data1
  ;
  ; OUTPUT PARAMETERS:
  ;       D - floating scalar giving the Kolmogorov-Smirnov statistic.   It
  ;               specifies the maximum deviation between the cumulative
  ;               distribution of the data and the supplied function
  ;       prob - floating scalar between 0 and 1 giving the significance level of
  ;               the K-S statistic.   Small values of PROB show that the
  ;               cumulative distribution function of DATA1 is significantly
  ;               different from DATA2
  ;
  ; EXAMPLE:
  ;       Test whether two vectors created by the RANDOMN function likely came
  ;       from the same distribution
  ;
  ;       IDL> data1 = randomn(seed,40)        ;Create data vectors to be
  ;       IDL> data2 = randomn(seed,70)        ;compared
  ;       IDL> kstwo, data1, data2, D, prob   & print,D,prob
  ;
  ; PROCEDURE CALLS
  ;       procedure PROB_KS - computes significance of K-S distribution
  ;
  ; REVISION HISTORY:
  ;       Written     W. Landsman                August, 1992
  ;       FP computation of N_eff      H. Ebeling/W. Landsman  March 1996
  ;       Fix for arrays containing equal values J. Ballet/W. Landsman Oct. 2001
  ;       Fix index when maximum difference is at array end Renbin Yan  Dec 2008
  ;       Handle large number when computing N_err  D. Schnitzeler/WL  Sep 2010
  ;-
  On_error, 2
  compile_opt idl2

  if ( N_params() LT 4 ) then begin
    print,'Syntax - KSTWO, data1, data2, d, prob'
    return
  endif

  n1 = N_elements( data1 )
  if ( N1 LE 3 ) then message, $
    'ERROR - Input data values (first param) must contain at least 4 values'

  n2 = N_elements( data2 )
  if ( n2 LE 3 ) then message, $
    'ERROR - Input data values (second param) must contain at least 4 values'

  sortdata1 = data1[ sort( data1 ) ]        ;Sort input arrays into
  sortdata2 = data2[ sort( data2 ) ]        ;ascending order

  fn1 = ( findgen( n1 +1 )  ) / n1          ;updated Dec 2008
  fn2 = ( findgen( n2 +1)  ) / n2

  j1 = 0l & j2 = 0l
  id1 = lonarr(n1+n2)  & id2 = id1
  i = 0l

  ; Form the two cumulative distribution functions, marking points where one
  ; must test their difference

  while ( j1 LT N1 ) and ( j2 LT n2 ) do begin

    d1 = sortdata1[j1]
    d2 = sortdata2[j2]
    if d1 LE d2 then j1 = j1 +1
    if d2 LE d1 then j2 = j2 +1

    id1[i] = j1   & id2[i] = j2
    i = i+1

  endwhile

  id1 = id1[0:i-1]   &  id2 = id2[0:i-1]

  ; The K-S statistic D is the maximum difference between the two distribution
  ; functions

  D = max( abs( fn1[id1] - fn2[id2] ) )
  N_eff =  long64(n1)*n2/ float(n1 + n2)           ;Effective # of data points
  PROB_KS, D, N_eff, prob                ;Compute significance of statistic

  return
end


; FUNCTION SEGMENTREGIONS
;   This code segments permanently shadowed regions from binary PSR map
;   Image has 0 values where not PSR, 1 = PSR.  Method is a region growing
;   traverse. Assigns a unique numeric identifier to all pixels in each
;   PSR.
Function SegmentRegions,image, minarea

  dimx = n_elements(image(*,0))-1
  dimy = n_elements(image(0,*))-1
  imstore = float(image)
  imstore(*,*) = 0.

  xchk = [-1, 0, 1, 0,-1, 1,-1, 1]
  ychk = [ 0, 1, 0,-1,-1, 1, 1,-1]

  id = long(1)

  for a = 1, dimx-2 do begin
    for b = 1, dimy-2 do begin
      if imstore(a,b) eq 0 and image(a,b) eq 1 then begin
        xlist = [long(a)]
        ylist = [long(b)]
        imstore(a,b) = id
        strt = long(0)
        xcnt = long(0)
        lend = long(1)
        while strt lt lend do begin
          xloc = xlist(strt)
          yloc = ylist(strt)
          for k = 0,7 do begin
            i = xloc+xchk(k)
            j = yloc+ychk(k)
            if i ge 0 and i le dimx and j ge 0 and j le dimy then begin
              if imstore(i,j) eq 0 and image(i,j) eq 1 then begin
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
        endif else imstore(xlist,ylist) = -1
      endif

    endfor
    if a mod 400 eq 0 then print, a
  endfor

  return, imstore
end




; PRO WEH_n_MaxTemp_vs_PSRPredict
; This code was produced to provide a predict capability based on the PSR areal density within CSETN's 30 km FOV.
; The code convolves a 30 km disk
; 
; Date:     October 30, 2024
; Author:   Tim McClanahan
; 
; 
; Inputs:  WEH Maps
;          Detection Significance Maps (WEH)
;          LRO Maps:
;             Illumination for PSR maps
;             Topo
;             Slope
;             Slope Azi angle
;             Max Temp
;
Pro WEH_n_MaxTemp_vs_PSRPredict

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
  loadct, 33
  r_curr(255) = 255
  g_curr(255) = 255
  b_curr(255) = 255
  r_curr(0) = 0
  tvlct, r_curr,g_curr,b_curr

  ; Read in LRO Maps
  ;
  illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
  latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
  longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
  dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
  TempMax = read_tiff(LRO_Maps+'SP_Diviner_TempMax.tiff')
  SlopeAziAngle = read_tiff(LRO_Maps+'Slope_Azi_Angle.tiff')
  Slope = read_tiff(LRO_Maps+'TopoSlope.tiff')

  ; Resize maps
  adim = 1440
  dem1 = congrid(dem,1440,1440,/center)
  longigrid = congrid(longigrid, adim,adim,/center)
  latitgrid = congrid(latitgrid, adim,adim,/center)
  latit = latitgrid
  longi = longigrid
  illum = congrid(illum,1440,1440,/center)
  illum_orig = illum

  ; Analyze all 10 map years if looped
  ;
  k = fltarr(3,3)
  k(*) = 1

  
; Read in WEH maps
  csetn_coll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')
  csetn_coll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_Suppress_sig220_UL4.tiff')
  

; This segment of the code derives Figure 7, PSR, non-PSR averaged WEH vs Latitude
;

ker = gaussian_function([2.5,2.5], width=15)
c = 7
for i = 0,14 do begin
  for j = 0,14 do begin
     if sqrt((i-c)^2. + (j-c)^2.) gt 7. then ker(i,j) = 0
  endfor
endfor
ker = ker / total(ker)

ill= illum
dt = where(illum eq 0)
ill(*) = 0
ill(dt) = 1
temp_avg = TempMax
psr_den = temp_avg

for i = 0,1439 do begin
  for j = 0,1439 do begin
    if latitgrid(i,j) le -70. then begin
       idat = ill(i-7:i+7,j-7:j+7)
       psrpix = where(idat eq 1 and ker gt 0., pcnt)
       psr_den(i,j) = float(pcnt)/149.
       tempdat = tempmax(i-7:i+7,j-7:j+7)
       dt = where(ker gt 0.)
       temp_avg(i,j) = total(tempdat*ker)
    endif
  endfor
endfor

latd = where(latitgrid le -70.)
psr_den1 = psr_den
psr_den1(latd) = psr_den(latd)*100.

lat1 = where(latitgrid le -82 and psr_den1 eq 0)
lat2 = where(latitgrid le -82 and psr_den1 gt 0 and psr_den1 le 10)
lat3 = where(latitgrid le -82 and psr_den1 gt 10 and psr_den1 le 20)
lat4 = where(latitgrid le -82 and psr_den1 gt 20 and psr_den1 le 30)
lat5 = where(latitgrid le -82 and psr_den1 gt 30 and psr_den1 le 40)
lat6 = where(latitgrid le -82 and psr_den1 gt 40 and psr_den1 le 50)
lat7 = where(latitgrid le -82 and psr_den1 gt 50 and psr_den1 le 60)
lat8 = where(latitgrid le -82 and psr_den1 gt 60 and psr_den1 le 70)
lat9 = where(latitgrid le -82 and psr_den1 gt 70 and psr_den1 le 80)
lat10 = where(latitgrid le -82 and psr_den1 gt 80 and psr_den1 le 90)
lat11 = where(latitgrid le -82 and psr_den1 gt 90)

psr_dens = illum
psr_dens(*) = 255
psr_dens(lat1) = 253
psr_dens(lat2) = 225
psr_dens(lat3) = 200
psr_dens(lat4) = 175
psr_dens(lat5) = 150
psr_dens(lat6) = 125
psr_dens(lat7) = 100
psr_dens(lat8) = 75
psr_dens(lat9) = 50
psr_dens(lat10) = 25
psr_dens(lat11) = 5

x = congrid(psr_dens(500:900,500:900), 600,600)
tv, x

b = bytarr(25, 256)
b(*,0:20) = 5
b(*,20:45) = 25
b(*,45:70) = 50
b(*,70:95) = 75
b(*,95:120) = 100
b(*,120:145) = 125
b(*,145:170) = 150
b(*,170:195) = 175
b(*,195:220) = 200
b(*,220:245) = 225
b(*,245:255) = 253
;
;
;  Latitude 82 S WITH PSR DENSITY =0 
;
;
;

latd = where(latitgrid le -80.5 and psr_den1 eq 0.)

pl = plot(temp_avg(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[170, 330], title='WEH Pixels(PSR=0%) vs FOV(AvgMaxTemp): -80.5$\deg$ to -90.0$\deg$', $
  xtitle='Max. Temp. ,K', ytitle='WEH pixels')

xvl = findgen(130)+180.
tempwehcoeff = linfit(temp_avg(latd),csetn_coll_weh(latd))
print, '> 86$\deg$: ',tempwehcoeff
wpall = plot(xvl, xvl*tempwehcoeff(1) + tempwehcoeff(0), color='black', thick=3,name='linfit',/overplot)

zro = fltarr(170)
pfz = plot(findgen(170)+160, zro, thick=2, color='black',linestyle=2,/overplot)

print, 'mean WEH, mean temp: =', mean(csetn_coll_weh(latd)), mean(temp_avg(latd))

print, "MaxTemp vs WEH:", correlate(temp_avg(latd),csetn_coll_weh(latd))

highweh = where(csetn_coll_weh gt 0.3 and psr_den1 eq 0 and latitgrid le -82)

print, "(Lat, Lon) Anomalously high WEH > 0.3 wt% and PSR density = 0% :", mean(latitgrid(highweh)), mean(longigrid(highweh))

read,xxx

;
;
;  Latitude 86 S
;
;
;

latd = where(latitgrid le -86)
latdpsr = where(latitgrid le -86 and ill eq 1.)

pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
          sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
          xrange=[-5, 105], title='FOV(PSR %) vs Pixel(WEH):  -86$\deg$ to -90$\deg$', $
          xtitle='FOV(PSR %)', ytitle='Pixel WEH')
          
xvl = findgen(100)
wehcoef86 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 86$\deg$ WEH Fit Coeffs: ',wehcoef86
wpall = plot(xvl, xvl*wehcoef86(1) + wehcoef86(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -86 and psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor

incavg = plot(ppvals, avgweh, color= 'red', thick=2, /overplot)

plt = plot(psr_den1(latd), temp_avg(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[30., 320.], $
  xrange=[-5, 105], title='FOV(PSR %) vs FOV(Avg. MaxTemp.):  -86$\deg$ to -90$\deg$', $
  xtitle='FOV(PSR %)', ytitle='FOV(Avg. MaxTemp)')

xvl = findgen(100)
tempcoef86 = linfit(psr_den1(latd),temp_avg(latd))
print, '> 86$\deg$: TempFit Coeffs',tempcoef86
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))

tpall = plot(xvl, xvl*tempcoef86(1) + tempcoef86(0), color='black', thick=3,/overplot)

;
;
;  Latitude 83 S
;
;
;

latd = where(latitgrid le -83 and latitgrid gt -86)
latdpsr = where(latitgrid le -83 and latitgrid gt -86 and ill eq 1.)


pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel(WEH): -83$\deg$ to -86$\deg$', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 83$\deg$: WEHCoeffs',wehcoef83
wpall = plot(xvl(0:60), xvl(0:60)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -83 and latitgrid gt -86 and psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor

ppvals(13) = 65
incavg = plot(ppvals(0:13), avgweh(0:13), color= 'red', thick=2, /overplot)

plt = plot(psr_den1(latd), temp_avg(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[30., 320.], $
  xrange=[-5, 105], title='FOV(PSR %) vs FOV(Avg. MaxTemp.): -83 $\leq$ Lat < -86$\deg$', $
  xtitle='FOV(PSR%)', ytitle='FOV(Avg. MaxTemp)')

xvl = findgen(100)
tempcoef83 = linfit(psr_den1(latd),temp_avg(latd))
print, '> 83$\deg$: TempCoeffs',tempcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))


tpall = plot(xvl(0:60), xvl(0:60)*tempcoef83(1) + tempcoef83(0), color='black', thick=3,/overplot)

;
;
;  Latitude 83:   No CABEUS S
;
;
;

latd = where(latitgrid le -83 and latitgrid gt -86 and (longigrid lt 290 or longigrid gt 335))
latdpsr = where(latitgrid le -83 and latitgrid gt -86 and (longigrid lt 290 or longigrid gt 335) and $
                ill eq 1.)


pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel WEH: -83$\deg$ to -86$\deg$: NoCabLons', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 83$\deg$ NO CAB: ',wehcoef83
wpall = plot(xvl(0:60), xvl(0:60)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -83 and latitgrid gt -86 and (longigrid lt 290 or longigrid gt 335) and $
                psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor

ppvals(13) = 65
incavg = plot(ppvals(0:13), avgweh(0:13), color= 'red', thick=2, /overplot)

plt = plot(psr_den1(latd), temp_avg(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[30., 320.], $
  xrange=[-5, 105], title='FOV(PSR %) vs FOV(Avg. MaxTemp.): -83$\deg$ to -86$\deg$: NoCabLons', $
  xtitle='FOV(PSR %)', ytitle='FOV(Avg. MaxTemp)')

xvl = findgen(100)
tempcoef83 = linfit(psr_den1(latd),temp_avg(latd))
print, '> 83$\deg$ NO CAB: ',tempcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))

tpall = plot(xvl(0:60), xvl(0:60)*tempcoef83(1) + tempcoef83(0), color='black', thick=3,/overplot)

;
;
;  Latitude 80 S NO STATION KEEPING
;
;
;



latd = where(latitgrid le -80.5 and latitgrid gt -83 and $
            (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315))
latdpsr = where(latitgrid le -80.5 and latitgrid gt -83 and $
            (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and ill eq 1.)

pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel WEH: -80.5$\deg$ to -83$\deg$, no SK', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 80$\deg$: WEH Coeffs',wehcoef83
print, 'n pixels = ',n_elements(latd)
wpall = plot(xvl(0:50), xvl(0:50)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -80.5 and latitgrid gt -83 and $
            (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and $
             psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor

incavg = plot(ppvals(0:9), avgweh(0:9), color= 'red', thick=2, /overplot)

plt = plot(psr_den1(latd), temp_avg(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[30., 320.], $
  xrange=[-5, 105], title='FOV(PSR %) vs FOV(Avg. MaxTemp.): -80.5$\deg$ to -83$\deg$, no SK', $
  xtitle='FOV(PSR %)', ytitle='FOV(Avg. MaxTemp)')

xvl = findgen(100)
tempcoef83 = linfit(psr_den1(latd),temp_avg(latd))
print, '> 80$\deg$: Temp Coeffs',tempcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))


tpall = plot(xvl(0:50), xvl(0:50)*tempcoef83(1) + tempcoef83(0), color='black', thick=3,/overplot)

;;;;
;;  Latitude 77 S to 80.5
;;;;;

latd = where(latitgrid le -77 and latitgrid gt -80.5 and $
    (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315))
latdpsr = where(latitgrid le -77 and latitgrid gt -80.5 and $
    (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and ill eq 1.)

pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel WEH: -77$\deg$ to -80.5$\deg$', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 77$\deg$: ',wehcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))


wpall = plot(xvl(0:35), xvl(0:35)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -77 and latitgrid gt -80.5 and $
            (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and $
             psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor

ppvals(7) = 35
incavg = plot(ppvals(0:7), avgweh(0:7), color= 'red', thick=2, /overplot)


latd = where(latitgrid le -74 and latitgrid gt -77 and $
  (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315))
latdpsr = where(latitgrid le -74 and latitgrid gt -77 and $
  (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and ill eq 1.)

pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel WEH: -74$\deg$ to -77$\deg$', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 74$\deg$: ',wehcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))


wpall = plot(xvl(0:25), xvl(0:25)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -74 and latitgrid gt -77 and $
            (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and $
             psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor


incavg = plot(ppvals(0:5), avgweh(0:5), color= 'red', thick=2, /overplot)

latd = where(latitgrid le -71 and latitgrid gt -74 and $
  (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315))
latdpsr = where(latitgrid le -71 and latitgrid gt -74 and $
  (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and ill eq 1.)

pl = plot(psr_den1(latd), csetn_coll_weh(latd), symbol='circle', font_size=15, Font_name='Times', $
  sym_size=0.3, sym_filled=1, linestyle=6, color='grey', yrange=[-0.25, 0.42], $
  xrange=[-5, 105], title='FOV(PSR %) vs Pixel WEH: -71$\deg$ to -74$\deg$', $
  xtitle='FOV(PSR %)', ytitle='Pixel WEH')

xvl = findgen(100)
wehcoef83 = linfit(psr_den1(latd),csetn_coll_weh(latd))
print, '> 74$\deg$: ',wehcoef83
print, 'n pixels = ', n_elements(latd)
print, 'n psr pixels = ',n_elements(latdpsr)
print, 'psr fraction of area', float(n_elements(latdpsr))/ float(n_elements(latd))


wpall = plot(xvl(0:18), xvl(0:18)*wehcoef83(1) + wehcoef83(0), color='black', thick=3,/overplot)

zro = fltarr(100)
pfz = plot(xvl, zro, thick=2, color='black',linestyle=2,/overplot)

avgweh = fltarr(20)
pvals = findgen(11)*10.
ppvals = findgen(20)
j = 0
for i = 0,9 do begin
  latdd = where(latitgrid le -71 and latitgrid gt -74 and $
    (longigrid lt 45 or (longigrid gt 135 and longigrid lt 225) or longigrid gt 315) and $
    psr_den1 ge pvals(i) and psr_den1 lt pvals(i+1))
  avgweh(j) = mean(csetn_coll_weh(latdd))
  avgweh(j+1) = avgweh(j)
  ppvals(j) = pvals(i)
  ppvals(j+1) = pvals(i+1)
  j = j + 2
endfor


incavg = plot(ppvals(0:3), avgweh(0:3), color= 'red', thick=2, /overplot)


read,xxx

return
end


