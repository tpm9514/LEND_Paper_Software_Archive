; Program:   DilatePSR_n_WEH_Props123.pro
; Comment:   The program erodes and dilates the PSR edges to demonstrate Properties 1,2,3.  
;            Enhanced WEH in PSRs, Instrumental blurring, Lower WEH with distance into non-PSR
; 
; Input:  Collimated WEH, Topography, Binary PSR map
; Output:  Generates Figs. 9a,b
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




pro prob_ks, D, N_eff, probks
  ;+
  ; NAME:
  ;       PROB_KS
  ; PURPOSE:
  ;       Return the significance of the Kolmogoroff-Smirnov statistic.  T
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
   print, term
   print, probks
  endfor

  probks = 1.          ;Sum did not converge after 100 iterations
  return

end


pro kstwo, data1, data2, D, prob
  ;+
  ; NAME:
  ;       KSTWO
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
; Description
;   This code segments permanently shadowed regions from binary PSR map
;   Image has 0 values where not PSR, 1 = PSR. Method is a region growing
;   traverse.  Assigns a unique numeric identifier to all pixel clusters, e.g. in each
;   PSR, minarea is the minimum area PSR to report in pixels.  Methods are described in 
;   Sonka et al., 1992.  PSR pixels are derived from average illumination conditions with 0% illumination
;   defined in Maps defined from Mazarico et al., 2011. 
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
          xloc = xlist[strt]
          yloc = ylist[strt]
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



; DilatePSR_WEH_Props123
; Description:  Observe the expected hydrogen concentration as a function of distance from PSR
; Show properties 1, 2, 3 as a function of latitude
;
Pro DilatePSR_n_WEH_Props123

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
  restore,SaveSets+'OceanHaline.sav'
  r_curr = r
  g_curr = g
  b_curr = b

  tvlct, r_curr,g_curr,b_curr

  ; Read in LRO maps
  ;
  illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
  latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
  longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
  dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
  SlopeAzi = read_tiff(LRO_Maps+'Slope_Azi_Angle.tiff')
  Slope = read_tiff(LRO_Maps+'TopoSlope.tiff')
  illum_orig = illum

  ; Diviner Max Temp
  ;
  adim = 1440

; Resize images to standard 1440 x 1440 pixels, 2 km wide pixels
  longigrid = congrid(longigrid, adim,adim,/center)
  latitgrid = congrid(latitgrid, adim,adim,/center)
  slope = congrid(slope,adim,adim,/center)
  slopeazi = congrid(slopeazi,adim,adim,/center)

  latitgrid_orig = latitgrid

  illum_orig = illum
  illum = congrid(illum,adim,adim)

; Read in Collimated WEH map for analysis
  csetn_coll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')

; Dilation and erosion kernel, 3x3 pixels
  kern = fltarr(3,3)
  kern(*) = 1

  psr = where(illum eq 0)
  psrimg = illum
  psrimg(*) = 0.
  psrimg(psr) = 1
  psrimg1 = psrimg

; Define latitude bands and range
  ns = 31
  lats = [-73.,-77.,-80.,-83.,-86.,-88.,-90.]
  latslist = ['73 S to 77 S', '77 S to 80 S', '80 S to 83 S', '83 S to 86 S', '86 S to 88 S', '88 S to 90 S']

; Results arrays, init = -999 for invalid
  WEHMeanfromPSR = fltarr(6,31)
  wehmeanfrompsr(*) = -999
  WEHSEMFromPSR = wehmeanfrompsr

  rimtrack = latitgrid
  rimtrack(*) = 0.

  ist = 5
  st = 420
  ed = 1020
  display = 1

  img = latitgrid
  img(*) = 0
  st = 420
  ed = 1020
  imgpr = latitgrid
  imgpr(*) = 0.
  imgpsr = img
  psr = where(illum eq 0.)
  illumpsr = illum
  illumpsr(*) = 0
  illumpsr(psr) = 1
  imgpsr = illumpsr
  imgpsr(*) = 0


; This loop evaluates within PSR areas, negative distance observations, iteratively erodes PSR edges, then avg's weh from new rims.
  for a = 0,5 do begin
    psrimg1 = psrimg
    for i = ist,0,-1 do begin
      oldimg = psrimg1
      psrimg1 = erode(psrimg1,kern)
      rim = oldimg-psrimg1
      data = where(rim eq 1 and latitgrid le lats(a) and latitgrid gt lats(a+1),nnn)
;  Analysis proceeds until the avail pixels are exhausted, i.e. <= 5.
      if n_elements(data) gt 5 then begin
        WEHMeanfromPSR(a,i) = mean(csetn_coll_weh(data))
        WEHSEMfromPSR(a,i) = stddev(csetn_coll_weh(data)) /sqrt(n_elements(data))
        print, lats(a),WEHMeanfromPSR(a,i), ' PSR depth from edge: ',(i*2)-(ist*2.)
        tvscl, psrimg(st:ed,st:ed)
      endif
    endfor
  endfor

  distf = findgen(ns)*2. - ist*2.
  rimtrack = psrimg
  rimtrack(*) = 0
  ict = 0

; This loop evaluates outside PSR areas, i.e. non-PSR, positive distance observations, iteratively dilates PSR edges, then avg's weh from new rims.
  for a = 0,5 do begin
    psrimg1 = psrimg
    for i = ist+1, ns-1 do begin
      if i lt ns-1 then begin
        oldimg = psrimg1
        psrimg1 = dilate(psrimg1,kern)
        rim = psrimg1-oldimg
        tvscl, rim(400:1000,400:1000)
        if a ge 3 then data = where(rim eq 1 and latitgrid le lats(a) and latitgrid gt lats(a+1),nnn)
        if a le 2 then data = where(rim eq 1 and latitgrid le lats(a) and latitgrid gt lats(a+1) and $
          (longigrid gt 315. or (longigrid ge 0 and longigrid lt 45) or $
          (longigrid gt 135 and longigrid lt 225)),nnn)
        rimtrack(data) = i
;  Analysis proceeds until the avail pixels are exhausted, i.e. <= 5.
        if n_elements(data) gt 5 then begin
          WEHMeanfromPSR(a,i) = mean(csetn_coll_weh(data))
          WEHSEMfromPSR(a,i) = stddev(csetn_coll_weh(data)) /sqrt(n_elements(data))
        endif
;        tvscl, psrimg1(400:1000,400:1000)
        print, a, i, distf(i), WEHMeanfromPSR(a,i), max(rimtrack)
;        read,xxx
      endif
      if i eq ns-1 then begin
        data = where(psrimg1 eq 0 and psrimg eq 0 and latitgrid le lats(a) and latitgrid gt lats(a+1),nnn)
        rimtrack(data) = i
        if n_elements(data) gt 5 then begin
          WEHMeanfromPSR(a,i) = mean(csetn_coll_weh(data))
          WEHSEMfromPSR(a,i) = stddev(csetn_coll_weh(data)) /sqrt(n_elements(data))
        endif
        ;        tvscl, psrimg1(400:1000,400:1000)
        print, WEHMeanfromPSR(a,i), max(rimtrack)
      endif
    endfor
  endfor
  print, ns
  tvscl, psrimg1

  psrpx = findgen(ns)*2 - ist*2.
  val = where(wehmeanfrompsr(0,*) ne -999.)
  p1 = errorplot(distf(val),WEHMeanFromPSR(0,val), WEHSEMfromPSR(0,val),title='Avg WEH wt% vs PSR Distance: -73$\deg$ to -83$\deg$',xtitle='Distance to PSR Edge [km]',$
    ytitle='Avg WEH [wt%]', xrange=[distf(0)-4,distf(30)+1],thick=2.,color='black',name='-73$\deg$:-77$\deg$',font_name='Times', font_size=14, yrange=[-0.03, 0.32])

  val = where(WEHMeanfromPSR(1,*) ne -999.)
  p83 = errorplot(distf(val),WEHMeanFromPSR(1,val), WEHSEMfromPSR(1,val),thick=2.,color='blue',name='-77$\deg$:-80$\deg$',/overplot)

  val = where(WEHMeanfromPSR(2,*) ne -999.)
  p80 = errorplot(distf(val),WEHMeanFromPSR(2,val), WEHSEMfromPSR(2,val),thick=2.,color='red',name='-80$\deg$:-83$\deg$',/overplot)

  val = where(wehmeanfrompsr(3,*) ne -999.)
  p11 = errorplot(distf(val),WEHMeanFromPSR(3,val), WEHSEMfromPSR(3,val),title='Avg WEH wt% vs PSR Distance: -83$\deg$ to -90$\deg$',xtitle='Distance to PSR Edge [km]',$
    ytitle='Avg WEH [wt%]', xrange=[distf(0)-4,distf(30)+1],thick=2.,color='purple',name='-83$\deg$:-86$\deg$',font_name='Times', font_size=14, yrange=[-0.03, 0.32])

  val = where(WEHMeanfromPSR(4,*) ne -999.)
  p75 = errorplot(distf(val),WEHMeanFromPSR(4,val), WEHSEMfromPSR(4,val),thick=2.,color='green',name='-86$\deg$:-88$\deg$',/overplot)

  val = where(WEHMeanfromPSR(5,*) ne -999.)
  p78 = errorplot(distf(val),WEHMeanFromPSR(5,val), WEHSEMfromPSR(5,val),thick=2.,color='dark orange',name='-88$\deg$:-90$\deg$',/overplot)

  head = ['Distance from PSR (km)', 'MeanWEH wt%','SEM WEH wt%']
  table_header = ['Fig 9b: -73 to -77 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9b_73S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(0,*),WEHSEMfromPSR(0,*), $
    header=head, table_header=table_header
    
  table_header = ['Fig 9b: -77 to -80 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9b_77S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(1,*),WEHSEMfromPSR(1,*), $
    header=head, table_header=table_header

  table_header = ['Fig 9b: -80 to -83 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9b_80S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(2,*),WEHSEMfromPSR(2,*), $
    header=head, table_header=table_header

  table_header = ['Fig 9a: -83 to -86 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9a_83S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(3,*),WEHSEMfromPSR(3,*), $
    header=head, table_header=table_header

  table_header = ['Fig 9a: -86 to -88 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9a_86S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(4,*),WEHSEMfromPSR(4,*), $
    header=head, table_header=table_header

  table_header = ['Fig 9a: -88 to -90 S: Mean WEH wt% vs Distance to PSR']
  write_csv, CSV + 'Fig_9a_88S_MeanWEH_vs_PSR_Distance.csv', distf,WEHMeanFromPSR(5,*),WEHSEMfromPSR(5,*), $
    header=head, table_header=table_header


  rho = fltarr(6,2)
  rho1 = fltarr(6,2)
  for i = 0,5 do begin
    val = where(wehmeanfrompsr(i,*) ne -999. and distf le 50.,nv)
    val1 = where(wehmeanfrompsr(i,*) ne -999. and distf le 50.,nv)
    x = fix(float(nv)/2.)
    rho(i,*) = r_correlate(distf(val),wehmeanfrompsr(i,val))
    rho1(i,*) = r_correlate(distf(val1),wehmeanfrompsr(i,val1))
;
;    plt1 = plot(distf(val(0:x+1)),wehmeanfrompsr(i,val(0:x+1)), title='Values from each Latitude Band: ' + strcompress(i,/remove_all), color='black')
;    plt2 = plot(distf(val(x+1:nv-1)),wehmeanfrompsr(i,val(x+1:nv-1)),color='red',/overplot)
;    
    kstwo,wehmeanfrompsr(i,val(0:x-1)),wehmeanfrompsr(i,val(x-1:nv-1)), D, prob
    
    Print, 'Table 3 Lat Band Entries'
    print, latslist(i), strcompress(n_elements(val)),',',strcompress(D), ',', strcompress(prob), ',', n_elements(val(0:x+1)), n_elements(val(x-1:nv-1)), ',',strcompress(rho(i)), ','
  endfor

  p1.save,Figures+'Fig9b_WEH_vs_PSR_Distance_73S.png'
  p11.save,Figures+'Fig9a_WEH_vs_PSR_Distance_83S.png'
  
  return
end
