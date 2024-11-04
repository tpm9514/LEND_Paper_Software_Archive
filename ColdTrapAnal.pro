; Program:   ColdTrapAnal.pro
; Comment:   The program 1) derives the PSR and Non-PSR WEH averages as a function of 2 deg latitude bins.
;            2) Correlates the PSRs observed WEH to their diameters.
;
; Input:  Collimated WEH, Topography, Binary PSR map
; Output:  Generates Figs. 7 and 8.
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




; PRO ColdTrapAnal
; This code is for Section 3.3
; The objective is perform a granulometry of all PSR's WEH to determine if the results are
; consistent with the hypothesis of instrumental blurring.  Linear fits are used to determine if true.
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
Pro ColdTrapAnal

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
  
  dt = where(illum eq 0)
  ill = illum
  ill(*) = 0
  ill(dt) = 1
  ill1 = segmentregions(ill,1)
  psravglat = fltarr(max(ill1)+1)
  psravglon = fltarr(max(ill1)+1)
  psravgdiam = fltarr(max(ill1)+1)
  psrarea = fltarr(max(ill1)+1)
  
  for k = 1, max(ill1) do begin
    dt = where(ill1 eq k,adat)
    numpsr = uniq(ill1(dt))
    psravglat(k) = mean(latitgrid(dt))
    psravgdiam(k) = sqrt(adat*4./!pi)*2.
    psrarea(k) = adat*4.
    lonpix = longigrid(dt)
    psravglon(k) = mean(lonpix)
    if max(lonpix) gt 345 and min(lonpix) lt 45. then begin
      l1 = where(lonpix gt 345)
      l2 = where(lonpix lt 45)
      psravglon(k) = mean([360.-l1,l2])
      if psravglon(k) lt 0. then psravglon(k) = 360.-psravglon(k)
    endif
  endfor  
   
; Read in WEH maps
  csetn_coll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')
  csetn_coll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_Suppress_sig220_UL4.tiff')
  

; This segment of the code derives Figure 7, PSR, non-PSR averaged WEH vs Latitude
;
ill2 = ill1

    lats = findgen(11)*2. -90
    coll_all = fltarr(10,2)
    coll_psr1 = fltarr(10,2)
    coll_psr2 = fltarr(10,2)
    avgpsrdiam = fltarr(10,2)
    coll_n_psr = fltarr(10,2)
    kstest = fltarr(10,2)
    kstest_supp = fltarr(10,2)
    
    psr_fract = fltarr(10)
    diff = psr_fract
    csetn_coll_weh_back_rem = csetn_coll_weh
;    
; Make Latitude Profiles
; i = 0 to 3 is for the high latitudes and all longitudes are used
; i > 3 is for the lower latitudes and eliminates 90 +/- 45 degs and 270 +/- 45 degs where 
; LRO station keeping and LEND power cycling was performed in mission years < 2.
;
    for i = 0,9 do begin
      if i le 3 then begin
; all is all pixels in the latitude band, psr1 is all PSR pixels in band, npsr is non-PSR pixels in band
        all = where(latitgrid ge lats(i) and latitgrid lt lats(i+1), act)
        psr1 = where(latitgrid ge lats(i) and latitgrid lt lats(i+1) and illum eq 0,pct1)
        npsr = where(latitgrid ge lats(i) and latitgrid lt lats(i+1) and illum gt 0.,npct)
      endif
      if i gt 3 then begin
      all = where(latitgrid ge lats(i) and latitgrid lt lats(i+1) and $
        ((longi ge 0 and longi le 45.) or $
        (longi gt 135 and longi lt 225.) or $
        (longi gt 315)),act)
      psr1 = where(latitgrid ge lats(i) and latitgrid lt lats(i+1) and $
        ((longi ge 0 and longi le 45.) or $
        (longi gt 135 and longi lt 225.) or $
        (longi gt 315)) and illum eq 0,pct1)
      npsr = where(latitgrid ge lats(i) and latitgrid lt lats(i+1) and $
        ((longi ge 0 and longi le 45.) or $
        (longi gt 135 and longi lt 225.) or $
        (longi gt 315)) and illum gt 0,npct)
      endif
            
; KS test of the psr, npsr pixel distributions for each latitude band
;
      kstwo, csetn_coll_weh(psr1),csetn_coll_weh(npsr),D, probks
      
      kstest(i,0) = D
      kstest(i,1) = probks
      coll_psr1(i,0) = mean(csetn_coll_weh(psr1))
      coll_n_psr(i,0) = mean(csetn_coll_weh(npsr))
      coll_all(i,0) = mean(csetn_coll_weh(all))
      coll_psr1(i,1) = stddev(csetn_coll_weh(psr1)) / sqrt(pct1)
      coll_n_psr(i,1) = stddev(csetn_coll_weh(npsr)) / sqrt(npct)
      coll_all(i,1) = stddev(csetn_coll_weh(all)) / sqrt(act)

      psr_fract(i) = float(pct1)/float(act)
      diff(i) = coll_psr1(i,0) - coll_all(i,0)
;
; Subtract the non-PSR WEH avg from the Collimated WEH map for each latitude bin, below.
; Used in Fig. 8.
;
      csetn_coll_weh_back_rem(all) = csetn_coll_weh_back_rem(all) - coll_n_psr(i,0)

      psr_fract(i) = float(pct1)/float(act)
      diff(i) = coll_psr1(i,0) - coll_all(i,0)
    endfor

; Make Fig 7 plot.   
;
    lats1 = lats+1.
    lats1 = lats1(0:9)
    flat = fltarr(10)
    flat(*) = 0.
    ;
    pcollb = errorplot(lats1,coll_psr1(*,0),coll_psr1(*,1),title='7)  PSR and Non-PSR WEH vs Latitude, >70$\deg$ S',xtitle='Latitude',ytitle='WEH wt%', $
      font_size=16,font_name='Times',name='PSR',xrange=[-90.5,-69.],thick=2, color='blue',yrange=[-0.03, 0.16])
    pcoll1 = errorplot(lats1,coll_n_psr(*,0),coll_n_psr(*,1),name='non-PSR',thick=2, color='red',/overplot)
    pcoll3 = errorplot(lats1,coll_all(*,0),coll_all(*,1),linestyle=3, name='Latit. Mean',color='black',thick=2,/overplot)
    sig = where(kstest(*,1) lt 0.01 and coll_psr1(*,0) gt coll_n_psr(*,0))
    nsig = where(kstest(*,1) ge 0.01 or coll_psr1(*,0) le coll_n_Psr(*,0))
    psig = plot(lats1(sig),coll_psr1(sig,0),linestyle=6, symbol='circle',color='blue',sym_size=1.1,name='Signif., p < 0.01',/sym_filled,/overplot)
    pnsig = plot(lats1(nsig),coll_psr1(nsig,0),linestyle=6,symbol='circle',color='grey',sym_size=1.1,name='Not Sigif., p $\geq$ 0.01',/sym_filled,/overplot)

    npss = plot(lats1,coll_n_psr(*,0),linestyle=6, symbol='square',color='red',sym_size=0.7,name='non-PSR',/sym_filled,/overplot)
    
    print, 'PCorr (PSR To PSRfract):',correlate(coll_psr2(*,0),psr_fract)

; Save plot as image and csv text files
; 
    pcollb.save,Figures+'Fig_7_PSR_nPSR_WEH_vs_Latitude.png'
    head = ['Latit','PSR Avg WEHwt%', 'PSR SEM WEHwt%', 'All Avg WEHwt%', 'All SEM WEHwt%', 'NonPSR Avg WEHwt%', 'NonPSR SEM WEHwt%']
    table_header = ['Fig. 7: PSR and nonPSR vs Latitude Band Avgs']
    write_csv, CSV + 'Fig_7_PSR_nonPSR_WEH_vs_LatBandAvg.csv',lats1,coll_psr1(*,0),coll_psr1(*,1),coll_all(*,0),coll_all(*,1), $
        coll_n_psr(*,0),coll_n_psr(*,1), header=head, table_header=table_header

; Check the area distribution of PSRs by latitude, reported in the paper, consistent with Mazarico et al., 2011
;
    psrlat1 = where(psravglat le -85,p1)
    psrlat2 = where(psravglat le -80 and psravglat ge -85., p2)
    psrlat3 = where(psravglat le -75 and psravglat ge -80., p3)
    psrlat4 = where(psravglat le -70 and psravglat ge -75., p4)
    
    print, mean(psravgdiam(psrlat1)), stddev(psravgdiam(psrlat1))
    print, mean(psravgdiam(psrlat2)), stddev(psravgdiam(psrlat2))
    print, mean(psravgdiam(psrlat3)), stddev(psravgdiam(psrlat3))
    print, mean(psravgdiam(psrlat4)), stddev(psravgdiam(psrlat4))


coll_weh1 = csetn_coll_weh

; Here we reassign csetn_coll_weh to its non-PSR background subtracted version.
  csetn_coll_weh = csetn_coll_weh_back_rem


; Region grow the binary PSR map for all PSR pixels > 75 S.
; Makes independent PSRs
;
    dt = where(illum eq 0 and latitgrid lt -75)

    ill = illum
    ill(*) = 0
    ill(dt) = 1
    rem = where(latitgrid gt -83 and ((longigrid gt 45 and longigrid lt 135) or (longigrid gt 225 and longigrid lt 315.)))
    ill(rem) = 0    
    illa = ill
    illum1 = segmentregions(ill,1)
    n_psr = max(illum1)
    print, 'N PSRs to Evaluate: ',n_psr
    
    areapsr = fltarr(n_psr+1)
    diameter = fltarr(n_psr+1)
    wehpsr = fltarr(n_psr+1)
    meanwehpsr = fltarr(n_psr+1)
    latval = fltarr(n_psr+1)
    lonval = fltarr(n_psr+1)
    allx = fltarr(n_psr+1)
    ally = fltarr(n_psr+1)

    for i = 1,n_psr do begin
      dat = where(illum1 eq i, tpix)
      areapsr(i) = float(tpix)*4.
      diameter(i) = sqrt((float(tpix)*4.)/!pi)*2.
      wehpsr(i) = max(csetn_coll_weh(dat))
      latval(i) = mean(latitgrid(dat))
      lonpix = longigrid(dat)
      lonval(i) = mean(lonpix)

      if max(lonpix) gt 330. and min(lonpix) lt 30 then begin
        poslon = where(lonpix gt 30.)
        lonpix(poslon) = lonpix(poslon) + 360.
        lonval(i) = mean(lonpix)
        if lonval(i) gt 360. then lonval(i) = lonval(i) - 360.
      endif

      allx(i) = mean(dat mod 1440.)
      ally(i) = mean(dat / 1440.)
    endfor

areapsr = areapsr(1:n_elements(AReapsr)-1)
wehpsr = wehpsr(1:n_elements(WEHpsr)-1)
diameter = diameter(1:n_elements(diameter)-1)
latval = latval(1:n_elements(latval)-1)
lonval = lonval(1:n_elements(lonval)-1)
allx = allx(1:n_elements(allx)-1)
ally = ally(1:n_elements(ally)-1)

latlo = where(latval le -75 and latval gt -83)
lathi = where(latval le -83)

p = plot(diameter(lathi),wehpsr(lathi),title='8)  PSR Diameter vs Collim. WEH wt%:  502 PSRs',xtitle='PSR Diameter, km',ytitle='WEH wt%',font_name='Times',font_size=15, $
         linestyle=6, symbol='circle',color='red',sym_filled=1,sym_size=0.75, xrange=[-1, 40], yrange=[-0.25,0.40], $
         name='A: 83$\deg$ to 90$\deg$ S')
pm = plot(diameter(latlo),wehpsr(latlo),linestyle=6, symbol='square',sym_size=0.7, color='grey',sym_filled=1, name='B: 75$\deg$ to 83$\deg$ S',/overplot)

x = linfit(diameter,wehpsr)
lin = findgen(38) 
line = x(0) + x(1)*lin

pp = plot(findgen(38),line,linestyle=3,color='black',name='Linfit(502 PSRs)',thick=2,/overplot)

print, ' '
print, 'All PSRs Fit: Linear model: ',x(0),x(1)
print, 'Shoemaker All: WEH Pred: ',x(0) + x(1)*37.
print, 'Haworth All: WEH Pred: ', x(0) + x(1)*35.
print, 'Cabeus All: WEH Pred: ',x(0) + x(1)*19.
print, 'Faustini All: WEH Pred: ', x(0) + x(1)*29.

xhi = linfit(diameter(lathi),wehpsr(lathi))
lin = findgen(38)
line1 = xhi(0) + xhi(1)*lin

xlo = linfit(diameter(latlo),wehpsr(latlo))
lin = findgen(38)
line2 = xlo(0) + xlo(1)*lin

pplo = plot(findgen(38),line1,linestyle=3,color='red',name='Linfit(A: 83$\deg$ to 90$\deg$ S)',thick=2,/overplot)
pphi = plot(findgen(38),line2,linestyle=3,color='grey',name='Linfit(B: 75$\deg$ to 83$\deg$ S)',thick=2,/overplot)

print, ' '
print, 'Shoemaker and Haworth Predictions for A, B linear fits'
print, 'Shoemaker A, Lo: WEH Pred: ',xlo(0) + xlo(1)*37.
print, 'Shoemaker B, Hi: WEH Pred: ',xhi(0) + xhi(1)*37.
print, 'Haworth A: WEH Pred: ', xlo(0) + xlo(1)*35.
print, 'Haworth B: WEH Pred: ', xhi(0) + xhi(1)*35.

highWEH = where(wehpsr gt 0.)
lowWEH = where(wehpsr le 0.)

kstwo, diameter(highweh), diameter(lowweh), D, probks
print, ' ' 
print, 'PSR Diameter tests:'
print, 'Diams WEH > 0: ', mean(diameter(highweh)), stddev(diameter(highweh))
print, 'Diams WEH < 0: ', mean(diameter(lowweh)), stddev(diameter(lowweh))
print, 'Test of PSR Diameters: 1) WEH > 0 vs 2) WEH < 0,  ', D, probks

;kstwo, highweh(highwehdiams), high

p.save,Figures+'Fig_8_PSRDiam_vs_CollimWEH.png'
print, ' '
print, 'Correlation: Rho values for A, B linear fits'
print, 'Latlo, A: rho:',r_correlate(diameter(latlo),wehpsr(latlo))
print, 'Lathi, B: rho:',r_correlate(diameter(lathi),wehpsr(lathi))
print, 'All PSR, rho:', r_correlate(diameter, wehpsr)

shoe_weh_pred = fltarr(100)
psr_run = fltarr(n_elements(areapsr))

head = ['img_x','img_y','PSR Mean Lat','PSR Mean Longi', 'Diameter km', 'MaxWEH wt%']
table_header = ['Fig 8: PSR Diameter vs Max WEH wt%']
write_csv, CSV + 'Fig_8_PSR_Diam_vs_MaxWEH.csv', allx,ally,latval,lonval,diameter,wehpsr, $
  header=head, table_header=table_header

for b = 0, 99 do begin
  psr_run(*) = 0.
  subset = randomu(seed, n_elements(areapsr))
  anal = where(subset gt 0.5)
  psr_run(anal) = 1
  x = linfit(diameter(anal),wehpsr(anal))
  shoe_weh_pred(b) = x(0) + x(1)*37.
endfor

print, 'Mean(Shoe_weh_pred), StdDev(shoe_weh_pred): ', mean(shoe_weh_pred), stddev(shoe_weh_pred)

PSRWEH_gt_0_77 = where(illum eq 0 and csetn_coll_weh gt 0. and latitgrid le -77., psrw77)
PSRWEH_le_0_77 = where(illum eq 0 and csetn_coll_weh le 0. and latitgrid le -77., psrnw77)
PSRALL_77 = where(illum eq 0 and latitgrid le -77., psrall77)

PSRWEH_gt_0_82 = where(illum eq 0 and csetn_coll_weh gt 0. and latitgrid le -82.,psrw82)
PSRWEH_le_0_82 = where(illum eq 0 and csetn_coll_weh le 0. and latitgrid le -82.,psrnw82)
PSRALL_82 = where(illum eq 0 and latitgrid le -82.,psrall82)

print, 'Percent of PSR w positive WEH poleward of 77 S:', float(psrw77)/float(psrall77)
print, 'Percent of PSR w positive WEH poleward of 82 S:', float(psrw82)/float(psrall82)




read,'End',xxx

return
end


