; Program:   NeutronEnhanced.pro
; Comment:   The program quantifies the neutron enhanced regions above 82 S.
; Input:  WEH, Statistics, Topo, Collimated, Uncollimated, Illum,  PSR maps
; Output:  Generates Figs. 4a,b
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
  ; 
  ; Source:   NASA Goddard Space Flight Center Astronomy Users Library
  ;    https://asd.gsfc.nasa.gov/archive/idlastro/ftp/pro/math/prob_ks.pro
  ;
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
  
;window,0,xsize=600,ysize=500
;  plot, sortdata1
;  oplot, sortdata2

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

;window,1,xsize=600,ysize=500
;plot, id1
;oplot, id2

  ; The K-S statistic D is the maximum difference between the two distribution
  ; functions

;window,2,xsize=600,ysize=500
;  plot, fn1[id1]
;  oplot, fn2[id2]


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




; PRO NeutronEnhanced
; This code is for Section 3.1
; The objective is show the neutron enhanced regions 
; consistent with the hypothesis of instrumental blurring. 
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
Pro NeutronEnhanced

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
  

  c = 56.
  sig1 = 2.20  ;  This one in the paper
  x1 = gaussian_function([sig1,sig1],width=111.)
  x1 = x1 / total(x1)
  coll_fov = x1

  tmax_smooth = convol(tempmax, coll_fov, /center, invalid=-999)

  ; Indices to Segment maps > -80 S with Pixel range
  st = 566
  ed = 874

  ; keep = pixels in 80 to S-pole band, kill = pixels > -80 S, used to make a circular S polar map below
  keep = where(latitgrid le -74.)
  kill = where(latitgrid gt -74.)

  years = indgen(10)+1

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


; Study of Non-PSR:  EFS vs PFS slopes for two longitude bands:
;efs1 = where(latitgrid le -82 and latitgrid ge -84 and slopeaziangle gt 150. and illum gt 0 and longigrid lt 180., e1)
;pfs1 = where(latitgrid le -82 and latitgrid ge -84 and slopeaziangle lt 30. and illum gt 0 and longigrid lt 180., p1)
;efs2 = where(latitgrid le -82 and latitgrid ge -84 and slopeaziangle gt 150. and illum gt 0 and longigrid ge 180., e2)
;pfs2 = where(latitgrid le -82 and latitgrid ge -84 and slopeaziangle lt 30. and illum gt 0 and longigrid ge 180., p2)

; Test 1:
;print, 'Test 1:  PFS mean, stderr: ',mean(csetn_coll_weh(pfs1)), stddev(csetn_coll_weh(pfs1))/sqrt(float(p1))
;print, 'Test 1:  EFS mean, stderr: ',mean(csetn_coll_weh(efs1)), stddev(csetn_coll_weh(efs1))/sqrt(float(e1))
;kstwo, csetn_coll_weh(pfs1), csetn_coll_weh(efs1), D, prob
;print, 'Test 1:  KS D, prob:', D, prob

; Test 2:
;print, 'Test 2:  PFS mean, stderr: ',mean(csetn_coll_weh(pfs2)), stddev(csetn_coll_weh(pfs2))/sqrt(float(p2))
;print, 'Test 2:  EFS mean, stderr: ',mean(csetn_coll_weh(efs2)), stddev(csetn_coll_weh(efs2))/sqrt(float(e2))
;kstwo, csetn_coll_weh(pfs2), csetn_coll_weh(efs2), D, prob
;print, 'Test 2:  KS D, prob:', D, prob


  Supp_start = csetn_coll_supp
  
  ; Validation of Table #1
  ps = where(illum eq 0)
  ill = illum
  ill(*)  =0
  ill(ps) = 1

  cab = where(latitgrid le -84.35 and latitgrid ge -84.45 and longigrid gt 311.7 and longigrid lt 311.9)
  Haw = where(latitgrid le -87.25 and latitgrid ge -87.35 and longigrid gt 4.2 and longigrid lt 4.4)
  Sho = where(latitgrid le -87.9 and latitgrid ge -88.1 and longigrid gt 43.6 and longigrid lt 43.8)
  Fau = where(latitgrid le -86.8 and latitgrid ge -87.0 and longigrid gt 81.1 and longigrid lt 81.3)

  ill(cab) = 2
  ill(Haw) = 2
  ill(sho) = 2
  ill(fau) = 2

  print, ' '
  print, 'PSR Max WEH Reports'
  print, 'Cabeus WEH: ',csetn_coll_weh(cab)
  print, 'Haw WEH: ',csetn_coll_weh(haw)
  print, 'SHoemaker WEH: ',csetn_coll_weh(sho)
  print, 'Faustini WEH: ',csetn_coll_weh(fau)

  print, 'Cabeus MxTemp: ',tempmax(cab)
  print, 'Haw MxTemp: ',tempmax(haw)
  print, 'SHoemaker MxTemp: ',tempmax(sho)
  print, 'Faustini MxTemp: ',tempmax(fau)

  tvscl, ill(500:900,500:900)

  ; Analysis of Neutron ENhanced area above 82 S

RegionPFS_nonPSR1 = where(latitgrid le -80 and latitgrid ge -82 and illum gt 0 and slopeaziangle le 30., p1)
RegionEFS_nonPSR1 = where(latitgrid le -80 and latitgrid ge -82 and illum gt 0 and slopeaziangle gt 150., e1)
RegionPFS_nonPSR2 = where(latitgrid le -80 and latitgrid ge -82 and longigrid le 180 and illum gt 0 and slopeaziangle le 30.,p2)
RegionEFS_nonPSR2 = where(latitgrid le -80 and latitgrid ge -82 and longigrid gt 180 and illum gt 0 and slopeaziangle gt 150.,e2)

print, ' '
print, 'Comparison of PFS vs EFS for low latitude band w few PSR in Discussion section'
kstwo, csetn_coll_weh(RegionPFS_nonPSR1), csetn_coll_WEH(RegionEFS_nonPSR1), D, Prob
print, 'EFS vs PFS, 82 to 84',D, Prob
kstwo, csetn_coll_weh(RegionPFS_nonPSR2), csetn_coll_WEH(RegionEFS_nonPSR2), D, Prob
print, 'EFS vs PFS, 82 to 84',D, Prob
print, ' '
print, 'PFS Coll WEH 1:', mean(csetn_coll_weh(RegionPFS_nonPSR1)), stddev(csetn_coll_weh(RegionPFS_nonPSR1))/sqrt(float(p1))
print, 'EFS Coll WEH 1:',mean(csetn_coll_weh(RegionEFS_nonPSR1)), stddev(csetn_coll_weh(RegionEFS_nonPSR1))/sqrt(float(e1))
print, 'PFS Coll WEH 2:', mean(csetn_coll_weh(RegionPFS_nonPSR2)), stddev(csetn_coll_weh(RegionPFS_nonPSR2))/sqrt(float(p2))
print, 'EFS Coll WEH 2:',mean(csetn_coll_weh(RegionEFS_nonPSR2)), stddev(csetn_coll_weh(RegionEFS_nonPSR2))/sqrt(float(e2))


NEpix = where(latitgrid le -82 and csetn_coll_supp ge 1.0, nepxct1)
PSRPix = where(latitgrid le -82 and illum eq 0., psrpxct1)
psrfract = where(latitgrid le -82 and csetn_coll_supp ge 1.0 and illum eq 0., pnect1)
Regionpix = where(latitgrid le -82, r82ct1)
SuppPix = where(latitgrid le -82. and csetn_coll_supp lt 1.0, Suppct1)

print, ' '
print, 'KS test of Slope AziAngle Distribs (All region) vs (NeutEnhanced)'
kstwo, SlopeAziangle(RegionPix),slopeaziangle(nepix), D, prob
print, 'KS test of Expected(Area 82) vs nenhance82 Pixel Distribs:', D, prob
kstwo, SlopeAziangle(PSRpix),slopeaziangle(nepix), D, prob
print, 'KS Test of Slope Azi Angle Distribs: All psr > 82S vs nenhance82:', D, prob

print, 'Stats for Fig 4'
print, 'Pixels, Area: Neutron Enhanced:', float(nepxct1), float(nepxct1*4.)
print, 'Pixels, Area: PSR: ', psrpxct1, float(psrpxct1)*4.
print, 'Pixels, ARea: All area >82 S:', r82ct1, float(r82ct1)*4.
print, 'Pixels that are NE and PSR:', pnect1

print, ' '
print, 'NE stats'
print, 'NE Fraction of Total > 82 S:', float(nepxct1)/float(r82ct1)
print, 'PSR Fraction of NE: ', float(pnect1)/float(nepxct1)
print, 'PSR fraction of Suppressed area:', float(psrpxct1) /float(suppct1) 

print, ' '
print, 'Temp Stats:  Expected'
print, 'Tempmax: NePix', mean(tempmax(nepix)), stddev(tempmax(nepix))
print, 'Tempmax: AllPix', mean(tempmax(RegionPix)), stddev(tempmax(RegionPix))
print, 'Tempmax: NePix', mean(tempmax(PSRPix)), stddev(tempmax(PSRpix))


NEPixels = fltarr(6)
PSRpixels = fltarr(6)
RegionPixels = fltarr(6)
slpdat = [0.,30.,60.,90.,120.,150.,180.]
slpinc = [15.,45.,75.,105.,135.,165.]

for i = 0,5 do begin
  NEpix82 = where(latitgrid le -82 and csetn_coll_supp gt 1.0 and slopeaziangle ge slpdat(i) and slopeaziangle lt slpdat(i+1), nepxct)
  PSRPix82 = where(latitgrid le -82 and illum eq 0. and slopeaziangle ge slpdat(i) and slopeaziangle lt slpdat(i+1), psrpxct)
  Region82 = where(latitgrid le -82 and slopeaziangle ge slpdat(i) and slopeaziangle lt slpdat(i+1), r82ct)
  NEPixels(i) = nepxct
  PSRPixels(i) = psrpxct
  RegionPixels(i) = r82ct
endfor
 

NEprcent = (NEPixels / total(Nepixels))*100.
PSRprcent = (PSRPixels / total(psrpixels))*100.
RegionPrcent = (RegionPixels/ total(Regionpixels))*100.

print, 'Slope Azi angle:  Neutron ENhanced %:', neprcent
print, 'Slope Azi angle:  PSR %:', PSRprcent
print, 'Slope Azi angle:  RegionPrcent %:', RegionPrCent

psrplot = plot(slpinc, neprcent, title='Slope Azimuth Angle Distributions, > 82$\deg$ S',xtitle='Slope Azi. Angle$\deg$', $
     ytitle='% of Distrib.', thick=2, font_name='Times',xrange=[0., 180.], color = 'red', font_size=14)
psrplot = plot(slpinc, psrprcent, thick=2, color='blue',/overplot)
regplot = plot(slpinc,regionprcent, thick=2,color='black',/overplot)

psrplot2 = plot(slpinc, neprcent, symbol='triangle',sym_size=1.5,sym_filled=1,color='red',name='Neut. Enhance.',/overplot)
psrplot2 = plot(slpinc, psrprcent, symbol='circle',sym_size=1.2,sym_filled=1,color='blue',name='PSR',/overplot)
regplot2 = plot(slpinc,regionprcent, thick=2, symbol='square',sym_size=1.1,sym_filled=1,name='Area > 82$\deg$ S',color='black',/overplot)

psrplot.save,Figures+'Fig_4b_PSR_nPSR_WEH_vs_Latitude.png'
head = ['Latit','PSR Avg WEHwt%', 'PSR SEM WEHwt%', 'All Avg WEHwt%', 'All SEM WEHwt%', 'NonPSR Avg WEHwt%', 'NonPSR SEM WEHwt%']
table_header = ['Fig 4b: PSR and nonPSR vs Latitude Band Avgs']


;write_csv, CSV+'Fig_4_NeutEnh_vs_PSR_SlopeAziAngle.csv',Header=string(slpinc) , d1,d2,d3,d4,d5,d6


print, 'Neutron Enhanced, MaxTemp: Mean, Std:', mean(TempMax(nepix)), stddev(TempMax(nepix))
print, 'PSR, MaxTemp: Mean, Std:', mean(TempMax(psrpix)), stddev(TempMax(psrpix))

  neareas = illum
  neareas(*) = 0
  neareas(nePix) = 1
  allneareas = segmentregions(neareas, 1)
  
  print, 'Number of Neutron Enhanced areas > 82 S:',max(allneareas)
  print, ' '

saa = fltarr(max(allneareas)+1)
ar = saa

for i = 1, max(allneareas) do begin
  dt = where(allneareas eq i, dct)
  saa(i) =  mean(slopeaziangle(dt))
  ar(i) = float(dct)
endfor

big = where(ar gt 25)
saa  = saa(big)
ar = ar(big)

  s = 590
  e = 850
  cw = csetn_coll_weh(s:e,s:e)
  cw_start  = cw
  lat = latitgrid(s:e,s:e)
  lat_start = lat
  lon = longigrid(s:e,s:e)
  lon_start = lon
  ill = illum(s:e,s:e)
  ill_start = ill
  saz = slopeaziangle(s:e,s:e)
  sup = csetn_coll_supp(s:e,s:e)

  ill = congrid(ill,500,500)
  lon = congrid(lon,500,500)
  lat = congrid(lat,500,500)
  cw  = congrid(cw,500,500)
  sup = congrid(sup, 500, 500)

  necw = fltarr(n_elements(cw(*,0)),n_elements(cw(0,*)))

  neuten = where(lat le -82. and sup gt 1.0)
  psr = where(lat le -82. and ill eq 0)
  Intersect = where(lat le -82 and ill eq 0 and sup gt 1.0)

  necw(*) = 0.
  necw(neuten) = 220
  necw(intersect) = 255
  window,1,xsize=500,ysize=500

  tvscl, necw
  necwpsr = fltarr(3,500,500)
  necwpsr(*) = 170

  illpsr = ill
  illpsr(*) = 0
  illpsr(psr) = 1
  k = fltarr(3,3)
  k(*) = 1
  illpsr1 = dilate(illpsr,k)
  illpsr2 = illpsr1-illpsr

  for i = 0,n_elements(necw(*,0))-1 do begin
    for j = 0,n_elements(necw(0,*))-1 do begin
      if necw(i,j) eq 220 then begin
        necwpsr(0,i,j) = 255
        necwpsr(1,i,j) = 70
        necwpsr(2,i,j) = 70
      endif
      if necw(i,j) eq 255 then begin
        necwpsr(0,i,j) = 255
        necwpsr(1,i,j) = 255
        necwpsr(2,i,j) = 100
      endif
      if lat(i,j) gt -82. then begin
        necwpsr(0,i,j) = 255
        necwpsr(1,i,j) = 255
        necwpsr(2,i,j) = 255
      endif
      if illpsr2(i,j) eq 1 then begin
        necwpsr(*,i,j) = 0
      endif
    endfor
  endfor


  write_tiff, DerivedMaps+'Fig_4a_LatDist82.tiff',lat,/float
  write_tiff, DerivedMaps+'Fig_4a_LonDist82.tiff',lon,/float
  write_tiff, DerivedMaps+'Fig_4a_NeutEnh_n_PSR.tiff',necw,/float
  
  tv, necwpsr,true=1
  img = tvrd(true=1)
  write_png,Figures+'Fig_4a_NegWEH.png', img
  all_PSR_spp = where(lat_start le -82. and cw_start gt 0. and ill_start eq 0.)
  psr_n_ne = where(ill_start eq 0 and cw_start le 0. and lat_start le -82.)
  psr_frags = cw_start
  psr_frags(*) = 0
  all_psr_supp = psr_frags
  psr_frags(psr_n_ne) = 1
  all_psr_supp(all_psr_spp) = 1
  psrfrag = segmentregions(psr_frags, 1)
  psrareas = segmentregions(all_psr_supp,1)
  
  psr_frag_area = fltarr(max(psrfrag)+1)
  for i = 1,n_elements(psr_frag_area)-1 do begin
    a = where(psrfrag eq i,act)
    psr_frag_area(i) = act*4.
  endfor
  
  psr_frag_area = psr_frag_area(1:n_elements(psr_frag_area)-1)
  psr_frag_diam = sqrt(psr_frag_area/ !pi)*2.
  
  print, 'Area and stddev in km: PSR areas that are NE:', mean(psr_frag_area), stddev(psr_frag_area)
  print, 'Diameter and stddev in km: PSR diams that are NE:', mean(psr_frag_diam), stddev(psr_frag_diam)

  psrs_a = fltarr(max(psrareas)+1)
  for i = 1,n_elements(psrs_a)-1 do begin
    p = where(psrareas eq i,pct)
    psrs_a(i) = pct*4.
  endfor
  psrs_a = psrs_a(1:n_elements(psrs_a)-1)
  psrs_diam = sqrt(psrs_a/!pi)*2.
  
  print, 'All Supp PSR and stddev in km: above 82:', mean(psrs_a), stddev(psrs_a)
  print, 'All Supp PSR Diameter and stddev in km: above 82 S:', mean(psrs_diam),stddev(psrs_diam)
   
  kstwo, psr_frag_diam,psrs_diam, D, prob
  print, 'KS Test:  Supp PSR diams vs NE PSR frags diam', D, prob


  Supp = supp_start
  Supp(*) = 0
  n = where(latitgrid le -82 and Supp_start ge 1.0)
  Supp(n) = 1
  Supp_seg = segmentregions(Supp, 1)

  ; Largest NE region near pole, NE1
  Ne1 = where(Supp_seg eq 79,ne1ct)   ; 79 is the region ID

  all88 = where(latitgrid le -88, all88ct)
  efs = where(latitgrid -88 and Supp_seg eq 79 and slopeaziangle gt 150, efsct79)
  epall = where(latitgrid le -88 and Supp_seg eq 79, epct79)
  pfs = where(latitgrid le -88 and Supp_seg eq 79 and slopeaziangle lt 30, pfsct79)

  print, 'Total pixels above 88 S:', all88ct
  print, ' '    
  print, 'Large NE region #79 km2', total(float(epct79))*4.
  print, 'Avg Illum for #79:', mean(illum(epall)), stddev(illum(epall))
  print, 'Avg Illum for 88 S', mean(illum(all88))
  print, 'TempMax for #79', mean(tempmax(epall)), stddev(tempmax(epall))
  print, 'EFS for region #79, km2', (float(efsct79)/float(epct79))*100.
  print, 'PFS for region #79, km2', (float(pfsct79)/float(epct79))*100.

read,'End Neutron Enhanced',xxxx


return
end
