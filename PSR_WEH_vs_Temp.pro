; Program:   PSR_WEH_vs_Temp.pro
; Comment:   The program analzes Haworth, Shoemaker, Faustini PSRs and correlates their internal WEH
; distributions to the PSRs internal maximum temperature distributions.
; 
; Input:  Collimated WEH, Topography, Latitude, Maximum Temperature maps
; Output:  Generates Figs. 6ad.
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



; Function CSETN_FOV_KERN
;    Function generates 2-D CSETN FOV kernels for map smoothing based on neutron transport modeling
;    Weights in collimated part of smoothing kernel = 0, outside 3-sigma = 0.
Function CSETN_FOV_Kern, sigmac, sigmau, km_per_pixel,wid_k,Comp

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets
; Field of View for kernels
;
restore,SaveSets+'CSETN_FOV2.sav'
  
  csetn_collimated = csetn_collimated2
  csetn_uncollimated = csetn_uncollimated2
  
  sigmac_max = sigmac*3.
  sigmau_max = sigmau*3.
  Coll_FOV = fltarr(wid_k, wid_k)
  UnColl_FOV = Coll_FOV

  center = wid_k/2.
  for i = 0,wid_k-1 do begin
    for j = 0,wid_k-1 do begin
      distance = sqrt((float(i)-center)^2+ (float(j)-center)^2.)*km_per_pixel
      if distance lt sigmac_max then Coll_FOV(i,j) = interpol(csetn_collimated,range,distance)
      if distance lt sigmau_max then UnColl_FOV(i,j) = interpol(csetn_uncollimated,range,distance)
    endfor
  endfor
  
  Coll_FOV = Coll_FOV / total(Coll_FOV)  
  UnColl_FOV = UnColl_FOV / total(UnColl_FOV)
  
  if comp eq 0 then return, Coll_FOV
  if comp eq 1 then return, UnColl_FOV
  return, -1
end

Function Dot, scx,scy, snx,sny
  v1 = [scx, scy]
  v2 = [snx, sny]
  theta = acos(total((v1 / sqrt(total(v1^2)))*(v2 / sqrt(total(v2^2)))))
  return,theta/!dtor
end


function mag,v1,v2
  return, sqrt(v1^2.+v2^2.)
end

function GaussVal, stdev, ptx
  return, exp(-(ptx)^2./(2.*stdev^2.))/(stdev*sqrt(2.*!pi))
end

function MakeKern, dim

  kern = fltarr(dim,dim)
  d2 = dim/2
  sig = dim / 6.
  for i = 0,dim-1 do begin
    for j = 0,dim-1 do begin
      d = sqrt(float((i-d2)^2+(j-d2)^2))
      kern(i,j) = gaussval(sig,d)
    endfor
  endfor

  kern = kern/total(kern)
  kern = kern/total(kern)

  return, kern
end


; FUNCTION SEGMENTREGIONS
; 
; Description:
;   This code segments areas including permanently shadowed regions from binary maps.
;   Assigns unique numerical identifier to each region
;   
;   Example for PSRs:
;   Image has 0 values where not PSR, 1 = PSR.  Method is a region growing
;   traverse. Assigns a unique numeric identifier to all pixels in each
;   PSR.
;
; Source:   T. McClanahan, as coded from description in Sonka et al., 1999
;
Function SegmentRegions,image, minarea

  dimx = n_elements(image(*,0))-1
  dimy = n_elements(image(0,*))-1

; image is the source binary image

; Imstore image contains the classified regions and their numerical identifiers
  imstore = float(image)
  imstore(*,*) = 0.

; At each pixel (a,b) the surrounding 8 pixels are evaluated
  xchk = [-1, 0, 1, 0,-1, 1,-1, 1]
  ychk = [ 0, 1, 0,-1,-1, 1, 1,-1]

  id = long(1)

; Traverse the image
  for a = 1, dimx-2 do begin
    for b = 1, dimy-2 do begin
      
      ; check if the pixel is in an ROI and has'nt been evaluated yet
      if imstore(a,b) eq 0 and image(a,b) eq 1 then begin
        xlist = [long(a)]
        ylist = [long(b)]
        imstore(a,b) = id
        strt = long(0)
        xcnt = long(0)
        lend = long(1)
        
        ; Grow the region by making a list of all connected pixels in image
        while strt lt lend do begin
          xloc = xlist(strt)
          yloc = ylist(strt)
          
          ; traverse the list and add pixels from the list if valid 
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



; Pro PSR_WEH_vs_TEMP
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
Pro PSR_WEH_vs_Temp

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
  sig1 = 2.25  ;  This one in the paper
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
  

hawc = where(ill1 eq 533.)
hawc1 = where(ill1 eq 533. and longigrid le 360. and longigrid gt 340.)
hawc2 = where(ill1 eq 533. and longigrid gt 0 and longigrid lt 20.)
shoc = where(ill1 eq 595.)
fauc = where(ill1 eq 664.)


haw = where(ill1 eq 533. and slope lt 3)
sho = where(ill1 eq 595. and slope lt 3)
fau = where(ill1 eq 664. and slope lt 3)

hawu = where(ill1 eq 533. and slope ge 3)
shou = where(ill1 eq 595. and slope ge 3)
fauu = where(ill1 eq 664. and slope ge 3)

ple = where(latit le -89.95)

Illcheck = ill1
illcheck(*) = 255
illcheck(hawc) = 130
illcheck(shoc) = 130
illcheck(fauc) = 130
illcheck(haw) = 0
illcheck(sho) = 0
illcheck(fau) = 0

illcheck(ple) = 10

loadct, 0
window,xsize=400,ysize=400
tvscl, congrid(illcheck(705:775,705:775), 400, 400)
psrs = congrid(illcheck(705:775,705:775),400,400)
write_png,Figures + 'Fig_6c_HawShoFau_LowHigh_Slope.png', psrs

hawcx = hawc mod 1440
hawcy = hawc / 1440
shocx = shoc mod 1440
shocy = shoc / 1440
faucx = fauc mod 1440
faucy = fauc / 1440


hawx = haw mod 1440
hawy = haw / 1440
shox = sho mod 1440
shoy = sho / 1440
faux = fau mod 1440
fauy = fau / 1440

latit = latitgrid
phlat = plot(latit(hawc),csetn_coll_weh(hawc),linestyle=6,symbol='circle',title='6a)  PSR Latitude vs WEH', ytitle='WEH wt%',xtitle='Latitude$\deg$',$
  sym_size=0.8, color='black',name='Haworth',sym_filled=1,xrange=[-88.8, -86.5],yrange=[-0.05,0.45],font_name='Times',font_size=17)
pslat = plot(latit(Shoc),csetn_coll_weh(shoc),linestyle=6,symbol='triangle',color='blue',sym_size=1.1,sym_filled=1,name='Shoemaker',/overplot)
pflat = plot(latit(fauc),csetn_coll_weh(fauc),linestyle=6,symbol='square',color='red',sym_size=0.75,sym_filled=1,name='Faustini',/overplot)

phlat.save,Figures+'Fig_6a_HawShoeFau_Latit_vs_WEH.png'

head = ['img_x','img_y','Latit deg','Longi East deg ', 'WEH wt%', 'MaxTemp K']
table_header = ['Fig 6a: Haworth, Latit_vs_WEH']
write_csv, CSV + 'Fig_6a_Haworth_Latit_vs_WEH_MaxTemp.csv', hawcx,hawcy,latitgrid(hawcx,hawcy),longigrid(hawcx,hawcy), $
                 csetn_coll_weh(hawcx, hawcy),tempmax(hawcx,hawcy), $
  header=head, table_header=table_header
table_header = ['Fig 6a: Shoemaker, Latit vs WEH']
write_csv, CSV + 'Fig_6a_Shoemaker_Latit_vs_WEH_MaxTemp.csv', shocx,shocy,latitgrid(shocx,shocy),longigrid(shocx,shocy), $
                 csetn_coll_weh(shocx,shocy),tempmax(shocx, shocy), $
  header=head, table_header=table_header
table_header = ['Fig 6a: Faustini, Latit vs WEH']
write_csv, CSV + 'Fig_6a_Faustini_Latit_vs_WEH_MaxTemp.csv', faucx,faucy,latitgrid(faucx,faucy),longigrid(faucx,faucy), $
                 csetn_coll_weh(faucx,faucy),tempmax(faucx,faucy), $
  header=head, table_header=table_header

ph = plot(tempmax(haw),csetn_coll_weh(haw),linestyle=6,symbol='circle',title='6d)  PSR Basins:  Max. Temp. vs WEH', ytitle='WEH wt%',xtitle='Max. Temp. K',$
          sym_size=0.8, color='black',name='Haworth Basin',sym_filled=1,xrange=[45.,95],yrange=[-0.05,0.45],font_name='Times',font_size=17)
ps = plot(tempmax(Sho),csetn_coll_weh(sho),linestyle=6,symbol='triangle',color='blue',sym_size=1.1,sym_filled=1,name='Shoemaker Basin',/overplot)
pf = plot(tempmax(fau),csetn_coll_weh(fau),linestyle=6,symbol='square',color='red',sym_size=0.75,sym_filled=1,name='Faustini Basin',/overplot)

h1 = linfit(tempmax(haw),csetn_coll_weh(haw))
s1 = linfit(tempmax(sho),csetn_coll_weh(sho))
f1 = linfit(tempmax(fau),csetn_coll_weh(fau))

xh = findgen(27) + 50.
hf = h1(0) + h1(1)*xh
xs = findgen(37) + 52.
sf = s1(0) + s1(1)*xs
xf = findgen(37) + 50.
ff = f1(0) + f1(1)*xf

h1pred = h1(0) + (h1(1)*tempmax(haw))
s1pred = s1(0) + (s1(1)*tempmax(sho))
f1pred = f1(0) + (f1(1)*tempmax(fau))

hr2 = 1. - total((csetn_coll_weh(haw)-h1pred)^2.) / total((csetn_coll_weh(haw)-mean(csetn_coll_weh(haw)))^2.)
sr2 = 1. - total((csetn_coll_weh(sho)-s1pred)^2.) / total((csetn_coll_weh(sho)-mean(csetn_coll_weh(sho)))^2.)
fr2 = 1. - total((csetn_coll_weh(fau)-f1pred)^2.) / total((csetn_coll_weh(fau)-mean(csetn_coll_weh(fau)))^2.)

print, 'R2 results:'
print, 'Haworth:  ', hr2
print, 'Shoemaker:  ', sr2
print, 'Faustini:  ', fr2


phh = plot(xh,hf,linestyle=3,color='black',name='Haw Fit',thick=3,/overplot)
phs = plot(xs,sf,linestyle=3,color='blue',name='Shoe Fit',thick=3,/overplot)
phf = plot(xf,ff,linestyle=3,color='red',name='Fau Fit',thick=3,/overplot)

ph.save,Figures+'Fig_6d_HawShoeFau_BasinWEH_vs_MaxTemp.png'

h1t = linfit(latit(haw),tempmax(haw))
s1t = linfit(latit(sho),tempmax(sho))
f1t = linfit(latit(fau),tempmax(fau))


xht = findgen(23)*(0.04) - 88.0
ht = h1t(0) + h1t(1)*xht
xst = findgen(28)*(0.035) - 88.6
st = s1t(0) + s1t(1)*xst
xft = findgen(23)*(0.035) - 87.6
ft = f1t(0) + f1t(1)*xft

h1tpred = h1t(0) + (h1t(1)*latit(haw))
s1tpred = s1t(0) + (s1t(1)*latit(sho))
f1tpred = f1t(0) + (f1t(1)*latit(fau))

htr2 = 1. - total((tempmax(haw)-h1tpred)^2.) / total((tempmax(haw)-mean(tempmax(haw)))^2.)
str2 = 1. - total((tempmax(sho)-s1tpred)^2.) / total((tempmax(sho)-mean(tempmax(sho)))^2.)
ftr2 = 1. - total((tempmax(fau)-f1tpred)^2.) / total((tempmax(fau)-mean(tempmax(fau)))^2.)

print, 'R2 Basin Temp results:'
print, 'Haworth:  ', htr2
print, 'Shoemaker:  ', str2
print, 'Faustini:  ', ftr2


phlatc = plot(latit(haw),tempmax(haw),linestyle=6,symbol='circle',title='6b)  PSR Basins:  Latitude vs Max. Temp.', ytitle='Max. Temp. K',xtitle='Latitude$\deg$',$
  sym_size=0.8, color='black',name='Haworth Basin',sym_filled=1,xrange=[-88.8, -86.5],yrange=[40,100],font_name='Times',font_size=17)
pslatc = plot(latit(Sho),tempmax(sho),linestyle=6,symbol='triangle',color='blue',sym_size=1.1,sym_filled=1,name='Shoemaker Basin',/overplot)
pflatc = plot(latit(fau),tempmax(fau),linestyle=6,symbol='square',color='red',sym_size=0.75,sym_filled=1,name='Faustini Basin',/overplot)

phht = plot(xht,ht,linestyle=3,color='black',name='Haw Fit',thick=3,/overplot)
phst = plot(xst,st,linestyle=3,color='blue',name='Shoe Fit',thick=3,/overplot)
phft = plot(xft,ft,linestyle=3,color='red',name='Fau Fit',thick=3,/overplot)


phlatc.save,Figures+'Fig_6b_BasinsOnly_WEH_vs_MaxTemp.png'
head = ['img_x','img_y','Latit deg','Longi deg E.', 'WEH wt%', 'MaxTemp. K']
table_header = ['Fig 6b: Haworth, Latit vs Max Temp']
write_csv, CSV + 'Fig_6b_Haworth_BasinOnly_Latit_vs_Max_Temp.csv', hawx,hawy,latitgrid(hawx,hawy),longigrid(hawx,hawy), $
                 csetn_coll_weh(hawx,hawy),tempmax(hawx,hawy), header=head, table_header=table_header
table_header = ['Fig 6b: Shoemaker, Latit vs Max Temp']
write_csv, CSV + 'Fig_6b_Shoemaker_BasinOnly_Latit_vs_Max_Temp.csv', shox,shoy,latitgrid(shox,shoy),longigrid(shox,shoy), $
                 csetn_coll_weh(shox,shoy),tempmax(shox,shoy), header=head, table_header=table_header
table_header = ['Fig 6b: Faustini, Latit vs Max Temp']
write_csv, CSV + 'Fig_6b_Faustini_BasinOnly_Latit_vs_Max_Temp.csv', faux,fauy,latitgrid(faux,fauy),longigrid(faux,fauy), $
                 csetn_coll_weh(faux, fauy),tempmax(faux,fauy), header=head, table_header=table_header
;

print, 'Haworth Correl Slope: ', h1
print, 'Shoemaker Correl Slope: ', s1
print, 'Faustini Correl Slope: ', f1
print, 'Haworth low slope area, km2: ', n_elements(Haw)*4.
print, 'Shoemaker low slope area, km2: ', n_elements(sho)*4.
print, 'Faustini low slope area, km2: ', n_elements(fau)*4.
print, 'Haworth Basin Thermal gradient: K / km', h1t(1)/30.3     ; h1t(1) is in degrees divide by 30.3 / deg
print, 'Shoemaker Basin Thermal gradient: K / km', s1t(1)/30.3   ; same
print, 'Faustini Basin Thermal gradient: K / km', f1t(1)/30.3    ; same
print, 'Haworth Max temp, mean +/- sigma:', mean(tempmax(haw)), stddev(tempmax(haw))
print, 'Shoemaker Max temp, mean +/- sigma:', mean(tempmax(sho)), stddev(tempmax(sho))
print, 'Faustini Max temp, mean +/- sigma:', mean(tempmax(fau)), stddev(tempmax(fau))

hawlon1 = 360.-longi(hawc1)
hawlon2 = longi(hawc2)
print, 'Haworth PSR mean Lat, Lon:', mean(latit(hawc)), mean([hawlon1, hawlon2])
print, 'Shoemaker PSR mean Lat, Lon:', mean(latit(shoc)), mean(longi(shoc))
print, 'Faustini PSRM mean Lat, Lon:', mean(latit(fauc)), mean(longi(fauc))

hawW = max(csetn_coll_weh(hawc), hwW)
HawT = tempmax(hawW)
SHoW = max(csetn_coll_weh(shoc), ShW)
ShoT = tempmax(ShoW)
FauW = max(csetn_coll_weh(fauc), FsW)
FauT = tempmax(FauW)


print, 'Haworth Max(PSR WEH), Lat, Lon, Temp:', HawW, Latit(hawc(HwW)), Longi(hawc(HwW)), tempmax(hawc(HwW))
print, 'Shoemaker Max(PSR WEH), Lat, Lon, Temp:', ShoW, Latit(fauc(ShW)), Longi(shoc(ShW)), tempmax(shoc(ShW))
print, 'Faustini Max(PSR WEH), Lat, Lon, Temp:', FauW, Latit(fauc(FsW)), Longi(fauc(FsW)), tempmax(fauc(FsW))


hawmean = mean(tempmax(haw))
shomean = mean(tempmax(sho))
faumean = mean(tempmax(fau))

hawhi = where(ill1 eq 533. and slope lt 3 and tempmax gt hawmean)
shohi = where(ill1 eq 595. and slope lt 3 and tempmax gt shomean)
fauhi = where(ill1 eq 664. and slope lt 3 and tempmax gt faumean)

hawlo = where(ill1 eq 533. and slope lt 3 and tempmax le hawmean)
sholo = where(ill1 eq 595. and slope lt 3 and tempmax le shomean)
faulo = where(ill1 eq 664. and slope lt 3 and tempmax le faumean)

print, 'Haworth: '
kstwo, hawhi, hawlo, d, probks
hd = d
hprobks = probks

print, d, probks
print, 'Shoemaker: '
kstwo, Shohi, Sholo, d, probks
sd = d
sprobks = probks

print, d, probks
print, 'Faustini: '
kstwo, fauhi, faulo, d, probks
fd = d
fprobks = probks

print, d, probks

hawcenterlat = median(latitgrid(hawc))
shocenterlat = median(latitgrid(shoc))
faucenterlat = median(latitgrid(fauc))
longi = longigrid
dt = where(longigrid ge 180)

hawcenterlon = mean(longi(hawc))
shocenterlon = mean(longigrid(shoc))
faucenterlon = mean(longigrid(fauc))

mxh = max(csetn_coll_weh(haw),hid)
mxs = max(csetn_coll_weh(sho),sid)
mxf = max(csetn_coll_weh(fau),fid)

print,'Correlate Haw: ', correlate(csetn_coll_weh(haw),tempmax(haw)), correlate(csetn_coll_weh(haw), tmax_smooth(haw))
print,'Correlate Shoe: ', correlate(csetn_coll_weh(sho),tempmax(sho)), correlate(csetn_coll_weh(sho), tmax_smooth(sho))
print,'Correlate Fau: ', correlate(csetn_coll_weh(fau),tempmax(fau)), correlate(csetn_coll_weh(fau),tmax_smooth(fau))
print, ' '

print,'Haworth Lat,Lon:', latitgrid(haw(hid)),longigrid(haw(hid))
print,'Haworth MaxWEH:',max(csetn_coll_weh(haw),hid)
print,'Haw center of mass:',hawcenterlat, 360+hawcenterlon
print,'Haw MaxTemp:', tempmax(haw(hid))
print, 'Haw Dist from center:',(latitgrid(haw(hid))-hawcenterlat)*30.3
xh = sort(tempmax(haw))
haws = tempmax(haw(xh))
fract = where(haws lt tempmax(haw(hid)))
print, 'Haw Percentile:',float(n_elements(fract)) / float(n_elements(haw))
print, ' '
print,'Shoe Lat,Lon:', latitgrid(sho(sid)),longigrid(sho(sid))
print,'SHoe MaxWEH:',max(csetn_coll_weh(sho),sid)
print,'Shoe center of mass:',shocenterlat,shocenterlon
print,'Shoe MaxTemp:',tempmax(sho(sid))
print, 'Shoe Dist from center:',(latitgrid(sho(sid))-shocenterlat)*30.3
xs = sort(tempmax(sho))
shos = tempmax(sho(xs))
fract = where(shos lt tempmax(sho(hid)))
print, 'Sho Percentile:',float(n_elements(fract)) / float(n_elements(sho))
print, ' '
print,'Fau Lat,Lon:', latitgrid(fau(fid)),longigrid(fau(fid))
print,'Fau MaxWEH:',max(csetn_coll_weh(fau),fid)
print,'Fau center of mass:',faucenterlat,faucenterlon
print,'Fau MaxTemp:',tempmax(fau(fid))
print, 'Fau Dist from center:',(latitgrid(fau(fid))-faucenterlat)*30.3
xf = sort(tempmax(fau))
faus = tempmax(fau(xf))
fract = where(faus lt tempmax(fau(fid)))
print, 'Fau Percentile:',float(n_elements(fract)) / float(n_elements(fau))


print, 'Table 2 Entries'
print, 'Haworth', n_elements(haw), n_elements(haw)*4., mean(tempmax(haw)), stddev(tempmax(haw)), h1t(1)/30.3, hd, hprobks
print, 'Shoemaker', n_elements(sho), n_elements(sho)*4., mean(tempmax(sho)), stddev(tempmax(sho)), s1t(1)/30.3, sd, sprobks
print, 'Faustini', n_elements(fau), n_elements(fau)*4., mean(tempmax(fau)), stddev(tempmax(fau)), f1t(1)/30.3, fd, fprobks


read,'END WEH vs TEMP', xxxx

return
end


