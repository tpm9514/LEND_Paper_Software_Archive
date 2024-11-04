; Program:   csetn_map82_Fig.pro
; Comment:   The program generates the Fig. 3ad.
; Input:  WEH, Statistics, Topo, Collimated, Uncollimated, Illum,  PSR maps
; Output:  Generates Fig. 3ad
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


; FUNCTION GAUSSVAL
;   Derive Gaussian weight given pixel position and sigma.
function GaussVal, stdev, ptx
  return, exp(-(ptx)^2./(2.*stdev^2.))/(stdev*sqrt(2.*!pi))
end

; FUNCTION MAKEKERN
;    Procedure is used to make 2-D gaussian smoothing kernel 
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
;   This code segments permanently shadowed regions from binary PSR map
;   Image has 0 values where not PSR, 1 = PSR. Method is a region growing
;   traverse.  Assigns a unique numeric identifier to all pixel clusters, e.g. in each
;   PSR, minarea is the minimum area PSR to report in pixels
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



; PRO READ_IMAGES
;    Read in maps from TIFF images
;    Tiff Maps are CSETN and SETN integral cps coverage. CSETN1-4 detectors observs are normalized to 1.275 cps
Pro Read_Images, Yr
  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
   
  case YR of
    1: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_1year_prisci_circle_2km_Sept2009_to_Sept2010.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_1year_prisci_circle_2km_Sept2009_to_Sept2010.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_1year_prisci_circle_2km_Sept2009_to_Sept2010.tiff')
      Yrs = '1'
      NameYr = 'Fig_S15ad_CSETN_Collim_80S_1Yr_2km.png'
      Title = 'CSETN: 1 Yr: Sept2009 to Sept2010'
    end
    2: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_2year_2km_July2009_to_July2011.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_2year_2km_July2009_to_July2011.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_2year_2km_July2009_to_July2011.tiff')
      Yrs = '2'
      NameYr = 'Fig_S16ad_CSETN_Collim_80S_2Yr_2km.png'
      Title = 'CSETN: 2 Yr: Jul2009 to Dec2011'
    end
    3: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_3year_2km_July2009_to_July2012.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_3year_2km_July2009_to_July2012.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_3year_2km_July2009_to_July2012.tiff')
      Yrs = '3'
      NameYr = 'Fig_S17ad_CSETN_Collim_80S_3Yr_2km.png'
      Title = 'CSETN: 3 Yr: Jul2009 to Dec2012'
    end
    4: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_4year_2km_July2009_to_July2013.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_4year_2km_July2009_to_July2013.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_4year_2km_July2009_to_July2013.tiff')
      Yrs = '4'
      NameYr = 'Fig_S18ad_CSETN_Collim_80S_4Yr_2km.png'
      Title = 'CSETN: 4 Yrs: July2009 to Dec2013'
    end
    5: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_5year_2km_July2009_to_July2014.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_5year_2km_July2009_to_July2014.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_5year_2km_July2009_to_July2014.tiff')
      Yrs = '5'
      NameYr = 'Fig_S19ad_CSETN_Collim_80S_5Yr_2km.png'
      Title = 'CSETN: 5 Yr: Jul2009 to Dec2014'
    end
    6: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_6year_2km_July2009_to_July2015.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_6year_2km_July2009_to_July2015.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_6year_2km_July2009_to_July2015.tiff')
      Yrs = '6'
      NameYr = 'Fig_S20ad_CSETN_Collim_80S_6Yr_2km.png'
      Title = 'CSETN: 6 Yr: Jul2009 to Dec2015'
    end
    7: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_7year_2km_July2009_to_July2016.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_7year_2km_July2009_to_July2016.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_7year_2km_July2009_to_July2016.tiff')
      Yrs = '7'
      NameYr = 'Fig_S21ad_CSETN_Collim_80S_7Yr_2km.png'
      Title = 'CSETN: 7.5 Yr: Jul2009 to Dec2016'
    end
    8: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_8year_2km_July2009_to_July2017.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_8year_2km_July2009_to_July2017.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_8year_2km_July2009_to_July2017.tiff')
      Yrs = '8'
      NameYr = 'Fig_S22ad_CSETN_Collim_80S_8Yr_2km.png'
      Title = 'CSETN: 8 Yr: Jul2009 to Dec2017'
    end
    9: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_9year_2km_July2009_to_July2018.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_9year_2km_July2009_to_July2018.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_9year_2km_July2009_to_July2018.tiff')
      Yrs = '9'
      NameYr = 'Fig_S23ad_CSETN_Collim_80S_9Yr_2km.png'
      Title = 'CSETN: 9 Yr: Jul2009 to Dec2018'
    end
    10: begin
      csetnmk = read_tiff(CountsMaps+'CSETNMK_25km_10+year_2km_July2009_to_Dec2019.tiff')
      csetnsk = read_tiff(CountsMaps+'CSETNSK_25km_10+year_2km_July2009_to_Dec2019.tiff')
      csetnct1 =read_tiff(CountsMaps+'CSETNCT1_25km_10+year_2km_July2009_to_Dec2019.tiff')
      Yrs = '10.5'
      NameYr = 'Fig_4ad_Fig_S13b_Fig_S25_CSETN_Collim_80S_10+Yr_2km.png'
      Title = 'CSETN: 10.5 Yr: Jul2009 to Dec2019'
    end
    else:  print, ' Error in yr ID???'
  endcase
return
end


; PRO CSETN_MAP82_FIG
; Task:  Generates pixel detection significance maps for the region above 80S.
;        Maps have 2 km x 2 km pixels. Process derives uncollimated_LN_backround
;        (LN = lunar neutron) by smoothing the CSETN LN neutron
;        suppression map. 
;        
;        Collimated_Map = CSETN_LN_Suppression - uncollimagted_LN_background
;
;        Pixel detection significance:
;
;        (Collimated_WEH_Map - Uncollimated_WEH_Map) / sqrt(TotStdErrMean^2. + TotStdErrMean^2)
;         
Pro CSETN_map82_FIG

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets
  
  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
;  restore,SaveSets + 'OceanHaline.sav'
  restore,SaveSets + 'CivitisCT.sav'
  r_curr = r
  g_curr = g
  b_curr = b
  tvlct, r_curr,g_curr,b_curr
  signif_pixels_yr = fltarr(10)
  signif_clusters_yr = fltarr(10)
  Fract_in_PSR = fltarr(10)
  Fract_near_PSR = fltarr(10)
  
  signif_Pix_counts = fltarr(1440,1440)
  
  psr_pixel_stats = fltarr(10,2)
  nopsr_pixel_stats = psr_pixel_stats
  all_pixel_stats = psr_pixel_stats
  
  yrs_3 = [10,7, 3]
  for yrids = 0,0 do begin
  yr = yrs_3(yrids) 
 
  
  illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
  latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
  longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
  dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
  
  Read_Images,yr
  
  ; Find where CSETN and SETN counts pixels eq 0, needed for year 1 map and normalize to make cps
  dt = where(csetnct1 gt 0)
  csetnsk(dt) = (csetnsk(dt)/csetnct1(dt))
  illum2 = illum
     
  adim = 1440
  
  ; Median smooth the topography for use in Detection Significance Figure (edge preserving, despeckling)
  dem1 = congrid(dem,1440,1440,/center)
  d = dem1
  for a = 3,1436 do begin
    for b = 3,1436 do begin
      d(a,b) = median(dem1(a-2:a+2,b-2:b+2))
    endfor
  endfor
  dem1 = d
  dem1(719-4:719+4,719) = max(dem1)
  dem1(719,719-4:719+4) = max(dem1)
  
  longigrid = congrid(longigrid, adim,adim,/center)
  latitgrid = congrid(latitgrid, adim,adim,/center)
  
  illum_orig = illum
  illum = congrid(illum,adim,adim)

; Indices to Segment maps > -80 S with Pixel range
  st = 599
  ed = 1440-599
  
; keep = pixels in 80 to S-pole band, kill = pixels > -80 S, used to make a circular S polar map below
  keep = where(latitgrid le -82.)
  kill = where(latitgrid gt -82.)
 
str_yr = '10_'

csetn_coll_WEH = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')
csetn_UnColl_WEH = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_UnColl_WEH_sig220_UL4.tiff')
csetn_coll_sigma = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_Statist_sig220_UL4.tiff')

csetn_coll_weh(kill) = max(csetn_coll_weh(keep))
csetn_uncoll_weh(kill) = max(csetn_uncoll_weh(keep))
csetn_coll_sigma(kill) = max(csetn_coll_sigma(keep))
dem1(kill) = max(dem1(keep))

csetn_coll_weh = csetn_coll_weh(st:ed,st:ed)
csetn_UnColl_Weh = csetn_uncoll_weh(st:ed,st:ed)
csetn_coll_sigma = csetn_coll_sigma(st:ed,st:ed)
latit = latitgrid(st:ed,st:ed)
illum = illum(st:ed,st:ed)
dem1 = dem1(st:ed,st:ed)

max_csetn_coll_weh = max(csetn_coll_weh)
min_csetn_coll_weh = min(csetn_coll_weh)

max_csetn_coll_sigma = max(csetn_coll_sigma)
min_csetn_coll_sigma = min(csetn_coll_sigma)

max_csetn_Uncoll_weh = max(csetn_UnColl_WEH)
min_csetn_Uncoll_weh = min(csetn_UnColl_WEH)

print, max_csetn_coll_weh
print, min_csetn_coll_weh
print, max_csetn_uncoll_weh
print, min_csetn_uncoll_weh

; The following lines regrid the segmented maps to a 400 x 400 pixel regions

latit = congrid(latit,400,400)
illum = congrid(illum,400,400)
dem1 = congrid(dem1,400,400)
csig= congrid(csetn_coll_sigma,400,400)
cweh = congrid(csetn_coll_weh,400,400)
ucweh = congrid(csetn_uncoll_weh,400,400)

ill = illum
noillum = where(illum eq 0. and latit lt -82.)
ill(*) = 0.
ill(noillum) = 1.
illpsr = segmentregions(ill,5.)
no = where(illpsr lt 0.)
illpsr(no) = 0.
id1 = where(illpsr gt 0.)
illpsr(id1) = 1.
k = fltarr(3,3)
k(*) = 1.
k(0,2) = 0.
k(2,0) = 0.
k(2,2) = 0.
k(0,0) = 0.
dilatepsr = dilate(illpsr,k)
psr_outlines = dilatepsr-illpsr

; Generate a color map to hold the detection significance map.
csetn_sig = bytarr(n_elements(cWEH(*,0)),n_elements(cweh(0,*)))
csetn_sig(*) = 255


; Rebin significant pixels as a 1-Tailed test of increasing significance above 3 sigma
; Colors are assigned in detection significance binsizes of 1.0
psr = where(psr_outlines gt 0)

kill = where(latit gt -82)

window,xsize=400,ysize=400
tvscl, csig
csig1 = tvrd(true=1)
tvscl, cweh
cweh1 = tvrd(true=1)
tvscl, ucweh
ucweh1 = tvrd(true=1)

    ; Assign color values to each significance category, PSR outlines are black
    ; Defined to start at green on color bar
    csetn_sig(kill) = 255
    csetn_sig(psr) = 0

    loadct, 0
    tvscl, dem1
    psrimg = dem1
    psr = where(psr_outlines eq 1)
    psrimg(psr) = 1

    window,xsize=400,ysize=400
    ;dem1 = sqrt(dem1)
    dem1 = bytscl(dem1)
    tvscl, dem1
    dm1 = tvrd(true=1)

    for a = 0,399 do begin
      for b = 0,399 do begin
        if psrimg(a,b) eq 1 then begin
          dm1(0,a,b) = 150
          dm1(1,a,b) = 150
          dm1(2,a,b) = 30
        endif
        if latit(a,b) gt -80 then begin
          dm1(0,a,b) = 255
          dm1(1,a,b) = 255
          dm1(2,a,b) = 255
        endif        
      endfor
    endfor

;    dm1(*,199-3:199+3) = 255
;    dm1(*,199,199-3:199) = 255

restore,SaveSets + 'CivitisCT.sav'
tvlct, r,g,b

; Define the 80S latitude line subtracting a dilated and non-dilated latitude map.
keep = where(latit ge -82)
lt80 = csetn_sig
lt80(*) = 0
lt80(keep) = 1
plt80 = dilate(lt80,k)
line80 = plt80-lt80
outline80 = where(line80 eq 1)
csetn_sig(outline80) = 0

; Define the color image array to store the South pole CSETN maps
big_img = bytarr(3,965,965)
big_img(*) = 255

; Define the color image array for the detection significance maps
sig_img = bytarr(965,480,3)
sig_img(*) = 255

big_img2 = bytarr(460,460,3)

; Display and capture the maps for display as RGB images

for a = 0,399 do begin
  for b = 0,399 do begin
    if psrimg(a,b) eq 1 then begin
      cweh1(0,a,b) = 40
      cweh1(1,a,b) = 40
      cweh1(2,a,b) = 40
      csig1(0,a,b) = 40
      csig1(1,a,b) = 40
      csig1(2,a,b) = 40
      ucweh1(0,a,b) = 40
      ucweh1(1,a,b) = 40
      ucweh1(2,a,b) = 40
    endif
    if latit(a,b) gt -82 then begin
      cweh1(0,a,b) = 255
      cweh1(1,a,b) = 255
      cweh1(2,a,b) = 255
      csig1(0,a,b) = 255
      csig1(1,a,b) = 255
      csig1(2,a,b) = 255
      ucweh1(0,a,b) = 255
      ucweh1(1,a,b) = 255
      ucweh1(2,a,b) = 255
    endif
  endfor
endfor

; Define the color bars for each map
b = bytarr(30,256)
for i = 0,255 do b(*,i) = i
window,xsize=30,ysize=256
tvscl, b
bb = tvrd(true=1)
tv, bb, true=1
write_png,Figures+'Fig_3_Civitis_ColorMap.png', bb

; Load in the CSETN maps to big_img
big_img(*,60:459,60:459) = ucweh1
big_img(*,60:459,520:919) = cweh1
big_img(*,526:925,60:459) = dm1
big_img(*,526:925,520:919) = csig1

big_img2(*,*,*) = 255
window,xsize=965, ysize=965
tv, big_img,true=1

write_png,Figures+'Fig_3ad_maps.png',big_img

read,xxx
endfor


return
end 


