; Program:   CSETN_map70.pro
; Comment:   The program generates Collim, UnCollim. Neut. Supp., WEH, Statistic maps as floating point TIFFs
; Input:  CSETN counts, CSETN time, variance maps
; Output:  Generates maps used for evaluations in Section 3.
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
  
;  sigmac_max = sigmac*3.
  sigmau_max = sigmau*4.
  Coll_FOV = fltarr(wid_k, wid_k)
  UnColl_FOV = Coll_FOV

  center = wid_k/2.
  for i = 0,wid_k-1 do begin
    for j = 0,wid_k-1 do begin
      distance = sqrt((float(i)-center)^2+ (float(j)-center)^2.)*km_per_pixel
      if distance lt sigmau_max then UnColl_FOV(i,j) = interpol(csetn_uncollimated,range,distance)
    endfor
  endfor
  
  UnColl_FOV = UnColl_FOV / total(UnColl_FOV)
  
  return, UnColl_FOV
end


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
;    Read in maps from TIFF images based on cumulative coverage from 7/2/2009 to mapping year.
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




; PRO CSETN_map70
; Task:  Generates CSETN COllimated, Uncollimated, Neut. Supp, WEH, Statistics south polar maps.
;        Maps are polar stereographic, 2 km x 2 km pixels, true at the south pole.   Maps are
;        Floating point tiffs.
;        
Pro CSETN_map70

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets
  
  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
  restore,SaveSets + 'OceanHaline.sav'
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
  
  yrs_3 = [10, 7, 4]
  for yrids = 0,0 do begin
  yr = yrs_3(yrids) 
 
  
  illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
  latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
  longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
  dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
  slope_azi = read_tiff(LRO_Maps+'Slope_Azi_angle.tiff')
  slope = read_tiff(LRO_Maps+'TopoSlope.tiff')
  TempMax = read_tiff(LRO_Maps+'SP_Diviner_TempMax.tiff')
  
  adim = 1440
  longigrid = congrid(longigrid, adim,adim,/center)
  latitgrid = congrid(latitgrid, adim,adim,/center)

  Read_Images,yr
  
  ; Find where CSETN and SETN counts pixels eq 0, needed for year 1 map and normalize to make cps
  dt = where(csetnct1 gt 0)
  csetnsk(dt) = (csetnsk(dt)/(csetnct1(dt)-1.))
  illum2 = illum
      
  latitgrid_orig = latitgrid
  illum_orig = illum
  illum = congrid(illum,adim,adim)

; Indices to Segment maps > -80 S with Pixel range
  st = 599
  ed = 1440-599
  
  ; keep = pixels in 80 to S-pole band, kill = pixels > -80 S, used to make a circular S polar map below
  keep = where(latitgrid le -70)
  kill = where(latitgrid gt -70)

  badpix = where(csetnct1 eq 0 and latitgrid le -65,bad_n)
  print, 'Invalid pixels = ',bad_n
  for i = 0,bad_n-1 do begin
   if latitgrid(badpix(i)) le -70 then print, 'Lat, Lon:',latitgrid(badpix(i)),longigrid(badpix(i))
  endfor
     
  csetn_subtract_back = (csetnmk - (0.685+0.339))    ; Isolate Collimated rate: Subtract sum of GCR and Uncollimated (0.685 + 0.339 = 1.024)  
;  csetn_ln_subtract_gcr = (csetnmk - (0.685))  ; Isolate LN rate: Subtracting GCR cps (Sanin et al., 2016)
  
  ; Normalize CSETN background subtracted LN map to background region rate (cps) in Lat band = -65 to -70 latitude
  backgd_cps = where(latitgrid le -65 and latitgrid gt -70 and csetnct1 gt 0.)
  print, 'Back_65_70_Mean = ',mean(csetn_subtract_back(backgd_cps))
  
  ; Calculate n suppression map, normalized to the average backgd cps rate between 65 and 70 S latitude.
  ; This step produces suppression = 1.0 at the background region -65S to -70S.
  csetn_suppress = csetn_subtract_back / mean(csetn_subtract_back(backgd_cps))          ;  n suppression as f(collimated rate)
    
  ; Derive smoothing kernels for CSETN (Uncollimated, Collimated) FOV's
  ; Produces 2-D smoothing kernels with total weights = 1.0.
  UnColl_Fov = CSETN_FOV_kern(5.01,36.03, 2., 150.,1)
  
  ; To generate Fig 2 uncomment these lines below
  ;unc = plot(findgen(150)*2.-150, uncoll_fov(*,75),title='CSETN UnCollimated Smoothing Filter',xtitle='Distance from Nadir, km',ytitle='Pixel Weight', $
  ;           font_name='Times', Font_size=14, thick=1.5,xrange=[-150, 150])
  ;unc.save,Figures+'Fig2_CSETN_UnColl_Filter.png'  
  
  ; Uncomment this section to substitute Gaussians for CSETN's FOV, used in SOM Appendix D.
  c = 56.
;  sig1 = 2.25  ;  This one in the paper
  sig1 = 2.2
  x1 = gaussian_function([sig1,sig1],width=111.)
  x1 = x1 / total(x1)
  coll_fov = x1

  ; For the COLLIMATED:  Isolate collimated and uncollimated n suppression:
  ; 1. Derive n suppression as f(collimated rate)
  ;    Derive the uncollimated neutron suppression at the SP, > -80 S as F(collimated rate)
  
; invalidate pixels if 0 obs
  csetn_supp = csetn_suppress
  inval = where(csetnct1 eq 0)
  csetn_supp(inval) = -999.
  
  csetn_uncoll_suppress1 = convol(csetn_supp, UnColl_FOV, /center,invalid=-999, /normalize)
  ;    Derive the collimated neutron suppression, this includes the uncollimated response as f(collimated rate)
  csetn_coll_suppress1 = convol(csetn_suppress, Coll_FOV, /center,invalid=-999)

; Subtract the scaled Uncollimated 
  csetn_uncoll_suppress1 = ((csetn_uncoll_suppress1 - 1.) * 0.5) + 1.
    
  ; Isolate the collimated only n suppression using a band pass filter.
  ;  The step creates the collimated only neutron suppression , background region suppression = 1.0  
  ;  A) or B) control if ULN is included in the map or subtracted
  ; A.  Subtract the ULN
   csetn_coll_suppress = (csetn_coll_suppress1 - csetn_uncoll_suppress1) + 1.
;   csetn_coll_suppress = csetn_coll_suppress * 0.98

   window,xsize=600,ysize=600

  ; B.  ULN not subtracted
   csetn_coll_w_uln_suppress = csetn_coll_suppress1
  
  ; Suppression to WEH conversion parms for CSETN collimated neutron suppression
  
  ; Sanin et al. 2016 parms
  at = double(1.2)
  bt = double(0.06)
  ct = double(-0.51)
    
  print, 'Neutron Suppress to WEH wt% Parms: a = ', at,'b = ', bt, 'c = ',ct

  ; To include regolith background WEH using Sanin et al. neutron suppression, subtract neutron suppression = 0.03 or 0.045 wt%
  ;csetn_coll_suppress = csetn_coll_suppress - 0.03
       
; Generate the Standard error map as a function of the Collimated and Uncollimated FOVs.
; Keep: is the set of pixels in the map region > 80deg S.
; These statistical derivations are relevant to the subtraction of the ULN (ONLY)

csetn_coll_stderr = csetn_coll_suppress
csetn_quad_stderr = csetn_coll_suppress
csetn_coll_stderr(*) = 0.
csetn_quad_stderr(*) = 0.

csetn_coll_stderr(keep) = sqrt(csetnsk(keep) / csetnct1(keep))

;csetn_uncoll_stderr(keep) = csetnsk_uncoll_pooled(keep) / sqrt(csetnct_uncoll_pooled(keep))

; Generate csetn std_error in quadrature to include the collimated and collimated as proxy for uncollimated Std_errors)
; A. IF ULN is subracted this operation is required, then add in quadrature
csetn_quad_stderr(keep) = sqrt(csetn_coll_stderr(keep)^2 + csetn_coll_stderr(keep)^2.)

; B. IF ULN is not subtracted this operation is used (does not add uncerts in quadrature)
;csetn_quad_stderr(keep) = csetn_coll_stderr(keep)

; Convert Neutron Suppression to WEH abundance, declare several maps
csetn_coll_weh = csetn_quad_stderr
csetn_coll_weh(*) = 0
csetn_uncoll_weh = csetn_coll_weh
csetn_coll_w_uln_WEH = csetn_coll_weh
csetn_ln_weh = csetn_coll_weh
csetn_coll_weh_unc = csetn_coll_weh
csetn_coll_weh_unc1 = csetn_coll_weh
csetn_coll_weh_unc2 = csetn_coll_weh
csetn_coll_weh_unc_avg = csetn_coll_weh
csetn_coll_supp_unc = csetn_coll_weh
csetn_coll_supp_unc1 = csetn_coll_weh
csetn_coll_sigma = csetn_coll_weh

; Convert collimated n suppression map to WEH, using collimated conversion parms from SOM.
CSETN_coll_weh(keep) = (-1.*(at)+sqrt((at^2.)+4*bt*(((csetn_coll_suppress(keep))^(1/ct))-1.)))/(2.*bt)

; Convert collimated n suppression map to WEH, using collimated conversion parms from SOM.
CSETN_Uncoll_WEH(keep) = (-1.*(at)+sqrt((at^2.)+4*bt*(((csetn_Uncoll_suppress1(keep))^(1/ct))-1.)))/(2.*bt)

; Convert collimated with uln suppression map to WEH, using collimated conversion parms from SOM.
CSETN_coll_w_uln_weh(keep) = (-1.*(at)+sqrt((at^2.)+4*bt*(((csetn_Uncoll_suppress1(keep))^(1/ct))-1.)))/(2.*bt)


; Generate output maps as floating point tiffs.  Maps are made in the DerivedMaps directory

window,1,xsize=600,ysize=600
tvscl, csetn_coll_weh(400:1000,400:1000)
window,2,xsize=600,ysize=600
tvscl, csetn_uncoll_weh(400:1000,400:1000)
str_yr = '10_'

nme = 'Test_0.50_MAPS70_CSETN_'+str_yr+'New_Coll_WEH_sig220_UL4.tiff'
descrip = 'File Name: '+nme
write_tiff,DerivedMaps+nme,csetn_coll_WEH,compression=0,/float ; Sigma above background

nme = 'Test_0.50_MAPS70_CSETN_'+str_yr+'New_Coll_Suppress_sig220_UL4.tiff'
descrip = 'File Name: '+nme
write_tiff,DerivedMaps+nme,csetn_coll_suppress,compression=0,/float ; Sigma above background

nme = 'Test_0.50_MAPS70_CSETN_'+str_yr+'New_UnColl_WEH_sig220_UL4.tiff'
descrip = 'File Name: '+nme
write_tiff,DerivedMaps+nme,csetn_uncoll_weh,compression=0,/float; Sigma above background

nme = 'Test_0.50_MAPS70_CSETN_'+str_yr+'New_UnColl_Suppress_sig220_UL4.tiff'
descrip = 'File Name: '+nme
write_tiff,DerivedMaps+nme,csetn_uncoll_suppress1,compression=0,/float ; Sigma above background;

nme = 'Test_0.50_MAPS70_CSETN_'+str_yr+'New_Coll_Statist_sig220_UL4.tiff'
descrip = 'File Name: '+nme
write_tiff,DerivedMaps+nme,csetn_quad_stderr,compression=0,/float ; Sigma above background

read,'PRogram End',xxx

endfor

return
end 


