
; Program:   Shoemaker_Faustini_PSR_Profiles.pro
; Comment:   The program generates the Fig. 5a1-4 and 5b1-4.   Study of Shoemaker, Faustini
; PSR longitude profiles of neutron suppression, WEH conversion, Topo, Max Temp
;  
; Input:  WEH, Statistics, Topo, Collimated, Uncollimated, Illum,  PSR maps
; Output:  Generates Fig. 5a1-4, 5b1-4
;
; Written by:   Timothy P. McClanahan, NASA GSFC              March 14, 2024
; Contact:   timothy.p.mcclanahan@nasa.gov
;
; The software is provided as is.   It is part of a data and software repository supporting a
; Planetary Science Journal publication:
;“Evidence for Widespread Hydrogen Sequestration within the Moon’s South Polar Cold Traps”

; The software is designated by NASA as a Creative Commons Zero repository.   NASA’s Software
; Release #:  GSC-19191-1.

; Repository  DOI:  10.5281/zenodo.10027812

; The software written in the Interactive Data Language 8.8.3, which requires a
; commercial license to run available at:   https://www.nv5geospatialsoftware.com/Products/IDL
;

; Function CSETN_FOV_KERN
;    Function generates 2-D CSETN FOV kernels for map smoothing based on neutron transport modeling
;    Weights in collimated part of smoothing kernel = 0, outside 3-sigma = 0.
Function CSETN_FOV_Kern, sigmac, sigmau, km_per_pixel,wid_k,Comp

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets
  ; Field of View for kernels
  ;
  restore,SaveSets+'CSETN_FOV2.sav'

;  csetn_collimated = csetn_collimated2
  csetn_uncollimated = csetn_uncollimated2

;  sigmac_max = sigmac*3.
  sigmau_max = sigmau*4.
;  Coll_FOV = fltarr(wid_k, wid_k)
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


; Function GENPooledStdDev_1
;  Generate the weighted average of variances for a pixel as F(FOV weight and N observs of pixels in FOV)
;  n_samps: set of # of pixel observations
;  s2_samps: set of pixel variances
;  coll_fov_weights: pixel FOV weights
;  return, the weighted pooled standard deviation of pixel
Function GenPooledStdDev_1, n_samps, s2_samps, coll_fov_weights
  ; Variance = Weighted average as a function of pixels FOV weight * n_samples / pixel
  V_c = total(((n_samps-1.)*coll_fov_weights)*s2_samps)/total((n_samps-1)*coll_fov_weights)
  ; StdDev
  S_c = sqrt(V_c)
  return,s_c
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


; PRO Shoemaker_Faustini_PSR_Profiles
; Task:  Generates pixel profiles through the largest PSRs.  
;
Pro Shoemaker_Faustini_PSR_Profiles

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

  yrs_3 = [10,7, 3]
  for yrids = 0,0 do begin
    yr = yrs_3(yrids)


    illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
    latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
    longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
    dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
    
; Used in paper

    csetn_coll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')
    csetn_uncoll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_UnColl_WEH_sig220_UL4.tiff')
    csetn_coll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_Suppress_sig220_UL4.tiff')
    csetn_uncoll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_UnColl_Suppress_sig220_UL4.tiff')
    TempMax = read_tiff(LRO_Maps+'SP_Diviner_TempMax.tiff')


;    PSR_Density = read_tiff(DerivedMaps+'PSR_Predict.tiff')
    
    Read_Images,yr

    ; Find where CSETN and SETN counts pixels eq 0, needed for year 1 map and normalize to make cps
    dt = where(csetnct1 gt 0)
    csetnsk(dt) = (csetnsk(dt)/csetnct1(dt))
    illum2 = illum

    adim = 1440

    ; Median smooth the topography for use in Detection Significance Figure (edge preserving, despeckling)

    longigrid = congrid(longigrid, adim,adim,/center)
    latitgrid = congrid(latitgrid, adim,adim,/center)
    dem = congrid(dem, adim,adim,/center)

    latitgrid_orig = latitgrid
    illum_orig = illum
    illum = congrid(illum,adim,adim)

    dt = where(illum eq 0)
    ill = illum
    ill(*) = 0
    ill(dt) = 1
    ill1 = segmentregions(ill,1)

    cabid = 375
    hawid = 533
    shoid = 595
    fauid = 664.

    cab = where(ill1 eq cabid)
    haw = where(ill1 eq hawid)
    haw1 = where(ill1 eq hawid and longigrid le 360. and longigrid gt 340.)
    haw2 = where(ill1 eq hawid and longigrid gt 0. and longigrid le 20.)
    sho = where(ill1 eq shoid)
    fau = where(ill1 eq fauid)

print, '--- Preliminary Check of PSRs Max WEH locations: '
print, 'Coll WEH, coll + uncoll suppression, coll suppression, tempmax, latit, longi'

    print, 'Shoemaker max WEH stats:'
    Shoewmax = max(csetn_coll_weh(sho),ids)
    shototsupp = (1.-csetn_coll_supp(sho(ids)))+(1.-csetn_uncoll_supp(sho(ids)))
    shototsupp = 1.-shototsupp

    print, shoewmax, shototsupp,csetn_coll_supp(sho(ids)), tempmax(sho(ids)), latitgrid(sho(ids)), longigrid(sho(ids))
    
    print, 'Faustini max WEH stats:'
    Fauwmax = max(csetn_coll_weh(fau),idf)
    Fautotsupp = (1.-csetn_coll_supp(fau(idf)))+(1.-csetn_uncoll_supp(fau(idf)))
    Fautotsupp = 1.-fautotsupp
    print, Fauwmax, fautotsupp, csetn_coll_supp(fau(idf)), tempmax(fau(idf)), latitgrid(fau(idf)), longigrid(fau(idf))

    ; Indices to Segment maps > -80 S with Pixel range
    st = 566
    ed = 1440-566

    ; keep = pixels in 80 to S-pole band, kill = pixels > -80 S, used to make a circular S polar map below
    keep = where(latitgrid le -80.)
    kill = where(latitgrid gt -80.)

; for paper
    lo = 620
    hi = 820
    
;    lo = 580
;    hi = 860
 
    c = 56.
    sig1 = 2.25  ;  This one in the paper
    x1 = gaussian_function([sig1,sig1],width=111.)
    x1 = x1 / total(x1)
    coll_fov = x1
    
    tmax_smooth = convol(tempmax, coll_fov, /center, invalid=-999)
  
    ltg = latitgrid(lo:hi,lo:hi)
    lng = longigrid(lo:hi,lo:hi)
    csc1w = csetn_coll_weh(lo:hi,lo:hi)
    csu1w = csetn_Uncoll_weh(lo:hi,lo:hi)
    tmax = tempmax(lo:hi,lo:hi)
    tmax_smooth = tmax_smooth(lo:hi,lo:hi)

    csc1s = csetn_coll_Supp(lo:hi,lo:hi)
    csu1s = csetn_Uncoll_Supp(lo:hi,lo:hi)
    ill = illum(lo:hi,lo:hi)
    dem1 = dem(lo:hi,lo:hi)

    cs_totw = csc1w + csu1w
    cs_tots = (1.-csc1s) + (1.-csu1s)
    
    cs_tots = 1.-cs_tots

;
;
; Shoemaker
;
;
cx = 100
cy = 100
nel = 55

dist = findgen(nel)*sqrt(2.)*2.
vecx = findgen(nel)+cx
vecy = findgen(nel)+cy
        
cc = plot(ltg(vecx,vecy),cs_tots(vecx,vecy),layout=[1,4,1],position=[0.08,0.72,0.50,0.96],xshowtext=0,dim = [850, 830], $
       font_name='Times',font_size=12, thick=2, ytitle='Tot. n Supp. $\epsilon$', $
       title='5c1 to 5c4)  Shoemaker PSR Profiles, $\it P$ to $\itC$',color='grey',yrange=[0.73, 1.02], xrange=[-90.2,-84.7],/current)
z = fltarr(n_elements(vecx))
z(*) = 1.
cz = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
shu = plot(ltg(vecx,vecy),csu1s(vecx,vecy), name='UL', thick=2, color='red',/overplot)
p = where(ill(vecx,vecy) eq 0)
psx = vecx(p)
psy = vecy(p)
sclp = dist(p)
spsr2 = plot(ltg(psx(0:12),psy(0:12)),cs_tots(psx(0:12),psy(0:12)), color = 'black',thick=5,/overplot)
spsr22 = plot(ltg(psx(13:14),psy(13:14)),cs_tots(psx(13:14),psy(13:14)), color = 'black',thick=5,/overplot)


CollimOnlys = csc1s
cc1 = plot(ltg(vecx,vecy),CollimOnlys(vecx,vecy),layout=[1,4,2],position=[0.08,0.46,0.50,0.7], $
          font_name='Times',font_size=12, xshowtext=0, thick=2, $
          ytitle='Coll. n Supp. $\epsilon$!ICL',color='grey',yrange=[0.73, 1.02], xrange=[-90.2,-84.7],/current)
cz2 = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
spsr2 = plot(ltg(psx(0:12),psy(0:12)),collimonlys(psx(0:12),psy(0:12)), color = 'black',thick=5,/overplot)
spsr22 = plot(ltg(psx(13:14),psy(13:14)),collimonlys(psx(13:14),psy(13:14)), color = 'black',thick=5,/overplot)
      
tm = plot(ltg(vecx,vecy),tmax(vecx,vecy),xrange=[-90.2,-84.7], color='grey', position =[0.08,0.26,0.50,0.44], $
        xshowtext=0,font_name='Times',font_size=12, thick=2, yrange=[30.,340.],/current)
spsr2 = plot(ltg(psx(0:12),psy(0:12)),tmax(psx(0:12),psy(0:12)), ytitle='Max. Temp. [K]',color = 'black',thick=5,/overplot)
spsr22 = plot(ltg(psx(13:14),psy(13:14)),tmax(psx(13:14),psy(13:14)), color = 'black',thick=5,/overplot)
      
tms = plot(ltg(vecx,vecy),tmax_smooth(vecx,vecy),xrange=[-90.2,-84.7], color='grey', position =[0.08,0.26,0.50,0.44], $
font_name='Times',font_size=12, thick=2, linestyle=2,/overplot)
        
topo = plot(ltg(vecx,vecy),dem1(vecx,vecy),xrange=[-90.2,-84.7], color='grey', position =[0.08,0.06,0.50,0.24], $
  font_name='Times', ytitle='Altit. [km]',xtitle='Latitude $\deg$, along 45.0 $\deg$ E.', $
  font_size=12, thick=2, yrange=[-5.0,3.6],/current)
spsr2 = plot(ltg(psx(0:12),psy(0:12)),dem1(psx(0:12),psy(0:12)), color = 'black',thick=5,/overplot)
spsr22 = plot(ltg(psx(13:14),psy(13:14)),dem1(psx(13:14),psy(13:14)), color = 'black',thick=5,/overplot)

profile = collimonlys(vecx,vecy)
minsupp = min(profile, mid)
print, 'Shoemaker Max Supp, Latitude, Longi: = ',Ltg(vecx(mid),vecy(mid)),lng(vecx(mid),vecy(mid))
print, 'Shoemaker Max Supp, Supp, MaxTemp: =', minsupp, tmax(vecx(mid),vecy(mid))


shoe_corr = r_correlate(collimonlys(vecx,vecy),tmax(vecx,vecy))
shoe_corr_sm = r_correlate(collimonlys(vecx,vecy),tmax_smooth(vecx,vecy))
nel = n_elements(vecx)
vcx = vecx(10:nel-10)
vcy = vecy(10:nel-10)
Shoe_corr1 = r_correlate(collimonlys(vcx,vcy),tmax(vcx,vcy))
Shoe_corr_sm1 = r_correlate(collimonlys(vcx,vcy),tmax_smooth(vcx,vcy))
      
      
sh = shoe_corr_sm
n = n_elements(dist)
print, 'Shoemaker temp vs collim_supp correlations, orig and smoothed:',shoe_corr, shoe_corr_sm
print, 'Shoemaker Suppression Only: temp vs collim_supp correlations, orig and smoothed:',shoe_corr1, shoe_corr_sm1
print, 'n  =', n
      
head = ['Figs 5c1 to 5c4: Shoemaker PSR Profiles:']
col_head = [ 'latit: deg', 'longi: deg E', 'total_Supp: eps', 'uncoll_Supp: eps', 'collim_only_Supp: eps', 'Altit: km', 'TempMax: K', 'TempMx_Sm: K']
write_csv, csv + 'Fig_5c14_Shoemaker_Coll_UnColl_Profiles.csv',ltg(vecx,vecy),lng(vecx,vecy),  $
   cs_tots(vecx,vecy),csu1s(vecx,vecy), collimonlys(vecx,vecy), dem1(vecx,vecy), tmax(vecx,vecy), tmax_smooth(vecx, vecy), $
   header = col_head, table_header=head

print, 'Validate the Max WEH location is in the profile and lat, lon'
print, ltg(vecx,vecy),lng(vecx,vecy), cs_tots(vecx,vecy)

print, 'Table 1 Entries:  Shoemaker Row', $
  latitgrid(sho(ids)), ',', longigrid(sho(ids)), ',', csetn_coll_supp(sho(ids)), ',', $
  csetn_coll_weh(sho(ids)), ',', tempmax(sho(ids)), ',', shoe_corr_sm(0), ',', shoe_corr_sm(1)


;
;
; Faustini PSR
;
;
stx = 146
sty = 107.
edx = 100.
edy = 100.

slope = (sty-edy)/(stx-edx)
x = findgen(78)+100.
y = (findgen(78)*slope) + 100.
vecx = x
vecy = y
dist = findgen(78)*(1.+slope)*2.


cc = plot(ltg(vecx,vecy),cs_tots(vecx,vecy),layout=[1,4,5],position=[0.56,0.72,0.98,0.96],xshowtext=0, font_name='Times', $
        font_size=12, thick=2, yrange=[0.78, 1.02], xrange=[-90.3,-84.6], $
        title='5d1 to 5d4)  Faustini PSR Profiles, $\it P$ to $\itD $',color='grey',/current)
z = fltarr(n_elements(vecx))
z(*) = 1.
cz = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
shu = plot(ltg(vecx,vecy),csu1s(vecx,vecy), name='UL', thick=2, color='red',/overplot)
p = where(ill(vecx,vecy) eq 0)
psx = vecx(p)
psy = vecy(p)
sclp = dist(p)
spsr1 = plot(ltg(psx(0:4),psy(0:4)),cs_tots(psx(0:4),psy(0:4)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(5:7),psy(5:7)),cs_tots(psx(5:7),psy(5:7)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(8:21),psy(8:21)),cs_tots(psx(8:21),psy(8:21)), color = 'black',thick=5,/overplot)

CollimOnlys = csc1s
cc1 = plot(ltg(vecx,vecy),CollimOnlys(vecx,vecy),layout=[1,4,6],position=[0.56,0.46,0.98,0.7], $
  font_name='Times',font_size=12, xshowtext=0, thick=2, $
  color='grey',yrange=[0.78, 1.02], xrange=[-90.3,-84.6],/current)
  cz2 = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
spsr1 = plot(ltg(psx(0:4),psy(0:4)),collimonlys(psx(0:4),psy(0:4)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(5:7),psy(5:7)),collimonlys(psx(5:7),psy(5:7)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(8:21),psy(8:21)),collimonlys(psx(8:21),psy(8:21)), color = 'black',thick=5,/overplot)

tm = plot(ltg(vecx,vecy),tmax(vecx,vecy),xrange=[-90.3,-84.6], layout=[1,4,8],color='grey', position =[0.56,0.26,0.98,0.44], $
  font_name='Times',xshowtext=0, font_size=12, thick=2, yrange=[30.,330.],/current)
spsr1 = plot(ltg(psx(0:4),psy(0:4)),tmax(psx(0:4),psy(0:4)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(5:7),psy(5:7)),tmax(psx(5:7),psy(5:7)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(8:21),psy(8:21)),tmax(psx(8:21),psy(8:21)), color = 'black',thick=5,/overplot)

tms = plot(ltg(vecx,vecy),tmax_smooth(vecx,vecy),xrange=[-90.3,-84.6], color='grey', $
  font_name='Times',font_size=12, thick=2, linestyle=2,/overplot)

topo = plot(ltg(vecx,vecy),dem1(vecx,vecy),xrange=[-90.3,-84.6], layout=[1,4,7],color='grey', position =[0.56,0.06,0.98,0.24], $
    font_name='Times',font_size=12, xtitle='Latitude $\deg$, along 81.2$\deg$ E.', $
    thick=2, yrange=[-4.5,2.0],/current)
spsr1 = plot(ltg(psx(0:4),psy(0:4)),dem1(psx(0:4),psy(0:4)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(5:7),psy(5:7)),dem1(psx(5:7),psy(5:7)), color = 'black',thick=5,/overplot)
spsr2 = plot(ltg(psx(8:21),psy(8:21)),dem1(psx(8:21),psy(8:21)), color = 'black',thick=5,/overplot)

profile = collimonlys(vecx,vecy)
minsupp = min(profile, mid)
print, 'Faustini Max Supp, Latitude, Longi: = ',Ltg(vecx(mid),vecy(mid)),lng(vecx(mid),vecy(mid))
print, 'Faustini Max Supp, Supp, MaxTemp: =', minsupp, tmax(vecx(mid),vecy(mid))

print, 'Validate the Max WEH location is in the Profile and lat, lon'
print, ltg(vecx,vecy),lng(vecx,vecy), cs_tots(vecx,vecy)


faust_corr = r_correlate(collimonlys(vecx,vecy),tmax(vecx,vecy))
faust_corr_sm = r_correlate(collimonlys(vecx,vecy),tmax_smooth(vecx,vecy))
nel = n_elements(vecx)
vcx = vecx(10:nel-10)
vcy = vecy(10:nel-10)
faust_corr1 = r_correlate(collimonlys(vcx,vcy),tmax(vcx,vcy))
faust_corr_sm1 = r_correlate(collimonlys(vcx,vcy),tmax_smooth(vcx,vcy))


n = n_elements(dist)
print, 'Faustini temp vs collim_supp correlations, orig and smoothed:',faust_corr, faust_corr_sm
print, 'Faustini Suppression Only: temp vs collim_supp correlations, orig and smoothed:',faust_corr1, faust_corr_sm1
print, 'n = ', n

cc.save,Figures+'Fig_5cd_Shoe_Faust_Profiles.png'


head = ['Figs 5d1 to 5d4: Faustini PSR Profiles:']
col_head = [ 'latit: deg', 'longi: deg E', 'total_Supp: eps', 'uncoll_Supp: eps', 'collim_only_Supp: eps', 'Altit: km', 'TempMax: K', 'TempMx_Sm: K']
write_csv, csv + 'Fig_5d14_Faustini_Coll_UnColl_Profiles.csv',ltg(vecx,vecy),lng(vecx,vecy),  $
  cs_tots(vecx,vecy),csu1s(vecx,vecy), collimonlys(vecx,vecy), dem1(vecx,vecy), tmax(vecx,vecy), tmax_smooth(vecx, vecy), $
  header = col_head, table_header=head

print, 'Table 1 Entries:  Faustini Row', $
    latitgrid(fau(idf)), ',', longigrid(fau(idf)), ',', csetn_coll_supp(fau(idf)), ',', $
    csetn_coll_weh(fau(idf)), ',', tempmax(fau(idf)), ',', faust_corr_sm(0), ',', faust_corr_sm(1)


read,'End of Shoe n Faust Run, 1 to end:',xxx

endfor

return
end