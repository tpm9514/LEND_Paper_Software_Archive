
; Program:   Cabeus_Haworth_PSR_Profiles.pro
; Comment:   The program generates the Fig. 5c1-4 and 5d1-4.   Study of Cabeus, Haworth
; PSR longitude profiles of neutron suppression, WEH conversion, Topo, Max Temp
;  
; Input:  WEH, Statistics, Topo, Collimated, Uncollimated, Illum,  PSR maps
; Output:  Generates Fig. 5c1-4,5d1-4
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
; 
; Description:
;   This code segments areas including permanently shadowed regions from binary maps.
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


; PRO Cabeus_ Haworth_PSR_Profiles
; Task:  Generates pixel profiles through the largest PSRs.  
;
Pro Cabeus_Haworth_PSR_Profiles

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


; Read supporting maps
    illum = read_tiff(LRO_Maps+'SP_Illum_45S_90S.tiff')
    latitgrid = read_tiff(LRO_Maps+'Latitude_45S_90S.tiff')
    longigrid = read_tiff(LRO_Maps+'Longitude_45S_90S.tiff')
    dem = read_tiff(LRO_Maps+'SP_LOLA_DEM_45S_90S.tiff')
    slope_azi = read_tiff(LRO_Maps+'Slope_Azi_angle.tiff')

    
; Used in paper

    csetn_coll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_WEH_sig220_UL4.tiff')
    csetn_uncoll_weh = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_UnColl_WEH_sig220_UL4.tiff')
    csetn_coll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_Coll_Suppress_sig220_UL4.tiff')
    csetn_uncoll_supp = read_tiff(DerivedMaps+'Test_0.50_MAPS70_CSETN_10_New_UnColl_Suppress_sig220_UL4.tiff')
    TempMax = read_tiff(LRO_Maps+'SP_Diviner_TempMax.tiff')
    
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

; PSR Identifiers in the Ill1 map.  The identifier is unique for all pixels in the PSR
    cabid = 375
    hawid = 533
    shoid = 595
    fauid = 664.

; select all PSR pixels for a given PSR;  Haworth is special case as its pixels span 0 deg lon.
    cab = where(ill1 eq cabid)
    haw = where(ill1 eq hawid)
    haw1 = where(ill1 eq hawid and longigrid le 360. and longigrid gt 340.)
    haw2 = where(ill1 eq hawid and longigrid gt 0. and longigrid le 20.)
    sho = where(ill1 eq shoid)
    fau = where(ill1 eq fauid)
    
print, '--- Preliminary Check of PSRs Max WEH locations: '
print, 'Coll WEH, coll + uncoll suppression, coll suppression, tempmax, latit, longi'
    print, 'Cabeus max WEH stats:'
    cabwmax = max(csetn_coll_weh(cab),idc)
    cabtotsupp = (1.-csetn_coll_supp(cab(idc)))+(1.-csetn_uncoll_supp(cab(idc)))
    cabtotsupp = 1.-cabtotsupp    
    print, cabwmax, cabtotsupp, csetn_coll_supp(cab(idc)), tempmax(cab(idc)), latitgrid(cab(idc)), longigrid(cab(idc))    


    print, 'Haworth max WEH stats:'
    Hawwmax = max(csetn_coll_weh(haw), idh)
    hawtotsupp = (1.-csetn_coll_supp(haw(idh)))+(1.-csetn_uncoll_supp(haw(idh)))
    hawtotsupp = 1.-hawtotsupp
    print,  hawwmax, hawtotsupp, csetn_coll_supp(haw(idh)),tempmax(haw(idh)), latitgrid(haw(idh)), longigrid(haw(idh))

 
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

    ; Cabeus peak
    stx = 36
    sty = 157.
    edx = 100.
    edy = 100.
    vx = [stx]
    vy = [sty]

    slope = (sty-edy)/(stx-edx)

    vecx = findgen(100)
    vecy = (findgen(100)*(-1.)*(57./64.))
    vecy = vecy + (sty-vecy(36))    
    dist = sqrt((vecy-vecy(36))^2. + (vecx-vecx(36))^2.)*2.

    d = where(dist lt 80.)
    dist = dist(d)
    dist(36:57) = dist(36:57)*(-1.)
    vecx = vecx(d)
    vecy = vecy(d)
        
    cc = plot(ltg(vecx,vecy),cs_tots(vecx,vecy),layout=[1,4,1],position=[0.08,0.72,0.50,0.96],xshowtext=0,dim = [850, 830], $
      font_name='Times',font_size=12, thick=2, $
      ytitle='Tot. n Supp. $\epsilon$',title='5a1 to 5a4)  Cabeus-1 PSR Profiles, $\it A$ to $\it A $',color='grey',yrange=[0.73, 1.02], xrange=[-87.2,-81.6],/current)
      z = fltarr(n_elements(vecx))
      z(*) = 1.
    cz = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
    shu = plot(ltg(vecx,vecy),csu1s(vecx,vecy), name='UL', thick=2, color='red',/overplot)
    p = where(ill(vecx,vecy) eq 0)
    spsr2 = plot(ltg(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)),cs_tots(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)), color = 'black',thick=5,/overplot)
    spsr22 = plot(ltg(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)),cs_tots(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)), color = 'black',thick=5,/overplot)
    spsr222 = plot(ltg(vecx(p(2:6)),vecy(p(2:6))),cs_tots(vecx(p(2:6)),vecy(p(2:6))), color = 'black',thick=5,/overplot)

    CollimOnlys = csc1s
    cc1 = plot(ltg(vecx,vecy),CollimOnlys(vecx,vecy),layout=[1,4,2],position=[0.08,0.46,0.50,0.7], $
          font_name='Times',font_size=12, xshowtext=0, thick=2, $
      ytitle='Collim. n Supp. $\epsilon$!ICL',color='grey',yrange=[0.73, 1.02], xrange=[-87.2,-81.6],/current)
      cz2 = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
      spsr3 = plot(ltg(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)),collimonlys(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)), color = 'black',thick=5,/overplot)
      spsr33 = plot(ltg(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)),collimonlys(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)), color = 'black',thick=5,/overplot)
      spsr333 = plot(ltg(vecx(p(2:6)),vecy(p(2:6))),collimonlys(vecx(p(2:6)),vecy(p(2:6))), color = 'black',thick=5,/overplot)
      
      tm = plot(ltg(vecx,vecy),tmax(vecx,vecy),xrange=[-87.2,-81.6], color='grey', position =[0.08,0.26,0.50,0.44], $
        font_name='Times',font_size=12, xshowtext=0, thick=2, ytitle='Max. Temp. [K]', yrange=[30.,320.],/current)
      spsr4 = plot(ltg(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)),tmax(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)), color = 'black',thick=5,/overplot)
      spsr44 = plot(ltg(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)),tmax(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)), color = 'black',thick=5,/overplot)
      spsr444 = plot(ltg(vecx(p(2:6)),vecy(p(2:6))),tmax(vecx(p(2:6)),vecy(p(2:6))), color = 'black',thick=5,/overplot)
      
      tms = plot(ltg(vecx,vecy),tmax_smooth(vecx,vecy),xrange=[-87.2,-81.6], color='grey', thick=2, linestyle=2,/overplot)
        
      topo = plot(ltg(vecx,vecy),dem1(vecx,vecy),xrange=[-87.2,-81.6], color='grey', position =[0.08,0.06,0.50,0.24], $
        font_name='Times',font_size=12, thick=2, ytitle='Altit. [km]',  $
        xtitle='Latitude $\deg$, along 311$\deg$ E.', yrange=[-4.3,0.7], /current)
      spsr5 = plot(ltg(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)),dem1(vecx(p(0):p(0)+1),vecy(p(0):p(0)+1)), color='black',thick=5,/overplot)
      spsr55 = plot(ltg(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)),dem1(vecx(p(1):p(1)+1),vecy(p(1):p(1)+1)), color = 'black',thick=5,/overplot)
      spsr555 = plot(ltg(vecx(p(2:6)),vecy(p(2:6))),dem1(vecx(p(2:6)),vecy(p(2:6))), color = 'black',thick=5,/overplot)

      profile = collimonlys(vecx,vecy)
      minsupp = min(profile, mid)

      tm = tmax(vecx,vecy)

print, 'Validate that the Max Supp location is in the profile'
print, ltg(vecx,vecy),lng(vecx,vecy), cs_tots(vecx,vecy)
      
; read,' End Cabeus:', xxx

      pdc = 0.
      zec = 0.
      pdcsm = 0.
      zecsm = 0.
      cab_corr = r_correlate(collimonlys(vecx,vecy),tmax(vecx,vecy), probd=pdc, zd = zec)
      cab_corr_sm = r_correlate(collimonlys(vecx,vecy),tmax_smooth(vecx,vecy), probd=pdcsm,zd=zecsm)

      vcx = vecx(10:49)
      vcy = vecy(10:49)
      cab_corr1 = r_correlate(collimonlys(vcx,vcy),tmax(vcx,vcy))
      cab_corr_sm1 = r_correlate(collimonlys(vcx,vcy),tmax_smooth(vcx,vcy))
            
      c = cab_corr_sm
      n = n_elements(dist)

      print, 'Cabeus temp vs collim_supp correlations, orig and smoothed:',cab_corr, cab_corr_sm, n
      print, 'Cabeus Suppression Only: temp vs collim_supp correlations, orig and smoothed, n:',cab_corr1, cab_corr_sm1, n
            
head = ['Figs 5a1 to 5a4: Cabeus PSR Profiles:']
col_head = [ 'latit: deg', 'longi: deg E', 'total_n_Supp: eps', 'uncoll_n_Supp: eps', 'collim_only_n_Supp: eps', 'Altit: km', 'TempMax: K', 'TempMx_Sm: K']
write_csv, csv + 'Fig_5a14_Cabeus_Coll_UnColl_Profiles.csv',ltg(vecx,vecy),lng(vecx,vecy),  $
   cs_tots(vecx,vecy),csu1s(vecx,vecy), collimonlys(vecx,vecy), dem1(vecx,vecy), tmax(vecx,vecy), tmax_smooth(vecx, vecy), $
   header = col_head, table_header=head

 print, 'Table 1 Entries:  Cabeus-1 Row', $
   latitgrid(cab(idc)), ',', longigrid(cab(idc)), ',', csetn_coll_supp(cab(idc)), ',', $
   csetn_coll_weh(cab(idc)), ',', tempmax(cab(idc)), ',', cab_corr_sm(0), ',', cab_corr_sm(1)


; Haworth peak
stx = 103
sty = 141.
edx = 100.
edy = 100.
print, ltg(stx,sty), lng(stx,sty),csc1w(stx,sty)

y = findgen(80) + 100.
x = findgen(80)/((sty-100.)/(stx-100.)) + 100.
vecx = x
vecy = y
ds = sqrt((x(0)-x(79))^2. + (y(0)-y(79))^2.)
dist = findgen(80)*(ds/79.)*2.


cch = plot(ltg(vecx,vecy),cs_tots(vecx,vecy),layout=[1,4,5],position=[0.56,0.72,0.98,0.96],xshowtext=0,font_name='Times',font_size=12, $
          thick=2,title='5b1 to 5b4)  Haworth PSR Profile, $\itP$ to $\itB$',color='grey',yrange=[0.78, 1.04], xrange=[-90.2,-84.6],/current)
z = fltarr(n_elements(dist))
z(*) = 1.
czh = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)
shuh = plot(ltg(vecx,vecy),csu1s(vecx,vecy), name='UL', thick=2, color='red',/overplot)
p = where(ill(vecx,vecy) eq 0)
psx = vecx(p)
psy = vecy(p)
sclp = dist(p)
spsr2h = plot(ltg(psx,psy),cs_tots(psx,psy), color = 'black',thick=5,/overplot)

CollimOnlys = csc1s
cc1 = plot(ltg(vecx,vecy),CollimOnlys(vecx,vecy), layout=[1,4,6], font_name='Times',font_size=12, xrange=[-90.2, -84.6], color='grey', $
  thick=2, position=[0.56,0.46,0.98,0.7], xshowtext=0, yrange=[0.78, 1.04],/current)
czch = plot(ltg(vecx,vecy),z,linestyle=2.,thick=2, /overplot)

spsr = plot(ltg(psx,psy),CollimOnlys(psx,psy), color = 'black',thick=5,/overplot)

tm = plot(ltg(vecx,vecy),tmax(vecx,vecy),xrange=[-90.2,-84.6], color='grey', position =[0.56,0.26,0.98,0.44], $
  font_name='Times',xshowtext=0, font_size=12, thick=2, yrange=[30.,330.],/current)
tm1 = plot(ltg(psx,psy),tmax(psx,psy),thick=4,color='black',/overplot);

tms = plot(ltg(vecx,vecy),tmax_smooth(vecx,vecy),xrange=[-90.2,-84.6], color='grey', position =[0.56,0.26,0.98,0.44], $
  font_name='Times',font_size=12, thick=2, linestyle=2, yrange=[30.,330.],/overplot)

topo = plot(ltg(vecx,vecy),dem1(vecx,vecy),xrange=[-90.2, -84.6], color='grey', position =[0.56,0.06,0.98,0.24], $
   xtitle='Latitude $\deg$, along 4.3$\deg$ E.',font_name='Times',font_size=12, thick=2, $
   yrange=[-4.5,5.5],/current)
topop1 = plot(ltg(psx,psy),dem1(psx,psy),thick=4,color='black',/overplot)

profile = collimonlys(vecx,vecy)
minsupp = min(profile, mid)
print, 'Haworth Max Supp, Latitude, Longi: = ',Ltg(vecx(mid),vecy(mid)),lng(vecx(mid),vecy(mid))
print, 'Haworth Max Supp, Supp, MaxTemp: =', minsupp, tmax(vecx(mid),vecy(mid))

print, 'Validate that the Max Supp location is in the profile'
print, ltg(vecx,vecy),lng(vecx,vecy), cs_tots(vecx,vecy)


pdh = 0.
zeh = 0.
pdhsm = 0.
zehsm = 0.
Haw_corr = r_correlate(collimonlys(vecx,vecy),tmax(vecx,vecy))
Haw_corr_sm = r_correlate(collimonlys(vecx,vecy),tmax_smooth(vecx,vecy))
vcx = vecx(10:49)
vcy = vecy(10:49)
Haw_corr1 = r_correlate(collimonlys(vcx,vcy),tmax(vcx,vcy))
Haw_corr_sm1 = r_correlate(collimonlys(vcx,vcy),tmax_smooth(vcx,vcy))

h = haw_corr_sm
n = n_elements(dist)
print, 'Haworth temp vs collim_supp correlations, orig and smoothed, n:',haw_corr, haw_corr_sm, n
print, 'Haworth1 Suppression Only: temp vs collim_supp correlations, orig and smoothed, n:',haw_corr1, haw_corr_sm1, n
print, '

cc.save,Figures+'Fig_5ab_Cabeus_Haw_Profiles.png'
    
head = ['Figs 5b1 to 5b4: Haworth PSR Profiles:']
col_head = [ 'latit: deg', 'longi: deg E', 'total_Supp: eps', 'uncoll_Supp: eps', 'collim_only_Supp: eps', 'Altit: km', 'TempMax: K', 'TempMx_Sm: K']
write_csv, csv + 'Fig_5b14_Haworth_Coll_UnColl_Profiles.csv',ltg(vecx,vecy),lng(vecx,vecy),  $
   cs_tots(vecx,vecy),csu1s(vecx,vecy), collimonlys(vecx,vecy), dem1(vecx,vecy), tmax(vecx,vecy), tmax_smooth(vecx, vecy), $
   header = col_head, table_header=head

 print, 'Table 1 Entries:  Haworth Row', $
   latitgrid(haw(idh)), ',', longigrid(haw(idh)), ',', csetn_coll_supp(haw(idh)), ',', $
   csetn_coll_weh(haw(idh)), ',', tempmax(haw(idh)), ',', haw_corr_sm(0), ',', haw_corr_sm(1)


read,'End of run, Cabeus and Haw Profiles, 1 to end: ',xxx

endfor

return
end
