  
; Name:  Pro Bandpass
; Description:  This software derives CSETNs collimated and uncollimated
;    response functions for the spatial bandpass filter. Profiles are derived from our GEANT4 modeling. We use the
;    collimated and uncollimated response to derive smoothing kernels.
;
;
Pro BandPass
  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
  restore,SaveSets + 'OceanHaline.sav'
  str = ''
  
  ; Read in GEANT4 modeled CSETN Collimated and Uncollimated, raw and fit
  openr,lun,'~/merrimac/LENDproc/LEND_Paper/Data/CSETN_FOV_Model/CSETN_FOV_data.txt',/get_lun
  range = fltarr(100)
  colraw = range    ; Collimated Raw
  colfit = range    ; Collimated Fit
  ucolraw = range   ; Uncollimated Raw
  ucolfit = range   ; Uncollimated Fit
  
  i = long(0)
  
  ; Read, in Text files and derive column data
  head = strarr(1)
  readf, lun, str   ;  Header line
  head(0) = str + ', ,'
  readf, lun, str   ;  Column descriptors
  cols = ['Range (km)', 'CSETN Coll Raw(cps)', 'CSETN Coll GaussFit(cps)', 'CSETN Uncoll Raw(cps)', 'CSETN UnColl GaussFit(cps)']
  while not eof(lun) do begin   
    readf, lun, str
    dat = str_sep(strcompress(str),' ')
    range(i) = float(dat(0))
    colraw(i) = float(dat(1))
    colfit(i) = float(dat(3))
    ucolraw(i) = float(dat(2))
    ucolfit(i) = float(dat(4))
    i = i + 1
  endwhile

range = range(0:i-1)
colraw = colraw(0:i-1)
colfit = colfit(0:i-1)
ucolraw = ucolraw(0:i-1)
ucolfit = ucolfit(0:i-1)
TotRaw = colraw + ucolraw
TotFit = colfit + ucolfit

write_csv,CSV + 'CSETN_FOV_Model.csv', range, colraw, colfit, ucolraw, ucolfit, table_header=head, header=cols

; Extract arrays using valid records, 0 to i-1

; Plot Fig showing CSETN Coll and UnColl spatial response as F(distance from nadir)
  p = plot(range,colraw+ucolraw, color='black',title='Spatial Response vs CSETN Total, Coll, UnColl',xtitle='Distance from nadir, km', $
    ytitle='Weight', xrange=[0,40], thick=2., name='Coll+UnColl',/ylog,yrange=[0.0001,0.1])
  p2 = plot(range,colraw, color='blue', thick=2, name='Coll',/overplot)

  p1 = plot(range,ucolraw, color='red',thick=2, name='UnColl',/overplot)


  head = ['Range km from nadir','Collim Resp. (wght)','Collim GFit','UnCollim Resp. (wght)','UnCollim GFit']
  str = 'File contains the CSETN Collimated and Uncollimated spatial response and Gaussian Fits from GEANT4 model, 1D profiles'
  csetndat = read_csv(CSV + 'CSETN_FOV_Model.csv', table_header=head, header= cols)

  pa = plot(range,colfit+ucolfit, color='black', linestyle=2, name='Fit Coll+UnColl',/overplot)
    
  pa2 = plot(range,colfit, color='blue',linestyle=2, name='Fit Coll',/overplot)

  pa1 = plot(range,ucolfit, color='red',linestyle=2, name='Fit UnColl',/overplot)


; Derive the UL and CL component fractions of the CSETN total response as F(distance from nadir)
  ColFitFract = colfit / Totfit    ; Collimated fraction
  UColFitFract = ucolfit / Totfit  ; UnCollimated fraction

; Derive fraction of CL and UL response
  ColResp = ColFitFract*colfit
  UColResp = uColFitFract*ucolfit

; Duplicate the range and the UL and CL responses for negative distances from nadir, for kernel 1D profile  
  r1 = [reverse(range)*(-1.), range]
  C1 = [reverse(ColResp),ColResp]
  U1 = [reverse(UColResp),UColResp]

  normfact = max(colraw+ucolraw)
  curaw = (colraw + ucolraw) / normfact
  crawnorm = colraw / normfact
  ucrawnorm = ucolraw / normfact
  pp = plot(range,curaw, color='black',title='Spatial Response vs CSETN Total, Coll, UnColl',xtitle='Distance from nadir, km', $
    ytitle='Weight', xrange=[0,40], thick=2., name='Coll+UnColl',yrange=[0.0,1.])
  pp2 = plot(range,crawnorm, color='blue', thick=2, name='Coll',/overplot)

  pp1 = plot(range,ucrawnorm, color='red',thick=2, name='UnColl',/overplot)

  head = ['Range km from nadir','Collim Resp. (wght)','Collim GFit','UnCollim Resp. (wght)','UnCollim GFit']
  str = 'File contains the CSETN Collimated and Uncollimated spatial response and Gaussian Fits from GEANT4 model, 1D profiles'
  csetndat = read_csv(CSV + 'CSETN_FOV_Model.csv', table_header=head, header= cols)

;  ppa = plot(range,(colfit+ucolfit)/normfact, color='black', linestyle=2, name='Fit Coll+UnColl',/overplot)

  ppa2 = plot(range,colfit/normfact, color='blue',linestyle=2, name='Fit Coll',/overplot)

  ppa1 = plot(range,ucolfit/normfact, color='red',linestyle=2, name='Fit UnColl',/overplot)
  
  pp.save,Figures + 'Fig_S7a_CSETN_Response_v_nadirdist.png'

  head = ['range km', 'total (cps)', 'collimated normalized(cps)', 'uncoll. normalized(cps)', 'collim gauss fit (cps)', 'uncollim gaussfit (cps)']
  write_csv,CSV+'Fig_S7a_Source_CollimUncollim_Profiles.csv',range, curaw, crawnorm, ucrawnorm, colfit/normfact, ucolfit/normfact,  $
     header=head, table_header=thead
  

  ; Derive the UL and CL component fractions of the CSETN total response as F(distance from nadir)
  ColFitFract = colfit / Totfit    ; Collimated fraction
  UColFitFract = ucolfit / Totfit  ; UnCollimated fraction

  ; Derive fraction of CL and UL response
  ColResp = ColFitFract*colfit
  UColResp = uColFitFract*ucolfit

  ; Duplicate the range and the UL and CL responses for negative distances from nadir, for kernel 1D profile
  r1 = [reverse(range)*(-1.), range]
  C1 = [reverse(ColResp),ColResp]
  U1 = [reverse(UColResp),UColResp]
 
; Save for making uncollimated response;  Used in CSETN_Map70.pro
  head = ['range km', 'collimated normalized(cps)', 'uncoll. normalized(cps)']
  thead = ['CSETN Collimated and Uncollimated response vs distance from nadir point']

; Normalize profiles to their AUC
  csetn_collimated2 = ColResp/total(ColResp)
  csetn_uncollimated2 = UColResp / total(UColResp)

; Write the Coll and UColl profiles out to CSV files and Make SaveSets
  save,file=SaveSets+'CSETN_FOV2.sav', csetn_collimated2, csetn_uncollimated2, range
  write_csv,CSV+'Fig_2_Source_CollimUncollim_Profiles.csv',range, $
    csetn_collimated2, csetn_uncollimated2, header=head, table_header=thead
  
  

; Display Profiles of Collim and Uncollim vs distance from Nadir
; Constrain the range of the UL and CL response
  collr = where(range lt 15.)   ; range of CL < 15 km from nadir
  collrng = range(collr)
  ColResp = ColResp(collr)
  ucollr = where(range lt 150.)  ; range of UL < 150 km from nadir
  ucollrng = range(ucollr)
  UCollResp = UColResp(ucollr)
  ColRespTot = [reverse(colresp),colresp]
  UcolRespTot = [reverse(ucollresp),ucollresp]
  ucr = [reverse(ucollrng)*(-1.), ucollrng]
  cr = [reverse(collrng)*(-1.), collrng]

; Plot 1D UL and CL kernel profiles
  pb1 = plot(cr,ColRespTot, color='black',title='Profile CSETN Collimated Kernel',xtitle='Distance from nadir, km', $
    ytitle='Rel. Weight', xrange=[-50,50], thick=2., name='Coll Resp')
;  pb1.save,Figures+'Fig_S9a_Coll_Spatial_Response.png'
  
  pb2 = plot(ucr,uColRespTot, color='black',title='Profile CSETN UnCollimated Kernel',xtitle='Distance from nadir, km', $
    ytitle='Rel. Weight', xrange=[-170,170], thick=2., name='UnColl Resp')
;  pb2.save,Figures+'Fig_S9b_UnColl_Spatial_Response.png'
    
  head = ['Range km from nadir','Collim Resp. (wght)','UnCollim Resp (wght)']
  str = 'File contains the CSETN Collimated and Uncollimated spatial response for spatial filtering, 1D profile'
;  write_csv,CSV+'Fig_S9ab_Coll_UnColl_Profiles.csv',head=head,table_header=str,r1,C1,U1

  read,'Check vals',xxx
  
  return
  end
  
  
; Name:  Pro SETNResp
; Description:  The objective of this software is to derive SETNs
;    response function. Profiles are derived from our GEANT4 modeling. 
;
Pro SETNResp
  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  common Maps,csetnmk,csetnsk,csetnct1,setnmk,setnsk,setnct1, yrs, NameYr, Title
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr

  device,retain=2
  device,decomposed=0
  restore,SaveSets + 'OceanHaline.sav'
  str = ''

  ; Read in modeled SETN, raw and fit
  openr,lun,CSETN_FOV + 'SETN_FOV_data.txt',/get_lun
  range = fltarr(500)
  ucolraw = range   ; SETN
  ucolfit = range   ; SETN Fit

  i = long(0)

  ; Read, in Text files and derive column data
  readf,lun,str
  readf,lun,str
  while not eof(lun) do begin
    readf, lun, str
    dat = str_sep(strcompress(str),' ')
    range(i) = dat(0)
    ucolraw(i) = dat(1)
    i = i + 1
  endwhile

  ; Extract arrays using valid records, 0 to i-1
  range = range(0:i-1)
  ucolraw = ucolraw(0:i-1)
  gfit = gaussfit(range,ucolraw)
  
  ; Plot Fig 5b showing CSETN Coll and UnColl spatial response as F(distance from nadir)
  p = plot(range,ucolraw, color='black',title='SETN spatial response relative to nadir=0.0 km',xtitle='Distance from nadir, km', $
    ytitle='Weight', xrange=[-160.,160.], thick=2., name='SETN')

  p1 = plot(range,gfit, color='red',thick=2, name='Gauss Fit',/overplot)
  p.save,Figures+'Fig_S7b_SETN_GEANT4_Profiles.png'

  head = ['Range km from nadir','SETN Resp. (wght)','SETN GFit']
  str = 'File contains the SETN Spatial Response w Gaussian Fit from GEANT4 model, 1D profiles'
  write_csv,CSV+'Fig_S7b_SETN_GEANT4_Profiles.csv',head=head,table_header=str,range,ucolraw,gfit

  read,'Check vals',xxx

  return
end