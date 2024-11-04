; Program:   lend_sp_pds_dld_mapping.pro
; Comment:   The program maps LEND PDS RDRDLD daily files in south polar stereographic format.   The mapping assumes
;            a single detector is mapping.   1 to 4 detectors may be independently mapped for each 1 Hz observation.
; 
; Input:  LEND RDRDLD4 daily files.   These files are not provided in the xenodo repository and must be downloaded
; from https://pds-geosciences.wustl.edu/missions/lro/default.htm
; 
; Output:  Generates Counts, neutron counts, variance maps.
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


; DLD_Head reads DLD record header 
Function DLD_Head

  Head = { DLDVers:     uint(0), $
    fileID:      uint(0), $
    TotRecs:     ulong(0), $
    StartSec:    ulong(0), $
    StartSubSec: ulong(0), $
    EndSec:      ulong(0), $
    EndSubSec:   ulong(0), $
    FileName:    bytarr(40)  }
  return, Head
end

; PDS_DLD reads DLD record 
; Records contain corrected DLD counts each detector
function PDS_DLD

  DLD = {  $
    lro_sec:            ulong64(0), $
    utc:                string(' ', format='(a23)'), $
    local_hour  :        byte(0), $
    local_min  :        byte(0), $
    boresight_latitude:  float(0.), $
    boresight_longitude: float(0.), $
    collect_dur:          ulong(0.), $
    stn1_bkd:           float(0.), $
    stn1_cnts:          float(0.), $
    setn_bkd:           float(0.), $
    setn_cnts:          float(0.), $
    stn2_bkd:           float(0.), $
    stn2_cnts:          float(0.), $
    stn3_bkd:           float(0.), $
    stn3_cnts:          float(0.), $
    csetn1_bkd:          float(0.), $
    csetn1_cnts:         float(0.), $
    csetn2_bkd:          float(0.), $
    csetn2_cnts:         float(0.), $
    csetn3_bkd:          float(0.), $
    csetn3_cnts:         float(0.), $
    csetn4_bkd:          float(0.), $
    csetn4_cnts:         float(0.), $
    shen_bkd:          fltarr(16), $
    shen_cnts:         fltarr(16), $
    sun_active:        byte(0), $
    nadir_ptg:         byte(0) }
  return, DLD
end


; XY_COORDS:  generates x,y map position, polar stereographic coordinates.
; parameters:  latitude, longitude of observation
;              position at center pixel = Lon0, lat0, R = radius of moon = 1737.4 km
; returns the x,y map position
;
function xy_coords,lat,lon,lon0,lat0,R
  k=(2*R)/((1.+sin(lat0*!dtor)*sin(lat*!dtor))+(cos(lat0*!dtor)*cos(lat*!dtor)*cos((lon-lon0)*!dtor)));
  x=k*cos(lat*!dtor)*sin((lon-lon0)*!dtor);
  y=k*(cos(lat0*!dtor)*sin(lat*!dtor)-sin(lat0*!dtor)*cos(lat*!dtor)*cos((lon-lon0)*!dtor));
  return, [x,y]
end



; LEND_SP_PDS_DLD_2km_Mapping Code
; Generates south polar maps of LEND CSETN and possible to map SETN observations, 
;
pro LEND_SP_PDS_DLD_2km_Mapping

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets
  device, retain=2
  device, decomposed=0
  loadct, 24
  
  r = 1737.4
  ; window,xsize=1100, ysize=1100

;  d1 = [ 63]
;  d2 = [ 31]
d1 = [77]
d2 = [38]

; Tiff File names:  
;   2km is pixel size for statistical counting maps
;   Dates of coverage are integral from July 2009 (except the 1 year mission map)
;   The maps are 1440 x 1440 south polar stereographic projection.   This creates 2 km pixels, true at the poles.
;   Map latitude range is 46.5 deg S to 90 deg S.
;
;   circle_mission.txt - this file contains only the low altitude circular mission, avg alt 50 km.   Sept 15 2009 to Dec 15, 2011
  savefiles = [ $
    '_25km_1year_prisci_circle_2km_Sept2009_to_Sept2010', $
    '_25km_2year_2km_July2009_to_July2011', $
    '_25km_3year_2km_July2009_to_July2012', $
    '_25km_4year_2km_July2009_to_July2013', $
    '_25km_5year_2km_July2009_to_July2014', $
    '_25km_6year_2km_July2009_to_July2015', $
    '_25km_7year_2km_July2009_to_July2016', $
    '_25km_8year_2km_July2009_to_July2017', $
    '_25km_9year_2km_July2009_to_July2018', $
    '_25km_10+year_2km_July2009_to_Dec2019' ]

  ; List of text files that contain DLD file references
  ; First two mappings study low and high solar potential
;   fnames = ['./MissionLists/circle_mission.txt']
  fnames = [ $
    './1yr_prisci_circle_mission.txt', $    ; 1 year of Science mapping phase: Sept 15, 2009 to Sept 15, 2010
    './2yr_mission.txt', $
    './3yr_mission.txt', $
    './4yr_mission.txt', $
    './5yr_mission.txt', $
    './6yr_mission.txt', $
    './7yr_mission.txt', $
    './8yr_mission.txt', $
    './9yr_mission.txt', $
    './10+yr_mission.txt' ]    ; full mission to date, over 10 years of observations, as above in SaveFiles list

; Directory where DLD data reside, Needs to be modified for other use
  dldfile = ''

; Loops through the mapping runs, accumulated mapping, set b = 0, generate all 10 maps,  set b = 9, is for just year 10+ mapping.
  for b = 9,9 do begin   ; year 10 mapping only

    print, 'In LEND PDS Mapping, CSETN and SETN: ', Savefiles(b), ' ',systime()

; 0.4 km x 0.4 km mappings, uses 25km mapping kernel = 63 x 63 uniform mapping kernel (disk) 
    CSETN = fltarr(7200,7200)
    CSETN(*,*) = 0.
    csetnct = csetn

    setn = csetn
    setnct = csetn

    line = ''

    head = DLD_Head()
    dld = PDS_DLD()
    cnt = long(0)

    lrorbitnum = long(0)
    
; Map = 2km x 2km pixels
    csetnMK = fltarr(1440,1440)
    csetnSK = csetnmk
    csetnct1 = csetnmk
    setnMK = csetnmk
    setnSK = csetnmk
    setnct1 = csetnmk

    ii = long(0)

    openr, filelun,MissionLists + fnames(b),/get_lun

    scilun = -1

    ii = long(0)

    ; Detector norm coeffs from LEND PDS processing doc
    ; https://pds-geosciences.wustl.edu/lro/lro-l-lend-2-edr-v1/lrolen_0xxx/document/lend_data_processing.pdf
    norms = [1.1106, 1.2942, 1.2320, 1.3167]
    
    ; normf = total(norms)/4.  Using normf = 1.275 cps / detector  1.275 * 4 detectors = 5.1 cps tot rate Sanin et al., 2016
    ; Normalizes observations detector standard
    normf = 1.275
    
    ; Loop through listing of DLD files in text files for duration of map
    while not eof(filelun) do begin

      readf,filelun,DLDfile
      if scilun ne -1 then free_lun,scilun
      openr, scilun, dlddir + DLDfile, $
        /get_lun,/swap_endian

      ; Loop through each DLD file
      while not eof(scilun) do begin

        readu, scilun, dld

        lrotime = dld.lro_sec / 256

        ; End of commissioning: LRO=274736160, UTC=2009/09/15 15:36  
        ; Start of elliptical: LRO = 345316653, UTC = 2011/12/11 17:17:33
        ; NOTE:  Do not use COMMISSIONING Phase in mapping
        ;if lrotime ge 345316653 then goto,finish

        ;Nadir and Sunactive flags in https://pds-geosciences.wustl.edu/lro/lro-l-lend-2-edr-v1/lrolen_0xxx/document/lend_rdr_sis.htm
        ; Use flags to assure nadir ptg and sun is quiet and SP latitude lt -46.5

        if dld.boresight_latitude le -46.5 and dld.nadir_ptg eq 1 and dld.sun_active eq 0 then begin

          detct = [dld.csetn1_cnts, dld.csetn2_cnts, dld.csetn3_cnts, dld.csetn4_cnts]
          
          ;
          ; lat, lon to x,y conversion for map with pixel resolution = 2km x 2km resolution
          xy1 = (xy_coords(dld.boresight_latitude,dld.boresight_longitude,0.,-90.,1737.4)/2.) + 719.5
          ; lat, lon to x,y conversion for map with pixel resolution = 10km x 2km resolution
          ; xy1 = (xy_coords(dld.boresight_latitude,dld.boresight_longitude,0.,-90.,1737.4)/10.) + 143.5
          cx1 = round(xy1(0))
          cy1 = round(xy1(1))

          ; Process mm = 0 to 3 CSETN detectors (4) if a valid observ then cnts ge 0
          ; normalize to normf as Function of their base rates (norms)
                             
          for mm = 0,3 do begin
            if detct(mm) ge 0. then begin
              
              ; make normalized observation
              csum = detct(mm)*(normf/norms(mm))
                            
              ; The following code calculates the running variance, counts and means as observations are mapped to each pixel.
              ; The method is based on Knuth et al. "The Art of Computer Programming", Vol 2, page 232, 3rd edition.
              ; The method was extracted to IDL for this study from John D. Cooks website: https://www.johndcook.com/blog/standard_deviation/
              ; Map 2km x 2km pixel resolution for statistical analysis
              csetnct1(cx1,cy1) = csetnct1(cx1,cy1) + 1.    ; Map with weighted observs / pixel
              mk1 = csetnmk(cx1,cy1)                      ; Mean rate at each pixel
              CSETNMK(cx1,cy1) = mk1 + (csum - mk1)/csetnct1(cx1,cy1)    ; Running average of mean rate
              CSETNSK(cx1,cy1) = CSETNSK(cx1,cy1) + (csum - mk1)*(csum - csetnmk(cx1,cy1))    ; Running variance / pixel
            endif
          endfor
        endif

        ; This if bloock is to create SETN Southern maps
        ; 
        ;if dld.setn_cnts ge 0 and dld.boresight_latitude lt -46.5 and dld.nadir_ptg eq 1 and dld.sun_active eq 0 then begin
          
          ; The following code calculates the SETN map running variance, counts and means as observations are mapped to each pixel.
          ; The method is based on Knuth et al. "The Art of Computer Programming", Vol 2, page 232, 3rd edition.
          ; The method was extracted to IDL from a method on John D. Cooks website: https://www.johndcook.com/blog/standard_deviation/
          ; Map 2km x 2km pixel resolution for statistical analysis
          ;setnct1(cx1,cy1) = setnct1(cx1,cy1) + 1.   ; Map with counts / pixel
          ;mk1 = setnmk(cx1,cy1)                      ; Mean CPS rate
          ;SETNMK(cx1,cy1) = mk1 + (csum - mk1)/setnct1(cx1,cy1)  ; Running average cps process
          ;SETNSK(cx1,cy1) = SETNSK(cx1,cy1) + (csum - mk1)*(csum - setnmk(cx1,cy1))  ; Running variance / pixel

        ; endif

        ii = ii + 1
        if ii mod 10000000 eq 0 then print, ii,' ', dld.utc, ' ',lrotime,' ',savefiles(b)

      endwhile

    endwhile

   ; Write tiffs only for the full mission 10+ year maps are large!
   ; if b eq 9 then begin
      ; write_tiff,Counts_Maps + 'CSETN'+savefiles(b)+'.tiff',csetn,compression=0,/float,description='CSETN'+savefiles(b)
      ; write_tiff,Counts_Maps + 'CSETNCT'+savefiles(b)+'.tiff',csetnct,compression=0,/float,description='CSETNCT'+savefiles(b)
      ; write_tiff,Counts_Maps + 'SETN'+savefiles(b)+'.tiff',setn,compression=0,/float,description='SETN'+savefiles(b)
      ; write_tiff,Counts_Maps + 'SETNCT'+savefiles(b)+'.tiff',setnct,compression=0,/float,description='SETNCT'+savefiles(b)
   ; endif

   ; Write tiffs for the running CSETN averages: mean = *mk, variance = *sk, observation counts = *ct maps
   ; All years
    write_tiff,CountsMaps + 'CSETNMK'+savefiles(b)+'.tiff',csetnmk,compression=0,/float,description='CSETNMK'+savefiles(b)
    write_tiff,CountsMaps + 'CSETNSK'+savefiles(b)+'.tiff',csetnsk,compression=0,/float,description='CSETNSK'+savefiles(b)
    write_tiff,CountsMaps + 'CSETNCT1'+savefiles(b)+'.tiff',csetnct1,compression=0,/float,description='CSETNCT1'+savefiles(b)

    ; Write tiffs for the running SETN averages: mean = *mk, variance = *sk, observation counts = *ct maps
    ; All Years
    ;write_tiff,CountsMaps + 'SETNMK'+savefiles(b)+'.tiff',setnmk,compression=0,/float,description='SETNMK'+savefiles(b)
    ;write_tiff,Counts_Maps + 'SETNSK'+savefiles(b)+'.tiff',setnsk,compression=0,/float,description='SETNSK'+savefiles(b)
    ;write_tiff,Counts_Maps + 'SETNCT1'+savefiles(b)+'.tiff',setnct1,compression=0,/float,description='SETNCT1'+savefiles(b)

  endfor

  read,xxx

  return
end


