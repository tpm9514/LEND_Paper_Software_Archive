; Program:   civitisct.pro
; Comment:   The program generates the civitis color map used in Fig. 3ad
; Input:  Cividis rgb color map from https://www.ncl.ucar.edu/Document/Graphics/ColorTables/cividis.shtml
; Output:  Color table and color bar
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


Pro CivitisCT
openr, lun, 'CividisCT.txt',/get_lun
str = ' '
r = fltarr(256)
g = fltarr(256)
b = fltarr(256)

readf, lun, str
readf, lun, str
for i = 0,255 do begin
  readf, lun, str
  vars = str_sep(str, ' ')
  r(i) = vars(0)
  g(i) = vars(1)
  b(i) = vars(2)
  print, str, r(i), g(i), b(i)
endfor

r = byte(r*255.)
g = byte(g*255.)
b = byte(b*255.)

save, file='~/merrimac/LENDproc/LEND_Paper/Data/SaveSets/CivitisCT.sav',r,g,b


read,xxx
return
end
