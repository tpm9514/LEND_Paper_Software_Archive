; Program:   gcr_correlate2.pro
; Comment:   The program is a theoretical study to determine if the detected Collimated and Uncollimated neutron flux is positively 
; correlated.  The study shows that GCR dependent variation during the solar cycle correlates CSETN's collimated and uncollimated neutrons
;
; Input:  A series of time profiles from Poisson distributed RVs that emulate the neutron emission flux of the last solar cycle.  
; 
; Output:  Shows that the detected Collimated and uncollimated lunar neutron emission flux is positively correlated.   We can't tell
; how much only to state that its positive.  Statistical uncertainties in the paper assume that the Collimated and Uncollimated 
; neutrons are from independent processes.   
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




; Description:
; This code is models the hypothesis that GCR variation induces LEND's ULN and CLN to be positively
; correlated.
; 
; Author:  Tim McClanahan     Date:  Aug 15, 2021
; 
Pro GCR_Correlate

; Trials report arrays means
ng = fltarr(1000)
con1 = ng
con2 = ng
con3 = ng
con4 = ng
con5 = ng

; Trials report arraysstddev 
sdcon1 = ng
sdcon2 = ng
sdcon3 = ng
sdcon4 = ng
sdcon5 = ng

; Trials
for j = 0,99 do begin

; C* = CLN time series, U* = ULN time series
dim = 2000000
c1 = fltarr(dim)
u1 = fltarr(dim)
c2 = fltarr(dim)
u2 = fltarr(dim)
c3 = fltarr(dim)
u3 = fltarr(dim)
c4 = fltarr(dim)
u4 = fltarr(dim)
c5 = fltarr(dim)
u5 = fltarr(dim)
c6 = fltarr(dim)
u6 = fltarr(dim)

;  Generate average neutron emission flux for a solar cycle
ap2 = sin((findgen(dim) / 320.)) * (-0.25)      ; Makes neutron emission at solar max 25% of solar min 
ap3 = sin((findgen(dim) / 320.)) * (-0.35)      ; Solar max 35 % of solar min
ap4 = sin((findgen(dim) / 320.)) * (-0.45)      ; Solar max is 45% of solar min
ap5 = sin((findgen(dim) / 320.)) * (-0.45) + sin(findgen(dim) / 30.) * 0.05   ; AP4 + added sinusoid of 0.05 amp
ap6 = sin((findgen(dim) / 320.)) * (-0.45) + sin(findgen(dim) / 30.) * 0.1    ; AP4 + added sinusoid of 0.1 amp

;  Case studies generate Poisson random variables as a function of mission time and avg neutron emmission flux

; No GCR variation:  Null Hypothesis
for i = 0,dim-1 do c1(i) = randomu(seed, poisson=1.)
for i = 0,dim-1 do u1(i) = randomu(seed, poisson=1.)
; GCR = -25%
for i = 0,dim-1 do c2(i) = randomu(seed, poisson= 1.+ ap2(i))
for i = 0,dim-1 do u2(i) = randomu(seed, poisson= 1.+ ap2(i))
; GCR = -35%
for i = 0,dim-1 do c3(i) = randomu(seed, poisson= 1.+ ap3(i))
for i = 0,dim-1 do u3(i) = randomu(seed, poisson= 1.+ ap3(i))
; GCR + -45%
for i = 0,dim-1 do c4(i) = randomu(seed, poisson= 1.+ ap4(i))
for i = 0,dim-1 do u4(i) = randomu(seed, poisson= 1.+ ap4(i))

; GCR + -45% + (0.05 amp sinusoid)
for i = 0,dim-1 do c5(i) = randomu(seed, poisson= 1.+ ap5(i))
for i = 0,dim-1 do u5(i) = randomu(seed, poisson= 1.+ ap5(i))

; GCR + -45% + (0.1 amp sinusoid)
for i = 0,dim-1 do c6(i) = randomu(seed, poisson= 1.+ ap6(i))
for i = 0,dim-1 do u6(i) = randomu(seed, poisson= 1.+ ap6(i))

ng(j) = correlate(c1,u1)
con1(j) = correlate(c2,u2)
con2(j) = correlate(c3,u3)
con3(j) = correlate(c4,u4)
con4(j) = correlate(c5,u5)
con5(j) = correlate(c6,u6)

print, 'J = ',j

endfor

nogcr = fltarr(dim)
nogcr(*) = 1
p = plot(nogcr,title='Modeled Lunar n vs Time (LEND era)',xtitle='Time',ytitle='GCR Rate',color='black', $
         name='NoGCR',thick=2, yrange=[0.,1.2])
p0 = plot(1.+ap2,color='purple',thick=2,name='-0.2',/overplot)
p1 = plot(1.+ap3,color='red',thick=2,name='-0.3',/overplot)
p3 = plot(1.+ap5,color='green',thick=2,name='-0.4 + low freq sinusoid',/overplot)
p4 = plot(1.+ap6,color='orange',thick=2,name='-0.4 + high freq sinusoid',/overplot)
p2 = plot(1.+ap4,color='blue',thick=2,name='-0.4',/overplot)

print, ' '
print, 'Means, Std for 100 Evaluations:'
print, 'No GCR Variation:',mean(ng), stddev(ng)
print, '-0.25 GCR:',mean(con1), stddev(con1)
print, '-0.35 GCR:',mean(con2), stddev(con2)
print, '-0.45 GCR:',mean(con3), stddev(con3)
print, '-0.45 + (0.05 Sin):',mean(con4), stddev(con4)
print, '-0.45 + (0.1 Sin):',mean(con5), stddev(con5)


read,xxx

return
end
