; Program:   HypothesisInstrumentalBlurring.pro
; Comment:   The program generates the Fig 1. profiles 
; 
; Input:  None - its a diagram
; Output:  Generates Fig 1.
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



Pro MakePlot
  device, decomposed = 0
  device, retain=2
  loadct,4

  @~/merrimac/LENDproc/LEND_Paper/Software/dirsets

  g1 = gaussian_function(15)
  g1 = g1 / total(g1)

  x = fltarr(1800)
  x(*) = 1.0
  x(280:290) = 0.9
  x(400:420) = 0.9
  x(520:560) = 0.9
  x(700:760) = 0.9
  x(920:1011) = 0.9
  x(1150:1290) = 0.9

  y = fltarr(1800)
  y(*) = 0.85
  y(300:340) = 0.8
  y(500:540) = 0.75
  y(700:740) = 0.7
  y(900:950) = 0.65

  xx1 = convol(x,g1,/center)
  yy1 = convol(y,g1,/center)

  gg = x
  gg(100:190) = g1*2.9+1.
  x = x(50:1400)
  xx1 = xx1(50:1400)
  gg = gg(50:1400)
  y = y(50:1400)
  yy1 = yy1(50:1400)

  p = plot(findgen(141)+20,(gg(20:160)+1)+0.03,thick = 2, font_name='Times', font_size='14',name='FOV',color='black', xrange=[-100,1400],yrange = [1.6, 2.15], $
    title='Convol FOV w PSR profiles, Variable (Width, Intensity)',ytitle='Intensity',xtitle=' Distance')
  p1 = plot(x+1, thick=1,color='blue', name='PSR profile 1',/overplot)
  p2 = plot(xx1+1, thick=3,color='grey', name='FOV Resp 1',/overplot)
  p3 = plot(y+1, thick=1,color='blue', name='PSR profile 2',/overplot)
  p4 = plot(yy1+1, thick=3,color='grey', name='FOV Resp 2',/overplot)

  p.save,Figures+'Fig_1_InstrumentalBlurHypothesis.png'
  head = ['Top Profile','Top Response','Bottom Profile', 'Bottom Response']
  write_csv,CSV+'Fig_1_InstrumentalBlurHypothesis.csv',head=head,x+1,xx1+1,y+1,yy1+1
  
  read,xxx
  return
end
