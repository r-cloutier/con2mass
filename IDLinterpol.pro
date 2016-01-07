; IDl procedure to interpolate between planet absolute magnitude and 
; mass. The IDL interpol function is better than similar python
; functions for this
;
; >>   gdl
; GDL> .compile IDLinterpol.pro
; GDL> IDLinterpol

pro IDLinterpol

; Get model data to interpolate over
readcol, 'foridl_mass_mag.dat', modelmass, modelmag
; Get planet magnitudes 
readcol, 'foridl_planetmag.dat', absmag
; Interpolate to get planet mass limits
massdet = interpol(modelmass, modelmag, absmag, /quadratic)
; Save data to file to be read back by python
openw, 1, 'forpython_planetmass.dat'
printf, 1, massdet, format='(F10.6)'
close, 1


stop
end
