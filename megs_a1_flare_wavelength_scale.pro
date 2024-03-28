pro megs_a1_flare_wavelength_scale

restore, 'Time_After_Gradual_Phase_Image.sav'

A1Locations=fltarr(12,2,1200)

for j=508,1023 do begin
  ;Fe XXII
  Ion=0
  ygauss=gradscale[1373:1382,j]
  nygauss=n_elements(ygauss)
  xgauss=make_array(nygauss,/integer,increment=1,/index)
  fitgauss=GAUSSFIT(xgauss,ygauss,a,nterms=3,chisq=q,sigma=s)
  A1Locations[ion,0,j]=1373.+a[1]
  A1Locations[ion,1,j]=13.58120

  ;Fe XXI
  Ion=1
  ygauss=gradscale[1418:1428,j]
  nygauss=n_elements(ygauss)
  xgauss=make_array(nygauss,/integer,increment=1,/index)
  fitgauss=GAUSSFIT(xgauss,ygauss,a,nterms=3,chisq=q,sigma=s)
  A1Locations[ion,0,j]=1418.+a[1]
  A1Locations[ion,1,j]=12.87520;nm

  ;Fe XX
  Ion=2
  ygauss=gradscale[1463:1472,j]
  nygauss=n_elements(ygauss)
  xgauss=make_array(nygauss,/integer,increment=1,/index)
  fitgauss=GAUSSFIT(xgauss,ygauss,a,nterms=3,chisq=q,sigma=s)
  A1Locations[ion,0,j]=1463.+a[1]
  A1Locations[ion,1,j]=12.18450;nm

  ;Fe XVIII
  Ion=3
  ygauss=gradscale[1649:1659,j]
  nygauss=n_elements(ygauss)
  xgauss=make_array(nygauss,/integer,increment=1,/index)
  fitgauss=GAUSSFIT(xgauss,ygauss,a,nterms=3,chisq=q,sigma=s)
  A1Locations[ion,0,j]=1649.+a[1]
  A1Locations[ion,1,j]=9.39320D;nm
endfor

;Fitting Peak Locations for All Rows
;index for each row
rowval=make_array(1200,/integer,increment=1,/index)

yfits=fltarr(12,1,1024)
x=make_array(1024,/float,increment=1,/index)

Ion=0
gdft=where(A1Locations[ion,0,*] gt 1376 and A1Locations[ion,0,*] lt 1378 and rowval gt 715 and rowval lt 824.8)
xcolfit=rowval[gdft]
ycolfit=A1Locations[ion,0,gdft]
fitcoeff=poly_fit(xcolfit,ycolfit,2,yfit=yfit,chisq=q,sigma=ps)
yfits[ion,0,*]=fitcoeff[2]*(x^2.)+fitcoeff[1]*x+fitcoeff[0]

ion=1
gdft=where(A1Locations[ion,0,*] gt 1421 and A1Locations[ion,0,*] lt 1424 and rowval gt 715.5 and rowval lt 823.2)
xcolfit=rowval[gdft]
ycolfit=A1Locations[ion,0,gdft]
fitcoeff=poly_fit(xcolfit,ycolfit,2,yfit=yfit,sigma=ps)
yfits[ion,0,*]=fitcoeff[2]*(x^2.)+fitcoeff[1]*x+fitcoeff[0]

ion=2
gdft=where(A1Locations[ion,0,*] gt 1466 and A1Locations[ion,0,*] lt 1469 and rowval gt 713.1 and rowval lt 824.3)
xcolfit=rowval[gdft]
ycolfit=A1Locations[ion,0,gdft]
fitcoeff=poly_fit(xcolfit,ycolfit,2,yfit=yfit,sigma=ps)
yfits[ion,0,*]=fitcoeff[2]*(x^2.)+fitcoeff[1]*x+fitcoeff[0]

Ion=3
gdft=where(A1Locations[ion,0,*] gt 1653 and A1Locations[ion,0,*] lt 1656 and rowval gt 719.2 and rowval lt 829.2)
xcolfit=rowval[gdft]
ycolfit=A1Locations[ion,0,gdft]
fitcoeff=poly_fit(xcolfit,ycolfit,2,yfit=yfit,sigma=ps)
yfits[ion,0,*]=fitcoeff[2]*(x^2.)+fitcoeff[1]*x+fitcoeff[0]

yccdA1=make_array(200,2048,/float)
xccdA1=make_array(2048,/float,increment=1,/index)
;for loop to fit the data to make the wavelength scale
saveyfit=fltarr(200,4)
scale_coeffA1=fltarr(200,4)

for j=700,899 do begin
  fitcoeff=poly_fit(yfits[0:3,0,j],A1Locations[0:3,1,j],3,yfit=yfit,chisq=c,sigma=ss)
  scale_coeffA1[j-700,*]=fitcoeff
  saveyfit[j-700,*]=yfit
  yccdA1[j-700,*]=fitcoeff[3]*(xccdA1^3.)+fitcoeff[2]*(xccdA1^2.)+fitcoeff[1]*xccdA1+fitcoeff[0]
endfor

save,A1Locations,file='MEGSA1_Time_After_Gradual_Phase.sav'
save,xccdA1,yccdA1,file='A1_EVE_wavelength_gradual_phase_scale.sav'
save,scale_coeffA1,file='MEGS_A1_Scale_Coeffiecients.sav'

end