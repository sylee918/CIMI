;--------------------------------------------------------------------------
;                              iono_pot.pro
;
; Procedure polts the ionospheric potential and/or Birkeland current, 
; conductance output from ram**.f, crc,**.f90 and cimi.f90
; 
; Modification history
; May 31, 2016
;  - Add plotting of potential at polar latitudes
; December 5, 2016
;  - Add calculation of Bob Robinson's conductance by using R2 current.
;--------------------------------------------------------------------------

Pro conductance_model_v2,fac0,facm1,facp1,mlt, sigmap,sigmah,enflux,avge
  ;  conductance provided by Bob Robinson
  ;  input the fac at any magnetic latitude (fac0). pos. -> upward current
  ;  the fac one degree equatorward (facm1)  
  ;  the fac one degree poleward (facp1)
  ;  and the magnetic local time as an integer from 0 to 23 (mlt)
  ; output:
  ;sigmap=Pedersen conductance
  ;hoverp=Hall to Pdersen ratio
  ;sigmah=Hall conductance
  ;enflux=precipitating particle energy flux in ergs/cm2-s

  ;  set the linear fit coefficients for the conductance model:
 spm0=[1.90641, 2.25115, 2.57007, 2.85271, 3.09056, 3.27706, 3.40765, 3.47969, $
      3.49254, 3.44751, 3.34787, 3.19886, 3.00767, 2.78347, 2.53740, 2.28253, $
      2.03393, 1.80862, 1.62557, 1.50574, 1.47203, 1.54932, 1.76444, 2.14620, 1.90641]
spm1=[-2.48804, -2.43332, -2.57739, -2.85514, -3.20905, -3.58924, -3.95343, -4.26700, $
       -4.50291, -4.64178, -4.67182, -4.58889, -4.39646, -4.10561, -3.73508, -3.31119, $
       -2.86790, -2.44681, -2.09712, -1.87565, -1.84686, -2.08283, -2.66325, -3.67544, -2.48804]
spp0=[0.321099, 0.447197, 0.517054, 0.542299, 0.533535, 0.500332, 0.451231, 0.393743, $
      0.334348, 0.278496, 0.230609, 0.194077, 0.171259, 0.163485, 0.171057, 0.193242, $
      0.228282, 0.273386, 0.324733, 0.377474, 0.425726, 0.462579, 0.480095, 0.469298, .321099]
spp1=[13.0462, 17.5440, 20.3049, 21.6573, 21.9016, 21.3101, 20.1273, 18.5696, $
      16.8254, 15.0550, 13.3910, 11.9378, 10.7718, 9.94146, 9.46723, 9.34159, $
       9.52900, 9.96593, 10.5609, 11.1943, 11.7188, 11.9588, 11.7108, 10.7433, 13.0462]
hpm=[1.00137, 1.01065, 1.01222, 1.00684, 0.995311, 0.978486, 0.957258, 0.932567, $
      0.905402, 0.876798, 0.847838, 0.819651, 0.793412, 0.770343, 0.751715, 0.738842, $ 
      0.733089, 0.735864, 0.748624, 0.772872, 0.810158, 0.862080, 0.930280, 1.01645, .0905402]
hpp=[0.944843, 1.04363, 1.11876, 1.17208, 1.20554, 1.22109, 1.22077, 1.20662, $
      1.18077, 1.14538, 1.10265, 1.05483, 1.00424, 0.953226, 0.904178, 0.859546, $ 
      0.821822, 0.793546, 0.777307, 0.775739, 0.791526, 0.827396, 0.886130, 0.970551, .944843]

    pedm1=0.
    pedp1=0.
    ;  if fac is changing sign, use the average values on either side
    if fac0*facm1 lt 0. then begin
      if facm1 lt 0. then pedm1=spm0(mlt)+spm1(mlt)*facm1
      if facm1 ge 0. then pedm1=spp0(mlt)+spp1(mlt)*facm1
      if facp1 lt 0. then pedp1=spm0(mlt)+spm1(mlt)*facp1
      if facp1 ge 0. then pedp1=spp0(mlt)+spp1(mlt)*facp1
      sigmap=.5*(pedm1+pedp1)
      hoverp=.5*(hpm(mlt)+hpp(mlt))  ;Hall conductance determined by average of ratios for pos. & neg. fac
      sigmah=sigmap*hoverp
    endif else begin
      if fac0 lt -.1 then begin
        sigmap=spm0(mlt)+spm1(mlt)*fac0
        sigmah=sigmap*hpm(mlt)
        hoverp=hpm(mlt)
      endif
      if fac0 gt .1 then begin
        sigmap=spp0(mlt)+spp1(mlt)*fac0
        sigmah=sigmap*hpp(mlt)
        hoverp=hpp(mlt)
      endif
      if fac0 ge -.1 and fac0 le .1 then begin
        pedm=spm0(mlt)+spm1(mlt)*fac0
        pedp=spp0(mlt)+spp1(mlt)*fac0
        sigmap=.5*(pedm+pedp)
        hoverp=.5*(hpm(mlt)+hpp(mlt))
        sigmah=sigmap*hoverp
      endif
    endelse
    avge=((1/.45)*hoverp)^(1./.85)
    if avge gt 0 then enflux=(sigmap/(40.*avge/(16.+avge*avge)))^2
    if avge le 0 then enflux=0
end   ; end of conductance_model_v2
 
;-------------  Main program -------------------------------------------------

close,/all

filein=' '  
read,' enter the file head of potential (i.e., 2013d151) => ', filein
print,' overlay? E-field(0), field-aligned current(1), cimi sigmaP(2) => '
read, '                 Robinson J||-driven sigmaP(3),      sigmaH(4) => ',ibir
read,' Sun direction? (1)Sun is left,  (2)Sun is up  => ',iSun
mtitle=' '
outfile=filein+'_poti'
if (ibir eq 0) then outfile=outfile+'_Efield'
if (ibir eq 1) then outfile=outfile+'_birk'
if (ibir eq 2) then outfile=outfile+'_cond'
if (ibir eq 3) then outfile=outfile+'_RcondP'
if (ibir eq 4) then outfile=outfile+'_RcondH'

; get colors for color table palette
if (ibir ne 1) then loadct,33      ;  blue - red 
if (ibir eq 1) then loadct,70      ;  red - white - blue
tvlct, red, green, blue, /get

; Read potential and conductance data 
openr,1,filein+'.pot'
readf,1,rc
readf,1,idim  
readf,1,irh   
readf,1,jdim 
readf,1,mdim           ; no. of longitudinal invariant
lath=fltarr(irh)   
lat2D=fltarr(idim,jdim)  &   BriN=lat2D   &   potentc=lat2D
readf,1,lat2D             ; latitude in degree
readf,1,lath              ; latitude in degree
readf,1,potentc         
readf,1,BriN            
itot=idim+irh     &   clat=fltarr(itot,jdim+1)    &  clatr=clat
for j=0,jdim do clat(idim:itot-1,j)=90.-lath(0:irh-1)   ; colatitude in degree
clat(0:idim-1,0:jdim-1)=90.-lat2D(0:idim-1,0:jdim-1)    ; colatitude in degree
clat(0:idim-1,jdim)=clat(0:idim-1,0)  
clatr(*,*)=clat(*,*)*!pi/180.                           ; colatitude in radian
mlt=fltarr(jdim)  &   ain=mlt 
phi=fltarr(jdim+1)   &   xb=phi   &    yb=phi
potent=fltarr(itot,jdim+1)   &  xx=potent    &    yy=potent   
pedpsi=potent    &   birk=potent   &   Efield=potent
potent1=fltarr(idim,jdim)    &  pedpsi1=potent1  &  dummy=potent1 
birk1=potent1    
potenth=fltarr(irh,jdim)
dummyy=fltarr(idim,jdim,mdim)
read,' enter number of time => ',ntime
hour=fltarr(ntime)    
potenta=fltarr(ntime,itot,jdim)   &   potenta(*,*,*)=0.
Efielda=fltarr(ntime,itot,jdim)   &   Efielda(*,*,*)=0.
pedpsia=potenta   &  birka=potenta
xxa=fltarr(ntime,itot,jdim+1)   &  yya=xxa      
iba=intarr(jdim)   &   ibamx=intarr(ntime)

; Constants needed for E field calculation when ibir=0
Rem=6.371e6     ; Earth radius in m
dphi2=2.*2.*!pi/jdim

for n=0,ntime-1 do begin    
    readf,1,hour1
    print,'n,hour ',n,hour1
    readf,1,iba
    readf,1,mlt                                       ; in hour, zero -> midnight
    readf,1,dummy        ; latS
    readf,1,dummy        ; phiS
    readf,1,dummy        ; ro
    readf,1,dummy        ; xmlto
    readf,1,dummy        ; bo
    readf,1,dummyy       ; y
    readf,1,potent1,format='(8f10.1)'
    readf,1,potenth,format='(8f10.1)'
    readf,1,pedpsi1      ; pedpsi
    readf,1,birk1        ; birk1 is current density into both hemispheres
    for j=0,jdim-1 do phi(j)=mlt(j)*!pi/12.           ; in radian,             
    phi(jdim)=phi(0)+2.*!pi
    if (iSun eq 2) then phi(*)=phi(*)-!pi/2.          ; zero -> dawn
    ibamx(n)=max(iba)
    hour(n)=hour1
    potenta(n,0:idim-1,*)=potent1(0:idim-1,*)/1000.     ; potential in kV
    potenta(n,idim:itot-1,*)=potenth(0:irh-1,*)/1000.     ; potential in kV
    if (ibir eq 1) then birka(n,0:idim-1,*)=0.5*birk1(0:idim-1,*) 
    if (ibir eq 2) then pedpsia(n,0:idim-1,*)=pedpsi1(0:idim-1,*)
    if (ibir ge 3) then begin
       for j=0,jdim-1 do begin
           imlt=round(mlt(j))
           if (imlt gt 23) then imlt=imlt-24
           for i=0,idim-1 do begin
               latm1=lat2D(i,j)-1.
               latp1=lat2D(i,j)+1.
               facInt=[0.,0.]
               if (latm1 ge lat(0) and latp1 le lat(idim-1)) then $
                  facInt=interpol(birk1(*,j),lat,[latm1,latp1])
               fac0=-birk1(i,j)        ; birk1 is pos. for downward current
               facInt(*)=-1.*facInt(*)
               conductance_model_v2,fac0,facInt(0),facInt(1),imlt,sigmap, $
                                    sigmah,enflux,avge
               if (ibir eq 3) then pedpsia(n,i,j)=sigmap  ;from both hemispheres
               if (ibir eq 4) then pedpsia(n,i,j)=sigmah  ;from both hemispheres
           endfor
       endfor
    endif
    ; Calculate E field if ibir eq 0
    if (ibir eq 0) then begin
       for j=0,jdim-1 do begin
           j0=j-1
           if (j0 lt 0) then j0=j0+jdim
           j1=j+1
           if (j1 ge jdim) then j1=j1-jdim
           for i=1,idim-1 do begin
               dsp=Rem*sin(clatr(i,j))*dphi2
               Eph=(potenta(n,i,j1)-potenta(n,i,j0))/dsp
               dsl=Rem*(clatr(i+1,j)-clatr(i-1,j))
               Eth=(potenta(n,i+1,j)-potenta(n,i-1,j))/dsl
               Efielda(n,i,j)=sqrt(Eph*Eph+Eth*Eth)*1e6     ; E in mV/m
           endfor
       endfor
    endif
    ; Convert to X-Y coordinates
    for j=0,jdim do begin
        for i=0,itot-1 do begin
            xxa(n,i,j)=clat(i,j)*cos(phi(j))
            yya(n,i,j)=clat(i,j)*sin(phi(j))
        endfor
    endfor
endfor
close,1

ytit='noon'
if (iSun eq 2) then ytit='dusk'

; Setup for plotting
halfl=50
ntick=2*halfl/10
if (ntick gt 10.) then ntick=ntick/2. 
dlat=2.*halfl/ntick
tickn=strarr(ntick+1)
for i=0,ntick do begin
    lat_1=90.-abs(-halfl+i*dlat)
    tickn(i)=string(fix(lat_1),'(i2)')
endfor
nlevel=59
blvls=fltarr(nlevel)  

; Determine potential levels 
dlevel=10.  
pothm=0.                ; middle level
lvlm=fix(nlevel*0.5)
lvls=intarr(nlevel)   &    clab=lvls
zero_or_1=1
for i=0,nlevel-1 do begin
    lvls(i)=pothm+(i-lvlm)*dlevel
    if (zero_or_1 eq 0) then zero_or_1=1 else zero_or_1=0
    clab(i)=zero_or_1     ; label every other contour
endfor
pc_color=intarr(nlevel)
pc_color(*)=0             ; black contours

; Determine color level
bcolor=intarr(nlevel)
colrmin=20.         &    colrmax=254.
if (ibir eq 0) then colrmin=120.
if (ibir eq 1) then begin
   colrmin=254.
   colrmax=0.
endif
dcolr=(colrmax-colrmin)/(nlevel-1)
for i=0,nlevel-1 do bcolor(i)=fix(colrmin+i*dcolr)

; Determine levels of E field
if (ibir eq 0) then begin
   Emin=0.   
   Emax=70.
   dlevel=(Emax-Emin)/(nlevel-1)
   for i=0,nlevel-1 do blvls(i)=Emin+i*dlevel
endif

; Determine levels of Birkeland current
if (ibir eq 1) then begin
   blmax=1.0 
   dlevel=2*blmax/(nlevel-1)
   for i=0,nlevel-1 do  blvls(i)=-blmax+i*dlevel
endif

; Determine levels of conductance
if (ibir ge 2) then begin
   pedpsimax=30.     
   pedpsimin=0.   
   dlevel=(pedpsimax-pedpsimin)/(nlevel-1)
   for i=0,nlevel-1 do blvls(i)=pedpsimin+i*dlevel
endif

!p.charthick=3
!p.charsize=1.8
set_plot,'ps'
  device,filename=outfile+'.ps',/inches,/color,yoffset=3.5,xoffset=0.3, $
         ysize=5.7,xsize=8.5,bits_per_pixel=24

; make ps file for converting to mpg
; device,filename=outfile+'.ps',/color,/landscape
; ps2pdf *.ps
; gm convert -delay 25 -quality 90 -rotate 180 *pdf poti.mpg

nl=-1
newplot:
nl=nl+1
if (nl ge ntime) then goto, endjob

xi=0.10 
xf=xi+0.48
yi=0.2 
yf=yi+0.72

xx(*,*)=xxa(nl,*,*)
yy(*,*)=yya(nl,*,*)
ibmx=ibamx(nl)

; Draw contours of E field
if (ibir eq 0) then begin
   Efield(*,0:jdim-1)=Efielda(nl,*,0:jdim-1)
   Efield(*,jdim)=Efield(*,0)
   if (min(Efield) lt blvls(0)) then Efield(where(Efield lt blvls(0)))=blvls(0)
   if (max(Efield) gt blvls(nlevel-1)) then $
      Efield(where(Efield gt blvls(nlevel-1)))=blvls(nlevel-1)
   contour,Efield,xx,yy,levels=blvls, $
           pos=[xi,yi,xf,yf],background=255, $
           c_colors=bcolor,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
           xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn,color=0, $
           xstyle=1,ystyle=1,/fill
endif

; Draw contours of Birkeland current
if (ibir eq 1) then begin
   birk(*,0:jdim-1)=birka(nl,*,0:jdim-1)
   birk(*,jdim)=birk(*,0)
   if (min(birk) lt blvls(0)) then birk(where(birk lt blvls(0)))=blvls(0)
   if (max(birk) gt blvls(nlevel-1)) then $
      birk(where(birk gt blvls(nlevel-1)))=blvls(nlevel-1)
   contour,birk,xx,yy,levels=blvls, $
           pos=[xi,yi,xf,yf],background=255, $
           c_colors=bcolor,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
           xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn,color=0, $
           xstyle=1,ystyle=1,/fill
endif 

; Draw contours of conductance
if (ibir ge 2) then begin
   pedpsi(*,0:jdim-1)=pedpsia(nl,*,0:jdim-1)     ; include both hemispheres
   pedpsi(*,jdim)=pedpsi(*,0)
   if (min(pedpsi) lt blvls(0)) then pedpsi(where(pedpsi lt blvls(0)))=blvls(0)
   if (max(pedpsi) gt blvls(nlevel-1)) then $
      pedpsi(where(pedpsi gt blvls(nlevel-1)))=blvls(nlevel-1)
   contour,pedpsi,xx,yy,levels=blvls,$
           pos=[xi,yi,xf,yf],background=255, $
           c_colors=bcolor,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
           xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn,color=0, $
           xstyle=1,ystyle=1,/fill
endif 

; plot latitude-longitude grid
npt=100
dang=2.*!pi/(npt-1)
ang=indgen(npt)*dang
radius=fltarr(npt)
for i=10,halfl,10 do begin                        ; plot latitudes
    radius(*)=float(i)
    oplot,/polar,radius,ang,linestyle=1,color=0
endfor                                             ; plot longitude
for i=0,23,6 do oplot,/polar,[0,halfl],[i*!pi/12.,i*!pi/12.],linestyle=1, $
    color=0

; Plot color bar
x0=xf+0.09
x1=x0+0.03
y0=yi   
dy=0.012
for i=0,nlevel-1 do begin
    y1=y0+i*dy
    y2=y1+dy*1.01
    polyfill,[x0,x1,x1,x0,x0],[y1,y1,y2,y2,y1],color=bcolor(i),/normal
endfor
xyouts,x1,y0,string(blvls(0),'(f4.1)'),color=0,size=1.6,/normal
xyouts,x1,y1-0.01,string(blvls(nlevel-1),'(f4.1)'),color=0,size=1.6,/normal
blabel='Birkeland current (microAmp/m2)'
if (ibir eq 0) then blabel='Electric field (mV/m)'
if (ibir ge 2) then blabel='Pedersen conductance (mho)'
if (ibir eq 4) then blabel='Hall conductance (mho)'
xyouts,x0-0.01,0.5*(y0+y2),blabel,color=0,orientation=90,alignment=0.5, $
       size=1.6,/normal
if (ibir eq 1) then begin
   xyouts,x0,y0-0.03,'upward current',size=1.3,color=0,/normal 
   xyouts,x0,y1+0.03,'downward current',size=1.3,color=0,/normal
endif

; Draw potential contours 
UT=hour(nl)
ihour=fix(UT)
imin=round((UT-float(ihour))*60.)
UT_lab=string(ihour,'(i3.3)')+':'+string(imin,'(i2.2)')
mtitle=filein+'  '+UT_lab
potent(*,0:jdim-1)=potenta(nl,*,0:jdim-1)
potent(*,jdim)=potent(*,0)
contour,potent,xx,yy,levels=lvls, $
        c_labels=clab,pos=[xi,yi,xf,yf], $
        title=mtitle,xtitle='magnetic latitude (degree)',ytitle=ytit,color=0,$
        background=225,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
        c_colors=pc_color,xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn, $
        xstyle=1,ystyle=1,thick=1,c_charsize=1.8,c_charthick=3,/noerase
oplot,xx(ibmx-1,*),yy(ibmx-1,*),color=255      ; draw the polar cap
;polyfill,xx(ibmx-1,*),yy(ibmx-1,*),color=255   ; paint the polor cap
xyouts,xi+0.02,yf-0.06,'MAX='+string(max(potent),'(f5.0)')+' kV',/normal
xyouts,xi+0.02,yf-0.11,'MIN='+string(min(potent),'(f5.0)')+' kV',/normal

erase
goto, newplot

endjob:
device,/close_file

end
