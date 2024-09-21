;-------------------------------------------------------------------------------
;
;                                plot_fls.pro
; 
; IDL procedure reads *.fls file and plot the pitch-angle averaged flux,
; pressure, density and pitch-angle distribution at the equatorial plane.
;
; Created on 18 September 2007 by Mei-Ching Fok, Code 673, NASA GSFC.
; 
; Modification History
; June 8, 2008
;   * Add the calculation of temperature
; June 10, 2008
;   * The pressure calculation is valid in relativitic limit.
; November 30, 2011
;   * Change calculation of temperature to mean energy
; May 12, 2014
;   * Add the option to plot flux mapped to ionosphere
; July 22, 2016
;   * Add the plotting of energy and pitch-angle distribution at maximum
; September 26, 2016
;   * Add the option to create gif files.
; February 6, 2017
;   * Add an option to make plot for the TWINS storm catalog.
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
pro plot_eqr,var,vmax,vmin,xi,yi,xf,yf,bar_lab,halfl,ro, $
                 mlto,bo,Lstar,Lstar_max,maxL,lat,mlt,x_wsize,y_wsize, $
                 dir,iovp,iopt,icr,blackc,ieq,piso,ii,jj,gthick,rxy
;-------------------------------------------------------------------------------

; Make levels
nlevel=59
if (xi gt 0.5 and iovp eq 4) then nlevel=8
mlvl=(nlevel-1)/2
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
vmid=0
if (piso eq 0) then vmid=0.5*(vmax+vmin)
dlvl1=(vmid-vmin)/mlvl         
dlvl2=(vmax-vmid)/mlvl         
colrmax=253          
colrmin=1
ncolor=colrmax-colrmin+1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
lvl(0)=vmin
for i=0,nlevel-1 do begin
    dlvl=dlvl1
    if (i gt mlvl) then dlvl=dlvl2
    if (i gt 0) then lvl(i)=lvl(i-1)+dlvl
    colr(i)=round(float(colrmin)+i*dcolr)
endfor

; force var equal or gt vmin 
if (min(var) lt vmin) then var(where(var le vmin))=vmin

; polyfill background in white and plot frome in black
if (xi lt 0.5) then polyfill,[0,1,1,0,0],[0,0,1,1,0],color=255,/normal
polyfill,[xi,xf,xf,xi,xi],[yi,yi,yf,yf,yi],color=blackc,/normal

; Setup plot grid and labels
dim=size(ro) 
xo=ro  &  yo=ro
if (ieq eq 1) then begin
   y_title='RE'
   xymax=halfl
   if (dir eq -1) then phi=mlto              ; Sun to the left
   if (dir eq 1) then phi=mlto+!pi           ; Sun to the right
   xo=ro*cos(phi)
   yo=ro*sin(phi)
endif else begin
   y_title='mlat (deg)'
   xymax=50.   
   if (dir eq -1) then phi=mlt               ; Sun to the left
   if (dir eq 1) then phi=mlt+!pi            ; Sun to the right
   for i=0,dim(1)-1 do begin
       xo(i,*)=(90.-lat(i))*cos(phi(*))
       yo(i,*)=(90.-lat(i))*sin(phi(*))
   endfor
endelse
ntick=4
dxy=2.*xymax/ntick
tickn=strarr(ntick+1)
for i=0,ntick do begin
    xy_1=-xymax+i*dxy
    xy_a=abs(xy_1)
    if (ieq ne 1) then xy_a=90.-xy_a
    tickn(i)=string(fix(xy_a),'(i2)')
    if (ieq eq 1 and xy_1*dir lt 0.) then tickn(i)='-'+tickn(i)
endfor

; plot flux or PA anisotropy
color1=blackc
if (iovp eq 0) then color1=255
contour,var,xo,yo,xrange=[-xymax,xymax],yrange=[-xymax,xymax], $
      levels=lvl,c_colors=colr,xstyle=1,ystyle=1,xticks=ntick,yticks=ntick, $
      xtickn=tickn,ytickn=tickn, $
      ytitle=y_title,pos=[xi,yi,xf,yf],charsize=1.3,color=color1,/fill,/noerase
if (ii ge 0) then oplot,[xo(ii,jj),xo(ii,jj)],[yo(ii,jj),yo(ii,jj)],psym=2, $
                  symsize=3,thick=3

; Add color bar and label
xb=xi
yb=yi-0.1*rxy
dx=(xf-xi)/nlevel
y1=yb+0.04
for i=0,nlevel-1 do begin
    x1=xb+i*dx
    x2=x1+dx*1.02
    polyfill,[x1,x2,x2,x1,x1],[yb,yb,y1,y1,yb],color=colr(i),/normal
endfor
xm=0.5*(xb+x2)
xyouts,xm,y1+0.01,bar_lab,alignment=0.5,size=1.7,color=blackc,/normal
if (xi gt 0.5 and iovp eq 4) then begin
   for i=2,nlevel do begin
       x1=xb+(i-1)*dx
       xyouts,x1,yb-0.04,string(i,'(i1)'),alignment=0.5,size=1.7, $
                 color=blackc,/normal
   endfor
endif else begin
   v_min=string(vmin,'(f4.1)')
   v_mid=string(fix(vmid),'(i1)')
   if (vmax lt 100.) then v_max=string(vmax,'(f4.1)')
   if (vmax ge 100.) then v_max=string(vmax,'(e7.1)')
   xyouts,xb,yb-0.04,v_min,alignment=0.5,size=1.7,color=blackc,/normal
   if (piso gt 0) then xyouts,xm,yb-0.03,v_mid,alignment=0.5,size=1.7, $
                              color=blackc,/normal
   xyouts,x2,yb-0.04,v_max,alignment=0.5,size=1.7,color=blackc,/normal
   if (xi gt 0.5 and iopt lt 5 and ieq eq 1) then begin
      xyouts,xb-0.05,yb-0.07,'field aligned',size=1.4,color=blackc,/normal
      xyouts,x2+0.05,yb-0.07,'perpendicular',alignment=1,size=1.4, $
             color=blackc,/normal
   endif
endelse

; fill cirle of radius ro(0,0) centered the earth with color2
color2=colrmin
if (piso gt 0) then color2=colr(mlvl)

; Draw earth 
npt=400                             ; no. of points in a circle
del = 2.0 * !pi / npt
if (ieq eq 1) then begin
   npt2=npt/2
   night_x=fltarr(npt2+2)     &    day_x=night_x
   night_y=fltarr(npt2+2)     &    day_y=night_y
   for i = 0,npt do begin
       i1 = i - npt2
       if (dir eq -1) then ang = float(i) * del + !pi/2. 
       if (dir eq 1) then ang = float(i) * del - !pi/2.
       cosa = cos(ang)
       sin_a = sin(ang)
       if (i le npt2) then day_x(i) = cosa
       if (i le npt2) then day_y(i) = sin_a
       if (i ge npt2) then night_x(i1) = cosa
       if (i ge npt2) then night_y(i1) = sin_a
   endfor
   day_x(npt2+1) = day_x(0)
   day_y(npt2+1) = day_y(0)
   night_x(npt2+1) = night_x(0)
   night_y(npt2+1) = night_y(0)
   polyfill,ro(0,0)*night_x,ro(0,0)*night_y,color=color2
   polyfill,ro(0,0)*day_x,ro(0,0)*day_y,color=color2
   polyfill,night_x,night_y,color=blackc
   polyfill,day_x,day_y,color=255
endif

; plot geosynchronous orbit, L* or B contours if ieq=1
if (ieq eq 1) then begin
   rad1=6.6  
   if (iovp eq 0) then rad1=halfl
   e1x=fltarr(npt+1)    &   e1y=e1x
   for i=0,npt do begin
       ang=i*del
       cosa=cos(ang)
       sin_a=sin(ang)
       e1x(i)=rad1*cosa
       e1y(i)=rad1*sin_a
   endfor
   if (iovp eq 0) then begin
      e2x=fltarr(npt+7)    &   e2y=e2x
      rad2=1.05*halfl
      e2x(0:npt)=e1x(0:npt)
      e2y(0:npt)=e1y(0:npt)
      e2x(npt+1:npt+6)=[rad2,rad2,-rad2,-rad2,rad2,rad2]
      e2y(npt+1:npt+6)=[0.,-rad2,-rad2,rad2,rad2,0.]   
   endif
   if (iovp eq 0) then polyfill,e2x,e2y,color=255
   if (iovp eq 1) then oplot,e1x,e1y,color=255,thick=gthick ;plot geosynchronous
   if (iovp eq 3) then begin
      blevel=[5e-9,1e-8,2e-8,5e-8,1e-7,2e-7,5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5]
    ; bcolor=[255,255,255,255,255,255,255,255,255,255,255,255,255]
      bcolor=blevel
      for i=0,12 do bcolor(i)=(i+1)*19
      contour,bo,xo,yo,levels=blevel,c_colors=bcolor,charsize=1.3, $
              thick=gthick,/overplot,/fill
   endif
   if (iovp eq 4 and xi ge 0.5) then begin
      Llevel=fltarr(nlevel)
      for i=0,nlevel-1 do begin
          Llevel(i)=1.+i*1.
          if (Llevel(i) eq ceil(Lstar_max)) then Llevel(i)=Lstar_max
          if (Llevel(i) ge Lstar_max) then colr(i)=255
      endfor
      contour,Lstar,xo,yo,levels=Llevel,c_colors=colr,/overplot,/fill
      colr(*)=0
      contour,Lstar,xo,yo,levels=Llevel(0:fix(Lstar_max)), $
              c_colors=colr(0:fix(Lstar_max)),/overplot
   endif
endif

; plot grids if iovp=2     
if (iovp eq 2) then begin                               ; plot grids
   for i=0,dim(1)-1 do oplot,xo(i,*),yo(i,*),color=255
   for j=0,dim(2)-1 do oplot,xo(*,j),yo(*,j),color=255
endif

return    ; end of plot_eqr
end


;-----------------------------------------------------------------------------
; main routine
;-----------------------------------------------------------------------------

close,/all

; get colors for color table palette
read,'colar table, 0=TWINS-webpage, 1=rainbow => ',icr
loadct,39       ;  rainbow with white
tvlct, red, green, blue, /get
blackc=0
gcolor=141
if (icr eq 1) then gcolor=93  
rxy=1.       ; the ratio of window size in x and y
if (icr ne 0) then rxy=1.5

; Setup plot size and format
iout=2
if (icr ne 0) then read,' plot format? (1) gif,  (2) ps => ',iout
y_wsize=360
x_wsize=y_wsize*rxy
if (iout eq 1) then begin
   yi=0.25
   plt_size=y_wsize*0.7/rxy  
endif else begin
   yi=0.15
   plt_size=y_wsize*(1.-2.*yi)/rxy  
endelse
yf=yi+plt_size/y_wsize
gthick=1.
if (iout eq 2) then gthick=4.
!p.color=blackc
!p.charthick=gthick
if (iout eq 2) then set_plot,'ps'

; read file name and energy, sina, lat information 
fhead=' '
read,'enter the file head (e.g., 2000a225_h) => ',fhead
openr,2,fhead+'.fls'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)   &   lat=fltarr(ir)   &  ece=energy
Ebound=fltarr(je+1)  &   mlt=fltarr(nmlt+1)   &  VarL=lat    &  mphi=fltarr(ip)
readf,2,varL   
readf,2,mphi   
readf,2,energy
readf,2,sina
read,'number of data set in time => ',ntime
read,'Sun to the?   -1 = left,  1 = right => ',dir
read,'half length of the plot => ', halfl
isav=2  &  lmax=1  &   jiso=1
if (icr ne 0) then begin
   read,'Save plots in (1)1 file, (2)individual files => ',isav
   read,'label flux maximum? 0=no, 1=yes => ',lmax
   read,'aniso = (1):(fper-fpar)/ftot, (2):(fper-fpar)/fpar => ',jiso
endif

; Calculate pitA
ig2=ig*2   &   pitA=fltarr(ig2)   &   pad=pitA
for m=0,ig2-1 do begin
    if (m lt ig) then pitA(m)=asin(sina(m))*180./!pi
    if (m ge ig) then pitA(m)=180.-pitA(ig2-1-m)
endfor

; Calculate Ebound    
for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
Ebound(0)=energy(0)*energy(0)/Ebound(1)
Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
 
; setup energy and species
nsp=strlen(fhead)-2
spe=strmid(fhead,nsp,2)                ; species
if (spe eq '_h') then species='H+'
if (spe eq 'he') then species='He+'
if (spe eq '_s') then species='SW H+'
if (spe eq 'Po') then species='PolarWind H+'
if (spe eq 'Pl') then species='PlasmasphericWind H+'
if (spe eq '_o') then species='O+'
IF (Spe eq '_e' or spe eq 'e1') then species='e-'

; setup arrays
lata=fltarr(ntime,nr)
ergflux=fltarr(ntime,nr,nmlt,je,ig)
plsfluxk=fltarr(ntime,nr,nmlt+1,je)    &  anisok=plsfluxk
plsfluxa=fltarr(ntime,nr,nmlt+1)       &  anisoa=plsfluxa
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa    &   boa=roa   &  plsDena=roa
denWPa=roa    &   Tpara=roa   &  Tpera=roa    &   Vpera=roa
Lstara=roa    &   CHpoa=roa   &  HIpoa=roa
iba=intarr(ntime,nmlt)     
ro=fltarr(nr,nmlt+1) & bo=ro & mlto=ro &  plsflux=ro & aniso=ro & Lstar=ro
Tpar=ro   &   Tper=ro   &   Vper=ro
houra=fltarr(ntime)     &   Bza=houra  &  Lstar_maxa=houra
fy=fltarr(ig)           &   ypar=fy    &   dmu=fy   &   cosa=fy  
parmod5=fltarr(6)
cosa(*)=cos(asin(sina(*)))
for m=0,ig-1 do begin
    if (m eq 0) then sina0=0.
    if (m gt 0) then sina0=0.5*(sina(m)+sina(m-1))
    if (m eq (ig-1)) then sina1=1.
    if (m lt (ig-1)) then sina1=0.5*(sina(m)+sina(m+1))
    dmu(m)=sqrt(1.-sina0*sina0)-sqrt(1.-sina1*sina1)
endfor
dE=fltarr(je)
for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; factor of calculating ion pressure, density and mean energy 
pf=8.*!pi/3.*1.e6*1.6e-16*1.e9    ; pressure in nPa
nf=4.*!pi                         ; density in cm-3

; Calculate ece, (E+Eo)/c/sqrt(E(E+2Eo)). c in cm/s
if (spe eq '_h' or spe eq '_s' or spe eq 'Po' or spe eq 'Pl') then mass=1.
if (spe eq 'he') then mass=4.  
if (spe eq '_o') then mass=16.
if (spe eq '_e' or spe eq 'e1') then mass=5.4462e-4
Eo=mass*1.673e-27*3.e8*3.e8/1.6e-16        ; rest-mass energy in keV
for k=0,je-1 do ece(k)=(energy(k)+Eo)/3.e10/sqrt(energy(k)*(energy(k)+2.*Eo))

; Read fluxes and calculate PA anisotropy
denWP1=0.
Tpar1=0.
Tper1=0.
HRPee=0.
HRPii=0.
Lstar_max=halfL
Lstar1=halfL
for n=0,ntime-1 do begin
    readf,2,hour,parmod1,parmod2,parmod3,Bz,parmod5,Lstar_max
    houra(n)=hour
    Bza(n)=Bz   
    Lstar_maxa(n)=Lstar_max   
    print,n,hour,Lstar_max
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
            readf,2,lat1,mlt1,xlatS,xmltS,ro1,mlto1,BriN,BriS,bo1,iba1,plsDen1,$
                    ompe,CHpower, $
                    HIpower,denWP1,Tpar1,Tper1,HRPee,HRPii,rpp,Lstar1
            lata(n,i)=lat1
            iba(n,j)=iba1
            roa(n,i,j)=ro1
            plsDena(n,i,j)=plsDen1
            denWPa(n,i,j)=denWP1
            mlt(j)=mlt1
            mltoa(n,i,j)=mlto1
            boa(n,i,j)=bo1
            Lstara(n,i,j)=Lstar1
            CHpoa(n,i,j)=CHpower
            HIpoa(n,i,j)=HIpower
            if (Lstar1 lt 0.) then Lstara(n,i,j)=999.   ; arbitrary large number
            Tpara(n,i,j)=Tpar1
            Tpera(n,i,j)=Tper1
            Vpera(n,i,j)=sqrt(2.*Tper1*1.6e-16/mass/1.67e-27)
            for k=0,je-1 do begin
                readf,2,fy
                if (i ge iba(n,j)) then fy(*)=0.
                ergflux(n,i,j,k,*)=fy(*)
                ; calculate pitch-angle (PA) averaged flux and PA anisotropy
                plsfluxk(n,i,j,k)=0.
             ;  anisok(n,i,j,k)=0.
                anisok(n,i,j,k)=-2.    ; beyond the range 
                fpara=0.
                fperp=0. 
                for m=0,ig-1 do begin
                    fydmu=fy(m)*dmu(m)
                    plsfluxk(n,i,j,k)=plsfluxk(n,i,j,k)+fydmu
                    fpara=fpara+fydmu*cosa(m)*cosa(m)
                    fperp=fperp+fydmu*sina(m)*sina(m)/2.
                endfor
                ft=fperp+fpara
                fpp=fperp-fpara
                if (ft gt 0. and jiso eq 1) then anisok(n,i,j,k)=fpp/ft
                if (fpara gt 0. and jiso eq 2) then anisok(n,i,j,k)=fpp/fpara
            endfor
        endfor
    endfor
    mltoa(n,*,*)=mltoa(n,*,*)*!pi/12.                 ; mlto in radian
    mlt(*)=mlt(*)*!pi/12.                             ; mlt in radian
    ; periodic boundary condition
    roa(n,*,nmlt)=roa(n,*,0)
    mltoa(n,*,nmlt)=mltoa(n,*,0)
    mlt(nmlt)=mlt(0)+2.*!pi
    boa(n,*,nmlt)=boa(n,*,0)
    Lstara(n,*,nmlt)=Lstara(n,*,0)
    CHpoa(n,*,nmlt)=CHpoa(n,*,0)
    HIpoa(n,*,nmlt)=HIpoa(n,*,0)
    plsDena(n,*,nmlt)=plsDena(n,*,0)
    denWPa(n,*,nmlt)=denWPa(n,*,0)
    Tpara(n,*,nmlt)=Tpara(n,*,0)
    Tpera(n,*,nmlt)=Tpera(n,*,0)
    Vpera(n,*,nmlt)=Vpera(n,*,0)
    plsfluxk(n,*,nmlt,*)=plsfluxk(n,*,0,*)
    anisok(n,*,nmlt,*)=anisok(n,*,0,*)
endfor
close,2

new_plot:

; choose plotting variable and which energy to be displaced
print,'plot? (1)flux,  (2)pressure,  (3)density,  (4)mean energy '
read, '      (5)Nhot/Ncold+Vr/Vper,  (6)Tpar+Tper (7)wave power => ',iopt
iplt=1
if (iopt le 4) then $
read,'plot? (0)L* (1)anisotropy (2)ionosphere map (3)spectra at peak => ',iplt
iovp=4
if (iplt gt 0 and iplt lt 2 and icr ne 0) then $
  read,'over plot? 0=white-background, 1=geosyn, 2=grid, 3=B contours => ',iovp
if (iopt ge 5) then begin
  ie1=0         ; entire energy range when plot total density and temperature
  ie2=je-1      ;
endif else begin
  print,'energy bin (keV): '
  for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f9.4," - ",f9.4)'
  read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2
endelse
iunit=0
if (iopt eq 1 and iplt ne 3) then $
   read,'unit? (0)#/(eVcm2ssr), (1)#/(keVcm2ssr), (2)#/(cm2ssr) => ',iunit
if (Ebound(ie1) gt 1.) then begin
   e0=Ebound(ie1)
   e1=Ebound(ie2+1)
   elab=string(round(e0),'(i4.4)')+'-'+string(round(e1),'(i4.4)')+'keV'
endif else begin
   e0=1000.*Ebound(ie1)
   e1=1000.*Ebound(ie2+1)
   elab=string(round(e0),'(i6.6)')+'-'+string(round(e1),'(i6.6)')+'eV'
endelse
if (iopt eq 1) then fmiddle=elab+'_flux' 
if (iopt eq 2) then fmiddle=elab+'_pressure'
if (iopt eq 3) then fmiddle=elab+'_density'
if (iopt eq 4) then fmiddle=elab+'_mean_energy'
if (iopt eq 6) then fmiddle='Temp'
if (iopt eq 7) then fmiddle='WavePower'
if (iplt eq 3) then fmiddle=elab+'_spec' 

; calculate energy-integrated flux or pressure or density or mean energy or Temp
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        for j=0,nmlt do begin
            plsfluxa(n,i,j)=0.
            anisoa(n,i,j)=0.
            weiSum=0.
            denSum=0.
            if (iopt lt 5) then begin
              for k=ie1,ie2 do begin    
                if (iopt eq 1) then wei=dE(k)
                if (iopt eq 2 or iopt eq 4) then wei=dE(k)*energy(k)*ece(k)
                if (iopt eq 3) then wei=dE(k)*ece(k)           
                weiSum=weiSum+wei
                if (iopt eq 4) then denSum=denSum+plsfluxk(n,i,j,k)*dE(k)*ece(k)
                plsfluxa(n,i,j)=plsfluxa(n,i,j)+plsfluxk(n,i,j,k)*wei  
                anisoa(n,i,j)=anisoa(n,i,j)+anisok(n,i,j,k)*wei  
              endfor
            endif
            if (iopt eq 1) then begin     ; differential flux
               if (iunit eq 0) then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum/1000.
               if (iunit eq 1) then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum
            endif
            if (iopt eq 2) then plsfluxa(n,i,j)=pf*plsfluxa(n,i,j)
            if (iopt eq 3) then plsfluxa(n,i,j)=nf*plsfluxa(n,i,j)
            if (iopt eq 4 and denSum gt 0.) then plsfluxa(n,i,j)= $
                                               plsfluxa(n,i,j)/denSum ; E in keV
            if (iopt eq 5 and plsDena(n,i,j) gt 0.) then $
                plsfluxa(n,i,j)=denWPa(n,i,j)/plsDena(n,i,j)
            if (iopt eq 6) then plsfluxa(n,i,j)=Tpara(n,i,j)
            if (iopt eq 7) then plsfluxa(n,i,j)=CHpoa(n,i,j)
            if (iopt lt 5) then anisoa(n,i,j)=anisoa(n,i,j)/weiSum
            if (iopt eq 6) then anisoa(n,i,j)=Tpera(n,i,j)
            if (iopt eq 7) then anisoa(n,i,j)=HIpoa(n,i,j)
        endfor
    endfor
endfor

; setup plot ranges
ilog=0
if (iplt le 2) then begin
   read,'scale?  0 = linear scale, or 1 = log scale => ',ilog
   if (ilog eq 1) then begin
      if(min(plsfluxa) lt 1.e-30)then plsfluxa(where(plsfluxa lt 1.e-30))=1.e-30
      plsfluxa=alog10(plsfluxa)
      if (iopt ge 5) then begin
         if (min(anisoa) lt 1.e-30) then anisoa(where(anisoa lt 1.e-30))=1.e-30
         anisoa=alog10(anisoa)
      endif
   endif
   yon=' '
   fmax=max(plsfluxa)    
   print,' fmax = ',fmax
   read,' Do you want to change it? (y/n) => ',yon
   if (yon eq 'y') then read,' enter new fmax => ',fmax
   fmin=0.   
   if (ilog eq 1) then fmin=fmax-3.   
   print,' fmin = ',fmin
   read,' Do you want to change it? (y/n) => ',yon
   if (yon eq 'y') then read,' enter new fmin => ',fmin
   amin=-1.
   amax=1.
   if (jiso eq 2) then begin
      amin=-0.5
      amax=1.5
   endif
   if (iopt ge 5) then begin  
      amax=max(anisoa)   
      print,' amax = ',amax
      read,' Do you want to change it? (y/n) => ',yon
      if (yon eq 'y') then read,' enter new amax => ',amax
      amin=0.
      if (ilog eq 1) then amin=amax-3.
      print,' amin = ',amin
      read,' Do you want to change it? (y/n) => ',yon
      if (yon eq 'y') then read,' enter new amin => ',amin
   endif
endif     ; end of if (iplt le 2)

; smooth data
yon='n'
if (iplt le 2 and icr ne 0) then read,' smooth the data? (y/n) => ',yon
nframe=0
if (yon eq 'y') then read,' enter number of intermediate frame => ',nframe
ntnew=(ntime-1)*(nframe+1)+1
ronew=fltarr(ntnew,nr,nmlt+1)  &  mltonew=ronew  
anisonew=ronew                 &  plsfluxnew=ronew
if (nframe eq 0) then begin
   hournew=houra        
   ronew=roa        
   mltonew=mltoa   
   anisonew=anisoa     
   plsfluxnew=plsfluxa
endif else begin
   hournew=interpol(houra,ntnew)
   for i=0,nr-1 do begin 
       for j=0,nmlt do begin
           ronew(*,i,j)=interpol(roa(*,i,j),ntnew)
           mltonew(*,i,j)=interpol(mltoa(*,i,j),ntnew)
           anisonew(*,i,j)=interpol(anisoa(*,i,j),ntnew)
           plsfluxnew(*,i,j)=interpol(plsfluxa(*,i,j),ntnew)
       endfor
   endfor
endelse

if (iout eq 1) then begin
   DEVICE, DECOMPOSED=-1
   window,1,xpos=200,ypos=200,xsize=x_wsize,ysize=y_wsize
endif
if (isav eq 1 and iout eq 2) then device,filename=fhead+'_'+fmiddle+'.ps', $
   /inches,/color,yoffset=3.0,xoffset=0.3,xsize=6.*x_wsize/y_wsize,ysize=6., $
   bits_per_pixel=24

; make ps file for converting to mpg
; device,filename=fhead+'_'+fmiddle+'.ps',/color,/landscape
; ps2pdf *.ps
; gm convert -delay 25 -quality 90 -rotate 180 *pdf Movie.mpg

; Plot plasma flux (or pressure, or density) and pitch angle anisotropy
maxL=max(Lstar_maxa)
for ipt=1,isav do begin
if (ipt gt 1) then begin
   fmiddle=elab+'_aniso'
   if (iplt eq 3) then fmiddle=elab+'_PAD' else begin
      loadct,18
      tvlct, red, green, blue, /get
      red(255)=255
      green(255)=255
      blue(255)=255
      red(0)=0
      green(0)=0
      blue(0)=0
      tvlct, red, green, blue
   endelse
endif
if (ipt eq 1 and icr eq 0) then print,'     hour       ro       mlto      Max '
for n=0,ntnew-1 do begin
    if (isav eq 2) then begin
       tlab=string(fix(n),'(i3.3)')
       device,filename=fhead+'_'+fmiddle+tlab+'.ps', $
              /inches,/color,xsize=6.*x_wsize/y_wsize, $
              ysize=6.,bits_per_pixel=24,/encapsulated
    endif
    hour=hournew(n)
    lat(*)=lata(n,*)
    ro(*,*)=ronew(n,*,*)
    mlto(*,*)=mltonew(n,*,*)
    bo(*,*)=boa(n,*,*)
    Lstar(*,*)=Lstara(n,*,*)
    Lstar_max=Lstar_maxa(n)
    aniso(*,*)=anisonew(n,*,*)
    plsflux(*,*)=plsfluxnew(n,*,*)
    ieq=1
    piso=0

    ; Local max value
    Atemp=where(ro lt halfl)
    maxflux=max(plsflux(Atemp))    ; max flux inside halfl
    ijmax=where(plsflux(0:nr-1,0:nmlt-1) eq maxflux)
    ijsize=size(ijmax)
    ii=-1   &   jj=-1      ; initial values
    if (lmax eq 1 and ijsize(1) eq 1) then begin
       jj=fix(ijmax/nr)
       ii=fix(ijmax-jj*nr)
    endif

    ; plot plasma distribution at equator or energy spectra at peak
    if (ipt eq 1) then begin
    x0=0.15/rxy
    xi=x0
    xf=xi+plt_size/x_wsize
    if (iopt eq 1) then begin
       if (iunit eq 0) then bar_lab='flux (/eV/cm2/sr/s)'
       if (iunit eq 1) then bar_lab='flux (/keV/cm2/sr/s)'
       if (iunit eq 2) then bar_lab='flux (/cm2/sr/s)'
    endif
    if (iopt eq 2) then bar_lab='pressure (nPa)'
    if (iopt eq 3) then bar_lab='density (cm^-3)'
    if (iopt eq 4) then bar_lab='mean energy (keV)'
    if (iopt eq 5) then bar_lab='hot/cold density ratio'
    if (iopt eq 6) then bar_lab='T_parallel (keV)'
    if (iopt eq 7) then bar_lab='Chorus Wave Power (pT)^2'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    if (iplt le 2) then begin
       plot_eqr,plsflux,fmax,fmin,xi,yi,xf,yf,bar_lab,halfl,ro, $
               mlto,bo,Lstar,Lstar_max,maxL,lat,mlt,x_wsize,y_wsize, $
               dir,iovp,iopt,icr,blackc,ieq,piso,ii,jj,gthick,rxy
       if (icr eq 0) then print,format='(4f10.4)', $
                       houra(n),ro(ii,jj),mlto(ii,jj)*12./!pi,plsflux(ii,jj)
    endif else begin
       plot,energy,plsfluxk(n,ii,jj,*)/1000,pos=[xi,yi,xf,yf], $
            xrange=[0,100],xstyle=1,xcharsize=1.5,ycharsize=1.5, $
            xtitle='[keV]',ytitle='[(s cm^2 eV sr)^-1]',/normal,/noerase
       oplot,energy,plsfluxk(n,ii,jj,*)/1000,color=gcolor,thick=6
       oplot,energy,plsfluxk(n,ii,jj,*)/1000,psym=1 
    endelse
    endif       ; end of if (ipt eq 1)

    ; plot anisotropy or ionosphere map or PAD at peak
    if (ipt eq isav) then begin
    if (isav eq 1) then xi=0.53  
    xf=xi+plt_size/x_wsize
    if (iplt eq 2) then ieq=-1
    if (ieq eq 1) then begin
       if (iplt le 1) then begin
          bar_lab='pitch angle anisotropy'
          if (iplt eq 0) then bar_lab='L*'
          if (iopt eq 2) then bar_lab='pressure anisotropy'
          if (iopt eq 4) then bar_lab='energy anisotropy'
          if (iopt eq 6) then bar_lab='T_perpendicular (keV)'
          if (iopt eq 7) then bar_lab='Hiss Wave Power (pT)^2'
          if (iopt ge 5 and ilog eq 1) then bar_lab='log '+bar_lab
          plot_eqr,aniso,amax,amin,xi,yi,xf,yf,bar_lab,halfl,ro, $
               mlto,bo,Lstar,Lstar_max,maxL,lat,mlt,x_wsize,y_wsize, $
               dir,iovp,iopt,icr,blackc,ieq,piso,ii,jj,gthick,rxy
       endif else begin
          for m=0,ig2-1 do begin
              pad(m)=0.
              if (m lt ig) then $
                 for k=ie1,ie2 do pad(m)=pad(m)+ergflux(n,ii,jj,k,m)*dE(k)
              if (m ge ig) then pad(m)=pad(ig2-1-m)
          endfor
          aveFlux=0.
          for m=0,ig-1 do aveFlux=aveFlux+pad(m)*dmu(m)
          if (aveFlux le 0) then pad(*)=1
          if (aveFlux gt 0) then pad(*)=pad(*)/aveFlux
          plot,pitA,pad,pos=[xi,yi,xf,yf],xrange=[0,180],yrange=[0,max(pad)], $
               xstyle=1,xtitle='Pitch Angle [Degrees]',xcharsize=1.5, $
               ycharsize=1.5,ytitle='Normalized PAD',/normal,/noerase
          oplot,pitA,pad,color=gcolor,thick=6      
          oplot,pitA,pad,psym=1
       endelse
    endif else begin
       plot_eqr,plsflux,fmax,fmin,xi,yi,xf,yf,bar_lab,halfl,ro, $
             mlto,bo,Lstar,Lstar_max,maxL,lat,mlt,x_wsize,y_wsize, $
             dir,iovp,iopt,icr,blackc,ieq,piso,ii,jj,gthick,rxy
    endelse
    endif      ; end of if (ipt eq isav)

    ;label the plot
    hours=string(fix(hour),'(i4)')
    minute=fix((hour-float(fix(hour)))*60.)
    minutes=string(minute,'(i2.2)')
    second=hour*3600.-fix(hour)*3600.-minute*60.
    sec=string(round(second),'(i2.2)')
    if (iplt le 2) then mlab=elab+' '+species
    if (iplt eq 3) then mlab='L = '+string(ro(ii,jj),'(f3.1)')+ $
                             ', MLT = '+string(mlto(ii,jj)*12./!pi,'(f4.1)')
    xyouts,0.50,yf+0.1,fhead+hours+':'+minutes+':'+sec,size=1.7, $
           alignment=0.5,color=blackc,/normal
    xyouts,x0,yf+0.03,mlab,size=1.8,color=blackc,/normal
    Llab='L*max = '+string(Lstar_max,'(f4.2)')
    if (iovp eq 4) then xyouts,xi,yf+0.03,Llab,size=1.8,color=blackc,/normal

    ; Make gif or ps file
    if (iout eq 1) then begin
       gimg = tvrd(0,0) 
       write_gif,fhead+'_'+fmiddle+'.gif',gimg,red(0:255),green(0:255), $
                 blue(0:255),/multiple
    endif
    if (iout eq 2) then erase
    if (isav eq 2) then device,/close_file
endfor    ; end of for n=0,ntnew-1
endfor    ; end of for ipt=1,isav

if (iout eq 1) then write_gif,fhead+'_'+fmiddle+'.gif',/close
if (isav eq 1 and iout eq 2) then device,/close_file
yon='n'
if (isav eq 1) then read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot

end

