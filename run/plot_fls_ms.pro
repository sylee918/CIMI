;-------------------------------------------------------------------------------
;
;                                plot_fls_ms.pro
; 
; IDL procedure modified from plot_fls.pro. This code is able to read and
; sum of flux, pressure, density, pitch-angle distribution and heating
; rate from multi-species.
;
; Created on 19 November 2015 by Mei-Ching Fok, Code 673, NASA GSFC.
;-------------------------------------------------------------------------------


;-------------------------------------------------------------------------------
pro plot_fls,var,vmax,vmin,xi,bar_lab,halfl,ro,mlto,bo,lat,mlt,rc, $
                 x_wsize,y_wsize,dir,iplt,iopt,blackc,ieq,jiso,iNS,ii,jj 
;-------------------------------------------------------------------------------

ionolab='At Northern Ionosphere'
if (iNS eq 2) then ionolab='At Southern Ionosphere'

; Make levels
nlevel=59
nlvl2=fix(nlevel/2)
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
vmid=0.5*(vmax+vmin)
;if (xi gt 0.5 and ieq eq 1) then vmid=0.
dlvl1=(vmid-vmin)/nlvl2
dlvl2=(vmax-vmid)/nlvl2
colrmax=253          
colrmin=1
ncolor=colrmax-colrmin+1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
lvl(0)=vmin
colr(0)=colrmin
for i=1,nlevel-1 do begin
    dlvl=dlvl1
    if (i gt nlvl2) then dlvl=dlvl2
    lvl(i)=lvl(i-1)+dlvl
    colr(i)=round(float(colrmin)+i*dcolr)
endfor
yi=0.40
plt_size=y_wsize*0.42
xf=xi+plt_size/x_wsize
yf=yi+plt_size/y_wsize

; force var equal or gt vmin 
if (min(var) lt vmin) then var(where(var le vmin))=vmin

; polyfill background color in black
polyfill,[xi,xf,xf,xi,xi],[yi,yi,yf,yf,yi],color=blackc,/normal

; Calculate altitude
Alt=(rc-1.)*6.375e3       ; altitdue in km
Alts=string(round(Alt),'(i3)')

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
   xymax=60.   
   if (dir eq -1) then phi=mlt               ; Sun to the left
   if (dir eq 1) then phi=mlt+!pi            ; Sun to the right
   for i=0,dim(1)-1 do begin
       for j=0,dim(2)-1 do begin
           xo(i,j)=(90.-lat(i,j))*cos(phi(i,j))
           yo(i,j)=(90.-lat(i,j))*sin(phi(i,j))
       endfor
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

; plot var
color1=blackc
if (iplt eq 0) then color1=255
contour,var,xo,yo,xrange=[-xymax,xymax],yrange=[-xymax,xymax], $
      levels=lvl,c_colors=colr,xstyle=1,ystyle=1,xticks=ntick,yticks=ntick, $
      xtickn=tickn,ytickn=tickn, $
      ytitle=y_title,pos=[xi,yi,xf,yf],charsize=1.3,color=color1,/fill,/noerase
if (ii ge 0) then oplot,[xo(ii,jj),xo(ii,jj)],[yo(ii,jj),yo(ii,jj)],psym=2

; Add color bar and label
xb=xi
yb=0.240
dx=0.0048
y1=yb+0.04
for i=0,nlevel-1 do begin
    x1=xb+i*dx
    x2=x1+dx*1.02
    polyfill,[x1,x2,x2,x1,x1],[yb,yb,y1,y1,yb],color=colr(i),/normal
endfor
v_min=string(vmin,'(f4.1)')
if (vmax lt 100.) then v_max=string(vmax,'(f4.1)')
if (vmax ge 100.) then v_max=string(vmax,'(e7.1)')
xyouts,0.5*(x2+xb),y1+0.01,bar_lab,alignment=0.5,size=1.4,color=blackc,/normal
xyouts,xb,yb-0.03,v_min,alignment=0.5,size=1.2,color=blackc,/normal
xyouts,x2,yb-0.03,v_max,alignment=0.5,size=1.2,color=blackc,/normal
; Add labels of field aligned and perpendicular when plotting anisotropy
if (xi gt 0.5 and ieq eq 1) then begin
   if (jiso eq 2) then begin
      x4=(x2-xb)/4.
      vmax2=string(vmax/2.,'(f4.1)')
   ;  xyouts,xb+x4,yb-0.03,'0.5',alignment=0.5,size=1.2,color=blackc,/normal
   ;  xyouts,xb+2.*x4,yb-0.03,'0',alignment=0.5,size=1.2,color=blackc,/normal
   ;  xyouts,xb+3.*x4,yb-0.03,vmax2,alignment=0.5,size=1.2,color=blackc,/normal
   endif
   xyouts,xb-0.05,yb-0.07,'field aligned',size=1.4,color=blackc,/normal
   xyouts,x2+0.05,yb-0.07,'perpendicular',alignment=1,size=1.4,color=blackc, $
          /normal
endif
; Add labels of equator and ionosphere when plotting heating rate
if (iopt ge 4 and ieq eq 1) then xyouts,0.5*(xb+x2),yb-0.07, $
   'At Equator',alignment=0.5,size=1.4,color=blackc,/normal
if (iopt ge 4 and ieq eq -1) then begin
   xyouts,0.5*(xb+x2),yb-0.07,ionolab,alignment=0.5,size=1.4,/normal 
   xyouts,0.5*(xb+x2),yb-0.11,'(Alt='+Alts+'km)',alignment=0.5,size=1.4,/normal
endif 

; fill cirle of radius ro(0,0) or 2Re centered the earth with color2
color2=colrmin
if (bar_lab eq 'pitch angle anisotropy') then begin
   m=(nlevel-1)/2
   color2=colr(m)
endif
ro00=ro(0,0)

; Draw earth 
npt=200                             ; no. of points in a circle
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
       sina = sin(ang)
       if (i le npt2) then day_x(i) = cosa
       if (i le npt2) then day_y(i) = sina
       if (i ge npt2) then night_x(i1) = cosa
       if (i ge npt2) then night_y(i1) = sina
   endfor
   day_x(npt2+1) = day_x(0)
   day_y(npt2+1) = day_y(0)
   night_x(npt2+1) = night_x(0)
   night_y(npt2+1) = night_y(0)
   polyfill,ro00*night_x,ro00*night_y,color=color2
   polyfill,ro00*day_x,ro00*day_y,color=color2
   polyfill,night_x,night_y,color=blackc
   polyfill,day_x,day_y,color=255
   oplot,day_x,day_y,color=blackc
endif

; Draw geosynchronous
if (ieq eq 1) then begin
   rad1=6.6  
   if (iplt eq 0) then rad1=halfl
   e1x=fltarr(npt+1)    &   e1y=e1x
   for i=0,npt do begin
       ang=i*del
       cosa=cos(ang)
       sina=sin(ang)
       e1x(i)=rad1*cosa
       e1y(i)=rad1*sina
   endfor
   if (iplt eq 0) then begin
      e2x=fltarr(npt+7)    &   e2y=e2x
      rad2=1.1*halfl
      e2x(0:npt)=e1x(0:npt)
      e2y(0:npt)=e1y(0:npt)
      e2x(npt+1:npt+6)=[rad2,rad2,-rad2,-rad2,rad2,rad2]
      e2y(npt+1:npt+6)=[0.,-rad2,-rad2,rad2,rad2,0.]   
   endif
endif

; plot geosynchronous orbit, or grids or B contours
if (ieq eq 1) then begin
   if (iplt eq 0) then polyfill,e2x,e2y,color=255
   if (iplt eq 1) then oplot,e1x,e1y,color=255,thick=3     ; plot geosynchronous
   if (iplt eq 2) then begin                               ; plot grids
      for i=0,dim(1)-1 do oplot,xo(i,*),yo(i,*),color=255
      for j=0,dim(2)-1 do oplot,xo(*,j),yo(*,j),color=255
   endif
   if (iplt eq 3) then begin
      blevel=[5e-9,1e-8,2e-8,5e-8,1e-7,2e-7,5e-7,1e-6,2e-6,5e-6,1e-5,2e-5,5e-5]
      bcolor=[255,255,255,255,255,255,255,255,255,255,255,255,255]
      contour,bo,xo,yo,levels=blevel,c_colors=bcolor,charsize=1.3,thick=3, $
              /overplot
   endif
endif

return
end


;-----------------------------------------------------------------------------
; main routine
;-----------------------------------------------------------------------------

close,/all

; get colors for color table palette
  loadct,39       ;  rainbow with white
  tvlct, red, green, blue, /get
  blackc=0

; Setup plot size and format
x_wsize=722
y_wsize=480
!p.color=blackc
!p.charthick=3.  
set_plot,'ps'

; read file name and energy, sina information 
fhead=' '
read,'enter the file head (e.g., 2000a225) => ',fhead
read,'number of species?  1=H+   2=H+,O+   3=H+,O+,He+ => ',isp
for is=2,isp+1 do begin
    if (is eq 2 ) then openr,is,fhead+'_h.fls'
    if (is eq 3 ) then openr,is,fhead+'_o.fls'
    if (is eq 4 ) then openr,is,fhead+'he.fls'
endfor
for is=2,isp+1 do readf,is,rc,ir,ip,je,ig 
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
energy=fltarr(je)    &   sina=fltarr(ig)      
varL=fltarr(ir)      &   mphi=fltarr(ip)
Ebound=fltarr(je+1)     
for is=2,isp+1 do readf,is,varL     
for is=2,isp+1 do readf,is,mphi     
for is=2,isp+1 do readf,is,energy
for is=2,isp+1 do readf,is,sina
read,'Sun to the?   -1 = left,  1 = right => ',dir
read,'half length of the plot => ', halfl
;read,'label flux maximum? 0=no, 1=yes => ',lmax
lmax=0
;read,'aniso = (1):(fper-fpar)/ftot, (2):(fper-fpar)/fpar => ',jiso
jiso=1
read,'Ionosphere map at (1)Northern (2)Southern ionosphere => ',iNS
read,'number of data set in time => ',ntime

; Calculate Ebound    
for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
Ebound(0)=energy(0)*energy(0)/Ebound(1)
Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
 
; setup species
species='H+'
splab='_H'
if (isp eq 2) then species='H+,O+'
if (isp eq 2) then splab='_H-O'
if (isp eq 3) then species='H+,O+,He+'
if (isp eq 3) then splab='_H-O-He'

; setup arrays
plsfluxk=fltarr(ntime,nr,nmlt+1,je,isp)    &  anisok=plsfluxk
ece=fltarr(je,isp)
plsfluxa=fltarr(ntime,nr,nmlt+1)           &  anisoa=plsfluxa
roa=fltarr(ntime,nr,nmlt+1)   &  mltoa=roa     &   boa=roa   &  plsDena=roa
denWPa=roa    &   Tpara=roa   &   Tpera=roa    &   mlta=roa
HRPea=roa     &   HRPia=roa   &   VHRPea=roa   &   VHRPia=roa
lata=roa      &   Bria=roa
iba=intarr(ntime,nmlt)     
ro=fltarr(nr,nmlt+1) & bo=ro & mlto=ro &  plsflux=ro & aniso=ro 
lat=ro               &   Bri=ro 
Tpar=ro   &   Tper=ro   &   Vper=ro    & ionoflux=ro  &  mlt=ro
houra=fltarr(ntime)     
fy=fltarr(ig)           &   ypar=fy      &    dmu=fy    &    cosa=fy
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

; factor of calculating ion pressure, density, mean energy and ionosphere Br
pf=8.*!pi/3.*1.e6*1.6e-16*1.e9    ; pressure in nPa
nf=4.*!pi                         ; density in cm-3
rcm=rc*6.375e6   
xme=8.0145e15                     ; Earth's magnetic moment in T m^3

; Calculate ece, (E+Eo)/c/sqrt(E(E+2Eo)). c in cm/s
for is=0,isp-1 do begin
    if (is eq 0) then mass=1.
    if (is eq 1) then mass=16.
    if (is eq 2) then mass=4.  
    Eo=mass*1.673e-27*3.e8*3.e8/1.6e-16        ; rest-mass energy in keV
    for k=0,je-1 do ece(k,is)= $
                         (energy(k)+Eo)/3.e10/sqrt(energy(k)*(energy(k)+2.*Eo))
endfor

; Read fluxes and calculate PA anisotropy
NorS='_N'
if (iNS eq 2) then NorS='_S'
mmax=ig-1
muSum=total(dmu(0:mmax))
for n=0,ntime-1 do begin
    for is=2,isp+1 do readf,is,hour,parmod1,parmod2,parmod3,Bz
    houra(n)=hour
    print,n,hour
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
            HRPea(n,i,j)=0.
            HRPia(n,i,j)=0.
            VHRPea(n,i,j)=0.
            VHRPia(n,i,j)=0.
            for is=2,isp+1 do begin
                readf,is,lat1,mlt1,lat2,mlt2,ro1,mlto1,bri1,bri2,bo1,iba1, $
                      plsDen1,ompe,CHpower,HIpower,denWP1,Tpar1,Tper1,HRPee, $
                      HRPii,rpp,Lstar,tvolume  ; FluxTube volume/magnetic flux
                HRPea(n,i,j)=HRPea(n,i,j)+HRPee
                HRPia(n,i,j)=HRPia(n,i,j)+HRPii
                VHRPea(n,i,j)=VHRPea(n,i,j)+HRPee*tvolume
                VHRPia(n,i,j)=VHRPia(n,i,j)+HRPii*tvolume
            endfor
            if (iNS eq 1) then begin
               lata(n,i,j)=lat1
               mlta(n,i,j)=mlt1
               Bria(n,i,j)=bri1
            endif else begin
               lata(n,i,j)=abs(lat2)
               mlta(n,i,j)=mlt2
               Bria(n,i,j)=bri2
            endelse
            iba(n,j)=iba1
            roa(n,i,j)=ro1
            plsDena(n,i,j)=plsDen1
            denWPa(n,i,j)=denWP1
            mltoa(n,i,j)=mlto1
            boa(n,i,j)=bo1
            Tpara(n,i,j)=Tpar1
            Tpera(n,i,j)=Tper1
            for k=0,je-1 do begin
                for is=0,isp-1 do begin
                  is2=is+2
                  readf,is2,fy
                  if (i ge iba(n,j)) then fy(*)=0.
                  ; calculate pitch-angle (PA) averaged flux and PA anisotropy
                  plsfluxk(n,i,j,k,is)=0.
                  anisok(n,i,j,k,is)=0.
                  fpara=0.
                  fperp=0. 
                  for m=0,mmax do begin
                       fydmu=fy(m)*dmu(m)
                       plsfluxk(n,i,j,k,is)=plsfluxk(n,i,j,k,is)+fydmu
                       fpara=fpara+fydmu*cosa(m)*cosa(m)
                       fperp=fperp+fydmu*sina(m)*sina(m)/2.
                  endfor
                  plsfluxk(n,i,j,k,is)=plsfluxk(n,i,j,k,is)/muSum 
                  ft=fperp+fpara
                  fpp=fperp-fpara
                  if (ft gt 0. and jiso eq 1) then anisok(n,i,j,k,is)=fpp/ft
                  if(fpara gt 0. and jiso eq 2)then anisok(n,i,j,k,is)=fpp/fpara
                endfor   ; end of is loop
            endfor       ; end of k loop
        endfor
    endfor
    mltoa(n,*,*)=mltoa(n,*,*)*!pi/12.                 ; mlto in radian
    mlta(n,*,*)=mlta(n,*,*)*!pi/12.                   ; mlt in radian
    ; periodic boundary condition
    lata(n,*,nmlt)=lata(n,*,0)
    roa(n,*,nmlt)=roa(n,*,0)
    mltoa(n,*,nmlt)=mltoa(n,*,0)
    mlta(n,*,nmlt)=mlta(n,*,0)
    Bria(n,*,nmlt)=Bria(n,*,0)
    boa(n,*,nmlt)=boa(n,*,0)
    plsDena(n,*,nmlt)=plsDena(n,*,0)
    denWPa(n,*,nmlt)=denWPa(n,*,0)
    Tpara(n,*,nmlt)=Tpara(n,*,0)
    Tpera(n,*,nmlt)=Tpera(n,*,0)
    plsfluxk(n,*,nmlt,*,*)=plsfluxk(n,*,0,*,*)
    anisok(n,*,nmlt,*,*)=anisok(n,*,0,*,*)
    HRPea(n,*,nmlt)=HRPea(n,*,0)
    HRPia(n,*,nmlt)=HRPia(n,*,0)
    VHRPea(n,*,nmlt)=VHRPea(n,*,0)
    VHRPia(n,*,nmlt)=VHRPia(n,*,0)
endfor
for is=2,isp+1 do close,is

new_plot:

; choose plotting variable and which energy to be displaced
read,'plot? 0=white-background, 1=geosyn, 2=grid, 3=B contours => ',iplt
print,'plot? (1)flux, (2)pressure, (3)density'
read,'      (4)HR-PlSp-e, (5)HR-PlSp-i => ',iopt
imap=2
iout=0
if (iopt lt 4) then read,'also plot? (1)anisotropy, (2)ionosphere map => ',imap
if (iopt ge 4) then begin
  ie1=0         ; entire energy range when plot energy deposition
  ie2=je-1      ;
; read,'Do you want to output heating rate in a file? (0)no, (1)yes => ',iout
endif else begin
  print,'energy bin (keV): '
  for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f9.4," - ",f9.4)'
  read, 'lower and upper energy bins (ie1,ie2) => ',ie1,ie2
endelse
iunit=0
if (iopt eq 1) then $
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
if (iopt eq 4) then fmiddle='HRPe'+NorS
if (iopt eq 5) then fmiddle='HRPi'+NorS

; calculate energy-integrated flux or pressure or density or mean energy or Temp
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        for j=0,nmlt do begin
            plsfluxa(n,i,j)=0.
            anisoa(n,i,j)=0.
            weiSum=0.
            weiSum0=0.
            if (iopt le 3) then begin
               for k=ie1,ie2 do begin    
                   for is=0,isp-1 do begin
                       wei=dE(k)
                       if (iopt eq 2) then wei=dE(k)*energy(k)*ece(k,is)
                       if (iopt eq 3) then wei=dE(k)*ece(k,is)           
                       weiSum=weiSum+wei
                       if (is eq 0) then weiSum0=weiSum0+wei
                       plsfluxa(n,i,j)=plsfluxa(n,i,j)+plsfluxk(n,i,j,k,is)*wei 
                       anisoa(n,i,j)=anisoa(n,i,j)+anisok(n,i,j,k,is)*wei  
                   endfor
               endfor
            endif
            if (iopt eq 1) then begin     ; differential flux
               if(iunit eq 0) then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum0/1000.
               if(iunit eq 1) then plsfluxa(n,i,j)=plsfluxa(n,i,j)/weiSum0
            endif
            if (iopt eq 2) then plsfluxa(n,i,j)=pf*plsfluxa(n,i,j)
            if (iopt eq 3) then plsfluxa(n,i,j)=nf*plsfluxa(n,i,j)
            if (iopt eq 4) then plsfluxa(n,i,j)=HRPea(n,i,j)/1.e3 ; keV/m^3/s to
            if (iopt eq 5) then plsfluxa(n,i,j)=HRPia(n,i,j)/1.e3 ; eV/cm^3/s
            if (iopt le 3) then anisoa(n,i,j)=anisoa(n,i,j)/weiSum
        endfor
    endfor
endfor

; setup plot ranges
read,'scale?  0 = linear scale, or 1 = log scale => ',ilog
if (ilog eq 1) then begin
   if (min(plsfluxa) lt 1.e-30) then plsfluxa(where(plsfluxa lt 1.e-30))=1.e-30
   plsfluxa=alog10(plsfluxa)
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
if (jiso eq 2 and iopt le 3) then begin
   amax=max(anisoa)
   print,' amax = ',amax
   read,' Do you want to change it? (y/n) => ',yon
   if (yon eq 'y') then read,' enter new amax => ',amax
   amin=-0.5
endif
if (iopt ge 4) then begin
   amin=0.
   if (ilog eq 1) then amin=7.     
   amax=5.
   if (ilog eq 1) then amax=10.    
endif

; smooth data
;read,' Do you want to smooth the data? (y/n) => ',yon
yon='n'
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

  device,filename=fhead+splab+'_'+fmiddle+'.ps',/inches,/color,yoffset=3.0, $
         xoffset=0.3,xsize=6.*x_wsize/y_wsize,ysize=6.,bits_per_pixel=24

; if (iout eq 1 ) then begin
;    openw,1,fhead+splab+'_'+fmiddle+'.dat'
;    printf,1,nr,nmlt+1,'      ; nlat,nmlt'
;    printf,1,' magnetic latitude (deg) at ionosphere'
;    printf,1,lat
;    printf,1,' magnetic local time (radian) at ionosphere'
;    printf,1,mlt
; endif

; make ps file for converting to mpg
; device,filename=fhead+splab+'_'+fmiddle+'.ps',/color,/landscape
; ps2pdf *.ps
; gm convert -delay 25 -quality 90 -rotate 180 *pdf Movie.mpg

; Plot plasma flux (or pressure, or density) and pitch angle anisotropy
for n=0,ntnew-1 do begin
    hour=hournew(n)
    lat(*,*)=lata(n,*,*)
    ro(*,*)=ronew(n,*,*)
    mlto(*,*)=mltonew(n,*,*)
    mlt(*,*)=mlta(n,*,*)
    Bri(*,*)=Bria(n,*,*)
    bo(*,*)=boa(n,*,*)
    aniso(*,*)=anisonew(n,*,*)
    plsflux(*,*)=plsfluxnew(n,*,*)

    ; plot plasma fluxes or ion pressure, or density, or mean E, or Tpar or HR
    xi=0.11
    ieq=1
    if (iopt eq 1) then begin
       if (iunit eq 0) then bar_lab='flux (/eV/cm2/sr/s)'
       if (iunit eq 1) then bar_lab='flux (/keV/cm2/sr/s)'
       if (iunit eq 2) then bar_lab='flux (/cm2/sr/s)'
    endif
    if (iopt eq 2) then bar_lab='pressure (nPa)'
    if (iopt eq 3) then bar_lab='density (cm^-3)'
    if (iopt eq 4) then bar_lab='PlSp e- heating (eV/cm^3/s)'
    if (iopt eq 5) then bar_lab='PlSp ion heating (eV/cm^3/s)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    Atemp=where(ro lt halfl)
    maxflux=max(plsflux(Atemp))    ; max flux inside halfl
    ijmax=where(plsflux(0:nr-1,0:nmlt-1) eq maxflux)
    ijsize=size(ijmax)
    ii=-1   &   jj=-1      ; initial values
    if (lmax eq 1 and ijsize(1) eq 1) then begin
       jj=fix(ijmax/nr)
       ii=ijmax-jj*nr
       print,'hour,ii,jj ',hour,ii,jj       
    endif
    plot_fls,plsflux,fmax,fmin,xi,bar_lab,halfl,ro,mlto,bo,lat,mlt,rc, $
                 x_wsize,y_wsize,dir,iplt,iopt,blackc,ieq,jiso,iNS,ii,jj 
    ; Output volume heating rate if iout=1
;   if (iout eq 1) then begin
;      if (iopt eq 4) then plsflux(*,*)=HRPea(n,*,*)/1.e3
;      if (iopt eq 5) then plsflux(*,*)=HRPia(n,*,*)/1.e3
;      printf,1,hour,'     ; hour'
;      printf,1,' ro in earth radius'
;      printf,1,ro
;      printf,1,' mlto in radian' 
;      printf,1,mlto
;      printf,1,' volume heating rate (eV/cm^3/s) at (ro,mlto)'
;      printf,1,plsflux
;   endif

    ; plot anisotropy or ionosphere map or hot/cold density ratio
    xi=xi+0.42
    if (iopt ge 4) then begin     ; calculate heating rate at ionosphere
       if (iopt eq 4) then begin
          plsflux(*,*)=VHRPea(n,*,*)
          bar_lab='PlSp e- heating (eV/cm^2/s)'
          if (ilog eq 0) then bar_lab='PlSp e- heating (10^9 eV/cm^2/s)'
       endif else begin
          plsflux(*,*)=VHRPia(n,*,*)
          bar_lab='PlSp ion heating (eV/cm^2/s)'
          if (ilog eq 0) then bar_lab='PlSp ion heating (10^9 eV/cm^2/s)'
       endelse
       for i=0,nr-1 do begin     ; 0.5*ionoflux becuase it is for one ionosphere
         for j=0,nmlt do ionoflux(i,j)=0.5*plsflux(i,j)*Bri(i,j)*1.e-1 ;eV/cm2/s
       endfor
       ; Output heat flux at ionosphere if iout=1
;      if (iout eq 1) then begin
;         printf,1,' heat flux (eV/cm^2/s) at ionosphere'
;         printf,1,ionoflux
;      endif
       if (ilog eq 0) then ionoflux=ionoflux/1.e9
       if (ilog eq 1) then begin
          if (min(ionoflux) lt 1.e-30) then $
             ionoflux(where(ionoflux lt 1.e-30))=1.e-30
          ionoflux=alog10(ionoflux)
          bar_lab='log '+bar_lab
       endif
    endif
    if (imap eq 2) then ieq=-1
    if (ieq eq 1) then begin
       bar_lab='pitch angle anisotropy'
       if (iopt eq 2) then bar_lab='pressure anisotropy'
       plot_fls,aniso,amax,amin,xi,bar_lab,halfl,ro,mlto,bo,lat,mlt,rc, $
                 x_wsize,y_wsize,dir,iplt,iopt,blackc,ieq,jiso,iNS,ii,jj 
    endif else begin
       if (iopt le 3) then plot_fls,plsflux,fmax,fmin,xi,bar_lab,halfl,ro, $
        mlto,bo,lat,mlt,rc,x_wsize,y_wsize,dir,iplt,iopt,blackc,ieq, $
        jiso,iNS,ii,jj 
       if (iopt ge 4) then plot_fls,ionoflux,amax,amin,xi,bar_lab,halfl,ro, $
        mlto,bo,lat,mlt,rc,x_wsize,y_wsize,dir,iplt,iopt,blackc,ieq, $
        jiso,iNS,ii,jj 
    endelse

    ;label the plot
    hours=string(fix(hour),'(i3)')
    minute=fix((hour-float(fix(hour)))*60.)
    minutes=string(minute,'(i2.2)')
    second=hour*3600.-fix(hour)*3600.-minute*60.
    sec=string(round(second),'(i2.2)')
    xyouts,0.07,0.92,fhead,size=1.8,color=blackc,/normal
    xyouts,0.5,0.92,hours+':'+minutes+':'+sec,size=1.8,alignment=0.5, $
           color=blackc,/normal 
    mtitle=elab+' '+species
    if (iopt ge 4) then mtitle='Heating from '+mtitle
    xyouts,0.45,0.86,mtitle,size=1.8,alignment=0.5,color=blackc,/normal

    erase
endfor

  device,/close_file
; if (iout eq 1) then close,1
  read,'Do you want to continue?  (y)yes, (n)no => ',yon
  if (yon eq 'y') then goto,new_plot

end
