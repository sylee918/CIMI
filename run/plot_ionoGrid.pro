;-----------------------------------------------------------------------------
; Program to read and plot precipitating energy flux and mean energy from file
; *.preci.
;
;  Created on 19 November 2011 by Mei-Ching Fok (Code 673/GSFC)
;
;  Modification history
;  November 9, 2011
;    * Add plotting Dst
;  December 19, 2011
;    * Adopt output with energy distribution
;  December 5, 2012
;    * Add calculation of Pedersen and Hall conductance.
;  January 5, 2021
;    * Add the option to plot total precipitation from N and S ionospheres.
; March 19, 2021
;    * Add the option to plot in magnetic longitude
; April 29, 2021
;    * Add the option to plot losscone asymmetry between N and S
;-----------------------------------------------------------------------------

close,/all

; get colors for color table palette
read,'colar table, 1=rainbow, 2=papco, 3=Pastels => ',icr
if (icr ne 2) then begin
   if (icr eq 1) then loadct,39            ;  rainbow with white
   if (icr eq 3) then loadct,18            ;  pastels
   tvlct, red, green, blue, /get
   red(0)=0  &  green(0)=0  &  blue(0)=0   ; make color=0 black
   tvlct, red, green, blue
endif else begin
   openr,1,'papco.tbl'
   readf,1,red,green,blue
   close,1
   tvlct, red, green, blue
endelse

; Read the grids and grid size
filehead=' '
read,'enter file name of precipitating flux (no extension) => ',filehead
read,' plot? (1)N_hemisphere, (2)S_hemisphere, (3)N+S_hemisphere => ',nf
j1=1
j2=2
if (nf eq 1) then j2=1
if (nf eq 2) then j1=2
if (nf eq 1 or nf eq 3) then openr,1,filehead+'_N.preci'
if (nf eq 2 or nf eq 3) then openr,2,filehead+'_S.preci'
Hhead='_N'
if (nf eq 2) then Hhead='_S'
if (nf eq 3) then Hhead='_NS'
iLC=0
if (nf eq 3) then read,'Plot ratio of N-S losscone? (0)no, (1)yes => ',iLC
read,' plot format? (1) gif,  (2) ps => ',iout
thk=iout*1.5
read,' azimuthal? (1)Sun is left, (2)Sun is up, (3)zero mlong is right => ',idir
dhead='_'+string(idir,'(i1)')
slab=' X mlong-0'
if (idir eq 3) then slab=' X Sun'
for m=j1,j2 do readf,m,rc,ir,ip,je,ik
for m=j1,j2 do readf,m,elon,ctp,stp   
dlong=2.*!pi/ip
re=6.3712e6      ; Earth radius in meter (re_m) defined in cimi.f90
rp0=re*rc        ; rc in m
altkm=(rp0-re)/1000. 
althead=string(fix(altkm),'(i4)')+'km altitude'
read,'enter number of data set in time => ',ntime
gride=fltarr(je)    &   ebound=fltarr(je+1)
Elab=strarr(je+1)   &   Elabi=Elab
Efluxk=fltarr(je+1)  &  Efluxk0=Efluxk    &   Efluxka=Efluxk
Ppre=Efluxk          &  Ppre0=Ppre
for m=j1,j2 do readf,m,gride

!p.background=255
!p.color=0
!p.charthick=thk
!p.charsize=2.-iout*0.5
;clab=255
clab=0   
if (icr eq 3) then clab=234

; Setup energy boundaries (ebound)
ebound(0)=0.
for k=1,je-1 do ebound(k)=sqrt(gride(k-1)*gride(k))
ebound(je)=gride(je-1)*gride(je-1)/ebound(je-1)
for k=0,je do begin
    Elab(k)=string(ebound(k),'(f6.1)')
    Elabi(k)=string(round(ebound(k)),'(i4.4)')
endfor
 
; Choose energy range
Ehead=' '
print,'plot energy(keV) : '
for k=0,je-1 do print,k,'  '+Elab(k),' - ',Elab(k+1)
read,' enter lower and upper energy bins => ',ie1,ie2
Ehead=Elabi(ie1)+'-'+Elabi(ie2+1)+'keV'

; Read precipitating energy flux [mW/m2] and mean eneryg [keV]
halfL=70.
rlab=halfL-5.
latmin=90.-halfL
Efluxa=fltarr(ntime,ir,ip+1)   &   meanEa=Efluxa    &    ALC=Efluxa
Eflux1=fltarr(ip)              &   Eflux1a=Eflux1   &    meanE1=Eflux1
plt_v=fltarr(ir,ip+1)
houra=fltarr(ntime)            &   Dsta=houra       &    mlon0=houra
xoc=fltarr(ntime,ir,ip+1)      &   yoc=xoc       ; grid at ionosphere
xo=fltarr(ir,ip+1)             &   yo=xo
mc=intarr(ntime,ir,ip)     ; number of point in losscone
total_flx=fltarr(ntime)
area=fltarr(3)
for n=0,ntime-1 do begin
   for m=j1,j2 do readf,m,hour,Dst
   houra(n)=hour
   Dsta(n)=Dst 
   total_flx(n)=0.
   for i=0,ir-1 do begin
       for j=0,ip-1 do begin
           for m=j1,j2 do begin
               readf,m,xlat1,mlt1,mlon1,area1,mc1,Bi 
               if (m eq j1) then begin
                  mlt=mlt1      ; use mlt at N if both hemispheres
                  mlong=mlon1   ; use mlon at N if both hemispheres
                  if (i eq 0 and j eq 0) then mlon0(n)=mlt1*!pi/12.
                  xlat=abs(xlat1) ; use lat @ N if both hemispher
                  BiN=Bi
               endif
               area(m)=area1
           endfor
           Efluxk(*)=0.
           Efluxka(*)=0.
           Ppre(*)=0.
           for m=j1,j2 do begin
               readf,m,Efluxk0
               Efluxk(*)=Efluxk(*)+Efluxk0(*)
               Efluxka(*)=Efluxka(*)+Efluxk0(*)*area(m)
               readf,m,Ppre0
               Ppre(*)=Ppre(*)+Ppre0(*)
           endfor
           mc(n,i,j)=ik-mc1
           ; find the N-S asymmetry in loss cone
           ALC(n,i,j)=0. 
           if (iLC eq 1 and BiN gt 0.) then begin
              BiB=sqrt(Bi/BiN)
              ALC(n,i,j)=(BiB-1.)/(BiB+1.)  ; (sinLCn-sinLCs)/(sinLCn+sinLCs)
           endif
           if (xlat lt latmin) then begin
              Efluxk(*)=0.
              Ppre(*)=0.
           endif
           Eflux1(j)=0.  
           Eflux1a(j)=0.  
           Ppretot=0. 
           Etot=0.
           for k=ie1,ie2 do begin
               if (Efluxk(k) gt 0.) then Eflux1(j)=Eflux1(j)+Efluxk(k)
               if (Efluxk(k) gt 0.) then Eflux1a(j)=Eflux1a(j)+Efluxka(k)
               if (Ppre(k) gt 0.) then begin
                  Etot=Etot+Ppre(k)*gride(k)
                  Ppretot=Ppretot+Ppre(k)
               endif
           endfor
           meanE1(j)=0.
           if (Ppretot gt 0.) then meanE1(j)=Etot/Ppretot
           total_flx(n)=total_flx(n)+Eflux1a(j)
           ; calculate grid 
           clat=90.-xlat 
           phi1=mlt*!pi/12.                        ; zero -> midnight
           if (idir eq 2) then phi1=phi1-!pi/2.    ; zero -> dawn
           if (idir eq 3) then phi1=mlong*!pi/180. ; zero -> zero mag longitude
           xoc(n,i,j)=clat*cos(phi1)
           yoc(n,i,j)=clat*sin(phi1)
       endfor  ; end of j loop
       xoc(n,i,ip)=xoc(n,i,0)
       yoc(n,i,ip)=yoc(n,i,0)
       Efluxa(n,i,0:ip-1)=Eflux1(0:ip-1)
       meanEa(n,i,0:ip-1)=meanE1(0:ip-1) 
   endfor      ; end of i loop
   total_flx(n)=total_flx(n)/1.e12         ; from mW to GW
endfor         ; end of time loop
for m=j1,j2 do close,m
meanEa(*,*,ip)=meanEa(*,*,0)
Efluxa(*,*,ip)=Efluxa(*,*,0)
if (nf eq 3) then ALC(*,*,ip)=ALC(*,*,0)
hour_max=ceil(max(houra)/24.)*24
hour_min=fix(min(houra)/24.)*24
x_tick=hour_max/24
Dst_max=ceil(max(Dsta)/10.)*10
Dst_min=ceil(min(Dsta)/10.)*10-10

; start plotting from the n0-th time
n0=1 
meanEa(0:n0-1,*,*)=0.
Efluxa(0:n0-1,*,*)=0.

; Setup color levels
nlevel=59
colr=intarr(nlevel)
colrmax=250.          &    colrmin=10.
dcolr=(colrmax-colrmin)/(nlevel-1)
for i=0,nlevel-1 do colr(i)=fix(colrmin+i*dcolr)
if (icr eq 3) then for i=0,nlevel-1 do colr(i)=fix(colrmax-i*dcolr)

; Setup contour levels
lvl=fltarr(nlevel,2)  &    maxV=strarr(2)   &   minV=maxV
icon=0
;read,'Do you want to plot conductance instead? (0)no, (1)yes => ',icon
; Set level max
lvlmax0=10.   ; max Pedersen conductance
lvlmax1=30.   ; max Hall conductance
yon=' '
if (icon eq 0) then begin   ; plot energy flux and mean energy or ALC
   lvlmax0=max(Efluxa)
   print,'max energy flux is => ',lvlmax0  
   read,'Do you want to change it ? (y)yes, (n)no => ',yon
   if (yon eq 'y') then read,'enter new max energy flux => ',lvlmax0
   lvlmax1=0.2
   if (iLC eq 0) then begin
      lvlmax1=max(meanEa)
      print,'max mean energy is => ',lvlmax1  
      read,'Do you want to change it ? (y)yes, (n)no => ',yon
      if (yon eq 'y') then read,'enter new max mean energy => ',lvlmax1
   endif
endif
; Set level min
lvlmin0=0.1   
lvlmin1=0.1  
if (iLC eq 1) then lvlmin1=-lvlmax1
; Set contour levels
for iplt=0,1 do begin
    isum=iplt+icon+1       ; log scale if isum=0
    if (iplt eq 0) then lvlmax=lvlmax0 else lvlmax=lvlmax1
    if (iplt eq 0) then lvlmin=lvlmin0 else lvlmin=lvlmin1
    if (isum eq 0) then lvlmin=lvlmax/1000.
    dlvl=(lvlmax-lvlmin)/(nlevel-1)
    if (isum eq 0) then dlvl=(alog10(lvlmax)-alog10(lvlmin))/(nlevel-1)
    for i=0,nlevel-1 do begin
        if (isum gt 0) then lvl(i,iplt)=lvlmin+i*dlvl
        if (isum eq 0) then lvl(i,iplt)=lvlmin*10.^(i*dlvl)
    endfor
    lvl2=nlevel-1
    if (iplt eq 1 and iLC eq 1) then lvl2=44
    if (isum gt 0) then maxV(iplt)=string(lvl(lvl2,iplt),'(f4.1)')
    if (isum eq 0) then maxV(iplt)=string(alog10(lvl(lvl2,iplt)),'(f4.1)')
    if (isum gt 0) then minV(iplt)=string(lvl(0,iplt),'(f4.1)')
    if (isum eq 0) then minV(iplt)=string(alog10(lvl(0,iplt)),'(f4.1)')
endfor

; Set up for contour plot
xb=fltarr(ip+1)         &      yb=xb
for j=0,ip-1 do begin
    phi1=j*2.*!pi/ip 
    xb(j)=halfL*cos(phi1)
    yb(j)=halfL*sin(phi1)
endfor
xb(ip)=xb(0)
yb(ip)=yb(0)
 ntick=2*halfL/20
;ntick=3          
dlat=2.*halfL/ntick
tickn=strarr(ntick+1)  &   ticknx=tickn    &   tickny=tickn
for i=0,ntick do begin
    lat_1=90.-abs(-halfL+i*dlat)
    tickn(i)=string(fix(lat_1),'(i2)')
endfor
for i=0,ntick do begin
    lat_1=75.-i*5.
    ticknx(i)=string(fix(lat_1),'(i2)')
    lat_1=75.+i*5.
    tickny(i)=string(fix(lat_1),'(i2)')
endfor
ptitle=strarr(2)
if (icon eq 0) then begin
   ftitle='_pflx_'
   if (iLc eq 1) then ftitle='_pflxLC_'
   ptitle(0)='Energy Flux'
   ptitle(1)='Mean Energy'
   if (iLC eq 1) then ptitle(1)='Loss Cone Asymmetry'
endif else begin
   ftitle='_cond_'
   ptitle(0)='Pedersen conductance'
   ptitle(1)='Hall conductance'
endelse

; Setup for color bars
y0=0.45
dy=0.007
y2=y0+nlevel*dy
yave=0.5*(y0+y2)
bar_lab=strarr(2)
bar_lab(0)='energy flux (mW/m2)'
bar_lab(1)='mean energy (keV)'
if (iLC eq 1) then bar_lab(1)='(sinLn-sinLs)/(sinLn+sinLs)'
if (icon eq 1) then bar_lab(*)='conductance (mho)'
; Set up for plotting latitude and longitude grid
npt=100
ang=fltarr(npt)       &    radius=ang
dang=2.*!pi/(npt-1)
for j=0,npt-1 do ang(j)=j*dang

; y position of line plots
yi=0.15
yf=yi+0.18

; make time labels
time_lab=strarr(ntime)
for i=0,ntime-1 do begin
    time=houra(i)
    ihour=fix(time)
    hour=string(ihour,'(i3)') 
    imin=round((time-float(ihour))*60.)
    minute=string(imin,'(i2.2)')
    time_lab(i)=hour+':'+minute+' '
endfor

; Setup scale for total energy flux
tot_flx_mx=max(total_flx)
print,'max total E flux (GW) => ',tot_flx_mx
read,'Do you want to change it? no(n), yes(y) => ',yon
if (yon eq 'y') then read,'enter max E flux => ',tot_flx_mx
tot_flx_mn=0.

; Calculate Pedersen and Hall conductance
if (icon eq 1) then begin
   sigmaP=Efluxa     &    sigmaH=Efluxa
   for n=0,ntime-1 do begin
       for i=0,ir-1 do begin
           for j=0,ip do begin
               Eave=meanEa(n,i,j)
               sigmaP1=40.*Eave*sqrt(0.5*Efluxa(n,i,j))
               sigmaP(n,i,j)=2.*sigmaP1/(16.+Eave*Eave)
               sigmaH(n,i,j)=0.45*sigmaP(n,i,j)*Eave^0.85
           endfor
       endfor
   endfor
endif

if (iout eq 2) then begin
   set_plot,'ps'
   device,filename=filehead+Hhead+ftitle+Ehead+dhead+'.ps',/inches,yoffset=3.,$
           xoffset=0.3,xsize=8.,ysize=5.36,/color,bits_per_pixel=24
endif else begin
   device,decomposed=-1
   window,1,xpos=200,ypos=200,xsize=640,ysize=430
endelse   

; Plot ionospheric energy flux and mean energy or conductance
y_title='magnetic latitude (deg)'
if (nf eq 2) then y_title='magnetic latitude south (deg)'
for i_img=n0,ntime-1 do begin
    xyouts,0.5,y2+0.06,filehead+'  '+Ehead+' '+althead,size=1.5, $
           alignment=0.5,/normal 
    hour=houra(i_img)
    xo(*,*)=xoc(i_img,*,*)
    yo(*,*)=yoc(i_img,*,*)
    philab=mlon0(i_img)        ; angle of mlong-0
    if (idir eq 2) then philab=philab-!pi/2.
    if (idir eq 3) then philab=-philab+!pi
            
    for iplt=0,1 do begin
        if (icon eq 0) then begin
           if (iplt eq 0) then begin
              plt_v(*,*)=Efluxa(i_img,*,*)
              ; locate max Eflux
              maxflux=max(plt_v,location)
              ind=array_indices(plt_v,location)
           endif else begin
                 if (iLC eq 0) then plt_v(*,*)=meanEa(i_img,*,*)
                 if (iLC eq 1) then plt_v(*,*)=ALC(i_img,*,*)
           endelse
           lvlmn=lvl(0,iplt)
           if (iplt eq 1 and iLC eq 1) then lvlmn=lvl(nlevel-1,iplt)
           if (iLC eq 0) then plt_v(where(Efluxa(i_img,*,*) lt lvl(0,0)))=lvlmn
        endif else begin
           if (iplt eq 0) then plt_v(*,*)=sigmaP(i_img,*,*)
           if (iplt eq 1) then plt_v(*,*)=sigmaH(i_img,*,*)
        endelse

        xi=0.07+iplt*0.5
        xf=xi+(y2-y0)*0.67    
        plot,[0,1],[0,0],xrange=[-halfL,halfL],yrange=[-halfL,halfL], $
    ;   plot,[0,1],[0,0],xrange=[15,30],yrange=[-15,0], $
             xstyle=1,ystyle=1,pos=[xi,y0,xf,y2],title=ptitle(iplt), $
             ytitle=y_title,xthick=thk,ythick=thk, $
    ;        xticks=ntick,yticks=ntick,xtickn=ticknx,ytickn=tickny,/noeras
             xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn,/noeras
    ;   polyfill,xb,yb                ; black background
    ;   contour,plt_v,xo,yo,levels=lvl(*,iplt),c_colors=colr,/fill,/overplot
    ;   xyouts,xi,y2-0.03,slab,color=clab,/normal 
        xyouts,xf-0.07,y2-0.03,time_lab(i_img),color=clab,/normal 

        ; Mark Sun and mlong-0 directions
        if (idir eq 1) then xyouts,xi+0.01,0.5*(y0+y2),'Sun',color=clab,/normal
        if (idir eq 2) then xyouts,0.5*(xi+xf),y2-0.03,'Sun',color=clab, $
                                   alignment=0.5,/normal
        if (idir eq 3) then xyouts,xf-0.01,0.5*(y0+y2),'mlong-0',color=clab, $
                                   alignment=1.,/normal
        xlab=rlab*cos(philab)
        ylab=rlab*sin(philab)
  ;     oplot,[xlab,xlab],[ylab,ylab],psym=7,color=clab
    
        ; Mark Eflux maximum
        if (icon eq 0 and iplt eq 0) then begin
           ii=ind(0)   &   jj=ind(1)
         ; mcS=string(mc(i_img,ii,jj),'(i2)')
         ; oplot,[xo(ii,jj),xo(ii,jj)],[yo(ii,jj),yo(ii,jj)],psym=2,color=255
         ; xyouts,xi+0.02,y0+0.1,mcS+' pt in LC at max Eflux',color=255,/normal
        endif

        ; plot grids on contours
        for i=0,ir-1 do oplot,xo(i,*),yo(i,*)     ; azimuthal grid
        for j=0,ip do oplot,xo(*,j),yo(*,j)       ; radial grid
        for i=0,ir-2 do oplot,xo(i:i+1,0),yo(i:i+1,0),color=253-i*3,thick=3
       ;    oplot,/polar,[0.,80.],[phi1,phi1],thick=iout,linestyle=1,color=255
        
        ; plot color bar
        x1=xf+0.04
        x2=x1+0.02
        lvl2=nlevel-1
        ymx=y2
        ymn=y0
        if (iplt eq 1 and iLC eq 1) then begin
           lvl2=44
           ymx=(yave+y2)/2.
           ymn=y0+(y2-y0)/8.
        endif
        for i=0,lvl2 do begin
            y1=ymn+i*dy
            polyfill,[x1,x2,x2,x1,x1],[y1,y1,ymx,ymx,y1],color=colr(i),/normal
        endfor
        xyouts,x2,ymn,minV(iplt),/normal
        if (iplt eq 1 and iLC eq 1) then xyouts,x2,yave+ymn-y0,' 0.',/normal
        xyouts,x2,y1-0.01,maxV(iplt),/normal
        xyouts,x1-0.01,yave,bar_lab(iplt),orientation=90,alignment=0.5,/normal
     
        ; plot total flux 
        if (iplt eq 0) then begin
           mtitle='total energy flux (GW)'
           plot,houra(n0:ntime-1),total_flx(n0:ntime-1),xthick=thk,ythick=thk,$
                xrange=[hour_min,hour_max],yrange=[tot_flx_mn,tot_flx_mx], $
                xstyle=1,ystyle=1,pos=[xi,yi,x2,yf],title=mtitle,thick=thk,$
                xtitle='hour',ytitle='GW',xticks=x_tick,/noerase
           oplot,[hour,hour],[tot_flx_mn,tot_flx_mx],thick=thk
        endif

        ; plot Dst
        if (iplt eq 1) then begin
           plot,houra(n0:ntime-1),Dsta(n0:ntime-1),xrange=[hour_min,hour_max], $
                yrange=[Dst_min,Dst_max],xstyle=1,ystyle=1,title='Dst', $
                pos=[xi,yi,x2,yf],thick=thk,xthick=thk,ythick=thk, $
                xtitle='hour',ytitle='Dst (nT)',xticks=x_tick,/noerase 
           oplot,[hour,hour],[Dst_min,Dst_max],thick=thk
        endif
    endfor

    ; Make gif file
    if (iout eq 1) then begin
       gimg = tvrd(0,0)
     ; gimg = tvrd(280,0,319,170)   ; Dst plot only
       write_gif,filehead+ftitle+Ehead+'.gif',gimg,red(0:255),green(0:255), $
                 blue(0:255),/multiple
    endif

    if (i_img lt ntime-1) then erase
endfor

if (iout eq 2) then device,/close_file
if (iout eq 1) then write_gif,filehead+ftitle+Ehead+'.gif',/close
end
