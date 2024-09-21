;-------------------------------------------------------------------------------
;
;                           plot_fls_L_Time.pro
;
; IDL procedure reads *.fls file and plot the pitch-angle local-time averaged
; flux in a L-time plot.
;
; Created on 9 November 2005 by Mei-Ching Fok, Code 612.2, NASA GSFC.
;
; Modification History
; March 18, 2010
;   * Overplot Dst is added.
; February 25, 2011
;   * New way to calculate fluxes at ionosphere
; October 28, 2011
;   * Read Dst from the *.fls file
; April 9, 2012
;   * Only plot flux at and below iba
; July 12, 2012
;   * add color bar option
; September 24, 2014
;   * add the option to plot invariant L or geocentric L (ro)
; November 5, 2015
;   * add the option to plot in L*
;-------------------------------------------------------------------------------

close,/all

; Setup plot  
x_wsize=600
y_wsize=450

; get colors for color table palette
  loadct,39       ;  rainbow with white
  tvlct, red, green, blue, /get
  !p.background=255
  !p.charsize=1.5

; read file name and energy, sina, Lshell0 information 
fhead=' '
read,'enter the file head (e.g., 2000a225_h) => ',fhead
openr,2,fhead+'.fls'
readf,2,rc,ir,ip,je,ig   ; je=number of energy grid, ig=number of pa grid
nr=ir                    ; number of radial grid
nmlt=ip                  ; number of local-time grid
Lshell0=fltarr(ir)  &   mphi=fltarr(ip)
energy=fltarr(je)   &   sina=fltarr(ig)   &   anor=sina
cosai=sina
Ebound=fltarr(je+1)
readf,2,Lshell0
readf,2,mphi
readf,2,energy
readf,2,sina
iplt=1
;read,'plot? 1=flux@equator,  2=flux@ionosphere => ',iplt
y_title='L shell'
iLs=0
if (iplt eq 1) then begin
   read,'plot at? 0=L*, 1=L-invariant, 2=L-geocentric, 3=L-McIlwain => ',iLs
   if (iLs eq 0) then y_title='L*' 
   if (iLs eq 1) then y_title='L-invariant' 
   if (iLs eq 2) then y_title='L-geocentric'
   if (iLs eq 3) then y_title='L-McIlwain'  
endif
alt=9999.
if (iplt eq 2) then read,'altitude (km) => ',alt
re_m=6.375e6                   ; Earth radius in m
xme=7.764748e15                ; magnetic dipole moment of the Earth
ri=(alt+re_m)/re_m             ; ionosphere distance in RE
rim=ri*re_m                    ; ionosphere distance in m
beq=xme/re_m^3                 ; B at equator at Earth surface (assume dipole)
read,'number of data set in time => ',ntime
read,'plot Dst? 0=no, 1=yes => ',iDst
read,'plot  1 = one frame,  2 = multiple frame w/ time bar  => ',imu   
read,' plot format? (1) gif,  (2) ps => ',iout
thk=iout*2.0
nmax=0
if (imu gt 1) then nmax=ntime-1
n2=ntime/2
parmod=fltarr(10)
!p.charthick=thk

; L-values for plotting
Lmin=1.
Lmax=7.

; define energy labels
for k=1,je-1 do Ebound(k)=sqrt(energy(k-1)*energy(k))
Ebound(0)=energy(0)*energy(0)/Ebound(1)
Ebound(je)=energy(je-1)*energy(je-1)/Ebound(je-1)
if (Ebound(je) lt 999.5) then ilab='(i3.3)'
if (Ebound(je) ge 999.5) then ilab='(i4.4)'
if (Ebound(je) lt 999.5) then flab='(f5.1)'
if (Ebound(je) ge 999.5) then flab='(f6.1)'
 
; setup energy and species
nsp=strlen(fhead)-1
spe=strmid(fhead,nsp,1)                ; species
if (spe eq 'h') then species='H+'
if (spe eq 'o') then species='O+'
if (spe eq 'e' or spe eq '1') then species='e-'
if (spe eq 's') then species='solar wind H+'

; setup arrays
ro=fltarr(nr,nmlt)   &   Lshell=ro   &   Lmc=ro   &   ompeij=ro   &   Lstar=ro
plsflux3=fltarr(nr,nmlt,je)  
plsfluxk=fltarr(ntime,nr,je)    
plsflux=fltarr(ntime,nr)     
denij=fltarr(ntime,nr,nmlt)   &   ompein=denij
hour=fltarr(ntime)        &     Dst=hour      &     Dsts=hour   &  day=hour
fy=fltarr(ig)      &      dmu=fy        &    cosa=fy
cosa(*)=cos(asin(sina(*)))
for m=0,ig-1 do begin
    if (m eq 0) then cosa0=1.
    if (m gt 0) then cosa0=0.5*(cosa(m)+cosa(m-1))
    if (m eq (ig-1)) then cosa1=0.
    if (m lt (ig-1)) then cosa1=0.5*(cosa(m)+cosa(m+1))
    dmu(m)=cosa0-cosa1
endfor
dE=fltarr(je)
for k=0,je-1 do dE(k)=Ebound(k+1)-Ebound(k)

; Define plasmapause
npp=40.    ; density in cm^-3 defining the plasmapause
print,'plasmapause is defined at density (cm-3) => ',npp
Lppm=fltarr(ntime)   &   Lppa=Lppm
Lppj=fltarr(nmlt)

; Read fluxes and calculate pa-mlt averaged flux
pi2=!pi/2.
anor(*)=asin(sina(*))/pi2
index=0.8
parmod8=fltarr(8)
Lstar_max=Lmax
Lstar1=Lmax
for n=0,ntime-1 do begin
    readf,2,hour1,Pdyn,Dst1,parmod8,Lstar_max
    hour(n)=hour1
    day(n)=hour1/24.
    print,n,hour1,Lstar_max
    Dst(n)=Dst1
    for i=0,nr-1 do begin
        for j=0,nmlt-1 do begin
            readf,2,lat1,xmlt,xlatS,xmltS,ro1,xmlto,briN,briS,bo,iba,den0, $
                    ompe0,CHp,HIp,denWP,Tpar,Tper,HRPe,HRPi,rppa,Lstar1
            cosl=cos(lat1*!pi/180.)
            Lshell(i,j)=rc/cosl/cosl
            bi=sqrt(4.-3.*ri/Lshell(i,j))*xme/rim^3       ; assume dipole
            ro(i,j)=ro1
            if (Lstar1 lt 0.) then Lstar1=Lstar_max
            Lstar(i,j)=Lstar1
            Lmc(i,j)=(Beq/bo)^0.33333
            denij(n,i,j)=den0*1.e-6     ; plasmasphere density in cm^-3
            if (ompe0 le 0.) then ompe0=ompeij(i-1,j)
            ompeij(i,j)=ompe0
            sqrtbb=sqrt(bi/bo)
            bibo=(bi/bo)^index
            if (iplt eq 2 and Lshell(i,j) gt ri) then begin     ; redo dmu
               for m=0,ig-1 do begin
                   sinai=sina(m)*sqrtbb 
                   if (sinai le 1.) then ai=asin(sinai)
                   if (sinai gt 1.) then ai=(1.-(1.-anor(m))^bibo)*pi2
                   cosai(m)=cos(ai)
               endfor
               for m=0,ig-1 do begin
                   if (m eq 0) then mu0=1.
                   if (m gt 0) then mu0=0.5*(cosai(m)+cosai(m-1))
                   if (m eq ig-1) then mu1=0.
                   if (m lt ig-1) then mu1=0.5*(cosai(m)+cosai(m+1))
                   dmu(m)=mu0-mu1
                   if (dmu(m) lt 0.) then dmu(m)=0.
               endfor
            endif
            for k=0,je-1 do begin
                readf,2,fy
                plsflux3(i,j,k)=0.    ; pitch-angle averaged flux
                if (i lt iba) then begin
                   for m=0,ig-1 do begin
                       plsflux3(i,j,k)=plsflux3(i,j,k)+fy(m)*dmu(m)
                   endfor
                endif
                if(iLs eq 0 and Lstar1 ge Lstar_max)then plsflux3(i,j,k)=0.
            endfor
        endfor  
    endfor     
    ; map flux, plasma density and ompe to fixed L grid: Lshell0 
    for j=0,nmlt-1 do begin
        if (Lshell(0,j) gt Lshell0(0)) then Lshell(0,j)=Lshell0(0) ;Limit Lshell
        if (Lmc(0,j) gt Lshell0(0)) then Lmc(0,j)=Lshell0(0)   ; Limit Lmc 
        if (iLs eq 0) then ii=interpol(findgen(nr),Lstar(*,j),Lshell0)
        if (iLs eq 1) then ii=interpol(findgen(nr),Lshell(*,j),Lshell0)
        if (iLs eq 2) then ii=interpol(findgen(nr),ro(*,j),Lshell0)
        if (iLs eq 3) then ii=interpol(findgen(nr),Lmc(*,j),Lshell0)
        den1=interpolate(denij(n,*,j),ii,missing=0.)
        denij(n,*,j)=den1(*)
        ompe1=interpolate(ompeij(*,j),ii,missing=0.)
        ompein(n,*,j)=ompe1(*)        ; interpolated ompe
        for k=0,je-1 do begin
            plsflux1=interpolate(plsflux3(*,j,k),ii,missing=0.)
            plsflux3(*,j,k)=plsflux1(*)
        endfor
    endfor
    ; find pa-mlt averaged flux
    for i=0,nr-1 do begin
        plsfluxk(n,i,0:je-1)=0.   
        if (Lshell0(i) gt ri or iplt eq 1) then begin
           for k=0,je-1 do plsfluxk(n,i,k)=total(plsflux3(i,*,k))/nmlt 
        endif
    endfor

    ; Locate plasmapause
    Lppm(n)=100.    ; arbitray a large value
    for j=0,nmlt-1 do begin
        Lpp1=interpol(Lshell0,denij(n,*,j),[npp,npp])
        Lppm(n)=min([Lppm(n),Lpp1(0)])  ; define plasmapause by min Lpp over mlt
        Lppj(j)=Lpp1(0)
    endfor
    Lppa(n)=total(Lppj)/nmlt        ; define plasmapause by average Lpp over mlt
endfor             ; n loop
close,2

; Find average, max and min ompe
ompe=fltarr(ntime,nr)  
ompe(*,*)=0.    &   ompeMX=ompe  &   ompeMN=ompe  ; initial values
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        j1=0
        for j=0,nmlt-1 do begin
            if (ompein(n,i,j) gt 0.) then begin
               if (j1 eq 0.) then ompeMN(n,i)=100.  ; initial value
               j1=j1+1
               ompe(n,i)=ompe(n,i)+ompein(n,i,j)
               ompeMX(n,i)=max([ompeMX(n,i),ompein(n,i,j)])
               ompeMN(n,i)=min([ompeMN(n,i),ompein(n,i,j)])
            endif
        endfor
        if (j1 gt 0) then ompe(n,i)=ompe(n,i)/j1  
    endfor
endfor

; Find average, max and min plasmasphere density
plsden=fltarr(ntime,nr)
plsden(*,*)=0.    &   plsdenMX=ompe  &   plsdenMN=ompe  ; initial values
for n=0,ntime-1 do begin
    for i=0,nr-1 do begin
        j1=0
        for j=0,nmlt-1 do begin
            if (denij(n,i,j) gt 0.) then begin
               if (j1 eq 0.) then plsdenMN(n,i)=9.e9  ; initial value
               j1=j1+1
               plsden(n,i)=plsden(n,i)+denij(n,i,j)
               plsdenMX(n,i)=max([plsdenMX(n,i),denij(n,i,j)])
               plsdenMN(n,i)=min([plsdenMN(n,i),denij(n,i,j)])
            endif
        endfor
        if (j1 gt 0) then plsden(n,i)=plsden(n,i)/j1
    endfor
endfor
    
; Scale Dst
if (iDst eq 1) then begin
   Lrange=Lmax-Lmin
   Drange10=ceil(max(Dst)/10.)+ceil(abs(min(Dst))/10.)
   Drange=10*Lrange*ceil(Drange10*1./Lrange)
   deltaD=Drange/Lrange
   if (deltaD gt 60.) then deltaD=60.
   Dmax=deltaD
   if (Dmax lt max(Dst)) then Dmax=ceil(max(Dst)/deltaD)*deltaD
   Dmin=Dmax-Drange
   if (Dmin gt min(Dst)) then begin
      Dmin=Dmin-deltaD
      Dmax=Dmin+Drange
   endif
   delta=Lmin-Dmin/deltaD
   for i=0,ntime-1 do Dsts(i)=Dst(i)/deltaD+delta
endif

if (iout eq 2) then begin
   set_plot,'ps'
   file_name=fhead+'_LTimePlot.ps'
   device,filename=file_name,/inches,/color,yoffset=3.0, $
          xoffset=0.3,xsize=7.,ysize=7.*y_wsize/x_wsize
endif

new_plot:
; choose which energy to be displaced
print,'energy bin (keV): '
for i=0,je-1 do print,i,Ebound(i:i+1),format='("  (",i2,") ",f9.4," - ",f9.4)'
print,fix(je),format='("  (",i2,") average fpe/fce")'
print,fix(je+1),format='("  (",i2,") maximum fpe/fce")'
print,fix(je+2),format='("  (",i2,") minimum fpe/fce")'
print,fix(je+3),format='("  (",i2,") average plasmasphere density")'
print,fix(je+4),format='("  (",i2,") maximum plasmasphere density")'
print,fix(je+5),format='("  (",i2,") minimum plasmasphere density")'
read, 'lower and upper bins (ie1,ie2) => ',ie1,ie2
iunit=0   
if (ie2 lt je) then begin
   read,'unit? (0)(MeVcm2ssr)^-1, (1)(keVcm2ssr)^-1, (2)(cm2ssr)^-1 => ',iunit
   e0=Ebound(ie1)
   e1=Ebound(ie2+1)
   fmiddle=string(round(e0),ilab)+'-'+string(round(e1),ilab)+'keV'
endif
if (ie2 eq je) then fmiddle='aveOmpe'
if (ie2 eq je+1) then fmiddle='maxOmpe' 
if (ie2 eq je+2) then fmiddle='minOmpe' 
if (ie2 eq je+3) then fmiddle='avePlsDen'
if (ie2 eq je+4) then fmiddle='maxPlsDen'
if (ie2 eq je+5) then fmiddle='minPlsDen'
ftail='.'
if (ie2 lt je) then ftail='_LTflux.'

; calculate total flux 
if (ie2 lt je) then begin
   for n=0,ntime-1 do begin
       for i=0,nr-1 do begin
           plsflux(n,i)=0.
           weiSum=0.
           for k=ie1,ie2 do begin    
               weiSum=weiSum+dE(k)
               plsflux(n,i)=plsflux(n,i)+plsfluxk(n,i,k)*dE(k)
           endfor
           if (iunit le 1) then plsflux(n,i)=plsflux(n,i)/weiSum  ; diff flux
           if (iunit eq 0) then plsflux(n,i)=1000.*plsflux(n,i)   ; diff flux
       endfor
   endfor
endif
if (ie2 eq je) then plsflux=Ompe
if (ie2 eq je+1) then plsflux=OmpeMX
if (ie2 eq je+2) then plsflux=OmpeMN
if (ie2 eq je+3) then plsflux=plsden
if (ie2 eq je+4) then plsflux=plsdenMX
if (ie2 eq je+5) then plsflux=plsdenMN

read,'flux in?  0 = linear scale, or 1 = log scale => ',ilog
; Find log value if ilog=1
if (ilog eq 1) then begin
   for n=0,ntime-1 do begin
       for i=0,nr-1 do begin
           if (plsflux(n,i) gt 1.e-20) then plsflux(n,i)=alog10(plsflux(n,i)) $
                                       else plsflux(n,i)=-20.
       endfor
   endfor
endif

; setup plot ranges
yon=' '
fmax=max(plsflux)    
print,' fmax = ',fmax
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmax => ',fmax
fmin=0.   
if (ilog eq 1) then fmin=fmax-5.
print,' fmin = ',fmin
read,' Do you want to change it? (y/n) => ',yon
if (yon eq 'y') then read,' enter new fmin => ',fmin
if(min(plsflux) lt fmin) then plsflux(where(plsflux lt fmin))=fmin
cDst=0   
;if (iDst eq 1) then read,'color of Dst curve? 0=dark, 255=white => ',cDst

; setup levels
nlevel=59
lvl=fltarr(nlevel)    &    colr=intarr(nlevel)
dlvl=(fmax-fmin)/(nlevel-1)
colrmax=254           
colrmin=1
dcolr=(float(colrmax)-float(colrmin))/(nlevel-1)
for i=0,nlevel-1 do begin
    lvl(i)=fmin+i*dlvl
    colr(i)=round(float(colrmin)+i*dcolr)
endfor

; smooth in time            
yon='n'
;read,'Do you want to smooth flux in time? y(yes), n(no) => ',yon
plsfluxt=plsflux
if (yon eq 'y') then begin
   read,'enter number of time sectors in smoothing => ',ntsec
   nts=fix(ntsec/2)
   for n=0,ntime-1 do begin
       n0=n-nts
       if (n0 lt 0) then n0=0
       n1=n+nts
       if (n1 ge ntime) then n1=ntime-1
       for i=0,nr-1 do plsfluxt(n,i)=(total(plsflux(n0:n1,i)))/(n1-n0+1)
   endfor
endif    

; Average in time            
yon='n'
;read,'Do you want to average flux in time? y(yes), n(no) => ',yon
plsfluxt=plsflux
if (yon eq 'y') then begin
   read,'enter number of time sectors in averaging => ',ntsec
   for n=0,ntime-1 do begin
       n0=fix(n/ntsec)*ntsec
       n1=n0+ntsec
       if (n1 ge ntime) then n1=ntime-1
       for i=0,nr-1 do plsfluxt(n,i)=(total(plsflux(n0:n1,i)))/(n1-n0+1)
   endfor
endif    

; smooth flux in L
yon='n'
;read,'Do you want to smooth flux in L? y(yes), n(no) => ',yon
plsfluxs=plsfluxt
if (yon eq 'y') then begin
   read,'enter number of L point before and after => ',nLL
   for i=0,nr-1 do begin
       i0=i-nLL
       if (i0 lt 0) then i0=0
       i1=i+nLL
       if (i1 ge nr) then i1=nr-1
       for n=0,ntime-1 do plsfluxs(n,i)=(total(plsfluxt(n,i0:i1)))/(i1-i0+1)
   endfor
endif

; determine plot plasmapause or not
iLpp=' '
read,'plot plasmapause? n=no, 1=min-Lpp, 2=averaged-Lpp => ',iLpp

; Plot plasma flux 
x0=0.12
xf=0.78
y0=0.16 
yf=0.81
ym=0.5*(y0+yf)
if (ie2 lt je) then begin
   elabel=string(e0,flab)+' - '+string(e1,flab)
   mtitle=elabel+' keV '+species
endif
if (ie2 eq je) then mtitle='MLT-averaged fpe/fce'
if (ie2 eq je+1) then mtitle='Max fpe/fce in MLT'
if (ie2 eq je+2) then mtitle='Min fpe/fce in MLT'
if (ie2 eq je+3) then mtitle='MLT-averaged plasmasphere density'
if (ie2 eq je+4) then mtitle='Max plasmasphere density in MLT'
if (ie2 eq je+5) then mtitle='Min plasmasphere density in MLT'
thr=12
hour_max=ceil(max(hour)/thr)*thr
hour_min=fix(min(hour)/thr)*thr
day_max=ceil(max(day))
day_min=fix(min(day))
for n=0,nmax do begin
  ; contour,plsfluxs,hour,Lshell0,yrange=[Lmin,Lmax],ystyle=1,xtitle='hour', $
  ;     ytitle=y_title,xrange=[hour_min,hour_max],xstyle=1, $
    contour,plsfluxs,day,Lshell0,yrange=[Lmin,Lmax],ystyle=1,xtitle='day', $
        ytitle=y_title,xrange=[day_min,day_max],xstyle=1, $
        title=mtitle,charsize=2.0,pos=[x0,y0,xf,yf],color=0,levels=lvl, $
        c_colors=colr,/fill

    ; plot plasmapause
    if (iLpp ne 'n') then begin
     ; if (iLpp eq '1') then oplot,hour,Lppm,thick=thk,max_value=Lmax,color=255
       if (iLpp eq '1') then oplot,day,Lppm,thick=thk,max_value=Lmax,color=255
       if (iLpp eq '1') then xyouts,x0+0.03,y0+0.2,'min-Lpp',color=255,/normal
     ; if (iLpp eq '2') then oplot,hour,Lppa,thick=thk,max_value=Lmax,color=255
       if (iLpp eq '2') then oplot,day,Lppa,thick=thk,max_value=Lmax,color=255
       if (iLpp eq '2') then xyouts,x0+0.03,y0+0.2,'ave-Lpp',color=255,/normal
    endif 

    ; Plot Dst
    if (iDst eq 1) then begin
     ; oplot,hour,Dsts,thick=3,color=cDst
       oplot,day,Dsts,thick=3,color=cDst
       for L=Lmin,Lmax-1 do begin
           D_st=string(Dmin+(L-Lmin)*deltaD,'(i4)')
        ;  xyouts,hour_max,L,D_st,size=2.0,color=0
           xyouts,day_max,L,D_st,size=2.0,color=0
       endfor
   ;   xyouts,hour_max,Lmax,' Dst',size=2.0,color=0
       xyouts,day_max,Lmax,' Dst',size=2.0,color=0
    endif

    ; draw color bar and label
    x1=xf+0.14
    x2=x1+0.03
    dy=(yf-y0)/(colrmax-colrmin+1)
    for i=colrmin,colrmax do begin
        y1=y0+i*dy
        y2=y1+dy*1.03
        polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=i,/normal  
    endfor
    if (iunit eq 0) then bar_lab='flux(/MeV/cm2/sr/s)'
    if (iunit eq 1) then bar_lab='flux(/keV/cm2/sr/s)'
    if (iunit eq 2) then bar_lab='flux(/cm2/sr/s)'
    if (ie2 ge je and ie2 le je+2) then bar_lab='fpe/fce'
    if (ie2 ge je+3) then bar_lab='density (/cm3)'
    if (ilog eq 1) then bar_lab='log '+bar_lab
    xyouts,x2,y0-0.02,string(fmin,'(f4.1)'),size=2.0,color=0,/normal
    xyouts,x2,yf,string(fmax,'(f4.1)'),size=2.0,color=0,/normal
    xyouts,x2+0.04,ym,bar_lab,alignment=0.5,orientation=90.,color=0,/normal
    xyouts,x0,0.93,strmid(fhead,0,nsp-1),size=2.,color=0,/normal
    iplt_lab=' '
    if (ie2 lt je) then begin
       if (iplt eq 1) then iplt_lab='equatorial flux'
       if(iplt eq 2) then iplt_lab=string(fix(alt),format='(i4)')+' km altitude'
       xyouts,xf,0.93,iplt_lab,alignment=1,color=0,/normal
    endif
    
    ; plot time bar when n gt 0
;   if (n gt 0) then begin
;      oplot,[hour(n),hour(n)],[Lshell0(0),Lshell0(ir-1)],color=cDst,thick=3
;      tlab=' t(hr)='+string(hour(n),'(f5.1)')
;      ali=0.
;      if (n gt n2) then ali=1.  
;      xyouts,hour(n),7.5,tlab,alignment=ali,color=cDst
;   endif

    ; Make gif file
    if (iout eq 1) then begin
       gimg = tvrd(0,0)
       write_gif,fhead+'_'+fmiddle+ftail+'gif',gimg,red(0:255),green(0:255), $
                 blue(0:255),/multiple
    endif
endfor

if (iout eq 1) then write_gif,fhead+'_'+fmiddle+ftail+'gif',/close
if (iout eq 2) then erase
read,'Do you want to continue?  (y)yes, (n)no => ',yon
if (yon eq 'y') then goto,new_plot
if (iout eq 2) then device,/close_file

end
