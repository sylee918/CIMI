;--------------------------------------------------------------------------
;                              equat_pot_h.pro
;
; Procedure polts the equatorial potential, and drift paths output from 
; ring current models.
;
; Input file: *.pot
;
; Modification history
; Septeber 23, 2016
;  - Add the capability to read file with high-latitude potential
; August 25, 2017
;  - Add the capability to plot electric field
;--------------------------------------------------------------------------

close,/all

; get colors for color table palette
read,'colar table, 1=rainbow, 2=blue-red => ',icr
if (icr eq 1) then loadct,39       ;  rainbow with white
if (icr eq 2) then loadct,33       ;  blue-red
tvlct, red, green, blue, /get

; constants
re_m=6.375e6         ; Earth radius in m
xme=7.750871e15      ; magnetic dipole moment of the Earth

; initial setup
filein=' '
read,' enter the file head of potential (i.e., 2013d151) => ', filein
openr,1,filein+'.pot'
read,' enter the half length of the plots => ', halfl
read,' Sun direction?  (1)left,  (-1)right  => ',isun
;read,' Do you want to highlight a contour? no(0), yes(1) => ', ihigh 
ihigh=0
read,' potential interval for convection in kV? => ',dlevelp
read,' add? corotation(0), DriftPath(1), EField(2), Ey(3), Vr(4), Vp(5) => ',imb
icharge=1
if (imb eq 1) then read,' drift paths for? ions(1), electrons(-1) => ', icharge
read,' plot format? (1) gif,  (2) ps => ',iout
read,' enter number of data set => ', ntime

; Read data
readf,1,rc         ; rc: ionospheric distance in RE
readf,1,idim
readf,1,irh
readf,1,jdim_1
readf,1,mdim           ; no. of longitudinal invariant
lat_1=fltarr(idim,jdim_1)   &   potentc_1=lat_1   &   BriN=lat_1
lat=fltarr(idim,jdim_1+1)
lath=fltarr(irh)
readf,1,lat_1           ; latitude in degree
readf,1,lath
readf,1,potentc_1       ; originally readf,1,alpha
readf,1,BriN 
jdim=jdim_1+1
hour_1=fltarr(jdim_1)
t_hour=fltarr(ntime)      
potent_1=fltarr(idim,jdim_1)   & ro_1=potent_1    &  mlto_1=ro_1    &  bo_1=ro_1
cond_1=potent_1   & birk_1=potent_1  &   bobi_1=ro_1   
xlatS=fltarr(idim,jdim_1)      & phiS=xlatS
potenth=fltarr(irh,jdim_1)
yo_1=fltarr(idim,jdim_1,mdim)
potent=fltarr(ntime,idim,jdim) & ro=potent      &  mlto=potent   &  bo=potent
bobi=potent  
yo=fltarr(ntime,idim,jdim,mdim)
potent1=fltarr(idim,jdim)   & ro1=potent1  & mlto1=ro1  &  bo1=ro1  
potent2=potent1   &    energy1=potent1     & ym=ro1     &  lat1=ro1 & bobi1=ro1
EorV=potent1      &    potentc=potent1    
yo1=fltarr(idim,jdim,mdim)
xe=fltarr(idim,jdim)           & ye=xe       
potcap=fltarr(jdim_1)
iba1=intarr(jdim_1)
iba=intarr(ntime,jdim)
cor1=2.*!pi/86400./1000.     ; constant factor in corotation
lat1(0:idim-1,0:jdim_1-1)=lat_1(0:idim-1,0:jdim_1-1)             ; in deg
lat1(0:idim-1,jdim-1)=lat1(0:idim-1,0)
potentc(0:idim-1,0:jdim_1-1)=potentc_1(0:idim-1,0:jdim_1-1)/1000. ; con. V to kV
potentc(0:idim-1,jdim-1)=potentc(0:idim-1,0)
for n=0,ntime-1 do begin
   readf,1,t_hr
   readf,1,iba1
   readf,1,hour_1
   readf,1,xlatS       
   readf,1,phiS       
   readf,1,ro_1
   readf,1,mlto_1
   readf,1,bo_1
   readf,1,yo_1
   readf,1,potent_1,format='(8f10.1)'
   readf,1,potenth,format='(8f10.1)'
   readf,1,cond_1
   readf,1,birk_1
   t_hour(n)=t_hr
   print,'n,t_hr ',n,t_hr
   iba(n,0:jdim_1-1)=iba1(0:jdim_1-1)
   iba(n,jdim-1)=iba(n,0)
   ro(n,0:idim-1,0:jdim_1-1)=ro_1(0:idim-1,0:jdim_1-1)               ; in RE
   ro(n,0:idim-1,jdim-1)=ro(n,0:idim-1,0)
   mlto(n,0:idim-1,0:jdim_1-1)=mlto_1(0:idim-1,0:jdim_1-1)           ; in hour
   mlto(n,0:idim-1,jdim-1)=mlto(n,0:idim-1,0)+24.
   bo(n,0:idim-1,0:jdim_1-1)=bo_1(0:idim-1,0:jdim_1-1)              
   bo(n,0:idim-1,jdim-1)=bo(n,0:idim-1,0)
   bobi(n,0:idim-1,0:jdim_1-1)=bobi_1(0:idim-1,0:jdim_1-1)              
   bobi(n,0:idim-1,jdim-1)=bobi(n,0:idim-1,0)
   yo(n,0:idim-1,0:jdim_1-1,0:mdim-1)=yo_1(0:idim-1,0:jdim_1-1,0:mdim-1)        
   yo(n,0:idim-1,jdim-1,0:mdim-1)=yo(n,0:idim-1,0,0:mdim-1)
   potent(n,0:idim-1,0:jdim_1-1)=potent_1(0:idim-1,0:jdim_1-1)/1000.  ; in kV
   potent(n,0:idim-1,jdim-1)=potent(n,0:idim-1,0)
endfor
close,1

; Setup for plotting
y_title='noon'
if (isun eq -1) then y_title='midnight'
ithick=1.5*iout
!p.charthick=ithick 
!p.charsize=1.5 
x1=0.1
x2=x1+0.38
x3=x2+0.1
x4=x3+0.38
y1=0.22
y2=y1+0.56
pos1=[x1,y1,x2,y2]
pos2=[x3,y1,x4,y2]
if (iout eq 2) then begin
   set_plot,'ps'
    device,filename=filein+'_poteq.ps',/inches,/color,yoffset=2.5,xoffset=0.3, $
           ysize=5.43,xsize=8.
endif
y_o_n=' '

; setup labels on axes
ntick=4
dxy=2.*halfl/ntick
tickn=strarr(ntick+1)
for i=0,ntick do begin
    xy_1=-halfl+i*dxy
    xy_a=abs(xy_1)
    tickn(i)=string(fix(xy_a),'(i2)')
    if (isun*xy_1 gt 0.) then tickn(i)='-'+tickn(i)
endfor

; Determine potential levels and colors
nlevel=100 
ip=nlevel-1      ; or smaller number as desired
plvl=fltarr(nlevel)  
plab1=intarr(nlevel)    &   plab0=plab1 
pline_style=intarr(nlevel)
pothm=0.             ; middle level
lvlm1=fix(nlevel*0.5)
zero_or_1=0
for i=0,nlevel-1 do begin
    plvl(i)=pothm+(i-lvlm1)*dlevelp
    if (zero_or_1 eq 0) then zero_or_1=1 else zero_or_1=0
    plab1(i)=zero_or_1     ; label every other contour
    plab0(i)=0
    pline_style(i)=0
endfor
hlvl=plvl
c_color=intarr(nlevel)
c_color(*)=0      ;  black colors
if (ihigh eq 1) then c_color(25)=225   ; highlight one contour

; Determine dphi2
dphi2=fltarr(jdim_1)
for j=0,jdim_1-1 do begin
    j0=j-1
    if (j0 lt 0) then j0=j0+jdim_1
    j1=j+1
    if (j1 ge jdim_1) then j1=j1-jdim_1
    dmlt2=hour_1(j1)-hour_1(j0)
    if (dmlt2 lt 0.) then dmlt2=dmlt2+24.
    dphi2(j)=dmlt2*!pi/12.
endfor

; Determine Efield or drift velocity level, color, etc if imb ge 2
rc_m=rc*re_m
xmer3=xme/rc_m^3
mlevel=26  
elvl=fltarr(mlevel)   &   colr=intarr(mlevel)
yb=y1-0.15
yi=yb-0.04
dx=(x4-x3)/mlevel
fmin=0.
fmax=3. 
if (imb ge 4 and halfl gt 8.) then fmax=50.     
if (imb ge 3) then fmin=-fmax
fmxlab=string(fix(fmax),'(i2)')
fmnlab=string(fix(fmin),'(i3)')
if (imb ge 2) then begin
   E_lab='E field'
   if (imb eq 3) then E_lab='Ey'
   if (imb ge 4) then E_lab='drift velocity (km/s)'
   dlvl=(fmax-fmin)/(mlevel-1)
   colrmax=254.
   colrmin=50.
   dcolr=(colrmax-colrmin)/(mlevel-1)
   for i=0,mlevel-1 do begin
       elvl(i)=fmin+i*dlvl
       colr(i)=round(colrmin+i*dcolr)
   endfor
endif

; Determine particle energy and reference point
if (imb eq 1) then begin
   read, ' enter hour index, nloc => ',nloc  
   choose_ij:
   read,' enter the i,j indices of the reference location => ',iloc,jloc
   radius=ro(nloc,iloc,jloc)
   mltl=mlto(nloc,iloc,jloc)        
   print,' reference ro, mlto => ', radius,mltl
   read, ' do you want a new set of i,j indices? yes(y), no(n) => ',y_o_n
   if (y_o_n eq 'y') then goto,choose_ij
   choose_mbin:
   read, ' enter m bin: (m=0)perpendicular, (m=mdim-1)parallel => ',mbin
   ym(*,*)=yo(nloc,*,*,mbin)
   pitchA=asin(ym(iloc,jloc))*180./!pi
   print,' pitch angle => ',pitchA
   read, ' do you want a new mbin? yes(y), no(n) => ',y_o_n
   if (y_o_n eq 'y') then goto,choose_mbin
   read, ' enter particle energy in keV at reference location => ', keV
   mu=keV*ym(iloc,jloc)*ym(iloc,jloc)/bo(nloc,iloc,jloc)
   print,' mu => ',mu
   mulab=string(mu,'(e8.1)')+' keV/T'
   for j=0,jdim-1 do begin
       potent2(*,j)=icharge*(potent(nloc,*,j)+potentc(*,j))
   endfor
   ; setup H level
   nhlvl=30
   hlvl=fltarr(nhlvl)
   keV1=keV*3
   hrange=max(potent2)-min(potent2)+keV1 
   hlvl0=min(potent2)-keV1/2.
   dlvl=hrange/(nhlvl-1)
   for i=0,nhlvl-1 do hlvl(i)=hlvl0+i*dlvl
   ; Determine energy levels
   lvlm=(mlevel-1)/2
   dlevel=max([2,round(keV/5.)])
   dlevel=dlevelp     
   eline_style=intarr(mlevel)  &  elab=eline_style
   eline_style(*)=1
   engm=round(keV)
   zero_or_1=0
   for i=0,mlevel-1 do begin
       elvl(i)=engm+(i-lvlm)*dlevel
       if (zero_or_1 eq 0) then zero_or_1=1 else zero_or_1=0
   ;   elab(i)=zero_or_1     ; label every other contour
       elab(i)=0             ; no label 
   endfor
endif

; start the time loop
it=-1
newplot:

  ; Calculate the x-y grid at the equator
  if (it lt ntime) then it=it+1  
  if (it ge ntime) then it=-1
  if (it lt 0) then goto, endjob
  potent1(0:idim-1,0:jdim-1)=potent(it,0:idim-1,0:jdim-1)
  ro1(0:idim-1,0:jdim-1)=ro(it,0:idim-1,0:jdim-1)
  mlto1(0:idim-1,0:jdim-1)=mlto(it,0:idim-1,0:jdim-1)*!pi/12.   ; mlto1 in rad
  bo1(0:idim-1,0:jdim-1)=bo(it,0:idim-1,0:jdim-1)
  bobi1(0:idim-1,0:jdim-1)=bobi(it,0:idim-1,0:jdim-1)
  yo1(0:idim-1,0:jdim-1,0:mdim-1)=yo(it,0:idim-1,0:jdim-1,0:mdim-1)
  for j=0,jdim-1 do begin
    if (iba(it,j) lt idim) then begin
       potent1(iba(it,j):idim-1,j)=potent1(iba(it,j)-1,j)
       for m=0,mdim-1 do yo1(iba(it,j):idim-1,j,m)=yo1(iba(it,j)-1,j,m)
    endif
  endfor
  for j=0,jdim-1 do begin
      for i=0,idim-1 do begin
          roh=ro1(i,j)
          if (roh gt halfl) then roh=halfl
          xe(i,j)=isun*roh*cos(mlto1(i,j))
          ye(i,j)=isun*roh*sin(mlto1(i,j))
      endfor
  endfor
  
  ; Add corotation 
  potent2=potent1
  if (imb le 1) then begin
     for j=0,jdim-1 do begin
       potent2(*,j)=icharge*(potent1(*,j)+potentc(*,j))  
     endfor
  endif
  
  ; Add energy1
  if (imb eq 1) then begin
     energy1(*,*)=mu*bo1(*,*)/yo1(*,*,mbin)^2 
     potent2(*,*)=potent2(*,*)+energy1(*,*)
  endif

  ; Populate potentials beyond iba
  for j=0,jdim-1 do begin
    if(iba(it,j) lt idim)then potent2(iba(it,j):idim-1,j)=potent2(iba(it,j)-1,j)
  endfor

  ; Calculate E field if imb ge 2
  if (imb ge 2) then begin
     for i=0,idim-1 do begin
         i0=i-1
         if (i0 lt 0) then i0=0
         i1=i+1
         if (i1 ge idim) then i1=idim-1
         for j=0,jdim-2 do begin
             j0=j-1
             if (j0 lt 0) then j0=j0+jdim_1
             j1=j+1
             sqbb=sqrt(bo1(i,j)/BriN(i,j))
             drr=rc_m*(lat1(i1,j)-lat1(i0,j))*!pi/180.
             rco=rc_m*cos(lat1(i,j)*!pi/180.)
             if (imb eq 3) then drr=re_m*(ro1(i1,j)-ro1(i0,j))
             if (imb eq 3) then rco=re_m*ro1(i,j)
             Er=0.
             if (drr gt 0.) then Er=(potent2(i1,j)-potent2(i0,j))/drr ;pot in kV
             Ep=(potent2(i,j1)-potent2(i,j0))/rco/dphi2(j)
             Eio=1.e3*sqrt(Er*Er+Ep*Ep)                   ; E in Volt/meter
             if(imb eq 2)then EorV(i,j)=1.e3*Eio*sqbb     ; E at equator (mV/m)
             if(imb eq 3)then EorV(i,j)=1e6*(Er*ye(i,j)+Ep*xe(i,j))/ro1(i,j) ;Ey
             if(imb eq 4)then EorV(i,j)=-Ep*sqbb/bo1(i,j) ; Vr(km/s) at eqautor
             if(imb eq 5)then EorV(i,j)=Er*sqbb/bo1(i,j)  ; Vp(km/s) at eqautor
             if (EorV(i,j) lt fmin) then EorV(i,j)=fmin
         endfor
         EorV(i,jdim-1)=EorV(i,0)
     endfor
  endif

  ; Draw magnetopause in both left and right panels  
  plot,xe(idim-1,0:jdim-1),ye(idim-1,0:jdim-1),pos=pos1, $
       xrange=[-halfl,halfl],yrange=[-halfl,halfl],xstyle=1,ystyle=1, $
       xtickformat='(a1)',ytickformat='(a1)',thick=2.0,background=255,color=0
  plot,xe(idim-1,0:jdim-1),ye(idim-1,0:jdim-1),pos=pos2, $
       xrange=[-halfl,halfl],yrange=[-halfl,halfl],xstyle=1,ystyle=1, $
       xtickformat='(a1)',ytickformat='(a1)', $
       thick=2.0,background=255,color=0,/noerase

  ; Draw potential contours
  contour,potent1,xe,ye,xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn, $
        levels=plvl(0:ip),c_labels=plab1(0:ip),pos=pos1,c_charsize=1.5, $
        xrange=[-halfl,halfl],yrange=[-halfl,halfl],xstyle=1,ystyle=1, $
        c_colors=c_color,c_linestyle=pline_style(0:ip),xtitle='RE', $
        ytitle=y_title,background=255,color=0,thick=ithick,c_charthick=ithick,$
        title='Potential (kV)',/noerase
        if (ihigh eq 1) then oplot,[xe(20,10),xe(24,36)],[ye(20,10),ye(24,36)],$
                              psym=2,color=225
  if (imb le 1) then contour,potent2,xe,ye,levels=hlvl, $
        xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn, $
        c_labels=plab0,pos=pos2,c_linestyle=pline_style,xtitle='RE', $
        xrange=[-halfl,halfl],yrange=[-halfl,halfl],xstyle=1,ystyle=1, $
        background=255,color=0,thick=ithick,c_charthick=ithick,/noerase

  ; Draw energy contours
  if (imb eq 1) then begin
     contour,energy1,xe,ye,levels=elvl,c_labels=elab,pos=pos2, $
             c_charthick=ithick,xticks=ntick,yticks=ntick,xtickn=tickn, $
             ytickn=tickn,xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
             xstyle=1,ystyle=1,c_linestyle=eline_style,background=255,color=0, $
             thick=ithick,c_charsize=1.5,/noerase
     ; renew iloc
     iloc=value_locate(ro1(*,jloc),radius)
     if (iloc lt idim-1) then begin
        iloc1=iloc+1
        if ((radius-ro1(iloc,jloc)) gt (ro1(iloc1,jloc)-radius)) then iloc=iloc1
     endif
     xx=xe(iloc,jloc)  
     yy=ye(iloc,jloc)  
     oplot,[xx,xx],[yy,yy],psym=2,color=254,symsize=2
  endif

  ; Draw electric field or drift velocity if (imb ge 2)
  EVtitle='Electric Field (mV/m)'
  if (imb eq 4) then EVtitle='Radial Drift'
  if (imb eq 5) then EVtitle='Azimuthal Drift'
  if (imb ge 2) then begin
     contour,EorV,xe,ye,levels=elvl,pos=pos2,xtitle='RE', $
     xrange=[-halfl,halfl],yrange=[-halfl,halfl],xstyle=1,ystyle=1, $
     title=EVtitle,c_colors=colr,/fill,/noerase
     for i=0,mlevel-1 do begin
         xi=x3+i*dx
         xf=xi+dx*1.02
         polyfill,[xi,xf,xf,xi,xi],[yi,yi,yb,yb,yi],color=colr(i),/normal
     endfor
     xyouts,x3,yi-0.04,fmnlab,alignment=0.5,/normal
     xyouts,x4,yi-0.04,fmxlab,alignment=0.5,/normal
     xyouts,0.5*(x3+x4),yi-0.05,E_lab,alignment=0.5,/normal
  endif

  ; label plots         
  UT=abs(t_hour(it))
  hour=string(fix(UT),'(i3.3)')
  imin=round((UT-float(fix(UT)))*60.)
  minute=string(imin,'(i2.2)')
  UT_lab=hour+':'+minute
  if (t_hour(it) lt 0.) then UT_lab='-'+UT_lab
  xyouts,0.1,0.96,filein+'  '+UT_lab+' hh:mm',size=2.0,color=0,/normal
  if (imb eq 1) then begin
     xyouts,-halfl,1.3*halfl,'solid lines: drift paths',color=0 
     xyouts,-halfl,1.1*halfl,'dotted lines: E (keV)',color=0
     xyouts,x3,yb,'mu ='+mulab,color=0,size=1.,/normal
     ekeV=energy1(iloc,jloc)
     pitchA=asin(ym(iloc,jloc))*180./!pi
     mlt1=mlto1(iloc,jloc)*12./!pi
     enlab=string(ekeV,'(f6.1)')
     pitlab=string(pitchA,'(f6.1)')
     rlab=string(ro1(iloc,jloc),'(f3.1)')
     mltlab=string(mlt1,'(f4.1)')
     xyouts,x3,yb-0.045,'*',size=2.,color=254,/normal
     xyouts,x3,yb-0.04,'   E = '+enlab+' keV, PA = '+pitlab+' deg', $
            size=1.1,color=0,/normal
     xyouts,x3,yb-0.08,'        at '+rlab+' RE, '+mltlab+' mlt', $
            size=1.1,color=0,/normal
  endif

  ; Make gif file
  if (iout eq 1) then begin
     gimg = tvrd(0,0)              ; same as tvrd(0,0,960,589)
   ; gimg = tvrd(30,150,450,439)   ; tvrd(x0,y0,dx,dy)
     write_gif,filein+'_poteq.gif',gimg,red(0:255),green(0:255), $
               blue(0:255),/multiple
  endif

erase
goto, newplot

endjob:
if (iout eq 1) then write_gif,filein+'_poteq.gif',/close
if (iout eq 2) then device,/close_file

end
