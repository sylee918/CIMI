;-------------------------------------------------------------------------------
;
;                           plot_fls_plDen.pro
;
;  Program and plots plasma density from fls file and output density at a 
;  boundary in the file *.plsbny
;
;  Created on 8 January 2007 by Mei-Ching Fok, Code 673, NASA GSFC
;-------------------------------------------------------------------------------
 
close,/all
;write_gif,/close
!p.charsize=2   
!p.charthick=3.   
!p.background=255

filehead=' '
read,'enter the file head (e.g., 2000a225_h) => ',filehead
openr,1,filehead+'.fls'  
read,'number of data set in time => ',ntime
read,'plot: (1) equatorial plasmasphere density,  (2) ionospheric map => ',iplt
rb=10.
iopt=-1
if (iplt eq 1) then begin
   read,'enter rb => ',rb
   read,'Sun on? (-1)left, (1)right => ',iopt
endif
halfl=ceil(rb)
readf,1,rc,ir,ip,je,ig
density=fltarr(ip+1,ir)  &   xo=density   &   yo=density    
ro=density  &  phi=density   
energy=fltarr(je)    &   sina=fltarr(ig)   &   fy=sina       
varL=fltarr(ir)  &    dnr=varL    
xp=fltarr(ip+1)   &   yp=xp    
rcimi=fltarr(ip)  &   dj=rcimi   &   mphi=rcimi
dummy=fltarr(8)
Ebound=fltarr(je+1)
readf,1,varL
readf,1,mphi
readf,1,energy
readf,1,sina

x_wsize=400
y_wsize=400

; get colors for color table palette
;loadct,39       ;  rainbow with white
 loadct,8        ;  green plasmasphere 
tvlct, red, green, blue, /get

!p.background=255
!p.color=0

; set up levels and colors
nlel=59
level=fltarr(nlel)
c_color=intarr(nlel)
 dmin=1.0          ; density
;dmin=30.          ; density
dmax=3000.
logmax=alog10(dmax)
logmin=alog10(dmin)
dlog=(logmax-logmin)/(nlel-1)
cmax=254.
cmin=1.
dc=(cmax-cmin)/(nlel-1)
for i=0,nlel-1 do begin
    logden=logmin+i*dlog
    level(i)=10.^logden
    c_color(i)=round(cmin+i*dc)
endfor

; Setup plasmapuase (PP)
read,'plasmapause? (0)no (1)cimi (2)byDensity (3)byGradient (4)2+3 => ',ipause
if (ipause ge 2) then read,'plasmapause density(cm-3) => ',densityP
jp=intarr(ip)   &   jx=jp
jmax=8 

; calculate arrays for earth drawing if iplt eq 1
if (iplt eq 1) then begin
  npt=400                             ; no. of points in a circle
  npt2=npt/2
  Night_x=fltarr(npt2+2)     &      day_x=Night_x
  Night_y=fltarr(npt2+2)     &      day_y=Night_y
  Re4_x=fltarr(npt+1)        &      Re4_y=Re4_x
  del = 2.0 * !pi / npt
  ang0=-iopt*!pi/2.
  for i = 0,npt do begin
      i1 = i - npt2
      ang = float(i) * del + ang0  
      cosa = cos(ang)
      sina = sin(ang)
      if (i le npt2) then day_x(i) = cosa
      if (i le npt2) then day_y(i) = sina
      if (i ge npt2) then Night_x(i1) = cosa
      if (i ge npt2) then Night_y(i1) = sina
  endfor
  day_x(npt2+1) = day_x(0)
  day_y(npt2+1) = day_y(0)
  Night_x(npt2+1) = Night_x(0)
  Night_y(npt2+1) = Night_y(0)
  ; 4-Re circle
  Re4_x(0:npt2)=4.*day_x(0:npt2)
  Re4_y(0:npt2)=4.*day_y(0:npt2)
  Re4_x(npt2:npt)=4.*Night_x(0:npt2)
  Re4_y(npt2:npt)=4.*Night_y(0:npt2)
endif

; parameters for plotting
x1=0.22
x2=0.88
y1=0.27
y2=0.93

; arrays and parameter for color bar
ncolor=cmax-cmin+1
temp=bytarr(ncolor,2)
for i=0,ncolor-1 do begin
   temp(i,0)=i+cmin
   temp(i,1)=i+cmin
endfor
bar_len=264
color_bar=congrid(temp,bar_len,20)

; setup fix mlt grid
nmlt=49
dmlt=24./(nmlt-1)
mltf=fltarr(nmlt)   &    rf=mltf    &    nf=mltf
for j=0,nmlt-1 do mltf(j)=j*dmlt     ; mltf in hour

; open file to write r and density at the boundary at fixed mlt
;openw,2,filehead+'.plsbnySM'
;printf,2,ntime,'     ; ntime'
;printf,2,nmlt,'    ; nmlt, mlt(nmlt) in hour, rb(nmlt) in RE, nb(nmlt) in cm-3'
;printf,2,mltf

; setup axes labels
ntick=4
dxy=2.*halfl/ntick
tickn=strarr(ntick+1)
for i=0,ntick do begin
    xy_1=-halfl+i*dxy
    xy_a=abs(xy_1)
    tickn(i)=string(fix(xy_a),'(i2)')
    if (xy_1*iopt lt 0.) then tickn(i)='-'+tickn(i)
endfor

; Setup plot
  if (iplt eq 1) then filetail='_plDen_eq.ps' 
  if (iplt eq 2) then filetail='_plDen_iono.ps' 
  set_plot,'ps'
  device,filename=filehead+filetail,/inches,/color,yoffset=3.5,xoffset=0.3, $
         ysize=6.0,xsize=6.0,bits_per_pixel=24

; Begin the time loop
for nt=1,ntime do begin
   ; read density 
   readf,1,hour1
   print,'nt,hour ',nt,hour1
   t=hour1*3600.
   hour=fix(t/3600.)
   minute=round((t-3600.*hour)/60.)
   tsec=hour*3600.+minute*60.
   for j=0,ir-1 do begin
       for i=0,ip-1 do begin
           readf,1,lat1,xmlt,xlatS,xmltS,ro1,xmlto,BriN,BriS,bo,ib,den1, $
                   dummy,rcim1
           for k=0,je-1 do readf,1,fy
           density(i,j)=den1*1.e-6      ; plasmasphere density in cm^-3
           phi(i,j)=xmlto*!pi/12.
           rcimi(i)=rcim1
           if (iplt eq 1) then begin
              xo(i,j)=-iopt*ro1*cos(phi(i,j))
              yo(i,j)=-iopt*ro1*sin(phi(i,j))
              ro(i,j)=ro1
           endif else begin
              xo(i,j)=xmlt
              yo(i,j)=lat1
           endelse
       endfor                         
   endfor
   ; Find plasmapause if ipause eq 1 (output from cimi.f90)
   if (ipause eq 1) then begin
      for i=0,ip-1 do begin
          jc=value_locate(ro(i,*),rcimi(i))
          xp(i)=-iopt*rcimi(i)*cos(phi(i,jc))
          yp(i)=-iopt*rcimi(i)*sin(phi(i,jc))
      endfor
   endif
   ; Find plasmapause if ipause ge 2
   if (ipause ge 2) then begin    ; define PP by densityP
      ; find PP by densityP
      for i=0,ip-1 do begin
          for j=2,ir-2 do begin
              if (densityP gt density(i,j)) then begin
                 jp(i)=j-1
                 goto,nextMLT
              endif
          endfor
          nextMLT:
      endfor
      jpp=jp
      ; find PP by max gradient 
      if (ipause ge 3) then begin  
         for i=0,ip-1 do begin
             im1=i-1
             if (im1 lt 0) then im1=im1+ip
             ip1=i+1
             if (ip1 ge ip) then ip1=ip1-ip
             dnr(*)=0.   
             for j=3,jp(i) do begin
                 dphi=abs(phi(ip1,j)-phi(im1,j))
                 if (dphi gt !pi) then dphi=2.*!pi-dphi
                 dn1=(density(i,j+1)-density(i,j-1))/(ro(i,j+1)-ro(i,j-1)) 
                 dn2=(density(ip1,j)-density(im1,j))/dphi/ro(i,j)
                 dnr(j)=sqrt(dn1*dn1+dn2*dn2)/density(i,j)  ; fractional change
             endfor
             dnmx=max(dnr,jmx)
             jx(i)=jmx
             jpp(i)=jx(i)                          ; define PP by max gradient
         endfor
      endif
      ; Smooth PP when ipause eq 4
      if (ipause eq 4) then begin
          ; smooth plasmapause in eastward direction
          i0=-1
          for i=1,ip-1 do begin
              jp1=abs(jpp(i)-jpp(i-1))
              if (i0 eq -1 and jp1 gt jmax) then begin
                 i0=i-1
                 goto,nexti1
              endif
              if (i0 ne -1 and jp1 gt jmax) then begin  ; redo jpp inside i0:i
                 for ii=i0+1,i-1 do begin
                     jnew=interpol([jpp(i0),jpp(i)],[i0,i],ii)
                     jpp(ii)=round(jnew)
                 endfor
                 i0=-1
              endif
              nexti1:  
          endfor
          ; test whether the plasmapause is smmoth after 1st smoothing
          for i=0,ip-1 do begin
              i1=i+1
              if (i1 ge ip) then i1=i1-ip
              dj(i)=abs(jpp(i)-jpp(i1))
              djmax1=max(dj)
          endfor
          ; smooth plasmapause in westward direction
          if (djmax1 gt jmax) then begin
             jpp0=jpp                           
             jpp=jx                         
             i0=-1
             for i=ip-2,0,-1 do begin
                 jp1=abs(jpp(i)-jpp(i+1))
                 if (i0 eq -1 and jp1 gt jmax) then begin
                    i0=i+1
                    goto,nexti2
                 endif
                 if (i0 ne -1 and jp1 gt jmax) then begin  ; redo jpp 
                    for ii=i+1,i0-1 do begin
                        jnew=interpol([jpp(i),jpp(i0)],[i,i0],ii)
                        jpp(ii)=round(jnew)
                    endfor
                    i0=-1
                 endif
                 nexti2:
             endfor
             ; test whether the plasmapause is smmoth after 2nd smoothing
             for i=0,ip-1 do begin
                 i1=i+1
                 if (i1 ge ip) then i1=i1-ip
                 dj(i)=abs(jpp(i)-jpp(i1))
             endfor
             if (max(dj) gt djmax1) then jpp=jpp0
          endif
      endif
      ; Calculate location of PP
      for i=0,ip-1 do begin
          jpp1=jpp(i)
          if (ipause eq 2) then jpp1=jpp1+1
          xp(i)=0.5*(xo(i,jpp(i))+xo(i,jpp1))
          yp(i)=0.5*(yo(i,jpp(i))+yo(i,jpp1))
      endfor
   endif
   
   ; periodic boundary condition
   density(ip,*)=density(0,*)
   xo(ip,*)=xo(0,*)
   if (iplt eq 2) then xo(ip,*)=xo(0,*)+24.
   yo(ip,*)=yo(0,*)
   xp(ip)=xp(0)
   yp(ip)=yp(0)
   
   ; plot density and time label
   polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=0,/normal
   if (iplt eq 1) then begin
      contour,density,xo,yo,position=[x1,y1,x2,y2],title=filehead, $
           xticks=ntick,yticks=ntick,xtickn=tickn,ytickn=tickn,xtitle='X (RE)',$
           ytitle='Y (RE)',xrange=[-halfl,halfl],yrange=[-halfl,halfl], $
           xstyle=1,ystyle=1,levels=level,c_colors=c_color,/fill,/noerase
      ; plot plasmapause and outer boundary
      contour,density,xo,yo,levels=level,c_colors=c_color,/overplot
      polyfill,Night_x,Night_y,color=1
      polyfill,day_x,day_y,color=255
   endif else begin
      contour,density,xo,yo,position=[x1,y1,x2,y2],title=filehead, $
           ytitle='mlat (deg)',xrange=[0,24],yrange=[20,70],xstyle=1, $
           xtitle='mlt (hour)',ystyle=1,levels=level,c_colors=c_color,/fill, $
           /noerase
      ; plot plasmapause and outer boundary
      contour,density,xo,yo,levels=level,c_colors=c_color,/overplot
   endelse
   hour_minute=string(hour,'(i3)')+':'+string(minute,'(i2.2)')
   xyouts,x1,0.87,hour_minute+' UT',size=1.5,color=255,/normal

   ; Draw plasmapause by gradient or by cimi output if ipause ge 1
   if (ipause ge 1 and max(density) gt 0.) then oplot,xp,yp,color=255,thick=3
   
   ; Add color bar and label
   dxx=(x2-x1)/nlel
   yi=y1-0.22
   yf=yi+0.05
   for i=0,nlel-1 do begin
       xi=x1+i*dxx
       polyfill,[xi,x2,x2,xi],[yi,yi,yf,yf],color=c_color(i),/normal
   endfor
   d10=(x2-x1)/alog10(dmax/dmin)
   for i=0,3 do begin
       if (i eq 0) then dlab=string(10^i,'(i1)')
       if (i eq 1) then dlab=string(10^i,'(i2)')
       if (i eq 2) then dlab=string(10^i,'(i3)')
       if (i eq 3) then dlab=string(10^i,'(i4)')
       xyouts,x1+i*d10,yi,'l',alignment=0.5,/normal
       xyouts,x1+i*d10,yi-0.04,dlab,alignment=0.5,/normal
   endfor
   xyouts,0.5*(x1+x2),yi-0.08,'density (cm^-3)',alignment=0.5,/normal
          
   ; get gif file
;  gimg=tvrd(0,0,x_wsize,y_wsize)
;  if (iplt eq 1) then filetail='_plDen_eq.gif'
;  if (iplt eq 2) then filetail='_plDen_iono.gif'
;  if (max(density) gt 100.) then write_gif,filehead+filetail,gimg,red(0:255), $
;            green(0:255),blue(0:255),/multiple   
   erase
endfor  
close,1
;close,2

;write_gif,filehead+'.gif',/close
device,/close_file

end
