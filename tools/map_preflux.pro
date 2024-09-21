;-----------------------------------------------------------------------------
; 
;                                 map_preflux.pro
;
; Program reads data from *.preci and plots precipitating energy flux on a
; geographic map.
;
; Created on 16 March 2022 by Mei-Ching Fok (Code 673/GSFC)
;-----------------------------------------------------------------------------

  close,/all

; get colors for color table palette
  loadct,39                               ;  rainbow with white
  tvlct, red, green, blue, /get
  iout=2
  thk=iout*1.5
  !p.background=255
  !p.color=0
  !p.charthick=thk
  !p.charsize=2.-iout*0.5

; Setup color levels
  nlevel=59
  colr=intarr(nlevel)
  colrmax=250.          &    colrmin=10.
  dcolr=(colrmax-colrmin)/(nlevel-1)
  for i=0,nlevel-1 do colr(i)=fix(colrmin+i*dcolr)
  
  re=6.3712e6      ; Earth radius in meter (re_m) defined in cimi.f90

; Read the grids and grid size
  filehead=' '
  read,' enter the file head of precipitation (i.e., 2013d151_e) => ',filehead
  openr,1,filehead+'_N.preci'
  openr,2,filehead+'_S.preci'
  for m=1,2 do readf,m,rc,ir,ip,je,ik
  for m=1,2 do readf,m,elon,ctp,stp  
  print,'elon,ctp,stp ',elon,ctp,stp
  rp0=re*rc        ; rc in m
  altkm=(rp0-re)/1000. 
  althead=string(round(altkm),'(i4)')+'km altitude'
  read,'plot precipitating (1)energy flux,  (2)mean energy => ',iplt
  read,'enter the min Eflux to plot  => ',Efluxmin
  read,'enter number of data set in time => ',ntime
  Mhead=' '
  if (iplt eq 1) then Mhead='_M_pflx_'
  if (iplt eq 2) then Mhead='_M_meanE_'
  ir2=ir*2
  gride=fltarr(je)    &   ebound=fltarr(je+1)
  Elab=strarr(je+1)   &   Elabi=Elab
  Efluxk=fltarr(je+1) &   Pprek=Efluxk
  Eflux=fltarr(ntime,ir,ip,2)    &   glat=Eflux       &  glon=Eflux
  Eflux2D=fltarr(ir2,ip)         &   glat2D=Eflux2D   &  glon2D=Eflux2D
  glon1D=fltarr(ip)              &   glat1D=glon1D    &  Eflux1D=glon1D
  Ppre1D=fltarr(ip)              &   PpreE1D=Ppre1D
  Ppre2D=fltarr(ir2,ip)          &   PpreE2D=Ppre2D
  Ppre=fltarr(ntime,ir,ip,2)     &   PpreE=Ppre
  glon0=fltarr(ntime,ir,2)       &   glat0=glon0      ; glon,glat at mlong=0
  mlata=fltarr(ntime,ir,ip)

; Setup energy boundaries (ebound)
  for m=1,2 do readf,m,gride
  ebound(0)=0.
  for k=1,je-1 do ebound(k)=sqrt(gride(k-1)*gride(k))
  ebound(je)=gride(je-1)*gride(je-1)/ebound(je-1)
  for k=0,je do begin
      Elab(k)=string(ebound(k),'(f6.1)')
      Elabi(k)=string(round(ebound(k)),'(i4.4)')
  endfor
  Ehead=' '
  print,'plot energy(keV) : '
  for k=0,je-1 do print,k,'  '+Elab(k),' - ',Elab(k+1)
  read,' enter lower and upper energy bins => ',ie1,ie2
  Ehead=Elabi(ie1)+'-'+Elabi(ie2+1)+'keV'
  read,' Output data in a fixed geographic grid? 0=no, 1=yes => ',iwrite

; Setup the fixed grid and open a file to write data
  if (iwrite eq 1) then begin
     ira=36  &   ipa=72
     glata=fltarr(ira)   &   Efluxa=fltarr(ipa,ira)      
     Pprea=Efluxa        &   PpreEa=Efluxa             &   meanEa=Efluxa
     glonp=fltarr(ipa+1) &   Efluxp=fltarr(ipa+1,ira)  &   meanEp=Efluxp
     glona=fltarr(ipa)  
     ddeg=360./ipa
     for j=0,ipa-1 do glona(j)=-180.+j*ddeg
     glonp(0:ipa-1)=glona(0:ipa-1)
     glonp(ipa)=glonp(0)+360.
     ddeg=180./ira
     for i=0,ira-1 do glata(i)=-90.+(i*1.+0.5)*ddeg
     openw,3,filehead+'.preci'
     printf,3,fix(ntime),ipa,ira,round(altkm),'    ; ntime,n_glon,n_glat,alt_km'
     printf,3,' glon in degree '
     printf,3,glona
     printf,3,' glat in degree '
     printf,3,glata
  endif

; Read precipitating energy flux [mW/m2] 
  houra=fltarr(ntime)            
  for n=0,ntime-1 do begin
     for m=1,2 do readf,m,hour
     houra(n)=hour
     for m=0,1 do begin
         m1=m+1
         for i=0,ir-1 do begin
             for j=0,ip-1 do begin
                 readf,m1,mlat1,mlt1,mlong1
                 mlata(n,i,j)=abs(mlat1)
                 mlatr=mlat1*!pi/180.           ; magnetic longitude in radian
                 mlongr=mlong1*!pi/180.         ; magnetic longitude in radian
                 ; convert mlatr and mlongr to geographic lat and long
                 stm=cos(mlatr)
                 ctm=sin(mlatr)
                 ctc=ctp*ctm-stp*stm*cos(mlongr)
                 glat1D(j)=90.-acos(ctc)*180./!pi
                 glon1=atan(stp*stm*sin(mlongr),ctm-ctp*ctc)
                 glon1=glon1*180./!pi+elon
                 glon2=glon1
                 if (glon1 gt 180.) then glon2=glon1-360.
                 if (glon1 lt -180.) then glon2=glon1+360.
                 glon1D(j)=glon2
                 if (j eq 0) then glat0(n,i,m)=glat1d(j)
                 if (j eq 0) then glon0(n,i,m)=glon1d(j)
                 ; read and calculate energy flux 
                 readf,m1,Efluxk
                 Eflux1D(j)=total(Efluxk(ie1:ie2))
                 ; read particle flux 
                 readf,m1,Pprek
                 Ppre1D(j)=total(Pprek(ie1:ie2))
                 PpreE1D(j)=0.
                 for k=ie1,ie2 do PpreE1D(j)=PpreE1D(j)+Pprek(k)*gride(k)
             endfor
             ; Sort arrays in accending glon order
             jj=sort(glon1D)
             for j=0,ip-1 do begin
                 glon(n,i,j,m)=glon1D(jj(j))
                 glat(n,i,j,m)=glat1D(jj(j))
                 Eflux(n,i,j,m)=Eflux1D(jj(j))
                 Ppre(n,i,j,m)=Ppre1D(jj(j))
                 PpreE(n,i,j,m)=PpreE1D(jj(j))
             endfor
         endfor  ; end of i loop  
     endfor      ; end of m loop
  endfor         ; end of time loop
  for m=1,2 do close,m

; Setup contour levels
  lvl=fltarr(nlevel)  
  yon=' '
  if (iplt eq 1) then begin
     lvlmax=max(Eflux)
     print,'max energy flux is => ',lvlmax  
     read,'Do you want to change it ? (y)yes, (n)no => ',yon
     if (yon eq 'y') then read,'enter new max energy flux => ',lvlmax
     lvlmin=Efluxmin      
  endif else begin
     read,'enter max mean energy => ',lvlmax  
     lvlmin=0.001
  endelse
  dlvl=(lvlmax-lvlmin)/(nlevel-1)
  for n=0,nlevel-1 do lvl(n)=lvlmin+n*dlvl
  maxV=string(lvlmax,'(f4.0)')
  minV=string(lvlmin,'(f4.2)')

; Setup for plotting Eflux and color bars
  xi=0.07
  xf=0.9
  y0=0.1 
  dy=0.01 
  y2=y0+nlevel*dy
  yave=0.5*(y0+y2)
  if (iplt eq 1) then bar_lab='energy flux (mW/m2)'
  if (iplt eq 2) then bar_lab='mean energy (keV)'

; make time labels
  time_lab=strarr(ntime)
  for i=0,ntime-1 do begin
      time=houra(i)
      ihour=fix(time)
      hours=string(ihour,'(i3)') 
      imin=round((time-float(ihour))*60.)
      minute=string(imin,'(i2.2)')
      time_lab(i)=hours+':'+minute+' '
  endfor

  if (iout eq 2) then begin
     set_plot,'ps'
     device,filename=filehead+Mhead+Ehead+'.ps',/inches,yoffset=3.,$
             xoffset=0.3,xsize=8.,ysize=5.36,/color,bits_per_pixel=24
  endif else begin
     device,decomposed=-1
     window,1,xpos=200,ypos=200,xsize=640,ysize=430
  endelse   

; Plot precipitating and output energy flux
  y_title='geographic latitude (deg)'
  x_title='geographic longitude (deg)'
  for i_img=0,ntime-1 do begin
      hour=houra(i_img)
      xyouts,0.5,y2+0.06,filehead+'  '+Ehead+' '+althead,size=1.5, $
             alignment=0.5,/normal 
      MAP_SET,/CYLINDRICAL,/GRID,/CONTINENTS,/LABEL,title='e- precipitation',$
              pos=[xi,y0,xf,y2],limit=[-90,-180,90,180],/noerase
      Eflux2D(ir:ir2-1,*)=Eflux(i_img,0:ir-1,*,0)
      Ppre2D(ir:ir2-1,*)=Ppre(i_img,0:ir-1,*,0)
      PpreE2D(ir:ir2-1,*)=PpreE(i_img,0:ir-1,*,0)
      glat2D(ir:ir2-1,*)=glat(i_img,0:ir-1,*,0)
      glon2D(ir:ir2-1,*)=glon(i_img,0:ir-1,*,0)
      for i=0,ir-1 do begin
          ir1=ir-1-i
          Eflux2D(ir1,*)=Eflux(i_img,i,*,1)
          Ppre2D(ir1,*)=Ppre(i_img,i,*,1)
          PpreE2D(ir1,*)=PpreE(i_img,i,*,1)
          glat2D(ir1,*)=glat(i_img,i,*,1)
          glon2D(ir1,*)=glon(i_img,i,*,1)
      endfor
      xyouts,xf-0.1,y2-0.03,time_lab(i_img),/normal 
      if (iwrite eq 0 and iplt eq 1) then contour,Eflux2D,glon2D,glat2D, $
                                          levels=lvl,c_colors=colr,/fill,/overplot
      if (iwrite eq 0 and iplt eq 2) then begin
         meanE2D=Eflux2D
         for i=0,ir2-1 do begin
             for j=0,ip-1 do begin
                 meanE2D(i,j)=0.
                 if (Eflux2D(i,j) ge Efluxmin) then $
                    meanE2D(i,j)=PpreE2D(i,j)/Ppre2D(i,j)
             endfor
         endfor
         contour,meanE2D,glon2D,glat2D,levels=lvl,c_colors=colr,/fill,/overplot
      endif

      ; find flux at (glata,glona) and output if iwrite eq 1
      if (iwrite eq 1) then begin
         Efluxi2D=fltarr(ir2,ipa)   &   glati2D=Efluxi2D    &   gloni2D=Efluxi2D
         Pprei2D=Efluxi2D           &   PpreEi2D=Efluxi2D
         ; map fluxes to glona grid
         for i=0,ir2-1 do begin
             gloni2D(i,*)=glona(*)
             loc=value_locate(glon2D(i,*),glona)
             djj=2
             for j=0,ipa-1 do begin
                jj1=loc(j)-djj+1
                if (jj1 lt 0) then jj1=0
                jj2=loc(j)+djj
                if (jj2 ge ip) then jj2=ip-1
                Efluxi2D(i,j)=total(Eflux2D(i,jj1:jj2))/(jj2-jj1+1)
                Pprei2D(i,j)=total(Ppre2D(i,jj1:jj2))/(jj2-jj1+1)
                PpreEi2D(i,j)=total(PpreE2D(i,jj1:jj2))/(jj2-jj1+1)
                glati2D(i,j)=total(glat2D(i,jj1:jj2))/(jj2-jj1+1)
             endfor
         endfor
         ; map fluxes to glata and glona grids
         for j=0,ipa-1 do begin
             loc=value_locate(glati2D(*,j),glata)
             dii=8 
             for i=0,ira-1 do begin
                 ii1=loc(i)-dii+1
                 if (ii1 lt 0) then ii1=0
                 ii2=loc(i)+dii
                 if (ii2 ge ir2) then ii2=ir2-1
                 Efluxa(j,i)=total(Efluxi2D(ii1:ii2,j))/(ii2-ii1+1)
                 Pprea(j,i)=total(Pprei2D(ii1:ii2,j))/(ii2-ii1+1)
                 PpreEa(j,i)=total(PpreEi2D(ii1:ii2,j))/(ii2-ii1+1)
             endfor
         endfor
         Efluxp(0:ipa-1,0:ira-1)=Efluxa(0:ipa-1,0:ira-1)  ; Efluxp for plotting
         ; find meanEa and meanEp
         for i=0,ira-1 do begin
             for j=0,ipa-1 do begin
                 meanEa(j,i)=0.
                 if (Pprea(j,i) gt 0.) then meanEa(j,i)=PpreEa(j,i)/Pprea(j,i)
                 meanEp(j,i)=0.     ; meanEp for plotting
                 if (Efluxa(j,i) ge Efluxmin) then meanEp(j,i)=meanEa(j,i)
             endfor
         endfor
         if (hour eq 9.) then Efluxp(0,*)=Efluxp(1,*)
         if (hour eq 9.) then meanEp(0,*)=meanEp(1,*)
         Efluxp(ipa,0:ira-1)=Efluxp(0,0:ira-1)
         meanEp(ipa,0:ira-1)=meanEp(0,0:ira-1)
         if (iplt eq 1) then contour,Efluxp,glonp,glata,levels=lvl, $
                                     c_colors=colr,/fill,/overplot 
         if (iplt eq 2) then contour,meanEp,glonp,glata,levels=lvl, $
                                     c_colors=colr,/fill,/overplot 
         ; print Eflux and meanE
         printf,3,hour, $
                  '  ; hour, Eflux(glon,glat) in mW/m2, meanE(glon,glat) in keV'
         printf,3,Efluxa
         printf,3,meanEa
      endif

      ; Mark LT-0 direction
      lon0=-hour*15.
      lon0=(lon0 mod 360.)
      if (lon0 lt -180.) then lon0=lon0+360.
      xyouts,lon0+2,0.,'--------------------------',orientation=90, $
             alignment=0.5,color=254
      lon2=lon0-2.
      if (lon2 lt -180.) then lon2=lon0+7.
      xyouts,lon2,0.,'midnight',orientation=90,alignment=0.5,color=254
    
      ; plot color bar
      x1=xf+0.04
      x2=x1+0.02
      lvl2=nlevel-1
      for i=0,lvl2 do begin
          y1=y0+i*dy
          polyfill,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color=colr(i),/normal
      endfor
      xyouts,x2+0.01,y0,minV,/normal
      xyouts,x2+0.01,y1-0.01,maxV,/normal
      xyouts,x1-0.01,yave,bar_lab,orientation=90,alignment=0.5,/normal
     
      ; Make gif file
      if (iout eq 1) then begin
         gimg = tvrd(0,0)
         write_gif,filehead+'_M_pflx_'+Ehead+'.gif',gimg,red(0:255), $
                   green(0:255),blue(0:255),/multiple
      endif
      if (i_img lt ntime-1) then erase

  endfor

  if (iwrite eq 1) then close,3
  if (iout eq 2) then device,/close_file
  if (iout eq 1) then write_gif,filehead+'_M_pflx_'+Ehead+'.gif',/close
  end
