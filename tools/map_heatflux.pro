;-----------------------------------------------------------------------------
;
;                                   map_heatflux.pro
;
; Program reads heat flux from *.fls and plots them on a geographic map.
;
; Created on 16 June 2023 by Mei-Ching Fok (Code 673/GSFC)
;-----------------------------------------------------------------------------

  close,/all

; get colors for color table palette
  loadct,39                               ;  rainbow with white
  tvlct, red, green, blue, /get
  !p.background=255
  !p.color=0
  !p.charthick=3.  

; Setup color levels
  nlevel=59
  colr=intarr(nlevel)
  colrmax=250.          &    colrmin=10.
  dcolr=(colrmax-colrmin)/(nlevel-1)
  for i=0,nlevel-1 do colr(i)=fix(colrmin+i*dcolr)
  
  re=6.3712e6      ; Earth radius in meter (re_m) defined in cimi.f90

; Read rc and ,elon,ctp,stp for conversion between MAG and GEO
  filehead=' '
  read,' enter the file head of the storm (i.e., 2013r151) => ',filehead
  openr,1,filehead+'_h_N.preci'
  readf,1,rc
  readf,1,elon,ctp,stp  
  close,1
  rp0=re*rc        ; rc in m
  altkm=(rp0-re)/1000. 
  althead=string(fix(altkm),'(i4)')+'km altitude'

; Read grid sizes and species label
  read,'number of species?  1=H+   2=H+,O+   3=H+,O+,He+ => ',isp
  for is=2,isp+1 do begin
      if (is eq 2 ) then openr,is,filehead+'_h.fls'
      if (is eq 3 ) then openr,is,filehead+'_o.fls'
      if (is eq 4 ) then openr,is,filehead+'he.fls'
  endfor
  species='H+'
  splab='_H'
  if (isp eq 2) then species='H+,O+'
  if (isp eq 2) then splab='_H-O'
  if (isp eq 3) then species='H+,O+,He+'
  if (isp eq 3) then splab='_H-O-He'
  for is=2,isp+1 do readf,is,rc,ir,ip,je,ig
  varL=fltarr(ir)    &   mphi=fltarr(ip)
  energy=fltarr(je)  &   sina=fltarr(ig)   &   fy=sina
  for is=2,isp+1 do readf,is,varL
  for is=2,isp+1 do readf,is,mphi
  for is=2,isp+1 do readf,is,energy
  for is=2,isp+1 do readf,is,sina

; Choose plasmaspheric electron or ion heating
  read,' plot? (1)PlSp-e heatflux, (2)PlSp-i heatflux => ',iopt
  if (iopt eq 1) then begin 
     fmiddle='_HFPe_M'
     PlSp_lab='PlSp e- heating (eV/cm^2/s)'
  endif else begin
     fmiddle='_HFPi_M'
     PlSp_lab='PlSp ion heating (eV/cm^2/s)'
  endelse
  map_title=PlSp_lab+' from Ring Current '+species
  read,'scale?  0 = linear scale, or 1 = log scale => ',ilog

; Read heat fluxes and map to geo coor.
  read,'enter number of data set in time => ',ntime
  houra=fltarr(ntime)        
  glata=fltarr(ntime,ir,ip,2)    &  glona=glata      &  hfluxa=glata
  glat2=fltarr(ip,2)             &  glon2=glat2      &  hflux2=glat2
  glat2D=fltarr(ir,ip)           &  glon2D=glat2D    &  hflux2D=glat2D
  for n=0,ntime-1 do begin
      for is=2,isp+1 do readf,is,hour
      houra(n)=hour
      print,n,hour
      for i=0,ir-1 do begin
          for j=0,ip-1 do begin
              HRP=0.     
              for is=2,isp+1 do begin
                  readf,is,lat1,mlt1,lat2,mlt2,ro1,mlto1,bri1,bri2,bo1,iba1, $
                        plsDen1,ompe,CHpower,HIpower,denWP1,Tpar1,Tper1,HRPee, $
                        HRPii,rpp,Lstar,tvolume  ; FluxTube volume/magnetic flux
                  if (iopt eq 1) then HRP=HRP+HRPee   ; heating rate (keV/m^3/s)
                  if (iopt eq 2) then HRP=HRP+HRPii   ;
              endfor
              VHRP=0.5*tvolume*HRP*0.1 ; 0.5 for one ionospher, 0.1 for eV/cm2/s
              for m=0,1 do begin       ; 2 hemispheres
                  if (m eq 0) then hflux2(j,m)=VHRP*bri1
                  if (m eq 1) then hflux2(j,m)=VHRP*bri2
                  ; find glat and glon
                  if (m eq 0) then mlatr=lat1*!pi/180.
                  if (m eq 1) then mlatr=lat2*!pi/180.
                  if (m eq 0) then mlongr=mphi(j)
                  if (m eq 1) then mlongr=mlongr+(mlt2-mlt1)*!pi/12.
                  stm=cos(mlatr)
                  ctm=sin(mlatr)
                  ctc=ctp*ctm-stp*stm*cos(mlongr)
                  glat2(j,m)=90.-acos(ctc)*180./!pi
                  glon1=atan(stp*stm*sin(mlongr),ctm-ctp*ctc)
                  glon1=glon1*180./!pi+elon
                  glon2(j,m)=glon1
                  if (glon1 gt 180.) then glon2(j,m)=glon1-360.
                  if (glon1 lt -180.) then glon2(j,m)=glon1+360.
              endfor
              for k=0,je-1 do begin
                  for is=2,isp+1 do begin
                      readf,is,fy
                  endfor
              endfor
          endfor
          ; Sort arrays in accending glon order
          for m=0,1 do begin
              jj=sort(glon2(*,m))
              for j=0,ip-1 do begin
                  glata(n,i,j,m)=glat2(jj(j),m)
                  glona(n,i,j,m)=glon2(jj(j),m)
                  hfluxa(n,i,j,m)=hflux2(jj(j),m)
              endfor
          endfor
      endfor
  endfor
  for is=2,isp+1 do close,is
  if (ilog eq 1) then begin
     if (min(hfluxa) lt 1.e-30) then hfluxa(where(hfluxa lt 1.e-30))=1.e-30
     hfluxa=alog10(hfluxa)
  endif

; Setup contour levels
  lvl=fltarr(nlevel)  
  lvlmax=max(hfluxa)
  yon=' '
  print,'max heat flux is => ',lvlmax  
  read,'Do you want to change it ? (y)yes, (n)no => ',yon
  if (yon eq 'y') then read,'enter new max heat flux => ',lvlmax
  lvlmin=0.      
  if (ilog eq 1) then lvlmin=lvlmax-3.
  dlvl=(lvlmax-lvlmin)/(nlevel-1)
  for n=0,nlevel-1 do lvl(n)=lvlmin+n*dlvl
  maxV=string(lvlmax,'(f4.1)')
  minV=string(lvlmin,'(f4.1)')

; Setup for plotting heat flux and color bars
  xi=0.07
  xf=0.9
  y0=0.1 
  dy=0.01 
  y2=y0+nlevel*dy
  yave=0.5*(y0+y2)
  bar_lab='heat flux (eV/cm^2/s)'

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

  set_plot,'ps'
  device,filename=filehead+splab+fmiddle+'.ps',/inches,yoffset=3.,$
             xoffset=0.3,xsize=8.,ysize=5.36,/color,bits_per_pixel=24

; Plot heat flux
  y_title='geographic latitude (deg)'
  x_title='geographic longitude (deg)'
  for i_img=0,ntime-1 do begin
      hour=houra(i_img)
      xyouts,0.5,y2+0.06,filehead+'  '+althead,alignment=0.5,/normal
      MAP_SET,/CYLINDRICAL,/GRID,/CONTINENTS,/LABEL,title=map_title, $
              pos=[xi,y0,xf,y2],limit=[-90,-180,90,180],/noerase
      for m=0,1 do begin
          hflux2D(*,*)=hfluxa(i_img,*,*,m)
          glat2D(*,*)=glata(i_img,*,*,m)
          glon2D(*,*)=glona(i_img,*,*,m)
          contour,hflux2D,glon2D,glat2D, $
                  levels=lvl,c_colors=colr,/fill,/overplot
      endfor
      xyouts,xf-0.1,y2-0.03,time_lab(i_img),/normal 

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
     
      if (i_img lt ntime-1) then erase
  endfor

  device,/close_file
  end
