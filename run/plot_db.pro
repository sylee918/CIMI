 close,/all

 loadct,39       ;  rainbow with white
 tvlct, red, green, blue, /get
 !p.color=0
 !p.background=255
 !p.charsize=2
 DEVICE, DECOMPOSED=0

 header=' '
 fhead=' ' 
 read,' enter the file head of .db (i.e., 2013d151) => ', fhead
 read,' enter the .db file to be plotted: 0=*.dbe, 1=*.dbi, 2=*.db => ',isp
 read,'Enter number of time => ',np    
 itime=1     ; 1 = hour, 2 = day
 if (isp eq 0) then openr,1,fhead+'.dbe'   
 if (isp eq 1) then openr,1,fhead+'.dbi'   
 if (isp eq 2) then openr,1,fhead+'.db'   
 sps='(i)'
 if (isp eq 0) then sps='(e)'
 if (isp eq 2) then sps='(i+e)'
 readf,1,header
 hour=fltarr(np)   &   Dst=hour  &  DstRC=hour  &  totE=hour  &  totEB=hour
 day=hour
 for m=0,np-1 do begin
     readf,1,hour1,Kp,AE,AL,Dst1,DstRC1,rcsum,totE1,totEB1
     hour(m)=hour1
     day(m)=hour1/24.
     Dst(m)=Dst1
     DstRC(m)=DstRC1
     totE(m)=totE1
     totEB(m)=totEB1
 endfor
 close,1

 DPS=-4e-30
 if (isp eq 0) then factor=3.0    ; bigger factor wider curve range, vice versa    
 if (isp eq 1) then factor=0.9    ; bigger factor wider curve range, vice versa    
 if (isp eq 2) then factor=1.2      
 if (isp eq 0) then Dst0=5.0      ; bigger Dst0 shifts curve up, vice versa
 if (isp eq 1) then Dst0=30.      ; bigger Dst0 shifts curve up, vice versa
 if (isp eq 2) then Dst0=30.
 DstEB=totEB    
 DstEB=totEB*DPS*factor+Dst0
 Dstmin=-450
 Dstmax=50.
 dn=1 
 DstEB0=DstEB
 for n=dn,np-1 do DstEB(n)=DstEB0(n-dn)
 print,' dn => ',dn

 tmin=0.
 daymax=ceil(max(day))
 hourmax=daymax*24
 if (itime eq 1) then begin
    time=hour
    time_lab='hour'
    tmax=hourmax
 endif else begin
    time=day 
    time_lab='day'
    tmax=daymax
 endelse
 plot,time,DstEB,title=fhead,xtitle=time_lab,ytitle='Dst', $
      xrange=[tmin,tmax],xstyle=1,pos=[0.1,0.12,0.88,0.92],ystyle=8, $
      yrange=[Dstmin,Dstmax],Ytick_get=Eticks
 oplot,time,Dst,linestyle=1     ; dotted line
 oplot,time,DstRC,linestyle=2   ; dashed line
 Dstmin=min(Eticks)    ; redefine Dstmin
 Dstmax=max(Eticks)    ; redefine Dstmax
 print,'Dstmin,Dstmax ',Dstmin,Dstmax
 Emin=(Dstmax-Dst0)/DPS/factor
 Emax=(Dstmin-Dst0)/DPS/factor
 print,'Emin,Emax ',Emin,Emax
 Axis,YAxis=1,yrange=[Emax,Emin],ystyle=1,ytitle='Ring Current Energy (keV)'

 Dst0s=string(Dst0,'(f3.0)')
 factors=string(factor,'(f4.2)')
 xyouts,0.70,0.85,'Dst0 = '+Dst0s,size=1.5,/normal
 xyouts,0.70,0.80,'factor = '+factors,size=1.5,/normal
 xyouts,0.5,0.85,'..... Dst',size=1.5,/normal
 xyouts,0.5,0.80,'-- Dst*',size=1.5,/normal
 xyouts,0.5,0.75,'___ Erc'+sps,size=1.5,/normal

 end
