; 
;                         plot_ece.pro
;
; Program reads the .ece file output from cimi.f90 and makes plots.
; 
; Created on 24 December 2013 by Mei-Ching Fok, Code 673, NASA GSFC.
;------------------------------------------------------------------------------

close,/all

; get colors for color table palette
loadct, 39        ;  rainbow with white
tvlct, red, green, blue, /get
blue(98)=0   &   green(98)=180        ; make cyan to dark green
red(254)=200  &  blue(254)=255   &   green(254)=0     ; make color 254 purple  
tvlct, red, green, blue

!p.color=0
!p.background=255
!p.charsize=2.5
!p.charthick=2.5
!p.thick=5.0

filehead=' '
read,' enter the file head of *.ece file (i.e., 2013d151_e) => ', filehead
read,' number of times in the file => ',ntime
read,' plot: (1)energy, (2)particle => ',iopt
openr,1,filehead+'.ece'
readf,1,je

; determine species
s1=strlen(filehead)-2
spe=strmid(filehead,s1,2)
if (spe eq '_h') then sp_lab=' H+'
if (spe eq 'he') then sp_lab='He+'
if (spe eq '_o') then sp_lab=' O+'
if (spe eq '_e' or spe eq 'e1') then sp_lab=' e-'

; setup arrayas
Lshell1=fltarr(4)      &   Lshell2=Lshell1
energy1=fltarr(je)     &   energy=fltarr(je+1)        
meanE=fltarr(je+2)
Elab=strarr(je+2)    
hour=fltarr(ntime)     &   Dst=hour
drift=fltarr(ntime,4,je+2)   &   driftIn=drift   &  driftOut=drift 
driftXB=drift    ; =driftIn+driftOut, net drift across the boundary
CEWave=drift  &   losscone=drift   &   Bchange=drift   &  esum=drift
Adrift=drift    &   CEWave=drift  &  esumB=drift

; Setup energy bins and mean energies
readf,1,energy1
energy(0)=0.
energy(1:je)=energy1(0:je-1)
for k=0,je do Elab(k)=string(energy(k),'(f9.4)')
Elab(je+1)='infinity'
meanE(*)=1.
if (iopt eq 2) then begin
   for k=1,je-1 do meanE(k)=sqrt(energy(k)*energy(k+1))
   meanE(0)=energy(1)*energy(1)/meanE(1)
   meanE(je)=energy(je)*energy(je)/meanE(je-1)
endif

; read data
header=' '
for n=0,ntime-1 do begin
    readf,1,hour1,Dst1
    hour(n)=hour1
    Dst(n)=Dst1
    for j=1,2 do readf,1,header
    for i=0,3 do begin
        if (i le 2) then begin
           readf,1,L_1,L_2
           Lshell1(i)=L_1
           Lshell2(i)=L_2
        endif 
        for k=0,je+1 do begin
            if (i le 2) then begin
              readf,1,driftIn1,driftOut1,drift1,CEWave1,losscone1,CoulCol, $
                      Bchange1,esum1
            endif else begin
              driftIn1=total(driftIn(n,0:2,k))
              driftOut1=total(driftOut(n,0:2,k))
              drift1=total(drift(n,0:2,k))
              CEWave1=total(CEWave(n,0:2,k))
              losscone1=total(losscone(n,0:2,k))
              Bchange1=total(Bchange(n,0:2,k))
              esum1=total(esum(n,0:2,k))
            endelse
            driftIn(n,i,k)=driftIn1/meanE(k)
            driftOut(n,i,k)=driftOut1/meanE(k)
            driftXB(n,i,k)=driftIn1/meanE(k)+driftOut1/meanE(k)
            drift(n,i,k)=drift1/meanE(k)
            CEWave(n,i,k)=CEWave1/meanE(k)  ;ChargeEx for ion, Wave gain for e-
            losscone(n,i,k)=losscone1/meanE(k)
            Bchange(n,i,k)=Bchange1/meanE(k)
            esum(n,i,k)=esum1/meanE(k)    
        endfor
    endfor
    readf,1,header
endfor
close,1    
Lshell1(3)=Lshell1(0)
Lshell2(3)=Lshell2(2)

; Setup plot size
x_wsize=800
y_wsize=600
DEVICE, DECOMPOSED=0
window,1,xpos=240,ypos=240,xsize=x_wsize,ysize=y_wsize

newplot:

; Make band label
read,' Plot band?  0(0) 1(1), 2(2), sum all bands(3) => ',iband
;Bandlab=string(Lshell1(iband),'(f3.1)')+' < L < '+string(Lshell2(iband),'(f3.1)')
Bandlab=string(round(Lshell1(iband)),'(i1)')+' < L < '+ $
        string(round(Lshell2(iband)),'(i2)')

; combine processes
Adrift=drift-driftXB      ; energy gain/loss due to adiabatic drift
esumB=esum-Bchange        ; total E without gain/loss from changing B

; Sum energy bin
print,'       k       energy(keV)'
for k=0,je do print,k,'    ',Elab(k),' - ',Elab(k+1)
read,'enter k1,k2 => ',k1,k2
ilog=' '
read,'log scale? (n)no, (y)yes => ',ilog
E1=Elab(k1)
E2=Elab(k2+1)
mtitle=E1+' < E(keV) < '+E2+sp_lab
if (k2 ge je) then mtitle='E(keV) > '+E1
sdriftIn=hour  &  sdriftOut=hour  &  sAdrift=hour  &  sCEWave=hour  
sdriftXB=hour  &  sLosscone=hour  &  sBchange=hour
sesumB=hour
for n=0,ntime-1 do begin
    sdriftIn(n)=total(driftIn(n,iband,k1:k2))
    sdriftOut(n)=total(driftOut(n,iband,k1:k2))
    sdriftXB(n)=total(driftXB(n,iband,k1:k2))
    sAdrift(n)=total(Adrift(n,iband,k1:k2))
    sCEWave(n)=total(CEWave(n,iband,k1:k2))
    sLosscone(n)=total(Losscone(n,iband,k1:k2))
    sBchange(n)=total(Bchange(n,iband,k1:k2))
    sesumB(n)=total(esumB(n,iband,k1:k2))
endfor

; determine plot scale
ymax=max([-min(sdriftOut),max(sesumB),max(abs(sAdrift)),max(abs(sLosscone))])
print,' ymax = ',ymax
read,' Do you want to change it? (y)yes, (n)no => ',header
if (header eq 'y') then read,' enter new ymax => ',ymax
ymin=-ymax
if (ilog eq 'y') then begin
   ymin=ymax/1.e3
   print,' ymin = ',ymin
   read,' Do you want to change it? (y)yes, (n)no => ',header
   if (header eq 'y') then read,' enter new ymin => ',ymin
endif

; Adjust plot data if log scale
pdriftIn=sdriftIn   &   pdriftOut=sdriftOut   &   pdriftXB=sdriftXB
pAdrift=sAdrift     &   pCEWave=sCEWave       &  pLosscone=sLosscone       
pBchange=sBchange   &   pesumB=sesumB
if (ilog eq 'y') then begin
   lymax=alog10(ymax)
   lymin=alog10(ymin)
   for n=0,ntime-1 do begin
       pdriftIn(n)=0.
       if (sdriftIn(n) gt ymin) then pdriftIn(n)=alog10(sdriftIn(n))-lymin
       pdriftOut(n)=0.           
       if (sdriftOut(n) lt -ymin) then pdriftOut(n)=lymin-alog10(abs(sdriftOut(n)))
       pdriftXB(n)=0.           
       if (sdriftXB(n) gt ymin) then pdriftXB(n)=alog10(sdriftXB(n))-lymin
       if (sdriftXB(n) lt -ymin) then pdriftXB(n)=lymin-alog10(abs(sdriftXB(n)))
       pAdrift(n)=0.           
       if (sAdrift(n) gt ymin) then pAdrift(n)=alog10(sAdrift(n))-lymin
       if (sAdrift(n) lt -ymin) then pAdrift(n)=lymin-alog10(abs(sAdrift(n)))
       pCEWave(n)=0.           
       if (sCEWave(n) gt ymin) then pCEWave(n)=alog10(sCEWave(n))-lymin
       if (sCEWave(n) lt -ymin) then pCEWave(n)=lymin-alog10(abs(sCEWave(n)))
       pLosscone(n)=0.           
       if (sLosscone(n) lt -ymin) then pLosscone(n)=lymin-alog10(abs(sLosscone(n)))
       pBchange(n)=0.           
       if (sBchange(n) gt ymin) then pBchange(n)=alog10(sBchange(n))-lymin
       if (sBchange(n) lt -ymin) then pBchange(n)=lymin-alog10(abs(sBchange(n)))
       pesumB(n)=0.           
       if (sesumB(n) gt ymin) then pesumB(n)=alog10(sesumB(n))-lymin
   endfor
   ; make plot label
   ytick1=2*(lymax-lymin)
   ylab=strarr(ytick1+1)
   for i=0,ytick1 do begin
       pval=lymin-lymax+i
       aval=abs(pval)
       sign=1.
       if (pval lt 0.) then sign=pval/aval
       yval=(ymin*10.^aval)*sign
       ylab(i)=string(yval,format='(e7.0)')
   endfor
   ; Reset ymin and ymax
   ymax=lymax-lymin  
   ymin=-ymax
endif

; Plot data
hourmn=fix(min(hour)/24.)*24.
hourmx=ceil(max(hour)/24.)*24.
;hourmx=168.
x_minor=6
x_tick=(hourmx-hourmn)/x_minor
for i=1,3 do if (x_tick gt 10) then x_tick=x_tick/2
y_title='Energy (keV)'
if (iopt eq 2) then y_title='# of particles'
if (ilog eq 'n') then begin
   plot,[hourmn,hourmx],[0,0],xrange=[hourmn,hourmx],yrange=[ymin,ymax], $
        xstyle=1,ystyle=1, ytitle=y_title,xtitle='hour',xticks=x_tick, $
        xminor=x_minor,title=mtitle,thick=1,pos=[0.25,0.15,0.95,0.9]
endif else begin
   plot,[hourmn,hourmx],[0,0],xrange=[hourmn,hourmx],yrange=[ymin,ymax], $
        xstyle=1,xminor=x_minor, $
        ystyle=1, ytitle=y_title,xtitle='hour',xticks=x_tick,title=mtitle, $
        yticks=ytick1,ytickname=ylab,thick=1,pos=[0.25,0.15,0.95,0.9]  
endelse
 if (iband eq 3) then oplot,hour,pesumB       ; plot total E only when all bands
 oplot,hour,pLosscone,color=254,thick=7                ; purple  
 oplot,hour,pAdrift,color=98                           ; green  
 oplot,hour,pCEWave,color=57                           ; blue
 oplot,hour,pBchange,color=180                         ; light green
 oplot,hour,pdriftIn,color=225,linestyle=2             ; dashed red
 oplot,hour,pdriftOut,color=225,linestyle=5            ; long-dashed red
 oplot,hour,pdriftXB,color=225,psym=6,nsum=3           ; red
 oplot,hour,pdriftXB,color=225                         ; red

; label the plot
xyouts,0.27,0.83,Bandlab,size=4,/normal
xyouts,0.27,0.17,filehead,size=1.5,charthick=1,/normal
if (ilog eq 'y') then xyouts,0.100,0.51,'!9+!X',font=-1,/normal

read,'Do you want to continue?  (y)yes, (n)no => ',header
if (header eq 'y') then goto,newplot

end
