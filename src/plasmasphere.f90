!-------------------------------------------------------------------------------
  subroutine plasmasphere_new(iyear,iday,ihour,rc,delt)
!-------------------------------------------------------------------------------
! A plasmasphere model calculates number of plasmaspheric ions per unit magnetic
! flux, Nion. This code is similar to Dan Ober's model (pbo_2.f) but this model
! takes user provided grid and does not do interpolation to its own grid.
! This code also considers the north-south asymmetry at the ionosphere.
!
! INPUT
! iyear, iday: year and day of the year
! rc: radial distance of ionosphere in Earth radius (RE)
! ir,ip: number of grids in mlat (xlati) and mlt (phi), respectively.
! xlati(ir,ip): mlat in radian at ionosphere in northern hemisphere
! phi(ip): mlt in radian (phi=0-->midnight) at ionosphere in northern hemisphere
! xlatiS(ir,ip), phiS(ip): mlat and mlt of the south conjugate points of 
!                         (xlati,phi)
! iba(ip): the last xlati index inside rb.
! ro(ir,ip): radial distance at the equator of the flux tube rooted at northern
!            ionosphere at xlati and phi.
! volume(ir,ip): flux tube volume per unit magnetic flux
! BriN(ir,ip): radial component of magnetic field at rc at (xlati,phi)
! vlE,vpE(ir,ip): plasmasphere drift velocity in latitude and longitude
! delt: simulation time in second of this call
!
! INPUT/OUTPUT
! Nion: number of ions per unit magnetic flux.    
!
! Created on 24 July 2007 by Mei-Ching Fok. Code 673, NASA GSFC.
!-------------------------------------------------------------------------------

  use constants, only: pi,re_m,i_one
  use cPlasmasphere_new, only: ir=>nlp,ip=>npp,xlati=>xlatp,phi=>phip, &
                         xlatiS=>xlatpS,phiS=>phipS,iba=>ibp,ro=>rp, &
                         volume=>volp,BriN=>Brip, &
                         varL=>varLp,vlE=>vlEp,vpE=>vpEp,ksai=>ksaip, &
                         Nion
  implicit none

  integer iyear,iday,ihour,imin,isec,nstep,nrun,i,j,n
  real dphi,dvarL,delt,rc,coslat,phi1,xsm,ysm,zsm,xgsm,ygsm,zgsm
  real xgse,ygse,zgse,dtmax,dt,NTro,denSat,Nsat(ir,ip)
  real dt1,cl(ir,ip),cp(ir,ip),vgseX,vgseY,vgseZ,cor,deltc
  character*1 Ndon(ir,ip),Sdon(ir,ip)

! Constants
  imin=0
  isec=0
  vgseX=-400.
  vgseY=0.
  vgseZ=0.

! Setup grid sizes
  dphi=2.*pi/ip
  dvarL=varL(2)-varL(1)

! Determine whether include corotation
  cor=2.*pi/86400.            ! corotation speed in rad/s
  deltc=delt*cor             ! corotation angle in delt
  if (deltc.lt.dphi) cor=0.  ! ignore corotation if deltc is smaller dphi

! find Ndon and Sdon
  call recalc_08(iyear,iday,ihour,imin,isec,vgseX,vgseY,vgseZ)
  do j=1,ip
     do i=1,iba(j)
        ! find Ndon
        coslat=cos(xlati(i,j))
        phi1=phi(j)+pi
        xsm=rc*coslat*cos(phi1)
        ysm=rc*coslat*sin(phi1)
        zsm=rc*sin(xlati(i,j))
        call smgsw_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,i_one)
        call gswgse_08(xgsm,ygsm,zgsm,xgse,ygse,zgse,i_one)
        if (xgse.gt.0.) Ndon(i,j)='d'
        if (xgse.le.0.) Ndon(i,j)='n'
        ! find Sdon
        coslat=cos(xlatiS(i,j))
        phi1=phiS(j)+pi
        xsm=rc*coslat*cos(phi1)
        ysm=rc*coslat*sin(phi1)
        zsm=rc*sin(xlatiS(i,j))
        call smgsw_08(xsm,ysm,zsm,xgsm,ygsm,zgsm,i_one)
        call gswgse_08(xgsm,ygsm,zgsm,xgse,ygse,zgse,i_one)
        if (xgse.gt.0.) Sdon(i,j)='d'
        if (xgse.le.0.) Sdon(i,j)='n'
     enddo
  enddo

! Setup the trough by trough(6.6)*dipolevolume(6.6) as in pbo_2.f
  NTro=9.40e20            ! number of particle per unit magnetic flux at trough 
  NTro=NTro/3.            ! MCF on Dec 04, 2019

! Calculate saturated Nion  
     do j=1,ip
        do i=1,iba(j)
           denSat=1.e6*10.0**(-0.3145*ro(i,j)+3.9043)   ! saturated density
           Nsat(i,j)=denSat*volume(i,j)
        enddo
     enddo

! determine time step 
  dtmax=30.                    ! maximum time step in second
  nstep=ceiling(delt/dtmax)
  dt=delt/nstep                ! time step in second
  
! calculate convection velocity and Courant numbers
  call Vconvect(ir,ip,cor,dvarL,dphi,dt,dt1,nrun,vlE,vpE,cl,cp)

! time loop to update Nion    
  do n=1,nstep    
     call trough(ir,ip,iba,NTro,Nion)
     call drift_pl(ir,ip,iba,nrun,cl,cp,ksai,Nsat,Ntro,Nion)
     call refill_loss(ir,ip,iba,dt,Nsat,BriN,Ndon,Sdon,Nion)
  enddo

  end subroutine plasmasphere_new


!-------------------------------------------------------------------------------
  subroutine Vconvect(ir,ip,cor,dvarL,dphi,dt,dt1,nrun,vlE,vpE,cl,cp)
!-------------------------------------------------------------------------------
! Routine calculates the convection velocities cl and cp
! Input: ir,ip,dvarL,dphi,dt,vlE,vpE
! Output: dt1,nrun,cl,cp

  use constants, only: pi
  implicit none

  integer ir,ip,i,j,nrun
  real dphi,vlE(ir,ip),vpE(ir,ip),cl(ir,ip),cp(ir,ip),dt,dt1,cmax,cmx,dvarL,cor

! Find nrun, new dt (dt1), cl and cp
  cmax=0.
  do i=1,ir
     do j=1,ip
        cl(i,j)=vlE(i,j)*dt/dvarL
        cp(i,j)=(vpE(i,j)+cor)*dt/dphi
        cmx=max(abs(cl(i,j)),abs(cp(i,j))) 
        cmax=max(cmx,cmax) 
     enddo
  enddo
  nrun=ifix(cmax/0.25)+1     ! nrun to limit the Courant number
  dt1=dt/nrun
  if (nrun.gt.1) then       ! reduce cl and cp if nrun > 1
     cl(1:ir,1:ip)=cl(1:ir,1:ip)/nrun
     cp(1:ir,1:ip)=cp(1:ir,1:ip)/nrun
  endif

  end subroutine Vconvect


!-------------------------------------------------------------------------------
  subroutine drift_pl(ir,ip,iba,nrun,cl,cp,ksai,Nsat,Ntro,Nion)
!-------------------------------------------------------------------------------
! Routine updates Nion due to drift (convection+corotation)
! Input: ir,ip,iba,nrun,vl,cl,cp,dt1,dlati,Nsat,Ntro
! Input/Output: Nion

  implicit none

  integer ir,ip,iba(ip),i,j,j_1,n,nrun
  real Nsat(ir,ip),Nion(ir,ip),f2(ir,ip),fb1(ip)
  real cl(ir,ip),cp(ir,ip),fal(ir,ip),fap(ir,ip),ksai(ir,ip),Ntro

! calculate f2
  Nion(1,1:ip)=Nsat(1,1:ip)    ! saturation density for first cell
  do i=1,ir
     do j=1,ip
        f2(i,j)=ksai(i,j)*Nion(i,j)
     enddo
  enddo

! flux at outer boundary
  do j=1,ip
     fb1(j)=NTro*ksai(ir,j)
  enddo

! Update f2 by drift
  do n=1,nrun
     call inter_flux(ir,ip,iba,cl,cp,fb1,f2,fal,fap)     
     do j=1,ip
        j_1=j-1
        if (j_1.lt.1) j_1=j_1+ip
        do i=2,iba(j)              
           f2(i,j)=f2(i,j)+cl(i-1,j)*fal(i-1,j)-cl(i,j)*fal(i,j) &
                          +cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j)  
           if (f2(i,j).lt.0.) then
              write(*,*) ' In plasmasphere, Error: f2(i,j).lt.0. '
              stop
           endif
        enddo
     enddo
  enddo 

! get tbe new Nion
  do j=1,ip
     Nion(1:ir,j)=f2(1:ir,j)/ksai(1:ir,j)
  enddo

  end subroutine drift_pl


!-------------------------------------------------------------------------------
  subroutine refill_loss(ir,ip,iba,dt,Nsat,BriN,Ndon,Sdon,Nion)
!-------------------------------------------------------------------------------
! Routine updates Nion due to refilling and loss.
! On the nightside, dN/dt = -N/(Bi*tau) and N = No*exp(-dt/tau)
! on the dayside, dN/dt = Fmax*(Nsat-N)/Nsat and 
!                 N = Nsat - (Nsat - No)*exp(-dt*Fmax/Nsat/Bi)
!
! Input: ir,ip,iba,dt,velume,Nsat,BriN,Ndon,Sdon
! Input/Output: Nion

  implicit none

  integer ir,ip,iba(ip),i,j
  real dt,Nsat(ir,ip),BriN(ir,ip),Nion(ir,ip)
  real tau,Fmax,factorN,factorD,tFNB
  character Ndon(ir,ip)*1,Sdon(ir,ip)*1

  Fmax=2.e12          ! limiting refilling flux in particles/m**2/sec
  tau=10.*86400.      ! nightside downward diffusion lifetime in second
  factorN=exp(-dt/tau)

  do j=1,ip
     do i=1,iba(j)
        ! at Northern ionosphere
        if (Ndon(i,j).eq.'d') then          ! dayside refilling
           tFNB=-dt*Fmax/Nsat(i,j)/BriN(i,j)
           factorD=exp(tFNB)
           Nion(i,j)=Nsat(i,j)-(Nsat(i,j)-Nion(i,j))*factorD
        else
           Nion(i,j)=Nion(i,j)*factorN      ! nightside diffusion
        endif

        ! at Southern ionosphere
        if (Sdon(i,j).eq.'d') then          ! dayside refilling
           tFNB=-dt*Fmax/Nsat(i,j)/BriN(i,j)
           factorD=exp(tFNB)
           Nion(i,j)=Nsat(i,j)-(Nsat(i,j)-Nion(i,j))*factorD
        else
           Nion(i,j)=Nion(i,j)*factorN      ! nightside diffusion
        endif

        ! limit Nion to Nsat
        if (Nion(i,j).gt.Nsat(i,j)) Nion(i,j)=Nsat(i,j)
     enddo
  enddo

  end subroutine refill_loss


!*******************************************************************************
  subroutine trough(ir,ip,iba,NTro,Nion)
!*******************************************************************************
! Routine fills up Nion at i>iba with NTro and makes sure the plasmasphere 
! density is not lower than NTro
! density.
!
! Input: ir,ip,NTro
! Input/Output: Nion

  implicit none

  integer ir,ip,i,j,iba(ip)
  real Nion(ir,ip),NTro

  do j=1,ip
     do i=1,iba(j)
        if (Nion(i,j).lt.NTro) Nion(i,j)=NTro
     enddo
     do i=iba(j)+1,ir
        Nion(i,j)=NTro
     enddo
  enddo

  end subroutine trough


!*******************************************************************************
      subroutine inter_flux(ir,ip,iba,cl,cp,fb1,f2,fal,fap)
!*******************************************************************************
!  Routine calculates the inter-flux, fal(i+0.5,j) and fap(i,j+0.5), using
!  2nd order flux limited scheme with super-bee flux limiter method.
!
!  Input: ir,ip,iba,cl,cp,fb1,f2
!  Output: fal,fap

   implicit none

   integer ir,ip,i,j,j_1,j1,j2
   integer iba(ip),ibm
   real cl(ir,ip),cp(ir,ip),f2(ir,ip),fal(ir,ip),fap(ir,ip),fb1(ip), &
        fwbc(0:ir+2,ip),xsign,fup,flw,x,r,xlimiter,corr

      fwbc(1:ir,1:ip)=f2(1:ir,1:ip)     ! fwbc is f2 with boundary condition

! Set up boundary condition
       fwbc(0,1:ip)=f2(1,1:ip)
       do j=1,ip
          fwbc(ir+1:ir+2,j)=fb1(j)
       enddo

! find fa*
      do j=1,ip
         j_1=j-1
         j1=j+1
         j2=j+2
         if (j_1.lt.1) j_1=j_1+ip
         if (j1.gt.ip) j1=j1-ip
         if (j2.gt.ip) j2=j2-ip
         ibm=max(iba(j),iba(j1))
         do i=1,ibm    
            ! find fal
            xsign=sign(1.,cl(i,j))
            fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j)   ! upwind
            flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)   ! LW
            x=fwbc(i+1,j)-fwbc(i,j)
            if (abs(x).le.1.e-27) fal(i,j)=fup
            if (abs(x).gt.1.e-27) then
               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
               if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
               if (r.le.0.) fal(i,j)=fup
               if (r.gt.0.) then
                  xlimiter=max(min(2.*r,1.),min(r,2.))
                  corr=flw-fup
                  fal(i,j)=fup+xlimiter*corr
                  if (fal(i,j).lt.0.) fal(i,j)=fup
               endif
            endif
            ! find fap
            xsign=sign(1.,cp(i,j))
            fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1)   ! upwind
            flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
            x=fwbc(i,j1)-fwbc(i,j)
            if (abs(x).le.1.e-27) fap(i,j)=fup
            if (abs(x).gt.1.e-27) then
               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j_1))/x
               if (xsign.eq.-1.) r=(fwbc(i,j2)-fwbc(i,j1))/x
               if (r.le.0.) fap(i,j)=fup
               if (r.gt.0.) then
                  xlimiter=max(min(2.*r,1.),min(r,2.))
                  corr=flw-fup
                  fap(i,j)=fup+xlimiter*corr
                  if (fap(i,j).lt.0.) fap(i,j)=fup
               endif
            endif
         enddo              ! end of do i=1,ir
      enddo                 ! end of do j=1,ip

      end subroutine inter_flux
