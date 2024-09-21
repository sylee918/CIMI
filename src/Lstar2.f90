!***********************************************************************
!                              calc_Lstar2  
!
! Routine modified from calc_Lstar. This routine calculates L * at given 
! (xlati,xmlt,K) in a given magnetic field configuration for the CIMI
! model of Mei-Ching H. Fok at (Code 673) NASA/GSFC.
! Calculation of magnetic flux and L* is made at ionophere
!
! 
! Input: xlati,xmlt,bo,iba,intB,ksai (field parameters)
! Output: L*, max L*
!
! Creaded by Suk-Bin Kang (Code 673) at NASA/GSFC on November 12, 2015
!
! Parameters
! - BEq: magnetic field strength at the Earth surface at magnetic equator
! - Bflux: magnetic flux enclosed by drift path at each grid.
! - ksai: output from cimi, magnetic flux per unit dphi*dvarL
! - Lstar: out L*
! - Lstar_max: L* of last closed drift shell.
! - mphi: magnetic longitude (in radian)
! - dphi: longitude grid size in radian at ionosphere.
! - xlat_i: interpolated MLAT.
! - b1,b2: variables to be used in xlat_i.
! - logbm_i: log10(bm) from i=1 to iba(j).
! - jdawn,jdusk,jnoon: MLT indices for dawn,dusk, and noon,
! - B_island: contains information if bm(i,j) is magnetic island
!           if 0, then normal, if 1, then magnetic island.
! - dBdMLAT: function that calculates integral B * dMLAT
!
! Modification history
!  * November 16, 2015 - correct interpolation of xlat_i ./. ii and ii+1
!                        to ii and ii1, skipping magnetic island.
!  * November 19, 2015 - L*(ir,ip) -> L*(ir,ip,ik)
!                        To avoid magnetic island
!                        locateB: ju=j(ii)+1 -> j(ii)
!                                 jl=j(1)    -> j(1)+1
!                        direct linear interpolation -> log linear "
!  * October 10, 2018  - add K=0 (bm=bo)                                        
!  * August 21, 2019 - Add the option of non-dipole internal field (intB)
!                      by Mei-Ching Fok
!  * September 17, 2019 - set ib0 to 1 -> calculate L* for all i (by MCFok)
!  * February 24, 2020 - Use magnetic longitude instead of SM phi to define
!                        the location at the ionosphere (by MCFok)
!**************************************************************************   

  subroutine calc_Lstar2(jnoon,Lstar,Lstar_max)

  use constants
  use cgrid,only: xlati,dlati,xlath,mphi,dphi,dvarL,ksai
  use cread1
  use cread2,only: rb,intB
  use cfield,only: bm,bo,ro,iba
  implicit none

  integer i,j,m,j0,ib0,i1,ii,ii1,jdawn,jdusk,jnoon,j1,j2,B_island(ir,ip),irb, &
          ji,jf
  real BEq,Bflux(ir,ip),dBdlat(ir),Lstar(ir,ip,0:ik),Lstar_max(0:ik),ffactor,&
       Bfluxg(ir,ip),Bflux1,Bflux1D(ir),xlat_i,dBdlat_i,b1,b2,logBm, &
       bmm(ir,ip,0:ik),xlatu(ir,ip)
  real,allocatable :: logbm_i(:)
  real,external :: dBdMLAT
 
  ib0=1                            ! inner boundary of L* calculation
  irb=ir-ib0+1

  BEq=xme/re_m**3          ! diople mangetic field on surface at equator

! Calculate Bfluxg from xlati upper bound (xlatu) to magnetic pole in dphi
! Bfluxg is magnetic flux per RE^2 
  do j=1,ip
     xlatu(1,j)=0.5*(xlath(1)*pi/180.+xlati(ir,j))  ! upper bound of xlati(ir)
     Bfluxg(1,j)=dBdMLAT(pi,BEq,rc,xlatu(1,j),mphi(j),intB)*dphi
  enddo
  ffactor=dphi*dvarL/(re_m*rc)**2     ! factor to convert ksai to Bflux
  do i=2,irb
     i1=ir+2-i
     do j=1,ip
        xlatu(i,j)=0.5*(xlati(i1,j)+xlati(i1-1,j))
        Bfluxg(i,j)=Bfluxg(i-1,j)+ksai(i1,j)*ffactor    
     enddo
  enddo

! Setup bmm
  bmm(1:ir,1:ip,0)=bo(1:ir,1:ip)
  bmm(1:ir,1:ip,1:ik)=bm(1:ir,1:ip,1:ik)

! Start m loop
  do m=0,ik

     ! Find magnetic island
     B_island(1:ir,1:ip)=0    ! 0=normal, 1=island
     do j=1,ip
        ii=1000
        do i=iba(j)-1,1,-1
           if (ii.lt.i) go to 20 
           if (bmm(i,j,m).lt.bmm(i+1,j,m)) then
              B_island(i,j)=1          
              ii=i
10            if(bmm(ii-1,j,m).lt.bmm(i+1,j,m)) then
                ii=ii-1
                B_island(ii,j)=1
                go to 10
              endif
        !     write(*,'(" bmm(i,j,m) has a island at ,m,j = ",2i3,"  i =
           endif
20         continue
        enddo
     enddo
   
     ! Caculates L* at each grid point
     Bflux=0.  ! Initialize magnetic fluxes
     Lstar(1:ib0-1,1:ip,m)=ro(1:ib0-1,1:ip)  ! set L*=ro at i<ib0
     do j0=1,ip                     ! MLT grid
        do i=ib0,iba(j0)
           Bflux(i,j0)=0.             ! initialization of magnetic flux  
           if (B_island(i,j0).eq.1) then 
              lstar(i,j0,m)=-1.*rb      
              go to 40                ! skipping magnetic island
           endif
           logBm=log10(bmm(i,j0,m))
           do j=1,ip
              ! locate drift shell (xlat_1) of constant bm at a mlt
              allocate(logbm_i(iba(j)))
              logbm_i(1:iba(j))=log10(bmm(1:iba(j),j,m))            !
              if (j.ne.j0.and.logBm.lt.logbm_i(1).and.logBm.ge.logbm_i(iba(j)))&
                 call locateB(logbm_i,iba(j),logBm,ii) ! Find the MLAT
              if (logBm.ge.logbm_i(1)) ii=1
              if (logBm.lt.logbm_i(iba(j))) ii=iba(j)
              if (j.eq.j0) ii=i
              if (ii.lt.iba(j)) then
                 ii1=ii
30               ii1=ii1+1
                 if(B_island(ii1,j).eq.1) go to 30
                 b2=logbm_i(ii1)-logBm
                 b1=logBm-logbm_i(ii)
                 xlat_i=(xlati(ii,j)*b2+xlati(ii1,j)*b1)/(b2+b1) !
              else
                 xlat_i=xlati(iba(j),j)                   
              endif
              ! Find Bflux from each local time sector
              Bflux1D(1:irb)=Bfluxg(1:irb,j)
              call lintp(xlatu(:,j),Bflux1D,irb,xlat_i,Bflux1) 
              Bflux(i,j0)=Bflux(i,j0)+Bflux1
              deallocate(logbm_i)
           enddo ! end of j
           lstar(i,j0,m)=2.*pi*BEq/Bflux(i,j0)
40         continue
       enddo    ! end of i loop
       if (iba(j0).lt.ir) lstar(iba(j0)+1:ir,j0,m)=rb     ! arbitrary
     enddo       ! end of j0 loop

     ! Determines L*max
     Lstar_max(m)=Lstar(iba(jnoon),jnoon,m) ! L* at noon magnetopause
     do j=1,ip  
        if(Lstar_max(m).gt.Lstar(iba(j),j,m))Lstar_max(m)=Lstar(iba(j),j,m)
     enddo
  enddo          ! end of m loop

  end subroutine calc_Lstar2


!**************************************************************************
!                           locateB
!  Routine modified from "locate1" in the CIMI model
!  Routine is return a value of j such that x is between xx(j) and xx(j+1),
!   excluding small humps and dips.
!  Routine is still valid if xx is a  monotonically increasing 
!   or decreasing function including small humps or dips by passing those.
!**************************************************************************
      subroutine locateB(xx,n,x,j)

      implicit none
      integer n,i,ii,iii,j,jl,ju,jm,j1,j2
      integer js(n)
      real xx(n),x

!  Make sure xx is increasing or decreasing monotonically
      js(1:n)=1000     ! arbitrary big value 
      ii=0        ! dummy index to count magnetic island
      i=n-1
10    i=i-1
      if (xx(i).le.xx(i+1)) then
         ii=ii+1
         js(ii)=i
20       if(xx(i-1).lt.xx(i+1)) then
           ii=ii+1
           i=i-1
           js(ii)=i
           go to 20
         endif
      endif
      if (i.gt.1) go to 10
       
   if (js(1).gt.999) then
      jl=0
      ju=n+1
40    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 40
      endif
      j=jl
   endif
   if (js(1).le.n) then
      jl=0
      ju=js(ii)
50    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 50
      endif
      j1=jl
      jl=js(1)
      ju=n+1
60    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 60
      endif
      j2=jl
      ! avoid a hump or dip
      j=j1
      if (x.le.xx(js(ii)-1).and.x.gt.xx(js(1)+1)) j=js(ii)-1
      if (x.le.xx(js(1)+1)) j=j2
   endif

   end subroutine locateB


!***************************************************************************
!                             dBdMLAT
!  Function calculates magnetic flux per unit PHI from mlat to magnetic pole
!  at the ionosphere. dBdMLAT is in magnetic flux per RE^2
!****************************************************************************
     function dBdMLAT(pi,BEq,rc,mlat,phij,intB)

     implicit none
     integer intB,m_one,ilat,i
     real BEq,rc,mlat,dBdMLAT,pi,pi2,rc2,phij,cosphi,sinphi,dlat,rcdlat, &
          mlat1,coslat,sinlat,xmag,ymag,zmag,xgsm,ygsm,zgsm,xgeo,ygeo,zgeo, &
          theta,xlong,Bri1,Btheta,Bphi


     m_one=-1
     ilat=20
     pi2=pi/2.
     rc2=rc*rc

     if (intB.eq.0) dBdMLAT=BEq*(1.-sin(mlat)**2)/rc   ! dipole

     if (intB.gt.0) then        ! IGRF
        cosphi=cos(phij)
        sinphi=sin(phij)
        dlat=(pi2-mlat)/ilat
        rcdlat=rc2*dlat
        dBdMLAT=0.
        do i=1,ilat
           mlat1=mlat+(i*1-0.5)*dlat
           coslat=cos(mlat1)
           sinlat=sin(mlat1)
           xmag=rc*coslat*cosphi
           ymag=rc*coslat*sinphi
           zmag=rc*sinlat
           call GEOGSW_08(xgeo,ygeo,zgeo,xmag,ymag,zmag,m_one)    ! mag to geo
           theta=acos(zgeo/rc)
           xlong=atan2(ygeo,xgeo)
           call IGRF_GEO_08(rc,theta,xlong,Bri1,Btheta,Bphi)
           Bri1=abs(Bri1)*1.e-9         ! Br in Tesla
           dBdMLAT=dBdMLAT+Bri1*rcdlat*coslat
        enddo
     endif

     end function dBdMLAT
