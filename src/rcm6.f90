!*****************************************************************************
!  
!                          rcm6.f90 
!
!  This is a modification of rcm5.f90.
!
!  This code use birk as input.
!
!  Created on November 23, 2018 by Mei-Ching Fok, Code 673, NASA GSFC.
!
!  Modification history
!  February 8, 2021
!    - Add index jdif to control the conversion between rcm and cimi 
!      azimuthal grids
!    - Change ri=re+100. to ri=re+hlosscone
!*****************************************************************************


!  rcm potential calculation machinery for
!  inclusion in Mei-Ching Fok ring current code
!
!  02-04-00    r w spiro
!
!
!**********************************************************
!  The following modules define and dimension
!  rcm-specific variables used to compute field aligned
!  currents and electric potential distribution.  They 
!  should be in Fok's main code
!**********************************************************
!
!***********************************************************
! Explanation of variables
!***********************************************************
!
! vfok   = output electric potential (volts) on fok grid
! xlati = fok latitude grid (radian)
! xmlt = fok local time grid (hours from midnight)
! volume = flux tube volume per unit magnetic flux at each spatial
!          grid point in meter/Tesla
! xmlto = local time of grid pt mapped to magnetic equatorial plane
!         (hours from midnight)
! ro = radial distance (Re) of grid pt mapped to magnetic equatorial
!      plane
!****************************************************************
!

      subroutine tube_volume(ro,xmlto,volume)
!
!********************************************************************
! Taken from the main of testdrive8.f
!
! Mei-Ching Fok          02-20-00
!********************************************************************
!
      use rcmgrid
      use conductance
      use plasma
      use potential

      dimension volume(idim,jdimf)
      dimension xmlto(idim,jdimf)
      dimension ro(idim,jdimf)
!
      pi=atan2(0.,-1.)

! loop over latitude index
      i=0
      do ic=idim,1,-1
         i=i+1
! loop over local time index
         do jc=1,jdimf
            j=jc+jdif
            if (j.gt.jdim) j=j-jdimf
!
! calculate vm from volume and put on rcm grid
! Note: fok grid extends from low latitude to high
!       rcm grid extends from high latitude to low
!       fok grid extends from midnight to midnight
!       rcm grid extends from noon to noon
!  ic and jc are fok grid coordinates
!  i and j are corresponding rcm grid coordinates.
!
            vm(i,j)=(volume(ic,jc)/6.370e6/1.e9)**(-2./3.)
            xxmlt=xmlto(ic,jc)+12.  !measure xxmlt from noon
            if(xxmlt.ge.24.) then
               xxmlt=xxmlt-24.  !ensure 0.< xxmlt <24.
            end if
            phi=xxmlt*pi/12.   !change from hours to radians
!
! calculate x and y mappings in magnetic equatorial plane
            xmin(i,j)=ro(ic,jc)*cos(phi)
            ymin(i,j)=ro(ic,jc)*sin(phi)
!           write(11,*) ic,jc,i,j,phi,xmin(i,j),ymin(i,j)
!
         end do
!
! extend to overlap region of rcm grid
         vm(i,jdim)=vm(i,3)
         vm(i,2)=vm(i,jdim-1)
         vm(i,1)=vm(i,jdim-2)
         xmin(i,jdim)=xmin(i,3)
         xmin(i,2)=xmin(i,jdim-1)
         xmin(i,1)=xmin(i,jdim-2)
         ymin(i,jdim)=ymin(i,3)
         ymin(i,2)=ymin(i,jdim-1)
         ymin(i,1)=ymin(i,jdim-2)
      end do
!     call outp(vm,idim,jdim,1,idim,1,1,jdim,1,0.,'vm',6,80)

      return
      end

     
      subroutine computev(vfok)
!
!*********************************************************************
! Control routine for rcm calculation to computed electric potential
!  distribution from (mu,K) representation of plasma distribution
!
! Bob Spiro               02-04-00
!
! modified                02-21-02 by r. w. spiro to include call to
!                         eqbndy to set equatorial bndy condition
!
!*********************************************************************
!
      use rcmgrid
      use conductance
      use plasma
      use current
      use potential
      use coefficients
!
      dimension vfok(idim,jdimf)

! Setup ibndy
      do j=1,jdim
         ibndy(j)=nint(ain(j))
      enddo

! specify conductance distribution
!     call cond(kp)
!
! calculate coefficients of elliptic pde
      call coeff
!
! calculate rhs of pde
      call cfive
!
! adjust rhs for bndy conditions
      call aphi
!
!*********************************************
! call to eqbndy added 02-21-02 by rws
!*********************************************
! calculate low-latitude bndy conditions
        call eqbndy
!
! solve elliptic magnetosphere/ionosphere coupling pde
      call gmresm (ain, c, vrcm, idim, jdim, jj1, jj2, ncoeff)
!     call outp(vrcm,idim,jdim,1,idim,1,1,jdim,1,0.,'v',6,80)
!
! output for creation of plot file
!     write(99,*) colat,aloct,xmin,ymin,vrcm,ain
!
! convert from v on rcm grid to v on fok grid
      call rcm2fok_v(vfok)
!     call outp(vfok,idim,jdimf,1,idim,1,1,jdimf,1,0.,'vfok',6,80)
!
      return   
      end               ! end of subroutine computev
!

!***************************************************************
      subroutine cond(kp)
!***************************************************************
! subroutine to compute conductance arrays for Mei-Ching Fok
!  electric potential calculation
!
!  06-11-99    r w spiro
!
      use rcmgrid
      use plasma
      use conductance
!
!     open(unit=18,file='rcmcon1',form='unformatted')
!     read(18) qtplam,qthall,qtped,ss
!     close(18)
!
      open(unit=19,file='rcmcon1a')  ! quiet time conductances depend on iday

      read(19,*) qtplam,qthall,qtped,ss
      close(19)

!     kp=3
      call condtot(kp,ExAC)
!
!     call outp(pedlam,idim,jdim,1,idim,1,1,jdim,1,0.,'pedlam',6,80)
!     call outp(pedpsi,idim,jdim,1,idim,1,1,jdim,1,0.,'pedpsi',6,80)
!     call outp(hall,idim,jdim,1,idim,1,1,jdim,1,0.,'hall',6,80)
!
      return
      end


! *******************************************************************
!                        c o n d t o t
! *******************************************************************
      subroutine condtot(kp,ExAC)

!  this subroutine calculates conductance enhancement due to auroral
!   electron precipitation and adds this to the quiet time conductances
!  last update: 11-06-86
!               02-07-96 frt - min_j array added
!               06-10-99 rws  mcfok version
!
!***********************************************************************
!
      use rcmgrid
      use conductance

      dimension awork(idim,jdim)

      pi=atan2(0.,-1.)
!
      do j=1,jdim
      do i=1,idim
! MCF 9 Mar 2017: ExAC is the factor to extend auroral conductance to low lat
         thelat=90.-colat(i,j)*180./pi
         thelat=90.-(colat(1,j)+ExAC*(colat(i,j)-colat(1,j)))*180./pi
! end
         themlt=12.+aloct(i,j)*12./pi
         if(themlt.gt.24.) then
            themlt=themlt-24.
         end if
!
! icase=3 ==> hall conductance
         icase=3
         call elemod(icase,kp,thelat,themlt,sigh)
!
! icase=4 ==> pedersen conductance
         icase=4
         call elemod(icase,kp,thelat,themlt,sigp)
!        write(23,*) 'i,j,lat,mlt,sigp,sigh ',i,j,thelat,themlt,
!    2                sigp,sigh


         pedpsi(i,j)=2.*sigp+qtped(i,j)
         pedlam(i,j)=2.*sigp/(sini(i,j)**2)+qtplam(i,j)
         hall(i,j)=2.*sigh/sini(i,j)+qthall(i,j)
! test case -- constant conductance
!        pedpsi(i,j)=10.
!        pedlam(i,j)=10./(sini(i,j)**2)
!        hall(i,j)=10./sini(i,j)
!
      enddo
      enddo

      njsmth=1
      nismth=1

      return
      end

      subroutine elemod(icase,ikp,glat,amlt,value)
!
!***********************************************************
!   subroutine to evaluate the hardy average auroral model
!     described in:  hardy, et al., "j. geophys. res.",
!     vol. 92, page 12,275-12,294 (1987).
!
!   please contact:   william mcneil
!                     radex, inc.
!                     three preston court
!                     bedford, massachusetts  01730
!                     (617)275-6767
!
!   with questions or comments.
!
!   inputs:
!
!     icase=1   energy flux
!           2   number flux
!           3   hall conductivity
!           4   pederson conductivity
!
!     ikp       kp division 0 to 6
!
!     glat      geomagnetic latitude
!
!     amlt      magnetic local time
!
!   outputs:
!
!     value     log10 energy flux in kev/cm**2-s-sr (icase=1)
!               log10 number flux in electrons/cm**2-s-sr (icase=2)
!               conductivity in mhos (icase=3,4)
!
!   internal variables
!
!     crd       coefficients for maximum function value
!     chat      coefficients for latitude of maximum value
!     cs1       coefficients for up-slope
!     cs2       coefficients for down-slope
!     cutl      low latitude cut-off value
!     cuth      high latitude cut-off value
!
!   files:
!
!     the file elecoef.dat must be present in the default
!     directory.
!
!   notes:
!
!     this version operates on vax/vms or ibm-pc machines.
!
!   21 june 1993 -- wjm
!***************************************************************
!
      common/savcom1/crd,chat,cs1,cs2,cutl,cuth
!
      dimension crd(13,7,4),chat(13,7,4),cs1(13,7,4),cs2(13,7,4)
      dimension cutl(4),cuth(4)
!
      character aline*80
!
!
      data ifirst/0/
!
      cutl(1)=6.
      cutl(2)=6.
      cutl(3)=0.
      cutl(4)=0.
!
      cuth(1)=7.
      cuth(2)=7.
      cuth(3)=.55
      cuth(4)=.55
!
      if(ifirst.eq.0)then
        ifirst=1
!
!   for plotting only, reset the poleward cutoffs
!
        open(71,file='elecoef.dat',form='formatted',&
                status='old')
        read(71,100)aline
  100   format(a80)
        do 10 jcase=1,4
        do  jkp=1,7
        do  jco=1,13
        read(71,101)crd(jco,jkp,jcase),chat(jco,jkp,jcase),&
                    cs1(jco,jkp,jcase),cs2(jco,jkp,jcase)
  101   format(26x,4f12.7)
        enddo
        enddo
   10   continue
        close(71)
      endif
!
!   return if glat < 50
!
      if(glat.lt.50) return
!
      kp=ikp+1
      xarg=amlt*3.14159265/12.
      rd=crd(1,kp,icase)
      hat=chat(1,kp,icase)
      s1=cs1(1,kp,icase)
      s2=cs2(1,kp,icase)
!
      do 11 j=1,6
      xa=j*xarg
      c=cos(xa)
      s=sin(xa)
      ipc=j+1
      ips=j+7
!
      rd=rd+c*crd(ipc,kp,icase)+s*crd(ips,kp,icase)
      hat=hat+c*chat(ipc,kp,icase)+s*chat(ips,kp,icase)
      s1=s1+c*cs1(ipc,kp,icase)+s*cs1(ips,kp,icase)
      s2=s2+c*cs2(ipc,kp,icase)+s*cs2(ips,kp,icase)
!
   11 continue
!
      value=epst(glat,rd,hat,s1,s2,cutl(icase),cuth(icase))
!
      return
      end

      function epst(clat,rd,hat,s1,s2,xmin,xmax)
!
      d=clat-hat
      ex=exp(d)
      xl=(1.-s1/s2*ex)/(1.-s1/s2)
      xl=alog(xl)
      ep=rd+s1*d+(s2-s1)*xl
!
      if (clat.lt.hat.and.ep.lt.xmin) ep=xmin
      if (clat.gt.hat.and.ep.lt.xmax) ep=xmax
!
      epst=ep
      return
      end
!
          
! *******************************************************************
!                      c o e f f
! *******************************************************************
      subroutine coeff
!
!
! code based on subroutine coeff in spiro.agu83.fort
!
!
!***********************************************************************
!
!
! this subroutine computes the coefficients c1,c2,c3 & c4.  these are
! coefficients of the elliptic magnetosphere-ionosphere coupling
!  equation that is solved in potent.  computed values
! of the coeffecients are stored in array c.
!
!  this subroutine called from subroutine computev
!
!
!***********************************************************************
!
!
!
      use rcmgrid
      use conductance
      use coefficients
      use potential
!
!
!
!
! work arrays
      dimension a(idim,jdim),b(idim,jdim),d(idim,jdim)
!
!
!
!
! initialize coefficient array
      do  j=1,jdim
      do  i=1,idim
      do  k=1,ncoeff
          c(k,i,j)=0.0
      enddo
      enddo
      enddo
!
!
!
      do  j=1,jdim
! MCF
   !  do  i=ibndy(j),idim
      do  i=1,idim
! MCF end
         a(i,j)=alpha(i,j)*pedpsi(i,j)/beta(i,j)
         b(i,j)=beta(i,j)*pedlam(i,j)/alpha(i,j)
         d(i,j)=2.*(b(i,j)/dlam**2+a(i,j)/dpsi**2)
      enddo   
      enddo   
!     call outp(a,idim,jdim,imin,imin+8,1,20,35,1,0.0,ilabel,'a',ntp,ncol)
!     call outp(b,idim,jdim,imin,imin+8,1,20,35,1,0.0,ilabel,'b',ntp,ncol)
!     call outp(d,idim,jdim,imin,imin+8,1,20,35,1,0.0,ilabel,'d',ntp,ncol)
!
!
      do 20 i=1,idim
    	  do 30 j=jj1,jj2,jint

! crude fix here frt 2/96
	  if(i.lt.ibndy(j))go to 30
!
            if (i.lt.idim.and.i.gt.ibndy(j)+1) then
               bb=b(i+iint,j)-b(i-iint,j)
! comment on 7/26/90 by yong, introduce signbe for hall

!              ee=hall(i+iint,j)-hall(i-iint,j)
               ee=signbe*hall(i+iint,j)-signbe*hall(i-iint,j)
! end of comment
            else if(i.eq.imax) then
               bb=3.*b(imax,j)-4.*b(imax-1,j)+b(imax-2,j)
! comment by yong on 7/26/90, introduce signbe for hall
!              ee=3.*hall(imax,j)-4.*hall(imax-1,j)+hall(imax-2,j)
               ee=3.*signbe*hall(imax,j)-4.*signbe*hall(imax-1,j)+&
                  signbe*hall(imax-2,j)
! end of comment by yong on 7/26/90
            else
               bmin=2.*b(i,j)-b(i+1,j)
               bc=0.5*b(i,j)
               if(bmin.lt.bc) bmin=bc
               bb=b(i+1,j)-bmin
!
! comment by yong on 7/26/90, add "signbe" to the following rcm lines
! where there is hall
               hmin=2.*signbe*hall(i,j)-signbe*hall(i+1,j)
               hc=0.5*signbe*hall(i,j)
               if(abs(hmin).lt.abs(hc)) hmin=hc
               ee=signbe*hall(i+1,j)-hmin
            end if
!
            cc=signbe*hall(i,j+jint)-signbe*hall(i,j-jint)
! end of change
            dd=(bb-cc*dlam/dpsi)*0.25
            aa=a(i,j+jint)-a(i,j-jint)
            ff=(aa+ee*dpsi/dlam)*0.25
            c(1,i,j)=(b(i,j)+dd)/(d(i,j)*dlam**2)
            c(2,i,j)=(b(i,j)-dd)/(d(i,j)*dlam**2)
            c(3,i,j)=(a(i,j)+ff)/(d(i,j)*dpsi**2)
            c(4,i,j)=(a(i,j)-ff)/(d(i,j)*dpsi**2)
!
!           if(i.le.imin+7.and.j.gt.20.and.j.lt.35) then
!              write(6,*) '****** i= ',i,'*****  j= ',j,'   ******'
!              write(6,*) 'bb= ',bb,'   ee= ',ee,'   cc= ',cc
!              write(6,*) 'dd= ',dd,'   aa= ',aa,'   ff= ',ff
!              write(6,*) 'c1,c2,c3,c4   ',c(1,i,j),c(2,i,j),c(3,i,j),
!    2                     c(4,i,j)
!           end if
!
!
   30    continue
!
!
         do 40 k=1,ncoeff
          if(jj1.gt.1) then
           do j=1,jj1-1
            c(k,i,j)=c(k,i,jdim-jwrap+j)
           end do
          end if
!
          if(jj2.lt.jdim) then
           do j=jj2+1,jdim
            c(k,i,j)=c(k,i,j-jdim+jwrap)
           end do
          end if
   40    continue
!
!
   20 continue
!
!
      return
      end
! *******************************************************************
!                   c f i v e
! *******************************************************************
      subroutine cfive
!***********************************************************************
!
!
!  program written by: r.w. spiro        last update: 06-26-85
!                                                     05-27-99
!                                streamlined for M.-C. Fok
!  algorithm by: r.a. wolf
!
! this subroutine computes c(5,i,j) from previously calculated birk(i,j)
!
!
!***********************************************************************
!
!
!
      use rcmgrid
      use coefficients
      use conductance
      use current

      dimension c5out(idim,jdim)


!
!
!
      ntp = 6
!
!  zero c(5,i,j) array
!
!
!
      do  i=1,imax
      do  j=1,jmax
         c(5,i,j)=0.0
      enddo
      enddo


      do  i=imin,imax
      do  j=1,jmax
         d=2.*(beta(i,j)*pedlam(i,j)/(alpha(i,j)*dlam**2)&
            +alpha(i,j)*pedpsi(i,j)/(beta(i,j)*dpsi**2))
         if(d.le.1.e-30) then
            c(5,i,j)=0.
         else
         c(5,i,j)=alpha(i,j)*beta(i,j)*(ri**2)*birk(i,j)/d
         end if
         c5out(i,j)=c(5,i,j)
      enddo
      enddo
!
! debug 08 feb 00
!     write(6,*) 'at end of subroutine cfive'
!     call outp(c5out,idim,jdim,1,idim,1,1,jdim,1,0.,'c5',6,80)

      return
      end
!
      subroutine hilat_bndy(vdrop,delmlt)
!
!*************************************************************
! subroutine to set rcm high-latitude boundary condition
!  and high-latitude boundary location.
!
! Modified 03-02-00
!   vdrop is one-half of the cross polar cap potential drop
!   delmlt is the number of hours to rotate the axis
!     of the potential pattern.  delmlt=2. rotates the
!     axis of the potential pattern from noon-midnight to
!     1400-0200.  delmlt=-2. rotates the axis to be from 1000
!     to 2200.  delmlt=0. makes the axis noon-midnight.
!*************************************************************
      use rcmgrid
      use potential
!
      pi=atan2(0.,-1.)
!     jdimf2=jdimf/2
!
      do j=1,jdim
!        ain(j)=2.                     ! high latitude boundary at xlati(ir-1)
! MCF
!        ain(j)=1.                     ! high latitude boundary at xlati(ir)
!        vpob(j)=-vdrop*sin(float(j-3)*pi/24. - delmlt*pi/12.)
         vpob(j)=-vdrop*sin(float(j-3)*2.*pi/float(jdim-3) - delmlt*pi/12.)
! MCF end
      end do                           ! MCF, March 14, 2002
!
      return
      end
!
      subroutine rcm2fok_v(vfok)
!
!*****************************************************************
! subroutine to rotate electric potential distribution from
!  rcm grid to fok grid
!
! 06/14/99      r w spiro
!
! 03/01/00      r w spiro
!  corrected error that returned potentials 12 hours out of phase
!*****************************************************************
!
      use rcmgrid
      use potential
      dimension vfok(idim,jdimf)

      do j=3,jdim-1
         jc=j-jdif
         if (jc.lt.1) jc=jc+jdimf
         if (jc.gt.jdimf) jc=jc-jdimf
         do i=1,idim
            ic=idim+1-i
            vfok(ic,jc)=vrcm(i,j)       ! both vfok and vrcm are in volt
         end do
      end do
!
      return
      end
!
!*************************************************************************
! subroutine to determine coefficients to implement polar  boundary cond
!
! acm --- june 30,1995
!
!************************************************************************
      subroutine aphi
!
      use rcmgrid
      use potential
      use coefficients
      use conductance
!
! please note that a,b,c,d,e now refer to the rcm code variables
! alpha, beta are from rcm, my alpha -> alp1, beta -> bet1
! will keep f as c5 (source)
! will keep gamma

!     implicit none

      parameter(isize=idim,jsize=jdim)
      parameter(maxdim=jsize)


      logical lsflag(maxdim)
!

      dimension a(idim,jdim),b(idim,jdim),d(idim,jdim)
!
      dimension alp1(maxdim),bet1(maxdim),gamma(maxdim)
      dimension vb(maxdim),vc(maxdim)
      dimension f(idim,jdim)
      dimension c5out(idim,jdim)
!
      real rlf,rlb

      external case3

      tol = 1e-06

      if(maxdim .ne. jdim) then
      write(*,*) 'stopping program --- arrays not dimensioned right'
      stop
      endif
!
      do j=1,jdim
         ibndy(j)=ibndy(j)
      end do

! zero f's, the source term

      do i=1,idim
       do j=1,jdim
          f(i,j) = 0.0
       enddo
      enddo

! alp1 = seperation between physical boundary and nearest grid line (above)
! alp1 = ibndy(j)-ain(j)
! ain = `grid' distance, not an integer
! need to determine ibndy --- first grid point

! next need to figure out betl and gamma
! first figure out if special case of point near grid line
! will reset min annd alpha here
!
! bet1 = separation between physical boundary above and phi_j
! also need to get potential at these points!

! next do loop replaced in call to read_vext frt 2/96
!     do j=2,jdim-1
!         ibndy(j) = int(ain(j)) + 1
!       enddo

!     ibndy(1) = ibndy(jdim-2)
!     ibndy(jdim)   = ibndy(3)

      do j=2,jdim-1
       diff = ain(j)
         idiff = int(diff)
         rmin = float(ibndy(j))
         alp1(j) = rmin - ain(j)
         icp = ibndy(j+1)-ibndy(j)
         icm = ibndy(j)-ibndy(j-1)

! first figure out if special case of point near grid line
! will reset min annd alpha here

       if( (diff - float(idiff)) .lt. tol) then
         lsflag(j) = .true.
         ibndy(j) = idiff
         alp1(j) = 0.0
         bet1(j) = 1.0
         gamma(j) = 1.0
       elseif ( (rmin-diff) .lt. tol) then
         lsflag(j) = .true.
         alp1(j) = 0.0
         bet1(j) = 1.0
         gamma(j) = 1.0
       endif

        if(.not. lsflag(j)) then
! get betl
         if( icp .le. 0) then
          bet1(j) = 1.0
          vb(j) = 1.0     ! will be v(i,j+1)
         else
           bet1(j) = (rmin-ain(j))/(ain(j+1)-ain(j))
! need vb here!

           rlf = sqrt( (1.0-bet1(j))**2 + (ain(j+1)-rmin)**2 )
           rlb = sqrt( bet1(j)**2 + (rmin-ain(j))**2 )
           vb(j) = vpob(j+1)*rlb + vpob(j)*rlf
           vb(j) = vb(j)/(rlf + rlb)

        endif
! get gamma
        if(icm .ge. 0) then
            gamma(j) = 1.0
            vc(j) = 1.0   ! will be v(i,j-1)
        else
           gamma(j) = (rmin-ain(j))/(ain(j-1)-ain(j))

           rlb = sqrt( (1.0-gamma(j))**2 + (ain(j-1)-rmin)**2 )
           rlf = sqrt( gamma(j)**2 + (rmin-ain(j))**2 )
           vc(j) = vpob(j)*rlb + vpob(j-1)*rlf
           vc(j) = vc(j)/(rlf + rlb)

        endif
      endif

      enddo
!20    format(1x,i4,2x,f6.4,2x,i4,2x,f8.4,2x,f8.4,2x,f8.4)
!22    format(1x,i4,2x,i4,2x,f6.4,2x,f6.4,2x,f6.4,2x,f6.4)

      lsflag(1) = lsflag(jdim-2)
      lsflag(jdim)   = lsflag(3)

      ibndy(1) = ibndy(jdim-2)
      ibndy(jdim)   = ibndy(3)

      bet1(1) = bet1(jdim-2)
      bet1(jdim)   = bet1(3)

      vb(1) = vb(jdim-2)
      vb(jdim)   = vb(3)

      gamma(1) = gamma(jdim-2)
      gamma(jdim)   = gamma(3)

      vc(1) = vc(jdim-2)
      vc(jdim)   = vc(3)
!
!
! calculate new coefficients for points near inner boundary
! 4 cases to consider: which reduces to 2
!  (1) 1 leg crossing, |
!  (2)                   | all three should be covered by case3
!  (3)                   |
!  (4) physical boundary pt close to grid point
!
! in the cases where gamma (or beta) is less than one, we put the contribution
! of these points into the source term and then set the coefficient to 0
! so that when sor runs that nearest neighbor is not used

      do j=2,jdim-1
        i = ibndy(j)

         if(lsflag(j)) then
             c(1,i,j) = 0.0
           c(2,i,j) = 0.0
           c(3,i,j) = 0.0
           c(4,i,j) = 0.0
           f(i,j) = vpob(j)
         else
       call case3(a,b,d,&
        alp1(j),bet1(j),gamma(j),i,j)
           if( j .gt. 13 .and. j .lt. 18 ) then
           endif
           if( gamma(j) .lt. 1.0) then
             f(i,j) = f(i,j) +c(4,i,j)*vc(j)
             c(4,i,j) = 0.0
           endif
           if( bet1(j) .lt. 1.0) then
             f(i,j) = f(i,j) +c(3,i,j)*vb(j)
             c(3,i,j) = 0.0
           endif
             f(i,j) = f(i,j) +c(2,i,j)*vpob(j)
           c(2,i,j) = 0.0
         endif

! c(1,i,j) -> v(i+1,j) or v4
! c(2,i,j) -> v(i-1,j) or va
! c(3,i,j) -> v(i,j+1) or v2
! c(4,i,j) -> v(i,j-1) or v3

      enddo
! periodicity --- use circle?
!
      do k=1,5
      c(k,ibndy(1),1) = c(k,ibndy(jdim-2),jdim-2)
        f(ibndy(1),1) = f(ibndy(jdim-2),jdim-2)
      c(k,ibndy(jdim),jdim) = c(k,ibndy(3),3)
        f(ibndy(jdim),jdim) = f(ibndy(3),3)
      enddo
!
! now need to determine if grid points above min need special attention
! loop through j(phi)
! for each j, determine what i(r) grid lines need nearest neighbors
! at these i(r) lines need to calc. bet1
! then need to call case3 and then get f
!

      do j=2,jdim-1
        i = ibndy(j)
        ipb = ibndy(j+1)
        im = ibndy(j-1)
! get the higher one
        ibig = ipb
        if(im .gt. ipb) ibig = im
        alp1(j) = 1.0
        do l=i+1,ibig-1
! commented out
           if ( ipb .gt. l) then
           rmin = float(l)
           bet1(j) = (rmin-ain(j))/(ain(j+1)-ain(j))

           rlf = sqrt( (1.0-bet1(j))**2 + (ain(j+1)-rmin)**2 )
           rlb = sqrt( bet1(j)**2 + (rmin-ain(j))**2 )
           vb(j) = vpob(j+1)*rlb + vpob(j)*rlf
           vb(j) = vb(j)/(rlf + rlb)
           else
           bet1(j) = 1.0
           vb(j)   = 1.0
        endif
        if( im  .gt. l) then
           rmin = float(l)
           gamma(j) = (rmin-ain(j))/(ain(j-1)-ain(j))

           rlb = sqrt( (1.0-gamma(j))**2 + (ain(j-1)-rmin)**2 )
           rlf = sqrt( gamma(j)**2 + (rmin-ain(j))**2 )
           vc(j) = vpob(j)*rlb + vpob(j-1)*rlf
           vc(j) = vc(j)/(rlf + rlb)
         else
           gamma(j) = 1.0
           vc(j)   = 1.0
        endif

      call case3(a,b,d,&
         alp1(j),bet1(j),gamma(j),l,j)

            if( gamma(j) .lt. 1.0) then
             f(l,j) = f(l,j) +c(4,l,j)*vc(j)
             c(4,l,j) = 0.0
            endif
            if( bet1(j) .lt. 1.0) then
             f(l,j) = f(l,j) +c(3,l,j)*vb(j)
             c(3,l,j) = 0.0
            endif
        c(1,l,1) = c(1,l,jdim-3)
        c(2,l,1) = c(2,l,jdim-3)
        c(3,l,1) = c(3,l,jdim-3)
        c(4,l,1) = c(4,l,jdim-3)
        f(l,1) = f(l,jdim-3)

        c(1,l,jdim) = c(1,l,3)
        c(2,l,jdim) = c(2,l,3)
        c(3,l,jdim) = c(3,l,3)
        c(4,l,jdim) = c(4,l,3)
        f(l,jdim) = f(l,3)
        enddo    ! l-loop
      enddo    ! j-loop

!     endif
! add source term, f, to c5
        do i=1,idim
          do j=1,jdim
            c(5,i,j) = f(i,j) + c(5,i,j)
            c5out(i,j)=c(5,i,j)
         enddo
      enddo
!
!     write(6,*) 'at end of aphi'
!     call outp(c5out,idim,jdim,1,idim,1,1,jdim,1,0.,
!    2          'c5',6,80)

        return
      end
! *******************************************************************
!                      c a s e 3
! *******************************************************************
      subroutine case3(a,b,d,&
           alp1,bet1,gamma,i,j)
!
!
! code based on subroutine coeff in rcm.f, copied from on 23-may-95
!
! acm --- june 30,1995
!
!***********************************************************************
!
!
! this subroutine computes the coefficients c1,c2,c3 & c4.  these are
! coefficients of the elliptic magnetosphere-ionosphere coupling
! equation that is solved in potent.  computed values
! of the coeffecients are stored in array c.
!
!  this subroutine called from subroutine aphi
!
! this computes the coefficients for the inner boundary with possible 3 legged
! crossing
!***********************************************************************
!
!
      use rcmgrid
      use coefficients
      use conductance



!
!
      dimension a(idim,jdim),b(idim,jdim),d(idim,jdim)
!
! modify the differentials for the cases where the crossings are not at grid
! lines

       k = dlam
       h = dlam*alp1
       dpf = dpsi*bet1
       dpb = dpsi*gamma

! added this exit case if i lt ibndy(i) - frt


       c(1,i,j) = 0.0
       c(2,i,j) = 0.0
       c(3,i,j) = 0.0
       c(4,i,j) = 0.0

!      if(i.lt.ibndy(j))return

!
         a(i,j)=alpha(i,j)*pedpsi(i,j)/beta(i,j)
         b(i,j)=beta(i,j)*pedlam(i,j)/alpha(i,j)

! d(i,j) depends on the differences which will not be the same as for
! interior points ---- use smp generated eij
!
       arcm = a(i,j)
       brcm = b(i,j)

! calculate the differences needed in the c's : aa,bb,cc,dd,ee,ff
! these need to be modified on the boundaries due to points being no
! longer in the physical region

! ii1 is imin+1
!            if (i.lt.imax.and.i.gt.ii1) then
! if alp1 = 1 then it is an interior point and we can use a central
! difference formula for [i]

            if (alp1 .eq. 1.0) then
               bb=b(i+iint,j)-b(i-iint,j)

! comment on 7/26/90 by yong, introduce signbe for hall
!              ee=hall(i+iint,j)-hall(i-iint,j)

               ee=signbe*hall(i+iint,j)-signbe*hall(i-iint,j)

! shouldn't need this case so erased
!            else if(i.eq.imax) then
! i'll just leave this as is for alp1<1 need forward difference

            else
               bmin=2.*b(i,j)-b(i+1,j)
               bc=0.5*b(i,j)
               if(bmin.lt.bc) bmin=bc
               bb=b(i+1,j)-bmin
               hmin=2.*signbe*hall(i,j)-signbe*hall(i+1,j)
               hc=0.5*signbe*hall(i,j)
               if(abs(hmin).lt.abs(hc)) hmin=hc
               ee=signbe*hall(i+1,j)-hmin
            end if
!
! now need to look at psi derivatives
! if bet1<1 then need backward difference --- must x by 2 because cent. diff.
! doesn't have the half
! if gamma<1 then need forward difference
! if both are < 1 then take derivative equal to 0
         if( bet1 .lt. 1.0 ) then
           if( gamma .lt. 1.0 ) then
            cc = 0.0
              aa = 0.0
           else
            cc = (signbe*hall(i,j) -signbe*hall(i,j-jint))*2.0
            aa = (a(i,j) - a(i,j-jint))*2
           endif
          elseif( gamma .lt. 1.0) then
            cc = (signbe*hall(i,j+jint) -signbe*hall(i,j))*2.0
            aa = (a(i,j+jint)-a(i,j))*2.0
          else
            cc=signbe*hall(i,j+jint)-signbe*hall(i,j-jint)
            aa=a(i,j+jint)-a(i,j-jint)
          endif

! put in dx parts, and factor of 2
       aa = aa/dpsi/2.
       bb = bb/dlam/2.
       cc = cc/dpsi/2.
       ee = ee/dlam/2.

!
! now calculate coefficients c1-c4 using smp generated formulas


      aij = (2 * brcm + bb * h + (-1) * cc * h) / (k * (h + k))

      bij = (-1) * (((-2) * brcm + bb * k + (-1) * cc * k) / (h * (h +&
          k)))
      cij = (2 * arcm + aa * dpb + dpb * ee) / (dpf * (dpb + dpf))

      dij = (-1) * (((-2) * arcm + aa * dpf + dpf * ee) / (dpb * (dpb +&
          dpf)))

! d(i,j) depends on the differences which will not be the same as for
! interior points ---- use smp generated eij
!
       arcm = a(i,j)
       brcm = b(i,j)

!     t1 = (-2) * ((brcm * h) / (h * k ** 2 + h ** 2 * k)) + (-2) * ((
!    $    brcm * k) / (h * k ** 2 + h ** 2 * k)) + (cc * h ** 2) / (h *
!    $    k ** 2 + h ** 2 * k) + (-1) * ((cc * k ** 2) / (h * k ** 2 +
!    $    h ** 2 * k)) + (-1) * ((dpb ** 2 * ee) / (dpb * dpf ** 2 +
!    $    dpb ** 2 * dpf)) + (dpf ** 2 * ee) / (dpb * dpf ** 2 + dpb **
!    $     2 * dpf)
!     eij = (-1) * ((aa * dpb ** 2) / (dpb * dpf ** 2 + dpb ** 2 * dpf)
!    $    ) + (aa * dpf ** 2) / (dpb * dpf ** 2 + dpb ** 2 * dpf) + (-2
!    $    ) * ((arcm * dpb) / (dpb * dpf ** 2 + dpb ** 2 * dpf)) + (-2)
!    $    * ((arcm * dpf) / (dpb * dpf ** 2 + dpb ** 2 * dpf)) + (-1) *
!    $    ((bb * h ** 2) / (h * k ** 2 + h ** 2 * k)) + (bb * k ** 2) /
!    $    (h * k ** 2 + h ** 2 * k) + t1
! rewritten 2/96 frt

      denom1 = (h + k ) * h * k
      denom2 = (dpb + dpf) * dpf * dpb

      t1 = (-2.) * ((brcm * h)      / denom1)&
         + (-2.) * ((brcm * k)      / denom1)&
         +         (cc * h**2)      / denom1&
         + (-1.) * ((cc * k**2)     / denom1)&
         + (-1.) * ((dpb**2 * ee)   / denom2)&
         +          (dpf**2 * ee)   / denom2

      eij = (-1.) * ((aa * dpb ** 2) / denom2)&
          +          (aa * dpf ** 2) / denom2&
          + (-2.) * ((arcm * dpb)    / denom2)&
          + (-2.) * ((arcm * dpf)    / denom2)&
          + (-1.) * ((bb * h ** 2)   / denom1)&
          +          (bb * k ** 2)   / denom1&
          + t1

      d(i,j) = -eij


! have to divide by d (the multiplier of v(i,j) [e] in n.r. sor)
      c(1,i,j) = aij/d(i,j)
      c(2,i,j) = bij/d(i,j)
      c(3,i,j) = cij/d(i,j)
      c(4,i,j) = dij/d(i,j)

      return
      end
!
! *******************************************************************
!                  c i r c l e
! *******************************************************************
!     subroutine circle(r,imax,jmax,jwrap)
      subroutine circular(r,imax,jmax,jwrap)  ! rename circle to circular
!                                             ! MCF March 13, 2002
!
!  author: r.w.spiro               last update: 2-1-84
!
!***********************************************************************
!  subroutine to circularize matrix r
!***********************************************************************
!
!
!  description of input parameters:
!         r=matrix to be circularized
!         imax= i dimension of r
!         jmax=j dimension of r
!         jwrap=wrap around size
!***********************************************************************
!
!
      dimension r(imax,jmax)
      jlast=jmax-jwrap
      do  i=1,imax
         do  j=1,jwrap-1
             r(i,j)=r(i,jlast+j)
         enddo
         r(i,jmax)=r(i,jwrap)
      enddo
      return
      end


!*****************************************************************************
      subroutine fok2rcm_grd(xlati,xmlt,hlosscone,sini_f)
!*****************************************************************************
!
! routine to convert from fok grid variables (xlati=latitude in radians,
!  xmlt=local time in hours) to rcm grid variables
!
!  rws      06-02-99
!
      use rcmgrid
!
      dimension xlati(idim,jdimf),sini_f(idim,jdimf)
      dimension xmlt(jdimf)
 
      pi=atan2(0.,-1.)
      a2=0.
      b2=0.
      aval=0.
      jwrap=3
      imax=idim
      jmax=jdim
      ii1=2
      ii2=imax-1
      jj1=jwrap
      jj2=jmax-1
      offset=0.
      dlam=1./float(idim-1)
      dpsi=(2.*pi)/float(jdim-jwrap)
      re=6375.
      ri=re+hlosscone
 
      ic=0
      do i=idim,1,-1
         ic=ic+1
!
         do j=1,jdim
            jc=j-jdif
            if (jc.lt.1) jc=jc+jdimf
            if (jc.gt.jdimf) jc=jc-jdimf
            theta=pi/2.-xlati(i,jc)
            stheta=sin(theta)
            sinit=2.*cos(theta)/sqrt(1.+3.*cos(theta)**2)
            beta(ic,j)=stheta
   !        sini(ic,j)=sinit
            sini(ic,j)=sini_f(i,jc)
            colat(ic,j)=theta
            aloct(ic,j)=xmlt(jc)*pi/12.-pi
            if(aloct(ic,j).lt.0.) then
               aloct(ic,j)=aloct(ic,j)+2.*pi
            else if(aloct(ic,j).ge.2.*pi) then
               aloct(ic,j)=aloct(ic,j)-2.*pi
            end if
            vcorot(ic,j)=-92400.*(1.-cos(theta)**2)
            bir(ic,j)=62000.*cos(theta)
         end do
      end do
!
      do i=1,idim
         do j=1,jdim
            if(i.eq.1) then
               alpha(i,j)=(2.*colat(2,j)-1.5*colat(1,j)-&
                           .5*colat(3,j))/dlam
            else if(i.eq.idim) then
               alpha(i,j)=(-2.*colat(idim-1,j)+1.5*colat(idim,j)+&
                            .5*colat(idim-2,j))/dlam
            else
               alpha(i,j)=(colat(i+1,j)-colat(i-1,j))/(2.*dlam)
            end if
         end do
      end do
!
      do j=1,jdim
         ibndy(j)=1
      end do
!
! write grid files
      open(unit=10,file='rcmcrd11')
      open(unit=11,file='rcmcrd21')
!
! MCF
      id=0
      write(11,*) idim,jdim
! MCF end
      write(10,*) id,a2,b2,aval,offset,dlam,dpsi,ri,re
! MCF
!     write(11,'(1x,4e15.8)') colat,aloct,sini
      write(11,'(1x,4e15.8)') colat
      write(11,'(1x,4e15.8)') aloct
      write(11,'(1x,4e15.8)') sini
! MCF end
      write(11,*) alpha
      write(11,*) beta
      write(11,*) vcorot
      write(11,*) bir
!
      close(unit=10)
      close(unit=11)
!
!     write(6,*) 'i,alpha,beta,colat,vcorot,bir,sini'
      do i=1,idim
!        write(6,*) i,alpha(i,27),beta(i,27),colat(i,27),
!    2    vcorot(i,27),bir(i,27),sini(i,27)
      end do
!
!     call outp(colat,idim,jdim,1,idim,1,1,jdim,1,0.,'colat',6,80)
!     call outp(aloct,idim,jdim,1,idim,1,1,jdim,1,0.,'aloct',6,80)
!
      return
      end
!
!
! *******************************************************************
!                 e q b n d y
! *******************************************************************
       subroutine eqbndy
!
!***********************************************************************
!
!  last update| 08-27-86                 written by| g.a.mantjoukis
!               02-21-02   r.w. spiro, mofified for inclusion in
!                          m.-c. fok version of rcm
!                      (1) information passing via modules
!                      (2) neutral wind arrays removed, neglect effect
!                          of neutral winds
!
!  subroutine to compute equatorial bndy condition.based on
!   mantjou.eqbndy.text
!
!  current conservation requires that at i=imax,
!  v(imax,j)=c(1,imax,j)*v(imax+1,j)+c(2,imax,j)*v(imax-1,j)+
!            c(3,imax,j)*v(imax,j+1)+c(4,imax,j)*v(imax,j-1)+c(5,imax,j)
!  where v(imax+1,j) is outside the modeling region.
!
!  the equatorial bndy condition gives an expression for v(imax+1,j) in
!  terms of quantities inside the modeling region and is of the form
!  v(imax+1,j)=ceq1*v(imax,j-1)+ceq2*v(imax,j)+ceq3*v(imax,j+1)+
!              ceq4*v(imax-1,j)+ceq5
!  where ceq1 through ceq5 are calculated below.
!
!  ss(j) is a cowling-type conductance (see mantjou.eqbndy.text)
!       integrated over the cross-section of the equatorial band,at
!       any given local time.
!
!
!  to set bnd cond to no current across imax (ie., no eq electrojet)
!   explicityly zero ss(j) for all j
!
!
!***********************************************************************
!
       use rcmgrid
       use conductance
       use coefficients
!
!
       bnorm=1.e-6
!
!  zero ss(j) for special condition of no current across imax
!     (ie., no equatorial electrojet)
!
!     do 10 j=1,jmax
!        ss(j)=0.0
!  10  continue

       do 20 j=jj1,jj2
          cf=alpha(imax,j)*dlam/beta(imax,j)/dpsi/pedlam(imax,j)
! comment by yong on 7/26/90. introduce 'signbe' for hall
          ceq1=cf*(signbe*hall(imax,j)-0.5*(ss(j+1)-4.0*ss(j)-ss(j-1))/&
          dpsi)
          ceq2=-4.*cf*ss(j)/dpsi
          ceq3=cf*(-signbe*hall(imax,j)+0.5*(ss(j+1)+4.0*ss(j)-ss(j-1))/&
          dpsi)
! end of change
          ceq4=1.0
          ceq5=0.
!
          den=1.-ceq2*c(1,imax,j)

          c(5,imax,j)=(c(5,imax,j)+ceq5*c(1,imax,j))/den
          c(4,imax,j)=(c(4,imax,j)+ceq1*c(1,imax,j))/den
          c(3,imax,j)=(c(3,imax,j)+ceq3*c(1,imax,j))/den
          c(2,imax,j)=(c(2,imax,j)+ceq4*c(1,imax,j))/den
          c(1,imax,j)=0.
   20  continue
!
       do  n=1,ncoeff
          do  j=1,jwrap-1
             c(n,imax,j)=c(n,imax,jmax-jwrap+j)
          enddo
          c(n,imax,jmax)=c(n,imax,jwrap)
       enddo   

       return
       end
