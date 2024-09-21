!
!                                sat_flux.f90
! 
! Program to extract energetic particle flux, cold plasma density, Lshell and
! electric field along satellite path from CIMI output in *.fls and *.rtp files
! 
! Input files: sat_flux.dat
!              satellite_sp.info, storm_satellite.pos
!              outname.rtp, outname_sp.fls, outname.pot
! Output file: outname_satellite_sp.flux
! 
! Created on August 2, 2022
!-------------------------------------------------------------------------------

  implicit none
  integer ntime,jday0,ir,ip,je,ig,nr,nth,nphi,ie,ipa,ipos,jyr,jday,jhour,jmin
  integer i,k,m,icall
  real,parameter :: pi=3.14159265
  real Lshell,tpos,xsm,ysm,zsm,plsDen,Efield,Amax,x1,x2,dmusum
  character storm*8,outname*8,satellite*5,species*2,header*80
  real,allocatable,dimension(:) :: EeV,energy,PitchA,sinPitchA,dmu,fluxPA
  real,allocatable,dimension(:,:) :: fluxk

! Read data in sat_flux.dat
  open(unit=1,file='sat_flux.dat',status='old')
  read(1,'(a8)') storm
  read(1,'(a8)') outname
  read(1,'(a5)') satellite
  read(1,'(a2)') species  
  read(1,*) jday0             ! day number corresponding to t=0
  read(1,*) ntime             ! no. of time in the *.fls file
  read(1,*) ir,ip,je,ig       ! dimension of data in *.fls file
  read(1,*) nr,nth,nphi       ! dimension of data in *.rtp file
  close(1)

! Read data in satellite_sp.info
  open(unit=2,file=satellite//species//'.info',status='old')
  read(2,*) ie,ipa       ! number of energy and pitch angle bins
  allocate (EeV(ie),energy(ie),PitchA(ipa),fluxPA(ie),fluxk(ie,ipa))
  allocate (sinPitchA(ipa),dmu(ipa))
  read(2,*) EeV    
  read(2,*) PitchA
  close(2)
  energy(:)=EeV(:)*1.e-3    ! energy in keV

! find sin(PitchA) and dmu of the satellite pitch angle grid
  Amax=pi
  if (pitchA(ipa).lt.90.) Amax=pi/2.
  do m=1,ipa
     sinPitchA(m)=sind(PitchA(m))
     if (m.eq.1) x1=0.
     if (m.gt.1) x1=0.5*(pitchA(m)+pitchA(m-1))*pi/180.
     if (m.eq.ipa) x2=Amax
     if (m.lt.ipa) x2=0.5*(pitchA(m)+pitchA(m+1))*pi/180.
     dmu(m)=sinPitchA(m)*(x2-x1)
  enddo
  dmuSum=sum(dmu)
 
! First call of find_CIMIpara to map cimi output on a fixed equatorial grid
  icall=1
  write(*,*) 'icall ',icall
  call find_CIMIpara(icall,outname,species,ntime,ir,ip,je,ig,nr,nth,nphi, &
         ie,ipa,tpos,xsm,ysm,zsm,energy,sinPitchA,Lshell,fluxk,plsDen,Efield)

! Read ipos from the storm_satellite.pos file
  open(unit=4,file=storm//'_'//satellite//'.pos',status='old')
  read(4,*) ipos 
  read(4,'(a80)') header

! Open file to write output
  open(unit=3,file=outname//'_'//satellite//species//'.flux')
  write(3,*) ipos,ie,'         ! ipos, ie '
  write(3,*) 'energy bins (keV)'
  write(3,'(7f11.3)') energy
  write(3,'(a)') ' iyr idy ihr min Lshell plsDen(m^-3) Efield(V/m)' &
                    //' flux (cm2 s sr keV)^-1'

! Start extracting flux and etc along satellite path
  do i=1,ipos
     read(4,*) jyr,jday,jhour,jmin,xsm,ysm,zsm       ! xsm,ysm,zsm in RE
     tpos=(jday-jday0)*86400.+jhour*3600.+jmin*60.   ! in second
     icall=icall+1
     write(*,*) 'icall ',icall
     call find_CIMIpara(icall,outname,species,ntime,ir,ip,je,ig,nr,nth,nphi, &
          ie,ipa,tpos,xsm,ysm,zsm,energy,sinPitchA,Lshell,fluxk,plsDen,Efield)
     ! find the pitch-angle flux, fluxPA
     do k=1,ie
        fluxPA(k)=0.
        do m=1,ipa
           fluxPA(k)=fluxPA(k)+fluxk(k,m)*dmu(m)
        enddo
        fluxPA(k)=fluxPA(k)/dmusum
     enddo
     ! write results
     write(3,'(4i4,f7.2,1p,2e12.3)') jyr,jday,jhour,jmin,Lshell,plsDen,Efield
     write(3,'(1p,7e11.3)') fluxPA
  enddo      ! end of do i=1,ipos
  close(4)
  close(3)

  end


!-------------------------------------------------------------------------------
  subroutine find_CIMIpara(icall,outname,species,ntime,ir,ip,je,ig,nr,nth,nphi,&
           ie,ipa,tpos,xsm,ysm,zsm,energy,sinPitchA,Lshell,fluxk,plsDen,Efield)
!-------------------------------------------------------------------------------
! Routine extracts CIMI output at satellite time (tpos) and 
! location (xsm,ysm,zsm)
! 
! inputs: icall,outname,species,ntime,ir,ip,je,ig,nr,nth,nphi,
!         ie,ipa,tpos,xsm,ysm,zsm,energy,sinPitchA
! output: Lshell,fluxk,plsDen,Efield

  implicit none
  real,parameter :: pi=3.14159265,re_m=6.3712e6
  integer i,j,k,m,ntime,ir,ip,je,ig,nr,nth,nphi,ie,ipa,icall,irh,indx(ip)
  integer je1,i0,i1,j0,j1,it1,it2,nn,ni,n,idummy(ip),ik
  real tpos,rsm,xsm,ysm,zsm,energy(ie),Lshell,fluxk(ie,ipa),plsDen,Efield
  real sinPitchA(ipa),dmu(ipa),tflux(ntime),fj(ip+1),var,fac1,fac2,dtf
  real Rrtp(nr),Trtp(nth),Prtp(nphi),theta,phi
  real BXrtp(nr,nth,nphi),BXrtpa(ntime,nr,nth,nphi)
  real BYrtp(nr,nth,nphi),BYrtpa(ntime,nr,nth,nphi)
  real BZrtp(nr,nth,nphi),BZrtpa(ntime,nr,nth,nphi)
  real Rortp(nr,nth,nphi),Rortpa(ntime,nr,nth,nphi)
  real MLTrtp(nr,nth,nphi),MLTrtpa(ntime,nr,nth,nphi)
  real ren(nr),xmltn(nphi),pot1(ir,ip),rn(ir,ip),rn1(ir)
  real gride(je),gridy(ig),ro(ir,ip),xmlto(ir,ip),bo1(ir,ip),den1(ir,ip)
  real flx1(ir,ip,je,ig),hour,dphi2,dummy1D(100),dummy2D(ir,ip),xmlt(ip)
  real flx2(ir,nphi,je,ig),den2(ir,nphi),bo2(ir,nphi),pot2(ir,nphi)
  real fr(ir),reni,flx3(ntime,nr,nphi,je,ig),den3(ntime,nr,nphi)
  real bo3(ntime,nr,nphi),pot3(nr,nphi),Er,Ep,Efd3(ntime,nr,nphi)
  real bo4(nr,nphi),bo5,BsatBeq,Lshell5(2),Efd4(nr,nphi),Efd5(2),den4(nr,nphi)
  real den5(2),flx4(nr,nphi,ig),flx5(2,je,ipa),sinAo,flxk(je)
  real xmlt1(ip+1),rj(ip+1),xmltj,Bsat,BXsat,BYsat,BZsat,Rosat,MLTosat
  character header*80,outname*8,species*2
  real,allocatable,dimension(:) :: xlath
  real,allocatable,dimension(:,:) :: potenth
  real,allocatable,dimension(:,:,:) :: dummy3D

! Read data in *.rtp and setup an equatorial grid when icall=1
  if (icall.eq.1) then
     ! read data from *.rtp
     open(unit=1,file=outname//'.rtp',status='old')
     read(1,'(a80)') header
     read(1,*) Rrtp            ! R in RE
     read(1,*) Trtp            ! theta in radian
     read(1,*) Prtp            ! phi in radian
     do n=1,ntime
        read(1,*) hour
        read(1,*) BXrtp
        read(1,*) BYrtp
        read(1,*) BZrtp
        read(1,*) Rortp
        read(1,*) MLTrtp
        tflux(n)=hour*3600.
        BXrtpa(n,1:nr,1:nth,1:nphi)=BXrtp(1:nr,1:nth,1:nphi)
        BYrtpa(n,1:nr,1:nth,1:nphi)=BYrtp(1:nr,1:nth,1:nphi)
        BZrtpa(n,1:nr,1:nth,1:nphi)=BZrtp(1:nr,1:nth,1:nphi)
        Rortpa(n,1:nr,1:nth,1:nphi)=Rortp(1:nr,1:nth,1:nphi)
        MLTrtpa(n,1:nr,1:nth,1:nphi)=MLTrtp(1:nr,1:nth,1:nphi)
     enddo
     close(1)
     ! Setup a fixed equatorial grid (ren,xmltn)
     ren(1:nr)=Rrtp(1:nr)
     do j=1,nphi
        xmltn(j)=Prtp(j)*12./pi
     enddo
     dphi2=Prtp(3)-Prtp(1)
  endif

! Read flux, density,Bo and potential and map to grid (ren,xmltn) if icall=1
  if (icall.eq.1) then   ! 2nd if (icall.eq.1)
     open(unit=2,file=outname//species//'.fls',status='old')
     read(2,'(a80)') header  
     read(2,*) dummy1D(1:ir)  ! make sure demension of dummy1D > ir and ip
     read(2,*) dummy1D(1:ip)  !
     read(2,*) gride(1:je)
     read(2,*) gridy(1:ig)
     open(unit=3,file=outname//'.pot',status='old')
     read(3,'(a80)') header  
     read(3,'(a80)') header  
     read(3,*) irh
     read(3,'(a80)') header  
     read(3,*) ik  
     allocate (potenth(irh,ip),xlath(irh),dummy3D(ir,ip,ik))
     read(3,*) dummy2D
     read(3,*) xlath
     read(3,*) dummy2D
     read(3,*) dummy2D

     do n=1,ntime
        ! read flux, density and Bo from *.fls
        read(2,*) hour
        write(*,*) 'in .fls, hour => ',hour
        do i=1,ir
           do j=1,ip
              read(2,*) dummy1D(1:4),ro(i,j),xmlto(i,j),dummy1D(1:2),bo1(i,j)
              read(2,*) den1(i,j)
              do k=1,je
                 read(2,*) flx1(i,j,k,1:ig)
              enddo
           enddo
        enddo

        ! read potentials from *.pot
        read(3,*) hour
        read(3,*) idummy
        read(3,*) xmlt
        do m=1,5
           read(3,*) dummy2D
        enddo
        read(3,*) dummy3D
        read(3,'(8f10.1)') pot1
        read(3,'(8f10.1)') potenth
        read(3,*) dummy2D
        read(3,'(1p,6e13.4)') dummy2D

        ! Setup grid (rn,xmltn) and map data on it
        do i=1,ir    ! map data on (rn,xmltn)     
           ! find new equatorial grid, rn(i,j) at fixed xmltn
           call indexx(ip,xmlto(i,:),indx)
           do j=1,ip                  ! rearrange xmlt in ascending order
              xmlt1(j)=xmlto(i,indx(j))
              rj(j)=ro(i,indx(j))
           enddo
           xmlt1(ip+1)=xmlt1(1)+24.
           rj(ip+1)=rj(1)
           do j=1,nphi
              xmltj=xmltn(j)
              if (xmltj.lt.xmlt1(1)) xmltj=xmltj+24.
              call lintp(xmlt1,rj,ip+1,xmltj,rn(i,j))
              if (i.gt.1.and.rn(i,j).lt.rn(i-1,j)) rn(i,j)=rn(i-1,j)
           enddo
           ! find energetic flux, plsp den, B & potential at rn(i,j),xmltn(j)
           do m=1,ig+3
              je1=je
              if (m.gt.ig) je1=1
              do k=1,je1
                 do j=1,ip
                    if (m.le.ig) fj(j)=flx1(i,indx(j),k,m)
                    if (m.eq.ig+1) fj(j)=den1(i,indx(j))
                    if (m.eq.ig+2) fj(j)=bo1(i,indx(j))
                    if (m.eq.ig+3) fj(j)=pot1(i,indx(j))
                 enddo
                 fj(ip+1)=fj(1)
                 do j=1,nphi
                    xmltj=xmltn(j)
                    if (xmltj.lt.xmlt1(1)) xmltj=xmltj+24.
                    call lintp(xmlt1,fj,ip+1,xmltj,var)
                    if (m.le.ig) flx2(i,j,k,m)=var
                    if (m.eq.ig+1) den2(i,j)=var
                    if (m.eq.ig+2) bo2(i,j)=var
                    if (m.eq.ig+3) pot2(i,j)=var
                 enddo
              enddo
           enddo   ! end of do m=1,ig+3
        enddo      ! end of do i=1,ir    ! map data on (rn,xmltn)

        ! map energetic flux, plsp den, B & E field on ren(i,j),xmltn(j)
        do j=1,nphi
           rn1(:)=rn(:,j)
           do m=1,ig+3
              je1=je
              if (m.gt.ig) je1=1
              do k=1,je1
                 do i=1,ir
                    if (m.le.ig) fr(i)=flx2(i,j,k,m)
                    if (m.eq.ig+1) fr(i)=den2(i,j)
                    if (m.eq.ig+2) fr(i)=bo2(i,j)
                    if (m.eq.ig+3) fr(i)=pot2(i,j)
                 enddo
                 do i=1,nr
                    reni=ren(i)
                    if (reni.lt.rn1(1)) reni=rn1(1)
                    var=0.
                    if (reni.le.rn1(ir)) call lintp(rn1,fr,ir,reni,var)
                    if (m.le.ig) flx3(n,i,j,k,m)=var
                    if (m.eq.ig+1) den3(n,i,j)=var
                    if (m.eq.ig+2) bo3(n,i,j)=var
                    if (m.eq.ig+3) pot3(i,j)=var
                 enddo
              enddo
           enddo   ! end of do m=1,ig+3
        enddo      ! end of do j=1,nphi
        ! find electric field at ren(1:nr),xmltn(1:nphi)
        do i=1,nr
           i0=i-1
           if (i0.lt.1) i0=1
           i1=i+1
           if (i1.gt.ir) i1=ir
           do j=1,nphi-1
              j0=j-1
              if (j0.lt.1) j0=nphi-1
              j1=j+1
              Er=-(pot3(i1,j)-pot3(i0,j))/(ren(i1)-ren(i0))
              Ep=-(pot3(i,j1)-pot3(i,j0))/dphi2/ren(i)
              Efd3(n,i,j)=sqrt(Er*Er+Ep*Ep)/re_m   ! Efield in Volt/m
           enddo
           Efd3(n,i,nphi)=Efd3(n,i,1)
        enddo
     enddo       ! end of do n=1,ntime
     close(3)
     close(2)
  endif          ! end of if (icall.eq.1) 2nd time

! Find L shell, energetic flux, plspDen, E field along satellite path if icall>1
  if (icall.gt.1) then
     rsm=sqrt(xsm*xsm+ysm*ysm+zsm*zsm)
     theta=acos(zsm/rsm)
     phi=atan2(ysm,xsm)
     if (phi.lt.0.) phi=phi+2.*pi     ! make phi from 0 to 2pi
     Lshell=99.       ! initial values
     Efield=0.        !
     plsDen=0.        !
     fluxk(:,:)=0.    !
     if (rsm.gt.ren(nr)) return                           ! outside CIMI domain
     if (theta.lt.Trtp(1).or.theta.gt.Trtp(nth)) return   ! outside CIMI domain
     call locate1(tflux,ntime,tpos,it1)
     if (it1.eq.0) it1=1
     it2=it1+1
     if (it2.gt.ntime) stop
     fac1=tflux(it2)-tpos
     fac2=tpos-tflux(it1)
     dtf=tflux(it2)-tflux(it1)
     do nn=1,2
        ni=it1+nn-1
        Lshell5(nn)=99.     ! initial values
        Efd5(nn)=0.         !
        den5(nn)=0.
        flx5(nn,:,:)=0.

        ! find satellite equatorial crossing point and B field there
        Rortp(1:nr,1:nth,1:nphi)=Rortpa(ni,1:nr,1:nth,1:nphi)
        MLTrtp(1:nr,1:nth,1:nphi)=MLTrtpa(ni,1:nr,1:nth,1:nphi)
        call lintp3(Rrtp,Trtp,Prtp,Rortp,nr,nth,nphi,rsm,theta,phi,Rosat)
        call lintp3(Rrtp,Trtp,Prtp,MLTrtp,nr,nth,nphi,rsm,theta,phi,MLTosat)
        if (Rosat.lt.ren(1)) Rosat=ren(1)
        if (Rosat.gt.ren(nr)) goto 999        ! outside CIMI domain
        Lshell5(nn)=Rosat      ! Lshell actually is the equatorial crossing
        bo4(1:nr,1:nphi)=bo3(ni,1:nr,1:nphi)
        call lintp2(ren,xmltn,bo4,nr,nphi,Rosat,MLTosat,bo5)
        
        ! find the magnetic feild at the satellite location
        BXrtp(1:nr,1:nth,1:nphi)=BXrtpa(ni,1:nr,1:nth,1:nphi)
        BYrtp(1:nr,1:nth,1:nphi)=BYrtpa(ni,1:nr,1:nth,1:nphi)
        BZrtp(1:nr,1:nth,1:nphi)=BZrtpa(ni,1:nr,1:nth,1:nphi)
        call lintp3(Rrtp,Trtp,Prtp,BXrtp,nr,nth,nphi,rsm,theta,phi,BXsat)
        call lintp3(Rrtp,Trtp,Prtp,BYrtp,nr,nth,nphi,rsm,theta,phi,BYsat)
        call lintp3(Rrtp,Trtp,Prtp,BZrtp,nr,nth,nphi,rsm,theta,phi,BZsat)
        Bsat=sqrt(BXsat*BXsat+BYsat*BYsat+BZsat*Bzsat)
        BsatBeq=sqrt(Bsat/bo5)

        ! Find electric field, plasmasphere density, flux at satellite location
        Efd4(1:nr,1:nphi)=Efd3(ni,1:nr,1:nphi)
        call lintp2(ren,xmltn,Efd4,nr,nphi,Rosat,MLTosat,var)
        Efd5(nn)=var*BsatBeq
        den4(1:nr,1:nphi)=den3(ni,1:nr,1:nphi)
        call lintp2(ren,xmltn,den4,nr,nphi,Rosat,MLTosat,var)
        den5(nn)=var
        do k=1,je
           flx4(1:nr,1:nphi,1:ig)=flx3(ni,1:nr,1:nphi,k,1:ig)
           do m=1,ipa
              flx5(nn,k,m)=0.
              sinAo=sinPitchA(m)/BsatBeq
              if (sinAo.gt.gridy(ig)) sinAo=gridy(ig)
              var=0.
              if (sinAo.ge.gridy(1)) call lintp3(ren,xmltn,gridy,flx4,nr,nphi, &
                                           ig,Rosat,MLTosat,sinAo,var)
              flx5(nn,k,m)=var   
           enddo
        enddo
 999    continue
     enddo          ! end of do nn=1,2

     ! Find Lshell, electric field, plasmasphere density, flux at tpos
     Lshell=(Lshell5(1)*fac1+Lshell5(2)*fac2)/dtf
     Efield=(Efd5(1)*fac1+Efd5(2)*fac2)/dtf
     plsDen=(den5(1)*fac1+den5(2)*fac2)/dtf
     fluxk(:,:)=0.    ! initial value
     do m=1,ipa
        do k=1,je
           flxk(k)=(flx5(1,k,m)*fac1+flx5(2,k,m)*fac2)/dtf
        enddo
        do k=1,ie
           if (energy(k).ge.gride(1).and.energy(k).le.gride(je)) & 
              call lintp(gride,flxk,je,energy(k),fluxk(k,m))
        enddo
     enddo
  endif             ! end of if (icall.gt.1)
  
  end subroutine find_CIMIpara


!-----------------------------------------------------------------------
      subroutine lintp(xx,yy,n,x,y)
!-----------------------------------------------------------------------
!  Routine does 1-D interpolation.  xx must be increasing or decreasing
!  monotonically. If x is beyound xx, it will be forced inside xx range.

      implicit none
      integer n,i,j,jl,ju,jm
      real xx(n),yy(n),x,x1,minxx,maxxx,y,d

!  Make sure xx is increasing or decreasing monotonically
      do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            write(*,*) ' lintp: xx is not increasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            write(*,*) ' lintp: xx is not decreasing monotonically '
            write(*,*) n,(xx(j),j=1,n)
            stop
          endif
      enddo

!  Make sure x is inside xx range
      minxx=minval(xx)
      maxxx=maxval(xx)
      x1=x
      if (x1.lt.minxx) x1=minxx
      if (x1.gt.maxxx) x1=maxxx

!    initialize lower and upper values
!
      jl=1
      ju=n
!
!    if not dne compute a midpoint
!
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
!
!    now replace lower or upper limit
!
        if((xx(n).gt.xx(1)).eqv.(x1.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
!
!    try again
!
      go to 10
      endif
!
!    this is j
!
      j=jl      ! if x.le.xx(1) then j=1
!                 if x.gt.xx(j).and.x.le.xx(j+1) then j=j
!                 if x.gt.xx(n) then j=n-1
      d=xx(j+1)-xx(j)
      y=(yy(j)*(xx(j+1)-x1)+yy(j+1)*(x1-xx(j)))/d

      end subroutine lintp


!-------------------------------------------------------------------------------
        subroutine lintp2(x,y,v,nx,ny,x1,y1,v1)
!-------------------------------------------------------------------------------
!  Routine does 2-D interpolation.  x and y must be increasing or decreasing
!  monotonically
!
        real x(nx),y(ny),v(nx,ny)

        call locate1(x,nx,x1,i)
        if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
        if (i.lt.1) i=1              ! extrapolation if out of range
        i1=i+1
        a=(x1-x(i))/(x(i1)-x(i))

        call locate1(y,ny,y1,j)
        if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
        if (j.lt.1) j=1              ! extrapolation if out of range
        j1=j+1
        b=(y1-y(j))/(y(j1)-y(j))

        q00=(1.-a)*(1.-b)
        q01=(1.-a)*b
        q10=a*(1.-b)
        q11=a*b
        v1=q00*v(i,j)+q01*v(i,j1)+q10*v(i1,j)+q11*v(i1,j1)

        end subroutine lintp2


!-------------------------------------------------------------------------------
        subroutine lintp3(x,y,z,v,nx,ny,nz,x1,y1,z1,v1)
!-------------------------------------------------------------------------------
!  This sub program takes 3-d interplation
!
      real x(nx),y(ny),z(nz),v(nx,ny,nz)

      call locate1(x,nx,x1,i)
      if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
      if (i.lt.1) i=1              ! extrapolation if out of range
      i1=i+1
      a=(x1-x(i))/(x(i1)-x(i))

      call locate1(y,ny,y1,j)
      if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
      if (j.lt.1) j=1              ! extrapolation if out of range
      j1=j+1
      b=(y1-y(j))/(y(j1)-y(j))

      call locate1(z,nz,z1,k)
      if (k.gt.(nz-1)) k=nz-1      ! extrapolation if out of range
      if (k.lt.1) k=1              ! extrapolation if out of range
      k1=k+1
      c=(z1-z(k))/(z(k1)-z(k))

      q000=(1.-a)*(1.-b)*(1.-c)*v(i,j,k)
      q001=(1.-a)*(1.-b)*c*v(i,j,k1)
      q010=(1.-a)*b*(1.-c)*v(i,j1,k)
      q011=(1.-a)*b*c*v(i,j1,k1)
      q100=a*(1.-b)*(1.-c)*v(i1,j,k)
      q101=a*(1.-b)*c*v(i1,j,k1)
      q110=a*b*(1.-c)*v(i1,j1,k)
      q111=a*b*c*v(i1,j1,k1)

      v1=q000+q001+q010+q011+q100+q101+q110+q111

      end subroutine lintp3 


!--------------------------------------------------------------------------
      subroutine locate1(xx,n,x,j)
!--------------------------------------------------------------------------
!  Routine return a value of j such that x is between xx(j) and xx(j+1).
!  xx must be increasing or decreasing monotonically. If not, the locate will
!  stop at the turning point.
!  If xx is increasing:
!     If x=xx(m), j=m-1 so if x=xx(1), j=0  and if x=xx(n), j=n-1
!     If x < xx(1), j=0  and if x > xx(n), j=n
!  If xx is decreasing:
!     If x=xx(m), j=m so if x=xx(1), j=1  and if x=xx(n), j=n
!     If x > xx(1), j=0  and if x < xx(n), j=n

      real xx(n)

!  Make sure xx is increasing or decreasing monotonically
      nn=n
      monoCheck: do i=2,n
         if (xx(n).gt.xx(1).and.xx(i).lt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
         if (xx(n).lt.xx(1).and.xx(i).gt.xx(i-1)) then
            nn=i-1
            exit monoCheck
         endif
      enddo monoCheck
      if (nn.ne.n) then
         write(*,*)'locate1: xx is not increasing or decreasing monotonically'
         write(*,*)'n,x ',n,x
         write(*,*)'xx ',xx
         stop
      endif

      jl=0
      ju=nn+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(nn).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
      j=jl

      end subroutine locate1


!--------------------------------------------------------------------------
      subroutine indexx(n,arrin,indx)
!--------------------------------------------------------------------------
!  Routine arranges an array in ascending order.

      real arrin(n)
      integer indx(n)
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      irj=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(irj)
          q=arrin(indxt)
          indx(irj)=indx(1)
          irj=irj-1
          if(irj.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.irj)then
          if(j.lt.irj)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=irj+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10

      end subroutine indexx
