module ModMate
 
 implicit none

 public :: MateFileName,&
           TimeMate,&
           iDoyStartMate,&
           read_mate_data,&
           unit_test_mate,&
           mate_geocorona


 private ! except

 character(len=80) :: MateFileName='MATE_geocorona.dat'

 real,allocatable,dimension(:) :: TimeMate   ! time in sec in MATE data
 !real,allocatable,dimension(:,:,:) :: &
 !                     xSMMate,ySMMate,zSMMate ! X,Y,Z in sm in RE
 real,allocatable,dimension(:) :: &
                      RMateGSE,LatMateGSE,LonMateGSE ! r,Lat,Lon in GSE in RE,deg,deg

 ! Mate grid info.
 integer nHourMate,nRMate,nLatMate,nLonMate,iDoyStartMate

 !real,allocatable :: HdenMate(:,:,:,:)  ! Mate exospheric H density (time,x,y,z)
 real,allocatable :: logHdenMate(:,:,:,:)  ! Mate exospheric log H density (time,x,y,z)

contains
!--------------------------------------------------------------------------
 subroutine read_mate_data
! This subroutine reads MATE exospheric H density
! 
! NOTE: geopack_2008 is used for coordinate transformation.
!       subroutine RECALC_08 needs to be called prior to this subroutine.
!--------------------------------------------------------------------------
 integer yr,doy,hr
 !real r1,LatGSE1,LonGSE1,xGSE,yGSE,zGSE,xGSW,yGSW,zGSW,xSM,ySM,zSM
 real r1,LatGSE1,LonGSE1,Hden1
 !external gswgse_08,smgsw_08  ! geopack_2008 subroutine

 real :: cDegToRad=0.01745329251  ! degree to rad
 real :: cDaytoSec=86400.
 real :: cHourToSec=3600.

 character(len=80) header 
 integer ihour,ir,iLat,iLon

 open(unit=4,file=trim(MateFileName),status='old')
 read(4,*) nRMate,nLatMate,nLonMate
 read(4,*) nHourMate
 read(4,*) header
 if (.not.allocated(TimeMate)) allocate(TimeMate(nHourMate))
 !if (.not.allocated(xSmMate)) allocate(xSmMate(nRMate,nLatMate,nLonMate))
 !if (.not.allocated(ySmMate)) allocate(ySmMate(nRMate,nLatMate,nLonMate))
 !if (.not.allocated(zSmMate)) allocate(zSmMate(nRMate,nLatMate,nLonMate))
 if (.not.allocated(RMateGSE)) allocate(RMateGSE(nRMate))
 if (.not.allocated(LatMateGSE)) allocate(LatMateGSE(nLatMate))
 if (.not.allocated(LonMateGSE)) allocate(LonMateGSE(nLonMate))
 !if (.not.allocated(HdenMate)) allocate(HdenMate(nHourMate,nRMate,nLatMate,nLonMate))
 if (.not.allocated(logHdenMate)) &
    allocate(logHdenMate(nHourMate,nRMate,nLatMate,nLonMate))
 do ihour=1,nHourMate
    do ir=1,nRMate
       do iLat=1,nLatMate
          do iLon=1,nLonMate
             !read(4,*) yr,doy,hr,r1,LatGSE1,LonGSE1,HdenMate(ihour,ir,iLat,iLon)
             read(4,*) yr,doy,hr,r1,LatGSE1,LonGSE1,Hden1
             logHdenMate(ihour,ir,iLat,iLon)=log10(Hden1)
             !if (ihour==1) then 
             !   if (iR==1.and.iLat==1.and.iLon==1) iDoyStartMate=doy
             !   xGSE=r1*cos(LatGSE1*cDegToRad)*cos(LonGSE1*cDegToRad)
             !   yGSE=r1*cos(LatGSE1*cDegToRad)*sin(LonGSE1*cDegToRad)
             !   zGSE=r1*sin(LatGSE1*cDegToRad)
             !   call gswgse_08(xGSW,yGSW,zGSW,xGSE,yGSE,zGSE,-1)
             !   call smgsw_08(xSM,ySM,zSM,xGSW,yGSW,zGSW,-1)
             !   xSMMate(iR,iLat,iLon)=xSM
             !   ySMMate(iR,iLat,iLon)=ySM
             !   zSMMate(iR,iLat,iLon)=zSM
             !endif
             if (ihour==1) then 
                if (iR==1.and.iLat==1.and.iLon==1) iDoyStartMate=doy
                RMateGSE(ir)=r1
                LatMateGSE(iLat)=LatGSE1
                LonMateGSE(iLon)=LonGSE1
             endif
             if (iR==1.and.iLat==1.and.iLon==1) &
                TimeMate(ihour)=(doy-iDoyStartMate)*cDayToSec + hr*cHourToSec
          enddo
       enddo
    enddo

 enddo
 close(4)


 end subroutine read_mate_data


!--------------------------------------------------------------------------
 subroutine mate_geocorona(t,xsm,ysm,zsm,hden)
!--------------------------------------------------------------------------
  real,intent(in) :: t,&         ! time in sec
                     xsm,ysm,zsm ! in RE
  real,intent(out) :: hden       ! 

  real :: xgse,ygse,zgse,xgsw,ygsw,zgsw
  real :: Rgse,Latgse,Longse,Thetagse,LonGseRad
  real :: logHden1,logHden2,logHden  ! Hden at TimeMate(iTme) and TimeMate(iTime+1)
  real :: cRadToDeg=57.2957795131

  integer :: iTime     ! index of TimeMate,
                       ! where TimeMate(iTime) is less than t

  external SPHCAR_08,gswgse_08,smgsw_08  ! geopack_2008 subroutine

  ! SM to GSE
  call smgsw_08(xsm,ysm,zsm,xGSW,yGSW,zGSW,1)
  call gswgse_08(xGSW,yGSW,zGSW,xGSE,yGSE,zGSE,1)

write(*,*) ' GST X,Y,Z test: ',xGSE,yGSE,zGSE

  ! X,Y,Z to R,Lat,Lon in GSE
  call SPHCAR_08(Rgse,ThetaGse,LongseRad,xGSE,yGSE,zGSE,-1)
  LatGSE=90.-ThetaGse*cRadToDeg
  LonGSE=LonGseRad*cRadToDeg

  ! interpolation
  !hden=1.e-3
  call locate1(TimeMate,nHourMate,t,iTime)
  write(*,*) 'Time =',iTime
  write(*,*) ' t =',t
  write(*,*) ' TimeMate(1) =',TimeMate(1)
  call lintp3(RMateGSE,LatMateGSE,LonMateGSE,logHdenMate(iTime,:,:,:),&
              nRMate,nLatMate,nLonMate,&
              Rgse,LatGSE,LonGSE,logHden1)
  if (iTime<nHourMate) then
     call lintp3(RMateGSE,LatMateGSE,LonMateGSE,logHdenMate(iTime+1,:,:,:),&
                 nRMate,nLatMate,nLonMate,&
                 Rgse,LatGSE,LonGSE,logHden2)
     ! 1D time interpolation
     logHden=( logHden1*(TimeMate(iTime+1)-t) + logHden2*(t-TimeMate(iTime)) ) &
          / (TimeMate(iTime+1) - TimeMate(iTime)) 
  else
     logHden=logHden1
  endif
  Hden=10.**logHden


  !call lintp4(TimeMate,RMateGSE,LatMateGSE,LonMateGSE,logHdenMate,&
  !            nHourMate,nRMate,nLatMate,nLonMate,&
  !            t,Rgse,LatGSE,LonGSE,Hden)

  end subroutine mate_geocorona

!-------------------------------------------------------------------------------
  subroutine lintp3(x,y,z,v,nx,ny,nz,x1,y1,z1,v1)
!-------------------------------------------------------------------------------
!  This sub program takes 3-d interplation.
!
  integer,intent(in) ::  nx,ny,nz
  real,intent(in)  ::  x(nx),y(ny),z(nz),v(nx,ny,nz),&
                       x1,y1,z1
  real,intent(out) :: v1

  real :: a,b,c,&
          q000,q001,q010,q011,&
          q100,q101,q110,q111

  integer i,j,k,i1,j1,k1

  call locate1(x,nx,x1,i)
  if (i.gt.(nx-1)) i=nx-1      ! extrapolation if out of range
  if (i.lt.1) i=1              ! extrapolation if out of range

  call locate1(y,ny,y1,j)
  if (j.gt.(ny-1)) j=ny-1      ! extrapolation if out of range
  if (j.lt.1) j=1              ! extrapolation if out of range

  call locate1(z,nz,z1,k)
  if (k.gt.(nz-1)) k=nz-1      ! extrapolation if out of range
  if (k.lt.1) k=1              ! extrapolation if out of range


  i1=i+1
  j1=j+1
  k1=k+1
  a=(x1-x(i))/(x(i1)-x(i))
  b=(y1-y(j))/(y(j1)-y(j))
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
 
 
!-------------------------------------------------------------------------------
        subroutine lintp4(t,x,y,z,v,nt,nx,ny,nz,t1,x1,y1,z1,v1)
!-------------------------------------------------------------------------------
!  Routine does 2-D interpolation.  x and y must be increasing or decreasing
!  monotonically
!
        integer,intent(in) :: nt,nx,ny,nz
        real,intent(in):: t(nt),x(nx),y(ny),z(nz),v(nt,nx,ny,nz)
        real,intent(out) :: v1
        integer it,ix,iy,iz,it1,ix1,iy1,iz1
        real t1,x1,y1,z1,t2,x2,y2,z2,dt,dx,dy,dz
        real q(0:1,0:1,0:1,0:1)  ! 4-D matrix

        !,q00,q01,q10,q11
        real mint,maxt,minx,maxx,miny,maxy,minz,maxz

        minx=minval(t)
        maxx=maxval(t)
        t2=t1
        if (t2.lt.mint) t2=mint      ! force x2 inside the x range
        if (t2.gt.maxt) t2=maxt      !
        call locate1(t,nt,t2,it)
        if (it.gt.(nt-1)) it=nt-1      
        if (it.lt.1) it=1               
        it1=it+1
        dt=(t2-t(it))/(t(it1)-t(it))

        mint=minval(x)
        maxt=maxval(x)
        x2=x1
        if (x2.lt.minx) x2=minx      ! force x2 inside the x range
        if (x2.gt.maxx) x2=maxx      !
        call locate1(x,nx,x2,ix)
        if (ix.gt.(nx-1)) ix=nx-1      
        if (ix.lt.1) ix=1               
        ix1=ix+1
        dx=(x2-x(ix))/(x(ix1)-x(ix))

        miny=minval(y)
        maxy=maxval(y)
        y2=y1
        if (y2.lt.miny) y2=miny      ! force y2 inside the y range
        if (y2.gt.maxy) y2=maxy      !
        call locate1(y,ny,y2,iy)
        if (iy.gt.(ny-1)) iy=ny-1      
        if (iy.lt.1) iy=1               
        iy1=iy+1
        dy=(y2-y(iy))/(y(iy1)-y(iy))

        minz=minval(z)
        maxz=maxval(z)
        z2=z1
        if (z2.lt.minz) z2=minz      ! force z2 inside the z range
        if (z2.gt.maxz) z2=maxz      !
        call locate1(z,nz,z2,iz)
        if (iz.gt.(nz-1)) iz=nz-1      
        if (iz.lt.1) iz=1               
        iz1=iz+1
        dz=(z2-z(iz))/(z(iz1)-z(iz))


        q(0,0,0,0)=(1.-dt)*(1.-dx)*(1.-dy)*(1.-dz)

        q(0,1,0,0)=(1.-dt)*    dx *(1.-dy)*(1.-dz)
        q(0,0,1,0)=(1.-dt)*(1.-dx)*    dy *(1.-dz)
        q(0,0,0,1)=(1.-dt)*(1.-dx)*(1.-dy)*    dz 
        q(1,0,0,0)=    dt *(1.-dx)*(1.-dy)*(1.-dz)

        q(1,1,0,0)=    dt *    dx *(1.-dy)*(1.-dz)
        q(1,0,1,0)=    dt *(1.-dx)*    dy *(1.-dz)
        q(1,0,0,1)=    dt *(1.-dx)*(1.-dy)*    dz 
        q(0,1,1,0)=(1.-dt)*    dx *    dy *(1.-dz)
        q(0,1,0,1)=(1.-dt)*    dx *(1.-dy)*    dz 
        q(0,0,1,1)=(1.-dt)*(1.-dx)*    dy *    dz 

        q(0,1,1,1)=(1.-dt)*    dx *    dy *    dz 
        q(1,0,1,1)=    dt *(1.-dx)*    dy *    dz 
        q(1,1,0,1)=    dt *    dx *(1.-dy)*    dz 
        q(1,1,1,0)=    dt *    dx *    dy *(1.-dz)

        q(1,1,1,1)=    dt *    dx *    dy *    dz 

        !function coeff(

        v1=q(0,0,0,0)*v(it ,ix ,iy ,iz ) &

          +q(1,0,0,0)*v(it1,ix ,iy ,iz ) &
          +q(0,1,0,0)*v(it ,ix1,iy ,iz ) &
          +q(0,0,1,0)*v(it ,ix ,iy1,iz ) &
          +q(0,0,0,1)*v(it ,ix ,iy ,iz1) &

          +q(1,1,0,0)*v(it1,ix1,iy ,iz ) &
          +q(1,0,1,0)*v(it1,ix ,iy1,iz ) &
          +q(1,0,0,1)*v(it1,ix ,iy ,iz1) &
          +q(0,1,1,0)*v(it ,ix1,iy1,iz ) &
          +q(0,1,0,1)*v(it ,ix1,iy ,iz1) &
          +q(0,0,1,1)*v(it ,ix ,iy1,iz1) &

          +q(1,1,1,0)*v(it1,ix1,iy1,iz ) &
          +q(1,1,0,1)*v(it1,ix1,iy ,iz1) &
          +q(1,0,1,1)*v(it1,ix ,iy1,iz1) &
          +q(0,1,1,1)*v(it ,ix1,iy1,iz1) &

          +q(1,1,1,1)*v(it1,ix1,iy1,iz1)

        end subroutine lintp4

!--------------------------------------------------------------------------
      subroutine locate1(xx,n,x,j)
!--------------------------------------------------------------------------
!  This subroutine is copied from cimi.f90

      real,intent(in) :: xx(n),x
      integer,intent(in) :: n
      integer,intent(out) :: j

      integer i,jl,ju,jm,nn

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

      if (j==0) j=1
      if (j==n) j=n

      end subroutine locate1


!-------------------------------------------------------------------------------
  subroutine unit_test_mate
!-------------------------------------------------------------------------------
 integer :: iRMate=1,&
            iLatMate=7,&
            iLonMate=1,&
            iTime
 real :: cSecToHour=1/3600.
 real :: xsm,ysm,zsm
 real :: xgse,ygse,zgse,xgsw,ygsw,zgsw
 real :: t0=0,Rgse0=5.1,Latgse0=0.,Longse0=15.,&
         Thetagse,LonGseRad
 real :: cRadToDeg=57.2957795131
 real :: Hden
 external RECALC_08,SPHCAR_08,gswgse_08,smgsw_08  ! geopack_2008 subroutine

 call read_mate_data

 call recalc_08(2008,164,0,0,0,-500.,0.,0.)

  ! (1) read_mate_data test
  write(*,*) ' R (RE), Lat (deg), Lon (deg)' 
  write(*,'(3f8.0)') RMateGSE(iRMate),LatMateGSE(iLatMate),LonMateGSE(iLonMate)
  write(*,*) 'Time, Hden' 
  write(*,'(5("(",0pf6.0,",",1pE10.3,")"))') &
       (TimeMate(iTime)*cSecToHour,&
        10.**logHdenMate(iTime,iRMate,iLatMate,iLonMate),iTime=1,nHourMate)

  ! (2) interpolation test 
  ThetaGSE=(90.-LatGSE0)/cRadToDeg
  LonGSERad=LonGSE0/cRadToDeg
  call SPHCAR_08(rgse0,ThetaGSE,LongseRad,xGSE,yGSE,zGSE,1)
  call gswgse_08(xGSW,yGSW,zGSW,xGSE,yGSE,zGSE,-1)
  call smgsw_08(xsm,ysm,zsm,xGSW,yGSW,zGSW,-1)
  call mate_geocorona(t0,xsm,ysm,zsm,hden)
  !call mate_geocorona(t0,Rgse0,Latgse0,Longse0,hden)
  
  write(*,*) 'Hden(interpol) Hden(Mate)'
  !write(*,'(4f9.1)') Hden, HdenMate(1,1,7,1),HdenMate(1,2,7,1),HdenMate(1,1,7,1)*0.8+HdenMate(1,2,7,1)*0.2
  !write(*,'(4f9.1)') Hden, 10.**logHdenMate(1,1,7,2),10.**logHdenMate(1,2,7,2),&
  !                   10.**(logHdenMate(1,1,7,2)*0.8 + logHdenMate(1,2,7,2)*0.2)
  write(*,'(4f9.1)') Hden, 10.**logHdenMate(1,7,7,2),10.**logHdenMate(1,8,7,2),&
                     10.**(logHdenMate(1,7,7,2)*0.8 + logHdenMate(1,8,7,2)*0.2)

  end subroutine unit_test_mate

end module ModMate
