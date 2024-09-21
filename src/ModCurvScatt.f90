!*******************************************************************************
!                           ModCurvScatt.f90
!
! Module field line curvature (FLC) scattering.
!
! Note: no variable has a ghost cell.
!*******************************************************************************
  module ModCurvScatt
    implicit none

    public :: &
      ! subroutines
      init_flc,&
      calc_flc_para,&
      calc_Daa_flc,&
      write_flc,&
      unit_test_Daa_flc,&
      unit_test_r_flc,&
      ! variables
      iflc,&
      tflc,&
      r_flc,&
      !zeta1_flc,&
      !zeta2_flc,&
      ReqMin,&
      Daa_flc

    private   ! except

    ! variables from cimi
    integer :: ir,ip,iw,ik,ns
    real :: q_e,re_m

    ! variables to cimi
    real,allocatable :: &
      r_flc(:,:),&       ! field line curvature radius
      p_flc(:,:),&       ! |e|*B*Rc
      zeta1_flc(:,:),&  ! Rc*d^(Rc)/ds^2
      zeta2_flc(:,:),&  ! Rc^2/B*d^(B)/ds^2
      Daa_flc(:,:,:,:,:) ! pitch-angle diffusion coeff. 
!    real tau_flc(ns,ir,ip,iw,ik)

    ! other variables
    integer iflc       ! iflc=0: no FLC,  iflc=1: FLC
    real zeta1,zeta2   ! variables in Young et al. 2008
    real :: ReqMin=4.  ! assume dipole curvature where <ReqMin (RE)
!   real :: tflc=0.    ! FLC calculation after tFLC in the simulation (sec)
    real tflc          ! FLC calculation after tFLC in the simulation (sec)

    ! variables for output files
    integer :: nout_flc=0   ! # of FLC output files
    real tstart,tstep

  contains
!***********************************************************************
!                        init_flc     
!
!  Routine receives parameters required to calculate pitch-angle
!   diffusion coefficient (Daa) due to field line curvature 
!   (curvature scattering) from cimi.
!***********************************************************************
  subroutine init_flc(ns_,ir_,ip_,iw_,ik_,echarge_,re_m_)
  integer,intent(in) :: ns_,ir_,ip_,iw_,ik_
  real,intent(in) :: echarge_,re_m_

  ns=ns_
  ir=ir_
  ip=ip_
  iw=iw_
  ik=ik_
  q_e=echarge_
  re_m=re_m_

  if (.not.allocated(r_flc)) &
     allocate(r_flc(ir,ip))

  if (.not.allocated(p_flc)) &
     allocate(p_flc(ir,ip))

  if (.not.allocated(zeta1_flc)) &
     allocate(zeta1_flc(ir,ip))

  if (.not.allocated(zeta2_flc)) &
     allocate(zeta2_flc(ir,ip))

  if (.not.allocated(Daa_flc)) &
     allocate(Daa_flc(ns,ir,ip,iw,ik))

  !! The followings show modules and variables in cimi main body
  !use cimigrid_dim, only: iw,ik
  !use constants, only: echarge, re_m

  end subroutine init_flc


!***********************************************************************
!                        calc_flc_para      
!  Routine calculates parameters required to calculate pitch-angle
!   diffusion coefficient (Daa) due to field line curvature 
!   (curvature scattering).
!***********************************************************************
  subroutine calc_flc_para(i,j,npf,ibmin,req,dss,xa1,ya1,za1,bba)
 
  !integer,parameter :: np=5
  integer,parameter :: np=7
  integer,intent(in) :: i,j,npf,ibmin  ! i=index of MLAT, 
                                       ! j=index of MLT,
                                       ! ibmin=index at minimum B, Bmin
  real,intent(in) :: req               ! distance from the earth center to Bmin
  real,dimension(npf),intent(in) :: dss,xa1,ya1,za1,bba
  real,dimension(np) :: xa,ya,za,ds,Rc
  real Rc0,ds0,ds1,ds2,dRc1,dRc2,B0,dB1,dB2
  integer i10,i01,is,is1,is2,ii,ii1

  ! (1) calculate curvature at bo
    i10=ibmin+(np-1)/2
    i01=ibmin-(np-1)/2
    xa(1:np)=xa1(i01:i10)
    ya(1:np)=ya1(i01:i10)
    za(1:np)=za1(i01:i10)
    ds(1:np)=dss(i01:i10)
    call calc_flc(np,xa,ya,za,ds,Rc)
    is=(np-1)/2+1
    is=minloc(Rc(2:np-1),dim=1)
    is=is+1
    ! the reason 9 points are needed here:
    !  For some cases, "ibmin" /= "is", and +- 3 points are enought in most cases
    !  in addition, it needs +- a point from "is" 
    !  to calculate second derivatives of Rc 
    Rc0=Rc(is)
    r_flc(i,j)=Rc0   ! field line curvature radius
    p_flc(i,j)=bba(ibmin)*q_e*Rc0*re_m  ! p corresponding to FLC in kg m/s

  ! (2) calculate zeta1,zeta2 at bo
    is1=is-1
    is2=is+1
    ds1=0.5*(ds(is1)+ds(is))
    ds2=0.5*(ds(is2)+ds(is))
    ds1=sqrt((xa(is1)-xa(is))**2 &
            +(ya(is1)-ya(is))**2 &
            +(za(is1)-za(is))**2 )
    ds2=sqrt((xa(is2)-xa(is))**2 &
            +(ya(is2)-ya(is))**2 &
            +(za(is2)-za(is))**2 )
    ds0=ds(is)
    if (req.gt.ReqMin) then   
       dRc1=Rc0-Rc(is1)
       dRc2=Rc(is2)-Rc0
       ! zeta1 as in Young et al. 2002; 2008
       zeta1=Rc0*(dRc2/ds2-dRc1/ds1)/ds0 
       !B0=bba(is)
       !dB1=B0-bba(is1)
       !dB2=bba(is2)-B0           
       B0=bba(ibmin)
       dB1=B0-bba(ibmin-1)
       dB2=bba(ibmin+1)-B0           
       ! Zeta2 as in Young et al. 2002; 2008
       zeta2=Rc0*Rc0/B0*(dB2/ds2-dB1/ds1)/ds0     
    else    ! assume a dipole < ReqMin
       zeta1=0.6666667
       zeta2=1.
    endif

    zeta1_flc(i,j)=zeta1
    zeta2_flc(i,j)=zeta2

  end subroutine calc_flc_para

  
!***********************************************************************
!                           calc_flc      
!  Routine calculates field line curvature 
!***********************************************************************
    subroutine calc_flc(n,xa,ya,za,ds,Rc)

    integer,intent(in) :: n
    real,dimension(n),intent(in) :: xa,ya,za,ds
    real,intent(out) :: Rc(n)
    integer i,i1,i_1
    real ds0,ds1,ds2,dx1,dx2,dy1,dy2,dz1,dz2,dx,dy,dz,ddx,ddy,ddz,&
         dxds1,dyds1,dzds1,dxds2,dyds2,dzds2

    !ds1=ds(1)
    dx1=xa(2)-xa(1)
    dy1=ya(2)-ya(1)
    dz1=za(2)-za(1)
    ds1=0.5*(ds(2)+ds(1))  
    !ds1=sqrt((xa(2)-xa(1))**2 &
    !        +(ya(2)-ya(1))**2 &
    !        +(za(2)-za(1))**2 )
    dxds1=dx1/ds1
    dyds1=dy1/ds1
    dzds1=dz1/ds1
    do i=2,n-1
       ds0=ds(i)
       i1=i+1
       i_1=i-1
       dx2=xa(i1)-xa(i)
       dy2=ya(i1)-ya(i)
       dz2=za(i1)-za(i)
       ds2=0.5*(ds(i1)+ds(i))  ! ds from i to i1
       !ds2=sqrt((xa(i1)-xa(i))**2 &
       !        +(ya(i1)-ya(i))**2 &
       !        +(za(i1)-za(i))**2 )
       dxds2=dx2/ds2
       dyds2=dy2/ds2
       dzds2=dz2/ds2
       dx=0.5*(xa(i1)-xa(i_1))
       dy=0.5*(ya(i1)-ya(i_1))
       dz=0.5*(za(i1)-za(i_1))
       ddx=(dxds2-dxds1)
       ddy=(dyds2-dyds1)
       ddz=(dzds2-dzds1)
       Rc(i)=(dx**2+dy**2+dz**2)**1.5&
            /sqrt((ddx*dy-dx*ddy)**2 &
                 +(ddy*dz-dy*ddz)**2 &
                 +(ddz*dx-dz*ddx)**2)/ds0
       dxds1=dxds2
       dyds1=dyds2
       dzds1=dzds2
    enddo
    Rc(1)=Rc(2)
    Rc(n)=Rc(n-1)

    end subroutine calc_flc


!***********************************************************************
!                        calc_Daa_flc      
!  Routine calculates parameters required to calculate pitch-angle
!   diffusion coefficient (Daa) due to field line curvature 
!   (curvature scattering).
!***********************************************************************
  subroutine calc_Daa_flc(n,i,j,p_I,Tbounce0,sina,Beq,BiN,BiS,mcN,mcS)

  !!use cimigrid_dim, only: iw,ik
  !!use ModCurvScatt, only: ie,pc_flc,zeta1,zeta2,Daa_flc,ekev_flc
  
  integer,intent(in) :: n,i,j,&
                        mcN,&            ! index of invariant K at BiN
                        mcS              ! index of invariant K at BiS
  real,intent(in) :: p_I(iw),&           ! momentum in mks unit
                     Tbounce0(iw,ik),&   ! bounce time in second
                     sina(ik),&          ! sin(pitch-angle)
                     Beq,&               ! B (nT) at the equator of minimum B
                     BiN,&               ! B (nT) at the northern footpoint
                     BiS                 ! B (nT) at the southern footpoint
  real sinaL,aL,Bi,a0sq(iK),aLsq
  ! variables as in Young et al. 2002; 2008
  real e, e2, e3,&  ! epsilon, epsilon^2, epsilon^3
       c,a1,a2,De,b0,w0,tb,D0,N2,Amx0,&
       a0(ik),&     ! pitch-angle in rad
       cosa(ik),&   ! cos(a0)
       sinacosa(ik),&  ! sina * cosa
       cosa_bar,&   ! cos(a0) at maximum of Amax
       Dmax         ! diffusion coefficient upper limit
  real,dimension(iw,ik) :: Amax,Aa0
  integer k,m,m_bar,m_max,m_L

  ! calculate loss cone
  Bi=min(BiN,BiS)   ! take the bigger loss cone
  if (Bi==BiN) then
     m_L=mcN
  else
     m_L=mcS
  endif
  sinaL=sqrt(Beq/Bi)
  aL=asin(sinaL)
  aLsq=aL**2

  ! calculate a0 and cosa
  do m=1,ik
     a0(m)=asin(sina(m))
     cosa(m)=cos(a0(m))
     a0sq(m)=a0(m)**2
     sinacosa(m)=sina(m)*cosa(m)
  enddo

  ! for coupled cimi p_I=p_I(k,m) and for standalone p_I=p_I(k)
  !do m=1,ik
  !enddo 
  do k=1,iw
     ! coefficients as in Young et al. 2002; 2008
     e=p_I(k)/p_flc(i,j)       ! epsilon
     !if (e.gt.0.584) e=0.584  ! slow diffusion limit
     if (e.gt.10.) e=10.       ! no slow diff limit
     if (e.gt.0.05) then
        e2=e*e
        e3=e2*e
        w0=1.051354+0.13513581*e-0.50787555*e2        ! omega
        if (w0>1.999999) w0=1.999999
        c=1.0663037-1.0944973/e+0.016679378/e2-0.00049938987/e3
        a1=-0.35533865+0.12800347/e+0.0017113113/e2    
        a2= 0.23156321+0.15561211/e-0.001860433/e2    
        b0=-0.51057275+0.93651781/e-0.0031690066/e2    
        De=-0.49667826-0.00819418/e+0.0013621659/e2  ! D(epsilon)
        Amax(k,:)=exp(c)*(zeta1_flc(i,j)**a1*zeta2_flc(i,j)**a2+De)
        Aa0(k,:)=sin(w0*a0(1:ik))*cosa(1:ik)**b0
     else 
        Amax(k,:)=0. 
     endif
  enddo

  do k=1,iw
     m_bar=maxloc(abs(Aa0(k,:)),dim=1)
     cosa_bar=cosa(m_bar)
     do m=1,ik
        Amx0=Amax(k,m)
        if (Amx0.ne.0.) then
           tb=Tbounce0(k,m)
           N2=1./Aa0(k,m_bar)/Aa0(k,m_bar)
           D0=0.5*Amx0*Amx0/tb
           !Dmax=4.9348/tb    ! da0=0.5pi during 0.5*tb
           if (m<m_L) then
              Dmax=2.*a0sq(m)/tb ! da0=a0 during 0.5*tb
           else
              Dmax=2.*aLsq/tb    ! da0=aL during 0.5*tb
           endif
           Daa_flc(n,i,j,k,m)=&
             D0*N2*(Aa0(k,m)/sinacosa(m))**2
           if (Daa_flc(n,i,j,k,m).gt.Dmax) Daa_flc(n,i,j,k,m)=Dmax
        else
           Daa_flc(n,i,j,k,m)=0.
        endif
     enddo
  enddo

  end subroutine calc_Daa_flc


!***********************************************************************
!                           write_flc     
! 
!  Routine prints curvature radius values 
!***********************************************************************
  subroutine write_flc(t,outname,ijs,js,ro,xmlto,bo,iw1,ik1,ener0,sina)

  integer,intent(in) :: ijs,iw1,ik1  ! # of species; index of energy and K grid
  integer,intent(in) :: js(ijs)      ! indices of species
  real,intent(in) :: t,ener0,sina    ! time; energy(keV); sina
  real,intent(in) :: ro(ir,ip),xmlto(ir,ip),bo(ir,ip) ! r, mlt, Beq
  character,intent(in) :: outname*8  ! output file name; species name
  character outname1*40

  real :: RadToDeg=57.29578049
  integer ispec

  ! output start time
  if (nout_flc==0) tstart=t
  if (nout_flc==1) tstep=t-tstart

  ! # of flc output files
  nout_flc=nout_flc+1

  if (t.eq.tstart) then 
     call system('mkdir -p FLC')
     open(unit=61,file='FLC/'//trim(outname)//'.parm',&
          form='unformatted',status='replace')
     write(61) tstart,tstep,nout_flc,ir,ip,ijs
     write(61) js(1:ijs)
     close(61)
  endif

  write(outname1,'(a,i8.8,a)') 'FLC/'//trim(outname)//&
        '_',int(t),'.flc'
  open(unit=60,file=trim(outname1),&
       form='unformatted',status='replace')
  write(60) t
  write(60) ener0,asin(sina)*RadToDeg
  write(60) ro(:,:)
  write(60) xmlto(:,:)
  write(60) bo(:,:)
  write(60) r_flc(:,:)
  do ispec=1,ijs
     write(60) Daa_flc(ispec,:,:,iw1,ik1)
  enddo
  close(60)

  open(unit=62,file='FLC/'//trim(outname)//'.end',&
       form='unformatted',status='replace')
  write(62) t
  close(62)
  
  end subroutine write_flc

!***********************************************************************
!                        unit_test_Daa_flc
!
!  This subroutine tests Daa_flc
!***********************************************************************
  subroutine unit_test_Daa_flc
  integer :: ir_=1,&
             ip_=1,&
             ns_=1,&
             iw_=40,&
             ik_=44

  real :: echarge_=1.60217663e-19,&  ! elementary charge in Coulumb
          re_m_=6.375e6,&            ! Earth's radius in m
          alt=1.e5,&                 ! altitude of loss cone 
          BE=3.12e-5,&               ! B at the surface at the equator in Tesla
          ReqMin_=4.                 ! lower limit for curvature calcuation

  real,allocatable :: &
    p_I(:),Tbounce_II(:,:),sina_I(:),v_I(:),E_I(:),a0_I(:)

  real :: L=4.,&          ! L-shell
          Emin=1000.,&    ! in eV
          Emax=100000.,&
          Beq,BiN,BiS,cosl2,sinaL,&
          rE,&
          Cmin=0.15,&
          Cmax=0.58,&
          dC,&
          a0min=2.,&
          a0max=88.,&
          da0,Ta0,Gamm,E0,&
          cDegToRad=0.0174533,&  ! deg to rad
          EvToKev=1.e-3,&        ! eV to keV
          mp=1.67262192369e-27,& ! proton rest mass in kg
          Ep=938272088.16,&      ! proton rest mass in eV
          EMspeed=2.99792458e8   ! speed of light in m/s

  integer k,m,mcN,mcS

  ! initialization
  ReqMin=ReqMin_
  call init_flc(ns_,ir_,ip_,iw_,ik_,echarge_,re_m_)

  ! variables from other code such as CIMI
  if (.not.allocated(p_I)) allocate(p_I(iw))
  if (.not.allocated(v_I)) allocate(v_I(iw))
  if (.not.allocated(E_I)) allocate(E_I(iw))
  if (.not.allocated(Tbounce_II)) allocate(Tbounce_II(iw,ik))
  if (.not.allocated(sina_I)) allocate(sina_I(ik))
  if (.not.allocated(a0_I)) allocate(a0_I(ik))

  ! loss cone at L=4
  Beq=BE/L**3
  cosl2=L*RE/(RE+alt)
  BiN=Beq*sqrt(4.+3*cosl2)/cosl2**3
  BiS=BiN
  sinaL=sqrt(Beq/BiN)
  mcN=ik

  ! pitch-angle grid
  da0=(a0max-a0min)/float(ik-1)
  do m=1,ik
     a0_I(m)=a0min+da0*(float(m-1))
     sina_I(m)=sin(a0_I(m)*cDegToRad)
     if (mcN==ik.and.sina_I(m)<sinaL) mcN=m
  enddo
  mcS=mcN

  ! gyroradius to curvature radius ratio
  zeta1_flc(1,1)=0.6666667  ! value in a dipole
  zeta2_flc(1,1)=1.
  r_flc(1,1)=1.333333
  p_flc(1,1)=r_flc(1,1)*BE/L**3*q_e*re_m
  dC=(Cmax-Cmin)/float(iw-1)

  do k=1,iw
     p_I(k)=p_flc(1,1)*(Cmin+dC*float(k-1))
     !v_I(k)=p_I(k)/mp
     !E_I(k)=0.5*p_I(k)**2/mp/q_e
     E0=mp*EMspeed**2
     E_I(k)=sqrt(E0**2 + (p_I(k)*EMspeed)**2)-E0
     Gamm=(1.+E_I(k)/E0)
     E_I(k)=E_I(k)/q_e
     v_I(k)=p_I(k)/mp/Gamm
  enddo


  ! bounce period
  do m=1,ik
     Ta0=1.38073-0.639693*sina_I(m)**0.75
     do k=1,iw
        Tbounce_II(k,m)=4.*L*re_m*Ta0/v_I(k)
        write(*,'(f11.1,2f10.3)') E_I(k),a0_I(m),Tbounce_II(k,m)
     enddo
  enddo


  ! calculate Daa_flc
  call calc_Daa_flc(1,1,1,p_I,Tbounce_II,sina_I,Beq,BiN,BiS,mcN,mcS)

  ! write Daa_flc
  open(unit=2,file='Daa_flc_test.dat')
  write(2,*) iw,ik,'   # of E and pitch-angle grid bins'
  write(2,*) ' energy in (keV)'
  write(2,'(10f11.1)') E_I(:)*EvTokeV
  write(2,*) ' gyro to curvature radius ratio'
  write(2,'(10f7.4)') p_I(:)/p_flc(1,1)
  write(2,*) ' pitch-angle (deg)'
  write(2,'(10f7.2)') a0_I(1:ik)
  do k=1,iw
     write(2,'(f11.1,f7.3,a)') E_I(k)*EvTokeV,p_I(k)/p_flc(1,1),' keV, epsilon'
     do m=1,ik
        write(2,'(f7.2,1pE11.4)') a0_I(m),Daa_flc(1,1,1,k,m)
     enddo
  enddo
  close(2)

  end subroutine unit_test_Daa_flc


!***********************************************************************
!                        unit_test_r_flc
!
!  This subroutine tests field line curvature radius r_flc
!***********************************************************************
  subroutine unit_test_r_flc
  integer :: ir_=1,&
             ip_=1,&
             ns_=1,&
             iw_=40,&
             ik_=44

  real :: echarge_=1.60217663e-19,&  ! elementary charge in Coulumb
          re_m_=6.375e6,&            ! Earth's radius in m 
          ReqMin_=5.                 ! lower limit for curvature calcuation

  ! variables for the subroutine calc_flc_para
  integer,parameter :: npf=9
  integer :: i=1,&
             j=1,&
             ibmin

  real :: req=4.
  real,dimension(npf) :: dss,xa1,ya1,za1,bba

  real :: dlambda=0.00174532925,& ! in rad, 0.1 deg
          L=4.,&                 ! L-shell
          BE=3.12e-5,&           ! B at the surface at the equator in Tesla
          r,lambda,cosl,sinl,dx,ddx,dz,ddz,dr,ddr,rc

  integer m

  ! initialization
  ReqMin=ReqMin_
  call init_flc(ns_,ir_,ip_,iw_,ik_,echarge_,re_m_)

  ! diploe magnetic field
  L=req
  ibmin=(npf-1)/2+1
  do m=1,npf
     lambda=dlambda*float(m-ibmin)
     cosl=cos(lambda)
     sinl=sin(lambda)
     r=L*cosl**2
     bba(m)=BE*sqrt(1.+3.*sinl**2)/r**3
     xa1(m)=-r*cosl
     ya1(m)=0.
     za1(m)=r*sinl
  enddo

  ! dss
  do m=2,npf-1
     dss(m)=0.5*sqrt((xa1(m-1)-xa1(m+1))**2 + (za1(m-1)-za1(m+1))**2)
  enddo
  dss(1)=dss(2)
  dss(npf)=dss(npf-1)

  ! curvature radius
  cosl=1.
  sinl=0.
  !cosl=cos(-dlambda*2.)
  !sinl=sin(-dlambda*2.)
  r=L*cosl**2
  dr=-2.*L*cosl*sinl
  ddr=2.*L*(sinl**2 - cosl**2)
  dx=dr*cosl - r*sinl
  dz=dr*sinl + r*cosl
  ddx=(ddr-r)*cosl - 2.*dr*sinl
  ddz=(ddr-r)*sinl + 2.*dr*cosl
  rc=(dx**2 + dz**2)**1.5/abs(ddx*dz-dx*ddz)

  call calc_flc_para(i,j,npf,ibmin,req,dss,xa1,ya1,za1,bba)

  write(*,'(a,f8.3,a)') 'curvature radius: ',rc,' RE, dipole'
  write(*,'(a,f8.3,a)') 'curvature radius: ',r_flc(i,j),' RE, simulation'
  write(*,'(2f8.3,a)') 0.66667,1.,' zeta1, zeta2 in dipole'
  write(*,'(2f8.3,a)') zeta1,zeta2,' zeta1, zeta2 in the simulation'

  end subroutine unit_test_r_flc


  end module           
