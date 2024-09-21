!*****************************************************************************
!  
!                         EMIC.f90 
!
!  This file contains diffusion coefficients accosiated EMIC waves
!
!  Diffusion coefficients are calculated according to the analytic
!   quasi-linear theory for field-aligned wave in Summer [2005].
!   Then, the coefficients are corrected by multiplying with 0.5 factor
!   according to Albert [2007].
!   The coefficients are bounce averaged using magnetic dipole field model,
!   according to Lyons et al. [1972].
!
!  We only consider pitch angle diffusion coefficients,
!   because energy and cross diffusion of electrons by EMIC waves are very small.
!
!  We use time-MLT-averaged intensities 
!    estimated from CRRES data [Kersten et al.,2014].
!
!  Created on November 07, 2014 by Suk-Bin Kang.
!
!  Modification History
!  * October 03 2015 - Add H band EMIC waves
!                    - Add calculation of plume location  
!  * October 22 2015 - irw = 21 -> 18 (in module D_EMIC in cimi4.f90)
!                    - wLshell = 1.25 - 11.25 -> 1.5 - 9.5
!
!*****************************************************************************
 module ModEmic
 use cimigrid_dim, only: ip
 implicit none

 integer,parameter :: irEMIC=16

 real LEmic(irEmic),MltEmic(ip),&
      EmicHIntens_I(irEMIC,ip,3),& 
      EmicHeIntens_I(irEMIC,ip,3) 

 contains
!*****************************************************************************
!                           readEMICCoef
!
!  Routine reads H and He band EMIC diffusion coeff.
!*****************************************************************************

  subroutine readEMICCoef
  end subroutine readEMICCoef

!*****************************************************************************'
!                          setEmicIntensity  
!
! Routine set the EMIC intensity base on the CRRES data
! Routine regards Kp=(2,4) as AE=(100,500) nT
! L = 3.5 - 7, every 0.5
! MLT = 12 - 18
!*****************************************************************************
  subroutine setEMICIntensity

 
  real EMIC_CRRES_H(24),EMIC_CRRES_He(24)
  
  data EMIC_CRRES_H /0.,0.,0.,5.e-5,6.25e-5,1.25e-5,0.,0.,&
                     2.5e-4,4.5e-4,4.75e-4,7.75e-4,1.e-3,1.3e-3,1.18e-3,4.25e-4,&
                     3.e-3,8.5e-3,7.e-3,7.e-3,2.05e-2,2.5e-2,2.e-2,1.25e-2/
  data EMIC_CRRES_He /2.e-4,3.e-5,2.5e-4,2.5e-4,5.e-4,0.001,5.e-4,5.e-4,&
                      3.5e-3,5.e-3,0.01,0.01,0.02,2.5e-2,2.5e-2,0.005,&
                      5.e-4,7.e-2,0.1,0.1,0.1,8.e-2,3.e-2,0.001/

! Set Lshell value L = 3.5 - 7
  do i=1,irEMIC
     LEmic(i)=i*0.5+1.       ! L = 1.5 - 9.5
  enddo
! Set MLT value
  do j=1,ip
     MltEmic(j)=(j-1)*0.5           ! mlt = 0.0 - 23.5
  enddo

! Set EMIC intensity
  do i=5,12
     EmicHIntens(i,25:37,1)=EMIC_CRRES_H(i-4)
     EmicHIntens(i,25:37,2)=EMIC_CRRES_H(i+4)
     EmicHIntens(i,25:37,3)=EMIC_CRRES_H(i+12)   
     EmicHeIntens(i,25:37,1)=EMIC_CRRES_He(i-4)
     EmicHeIntens(i,25:37,2)=EMIC_CRRES_He(i+4)
     EmicHeIntens(i,25:37,3)=EMIC_CRRES_He(i+12)  ! for kp>4, MAX ~ 0.1nT^2
  enddo

  end subroutine setEmicIntensity

!*****************************************************************************'
!                          get_emic_power  
!
! Routine calculates the EMIC intensity
! Routine regards Kp=(2,4) as AE=(100,300) nT
! Input ro,xmlt,iae
! Output Bpower1,Bpower2
!*****************************************************************************
  subroutine get_emic_power
  use cfield, only: ro,xmlto
  use convect,only: AE
  use cWpower, only: &
     !  inside plasmasphere
          EmicHPowerInPs ,& ! H band
          EmicHePowerInPs,& ! He band
          !EmicOPowerInPs,& ! O band
     !  outside plasmasphere
          EmicHPowerOutPs ,& ! H band
          EmicHePowerOutPs,& ! He band
          !EmicOPowerOutPs    ! He band

  integer  iae,i,j
  real Bpower1,Bpower2,mlt,ro,xmlt,EBpower(irw,ip)
 
! Determine AE level in 3 bands: 0-100, 100-300, >300
  if (AE<100.) iae=1
  if (AE>=100..and.AE<300.) iae=2
  if (AE>=300.) iae=3
 
  do i=1,ir
     do j=1,ip
        call lintp2(LEmic,MltEmic,EmicHIntens(1:irEmic,1:ip,iae),irEmic,ip,&
                    ro(i,j),xmlto(i,j),EmicHPowerInPs(i,j))
        call lintp2(LEmic,MltEmic,EmicHeIntens(1:irEmic,1:ip,iae),irEmic,ip,&
                    ro(i,j),xmlto(i,j),EmicHePowerInPs(i,j))
        EmicHPowerOutPs(i,j)=EmicHPowerInPs(i,j)
        EmicHePowerOutPs(i,j)=EmicHePowerInPs(i,j)
     enddo
  enddo

  end subroutine get_emic_power

 end module ModEmic
