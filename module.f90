!
!   module.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!
module parameters
  implicit none
  
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  integer :: npButterworth
  real(kind(0d0)) :: fmin,fmax
  real(kind(0d0)) :: eps=1.d-5 ! tolerance during MT inversion
  integer :: nmt
  integer :: np ! integer length of synthetics (to be considered)
  integer :: npData ! integer lenght of observed
  integer :: ntStep
  real(kind(0d0)) :: dt, tlen, tlenData
  real(kind(0d0)), allocatable :: twin(:,:)
  integer, allocatable :: itwin(:,:)
  integer :: ntwin
  character(200), allocatable :: filenames(:)
  real(kind(0d0)), allocatable :: obsRaw(:,:), obsFilt(:,:)
  integer :: calculMode ! 0=normal; 1=filter and stop
  character(200) :: workingDir
  character(200) :: resultDir
  character(200) :: inversionName
  integer :: nConfiguration ! The number of configurations (source location in 3D space, model)
  character(200) :: SGTinfo
  character(200) :: parentDir
  character(200) :: eventName
  character(200) :: stationName
  real(kind(0d0)) :: stla,stlo
  character(120) :: INFO_TSGT,INFO_RSGT,rsampletxt,modelcard
  character(120) :: Poutputdir,psvmodel,modelname,DSMconfFile
  integer :: imin,imax,maxlmax
  real(kind(0d0)) :: omegai
  real(kind(0d0)) :: r0delta,r0max,r0min,ratc,ratl,rdelta_,re,rmax_,rmin_
  integer :: rsgtswitch,synnswitch,tsgtswitch
  real(kind(0d0)) :: thetadelta,thetamax,thetamin
  real(kind(0d0)), parameter :: phimin=-180.d0
  real(kind(0d0)), parameter :: phimax= 180.d0
  real(kind(0d0)), parameter :: phidelta=1.d0
  real(kind(0d0)) :: tlenFull
  integer :: r_n, theta_n,phi_n,iWindowStart,iWindowEnd,lsmooth,np0,np1
  integer :: iConfTheta,iConfPhi,iConfR
  real(kind(0d0)), allocatable :: r_(:),thetaD(:),phiD(:),omega(:)
  real(kind(0d0)) :: dtn,start,end,samplingHz
end module parameters


module tmpSGTs
  implicit none
  ! in frequency domain
  complex(kind(0d0)), allocatable :: rsgtF(:,:)
  complex(kind(0d0)), allocatable :: h3(:,:),h4(:,:)
  complex(kind(0d0)), allocatable :: rsgtomega(:,:,:)
  complex(kind(0d0)), allocatable :: u_freq(:)
  ! in time domain
  real(kind(0d0)), allocatable :: t(:),u(:,:),u0(:,:),v(:),v0(:,:),hu(:),hu0(:,:)
  real(kind(0d0)), allocatable :: fwin(:,:)
  integer, allocatable :: nt1(:),nt2(:)
  real(kind(0d0)), allocatable :: denomv(:),denomu(:),coeff(:,:,:,:),coeffV(:,:)
  integer, parameter :: num_rsgtPSV=10
  integer, parameter :: num_rsgtSH=5
  integer, parameter :: num_rsgt=10
end module tmpSGTs


