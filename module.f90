!
!   module.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!
module parameters
    implicit none
  
    !real(kind(0d0)), parameter :: pi=3.1415926535897932d0
    integer :: npButterworth
    real(kind(0d0)) :: fmin,fmax
    real(kind(0d0)) :: eps=1.d-5 ! tolerance during MT inversion
    integer :: nmt
    !integer :: npDSM ! integer length of synthetics (to be considered)
    integer :: npData ! integer length of observed
    integer :: ntStep
    real(kind(0d0)) :: dt, tlenData ! tlenDSM
    real(kind(0d0)), allocatable :: twin(:,:),twinObs(:,:)
    real(kind(0d0)) :: movingWindowStep
    integer :: NmovingWindowDimension
    integer, allocatable :: itwin(:,:)
    integer, allocatable :: itwinObs(:,:)
    integer :: ntwin,ntwinObs
    real(kind(0d0)), allocatable :: fMovingWindowStart(:), fMovingWindowEnd(:), fEachShift(:,:)
    integer, allocatable :: iMovingWindowStart(:),iMovingWindowEnd(:)
    integer, allocatable :: totalNumberInWindowDimension(:)
    integer :: nTimeCombination
    
    integer, allocatable :: iEachWindowStart(:,:), iEachWindowEnd(:,:)



    character(200), allocatable :: filenames(:)
    real(kind(0d0)), allocatable :: obsRaw(:,:), obsFilt(:,:),obsFiltTapered(:,:)
    integer :: calculMode ! 0=normal; 1=filter and stop
    character(200) :: workingDir
    character(200) :: resultDir
    character(200) :: inversionName
    integer :: nConfiguration ! The number of configurations (source location in 3D space, model)
    real(kind(0d0)), allocatable :: conf_depth(:),conf_lat(:),conf_lon(:),conf_gcarc(:),conf_azimuth(:)

    character(200) :: SGTinfo
    character(200) :: parentDir
    character(200) :: eventName
    character(200) :: stationName
    real(kind(0d0)) :: stla,stlo
    character(200) :: INFO_TSGT,INFO_RSGT,rsampletxt,modelcard
    character(200) :: Poutputdir,psvmodel,modelname,DSMconfFile
    integer :: imin,imax,maxlmax
    real(kind(0d0)) :: omegai
    real(kind(0d0)) :: r0delta,r0max,r0min,ratc,ratl,rdelta_,re,rmax_,rmin_
    integer :: rsgtswitch,synnswitch,tsgtswitch
    real(kind(0d0)) :: thetadelta,thetamax,thetamin
    !real(kind(0d0)), parameter :: phimin=0.d0
    !real(kind(0d0)), parameter :: phimax=0.d0
    !real(kind(0d0)), parameter :: phidelta=1.d0
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
    real(kind(0d0)), allocatable :: t(:),rsgtTime(:,:),u0(:,:),v(:),v0(:,:),hu(:),hu0(:,:)
    real(kind(0d0)), allocatable :: fwin(:,:)
    integer, allocatable :: nt1(:),nt2(:)
    real(kind(0d0)), allocatable :: denomv(:),denomu(:),coeff(:,:,:,:),coeffV(:,:)
    integer, parameter :: num_rsgtPSV=10
    integer, parameter :: num_rsgtSH=5
    integer, parameter :: num_rsgt=10

    real(kind(0d0)), allocatable :: tmparray(:,:,:)
  
end module tmpSGTs


module angles
    implicit none
    !real(kind(0d0)), allocatable :: phi00(:),phi0(:),theta0(:)
    integer :: ntheta,nphi,nr
    real(kind(0d0)), allocatable :: phitheta(:,:),thetaphi(:,:)
    real(kind(0d0)), allocatable :: latgeo(:,:), longeo(:,:)
    real(kind(0d0)), allocatable :: phi(:,:), theta(:,:)
    real(kind(0d0)), allocatable :: crq(:,:),srq(:,:),crq2(:,:),srq2(:,:)
    !real(kind(0d0)), allocatable :: csq(:,:),ssq(:,:),csq2(:,:),ssq2(:,:)
    !real(kind(0d0)), allocatable :: cqs(:,:),sqs(:,:),cqs2(:,:),sqs2(:,:)
    real(kind(0d0)), allocatable :: deltar(:,:),deltas(:,:)
    real(kind(0d0)) :: slat,slon,sdep,rlat,rlon
    real(kind(0d0)) :: gcarcmin,gcarcmax,dgcarc,azimuthmin,azimuthmax,dazimuth
    real(kind(0d0)) :: radiusmin,radiusmax,dradius
    real(kind(0d0)), allocatable  :: gcarc(:), azimuth(:),radius(:)
    integer,allocatable :: ithetaD(:), iradiusD(:)
    !real(kind(0d0)) :: pi=4.d0*datan(1.d0)
    real(kind(0d0)), parameter :: pi=3.1415926535897932d0
    real(kind(0d0)) :: radian2degree=180.d0/pi
    real(kind(0d0)) :: degree2radian=pi/180.d0
end module angles
