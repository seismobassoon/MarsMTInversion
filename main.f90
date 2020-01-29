!
!   main.f90
!   MarsInversion
!
!   Created by fuji on 04/12/2019.
!   Copyright 2019 nfuji. All rights reserved.
!

program MarsInversion



  use parameters
  use tmpSGTs
  implicit none

  integer :: mtcomp,jmtcomp
  integer ::icomp,iWindow,it,jjj
  integer :: iMovingWindow,iConfiguration,iMovingWindowStep
  real(kind(0d0)), allocatable :: taperDSM(:),taperOBS(:)
  real(kind(0d0)), allocatable :: tmparray(:,:,:)
  real(kind(0d0)), allocatable :: GreenArray(:,:,:)
  real(kind(0d0)), allocatable :: obsArray(:,:),obsRawArray(:,:)
  real(kind(0d0)), allocatable :: modArray(:,:),modRawArray(:,:)
  real(kind(0d0)), allocatable :: filtbefore(:),filtafter(:)
  real(kind(0d0)) :: xfwin
  real(kind(0d0)), allocatable :: ata(:,:),atd(:),atainv(:,:)
  real(kind(0d0)), allocatable  :: mtInverted(:,:,:)
  real(kind(0d0)), allocatable :: misfitTaper(:,:,:)
  real(kind(0d0)), allocatable :: misfitRaw(:,:,:)
  character(200) :: synfile,tmpfile,list
  real(kind(0d0)) :: dummyFloat
  real(kind(0d0)), allocatable :: varZ(:,:),varN(:,:),varE(:,:)
  real(kind(0d0)), allocatable :: modZ(:,:),modN(:,:),modE(:,:)
  real(kind(0d0)), allocatable :: varRawZ(:,:),varRawN(:,:),varRawE(:,:)
  real(kind(0d0)), allocatable :: modRawZ(:,:),modRawN(:,:),modRawE(:,:)
  
  real(kind(0d0)) :: fakeMT(1:6)
110 format(a200)
  ! making taper function

  call pinput
  
  !print *, np,ntwin
  allocate(taperDSM(1:npDSM))
  allocate(taperOBS(1:npData))



  !print *, iMovingWindowStart(1), iMovingWindowStart(2)
  !print *, iMovingWindowEnd(1), iMovingWindowEnd(2)
  !print *, nTimeCombination, npData, npDSM, ntStep


  allocate(ata(1:nmt,1:nmt))
  allocate(atainv(1:nmt,1:nmt))
  allocate(atd(1:nmt))
  allocate(mtInverted(1:nmt,1:nTimeCombination,1:nConfiguration))
  allocate(misfitTaper(1:nmt,1:nTimeCombination,1:nConfiguration))
  allocate(misfitRaw(1:nmt,1:nTimeCombination,1:nConfiguration))

  allocate(varZ(1:nTimeCombination,1:nConfiguration))
  allocate(varN(1:nTimeCombination,1:nConfiguration))
  allocate(varE(1:nTimeCombination,1:nConfiguration))
  allocate(modZ(1:nTimeCombination,1:nConfiguration))
  allocate(modN(1:nTimeCombination,1:nConfiguration))
  allocate(modE(1:nTimeCombination,1:nConfiguration))

  allocate(varRawZ(1:nTimeCombination,1:nConfiguration))
  allocate(varRawN(1:nTimeCombination,1:nConfiguration))
  allocate(varRawE(1:nTimeCombination,1:nConfiguration))
  allocate(modRawZ(1:nTimeCombination,1:nConfiguration))
  allocate(modRawN(1:nTimeCombination,1:nConfiguration))
  allocate(modRawE(1:nTimeCombination,1:nConfiguration))

  ! making taper for synthetics
  taperDSM=0.d0
  do iWindow=1,ntwin
     do it=itwin(1,iWindow),itwin(4,iWindow)
        if((it.gt.itwin(1,iWindow)).and.(it.lt.itwin(2,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwin(1,iWindow))/dble(itwin(2,iWindow)-itwin(1,iWindow)))
           taperDSM(it)=xfwin*xfwin
        elseif((it.ge.itwin(2,iWindow)).and.(it.le.itwin(3,iWindow))) then
           taperDSM(it)=1.d0
        elseif((it.gt.itwin(3,iWindow)).and.(it.lt.itwin(4,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwin(4,iWindow))/dble(itwin(3,iWindow)-itwin(4,iWindow)))
           taperDSM(it)=xfwin*xfwin
        endif
     enddo
  enddo
  
  ! making taper for observed data
  taperOBS=0.d0
  do iWindow=1,ntwinObs
     do it=itwinObs(1,iWindow),itwinObs(4,iWindow)
        if((it.gt.itwinObs(1,iWindow)).and.(it.lt.itwinObs(2,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwinObs(1,iWindow))/dble(itwinObs(2,iWindow)-itwinObs(1,iWindow)))
           taperOBS(it)=xfwin*xfwin
        elseif((it.ge.itwinObs(2,iWindow)).and.(it.le.itwinObs(3,iWindow))) then
           taperOBS(it)=1.d0
        elseif((it.gt.itwinObs(3,iWindow)).and.(it.lt.itwinObs(4,iWindow))) then
           xfwin=dsin(0.5d0*pi*dble(it-itwinObs(4,iWindow))/dble(itwinObs(3,iWindow)-itwinObs(4,iWindow)))
           taperOBS(it)=xfwin*xfwin
        endif
     enddo
  enddo


  ! we apply the butterworth filter to the Green's functions


  allocate(tmparray(1:npDSM,1:3,1:nmt))
  allocate(GreenArray(1:npDSM,1:3,1:nmt))
  allocate(obsArray(1:npDSM,1:3),obsRawArray(1:npDSM,1:3))
  allocate(filtbefore(1:npDSM),filtafter(1:npDSM))

  allocate(modArray(1:npDSM,1:3))
  allocate(modRawArray(1:npDSM,1:3))


  mtInverted=0.d0
  ! Grande boucle pour chaque configuration


  if(calculMode.eq.2) then
     allocate(rsgtomega(1:num_rsgtPSV,imin:imax,1:theta_n))
     allocate(u(iWindowStart:iWindowEnd,1:num_rsgtPSV))
  endif


  do iConfiguration=1,nConfiguration
    
    iConfR=(iConfiguration-1)/(phi_n*theta_n)+1
    iConfTheta=mod((iConfiguration-1),(r_n*theta_n))/phi_n+1
    iConfPhi=mod(mod((iConfiguration-1),(r_n*theta_n)),phi_n)+1
    
    print *, iConfiguration, iConfR, iConfTheta, iConfPhi, "henlo"

    

    if(iConfPhi.eq.1) then
        ! SSGT reading
        call rdsgtomega(r_(iConfR),0.d0,num_rsgtPSV,num_rsgtPSV,20)
        call rdsgtomega(r_(iConfR),0.d0,num_rsgtSH,num_rsgtPSV,10)
    endif

        ! here we construct 10 SGTs
    call vectorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomega(1:num_rsgtPSV,imin:imax,iConfTheta),u(iWindowStart:iWindowEnd,1:num_rsgtPSV),omegai,tlenFull,iWindowStart,iWindowEnd)
    stop
        ! NF redefine u and reconsider the order of rsgtomega

        ! Hey here, we can even rotate with 360 phi!!!
        ! mais bon c'est just trop couteux

     !


     ! Make a fake series of observed waveforms for synthetic inversion
     ! nothing to do with the inversion, comment out when we are on the production mode
     !if(iConfiguration.eq.1) then
     !   fakeMT=(/1.d0,2.d0,-10.d0,1.d0,3.d0,-2.d0/)
     !   open(unit=4,file="fakeZ.data.fake",status='unknown')
     !   open(unit=5,file="fakeN.data.fake",status='unknown')
     !   open(unit=6,file="fakeE.data.fake",status='unknown')
     !   do it=1,np
     !      write(4,*) dot_product(fakeMT(1:6),(tmparray(it,1,1:6)))
     !      write(5,*) dot_product(fakeMT(1:6),(tmparray(it,2,1:6)))
     !      write(6,*) dot_product(fakeMT(1:6),(tmparray(it,3,1:6)))
     !   enddo
     !   close(4)
     !   close(5)
     !   close(6)
     !endif

     ! Here we first filter Green's function as a whole and taper them
     
     do mtcomp=1,nmt
        do icomp=1,3
           filtbefore(1:npDSM)=tmparray(1:npDSM,icomp,mtcomp)
           call bwfilt(filtbefore,filtafter,dt,npDSM,1,npButterworth,fmin,fmax)
           tmparray(1:npDSM,icomp,mtcomp)=filtafter(1:npDSM)
           GreenArray(1:npDSM,icomp,mtcomp)=filtafter(1:npDSM)*taperDSM(1:npDSM)
        enddo
     enddo
     
  


     
     ! Construct AtA
     
     
     ata=0.d0
     do mtcomp=1,nmt
        do jmtcomp=1,mtcomp
           ata(mtcomp,jmtcomp)=sum(GreenArray(1:npDSM,1:3,mtcomp)*GreenArray(1:npDSM,1:3,jmtcomp))
        enddo
     enddo
     
     ! AtA is symmetric
     
     do mtcomp=1,nmt
        do jmtcomp=mtcomp,nmt
           ata(mtcomp,jmtcomp)=ata(jmtcomp,mtcomp)
        enddo
     enddo
     
     
     ! Try the inverse of ata
     
     !call inverse(ata,atainv,mtcomp)

     !do mtcomp=1,nmt
     !   do jmtcomp=1,nmt
     !      print *, mtcomp,jmtcomp,ata(mtcomp,jmtcomp),atainv(mtcomp,jmtcomp)
     !  enddo
     !enddo
     !stop
     
!!!!!!! NOW WE START MT INVERSION FOR EACH TIME iMovingWindow
     
     
     ! For the observed data we taper the whole signal of a length of npData
     
     ! For the observed data we filter the whole signal of a length of npData
     
     do icomp=1,3
        call bwfilt(obsRaw(1:npData,icomp),obsFilt(1:npData,icomp),dt,npData,1,npButterworth,fmin,fmax)
     enddo
     

     ! Raw data and filtered data are written as fort.11-13 for references
     

     open(11,file="obsZtreated.txt",status='unknown')
     open(12,file="obsNtreated.txt",status='unknown')
     open(13,file="obsEtreated.txt",status='unknown')
     do it=1,npData
        write(11,*) dble(it)*dt,obsRaw(it,1),obsFilt(it,1)
        write(12,*) dble(it)*dt,obsRaw(it,2),obsFilt(it,2)
        write(13,*) dble(it)*dt,obsRaw(it,3),obsFilt(it,3)
     enddo
     close(11)
     close(12)
     close(13)

     if(calculMode.eq.1) stop
     
     ! Here we taper the observed function for each moving window of np
     
     do iMovingWindowStep=1,nTimeCombination
        iMovingWindow=(iMovingWindowStep-1)*ntStep+1
        do icomp=1,3
           obsRawArray(1:npDSM,icomp)=obsFilt(iMovingWindow:iMovingWindow+npDSM-1,icomp)
           obsArray(1:npDSM,icomp)=obsRawArray(1:npDSM,icomp)*taperDSM(1:npDSM)
        enddo
        atd=0.d0
        ! Atd construction
        
        do mtcomp=1,nmt
           atd(mtcomp)=sum(GreenArray(1:npDSM,1:3,mtcomp)*obsArray(1:npDSM,1:3))
        enddo
       
      
        ! MT inversion by CG
        call invbyCG(nmt,ata,atd,eps,mtInverted(1:nmt,iMovingWindowStep,iConfiguration))
    
        ! MT inversion by CG for 5 components (without Mpp)

        !call invbyCG(5,ata(1:5,1:5),atd(1:5),eps,mtInverted(1:5,iMovingWindowStep,iConfiguration))
        !mtInverted(6,iMovingWindowStep,iConfiguration)= &
        !     -(mtInverted(1,iMovingWindowStep,iConfiguration) &
        !     + mtInverted(4,iMovingWindowStep,iConfiguration))
        


        !mtInverted(1:nmt,iMovingWindowStep,iConfiguration)=matmul(atainv,atd)


        ! residual evaluation with/without tapering
        modRawArray=0.d0
        modArray=0.d0
        do mtcomp=1,nmt
           modRawArray(1:npDSM,1:3)=modRawArray(1:npDSM,1:3) &
                +tmparray(1:npDSM,1:3,mtcomp)*mtInverted(mtcomp,iMovingWindowStep,iConfiguration)
           modArray(1:npDSM,1:3)=modArray(1:npDSM,1:3) &
                +GreenArray(1:npDSM,1:3,mtcomp)*mtInverted(mtcomp,iMovingWindowStep,iConfiguration)
        enddo

        write(list,'(I7,".",I7)') iConfiguration,iMovingWindowStep
        do jjj=1,15
           if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
        enddo

        tmpfile=trim(resultDir)//'/'//trim(list)//"modRaw.dat"
        open(unit=21,file=tmpfile,status='unknown')

        tmpfile=trim(resultDir)//'/'//trim(list)//"mod.dat"
        open(unit=22,file=tmpfile,status='unknown')
        
        tmpfile=trim(resultDir)//'/'//trim(list)//"obsRaw.dat"
        open(unit=23,file=tmpfile,status='unknown')

        tmpfile=trim(resultDir)//'/'//trim(list)//"obs.dat"
        open(unit=24,file=tmpfile,status='unknown')

        varZ(iMovingWindowStep,iConfiguration)=0.d0
        varN(iMovingWindowStep,iConfiguration)=0.d0
        varE(iMovingWindowStep,iConfiguration)=0.d0
        modZ(iMovingWindowStep,iConfiguration)=0.d0
        modN(iMovingWindowStep,iConfiguration)=0.d0
        modE(iMovingWindowStep,iConfiguration)=0.d0
        varRawZ(iMovingWindowStep,iConfiguration)=0.d0
        varRawN(iMovingWindowStep,iConfiguration)=0.d0
        varRawE(iMovingWindowStep,iConfiguration)=0.d0
        modRawZ(iMovingWindowStep,iConfiguration)=0.d0
        modRawN(iMovingWindowStep,iConfiguration)=0.d0
        modRawE(iMovingWindowStep,iConfiguration)=0.d0

        do it=1,npDSM

           write(21,*) dt*dble(it), modRawArray(it,1), modRawArray(it,2), modRawArray(it,3)
           write(22,*) dt*dble(it), modArray(it,1), modArray(it,2), modArray(it,3)
           write(23,*) dt*dble(it), obsRawArray(it,1), obsRawArray(it,2), obsRawArray(it,3)
           write(24,*) dt*dble(it), obsArray(it,1), obsArray(it,2), obsArray(it,3)

           varZ(iMovingWindowStep,iConfiguration)= &
                varZ(iMovingWindowStep,iConfiguration)+obsArray(it,1)**2
           varN(iMovingWindowStep,iConfiguration)= &
                varN(iMovingWindowStep,iConfiguration)+obsArray(it,2)**2
           varE(iMovingWindowStep,iConfiguration)= &
                varE(iMovingWindowStep,iConfiguration)+obsArray(it,3)**2
           modZ(iMovingWindowStep,iConfiguration)= &
                modZ(iMovingWindowStep,iConfiguration)+(modArray(it,1)-obsArray(it,1))**2
           modN(iMovingWindowStep,iConfiguration)= &
                modN(iMovingWindowStep,iConfiguration)+(modArray(it,2)-obsArray(it,2))**2
           modE(iMovingWindowStep,iConfiguration)= &
                modE(iMovingWindowStep,iConfiguration)+(modArray(it,3)-obsArray(it,3))**2
           varRawZ(iMovingWindowStep,iConfiguration)= &
                varZ(iMovingWindowStep,iConfiguration)+obsRawArray(it,1)**2
           varRawN(iMovingWindowStep,iConfiguration)= &
                varN(iMovingWindowStep,iConfiguration)+obsRawArray(it,2)**2
           varRawE(iMovingWindowStep,iConfiguration)= &
                varE(iMovingWindowStep,iConfiguration)+obsRawArray(it,3)**2
           modRawZ(iMovingWindowStep,iConfiguration)= &
                modZ(iMovingWindowStep,iConfiguration)+(modRawArray(it,1)-obsRawArray(it,1))**2
           modRawN(iMovingWindowStep,iConfiguration)= &
                modN(iMovingWindowStep,iConfiguration)+(modRawArray(it,2)-obsRawArray(it,2))**2
           modRawE(iMovingWindowStep,iConfiguration)= &
                modE(iMovingWindowStep,iConfiguration)+(modRawArray(it,3)-obsRawArray(it,3))**2

        enddo

        close(21)
        close(22)

        
        
     enddo

     
  enddo


  open(unit=1,file=trim(inversionName)//".inv_result",status='unknown')
  open(unit=2,file=trim(inversionName)//".raw_var",status='unknown')
  open(unit=3,file=trim(inversionName)//".tap_var",status='unknown')
  do iConfiguration=1,nConfiguration
     do iMovingWindowStep=1,nTimeCombination
        iMovingWindow=iMovingWindowStep*ntStep
        iMovingWindow=(iMovingWindowStep-1)*ntStep+1
        write(1,*) iConfiguration, dble(iMovingWindow-1)*dt, &
             mtInverted(1:nmt,iMovingWindowStep,iConfiguration)
        write(2,*) iConfiguration, dble(iMovingWindow-1)*dt, &
             modRawZ(iMovingWindowStep,iConfiguration),modRawN(iMovingWindowStep,iConfiguration), &
             modRawE(iMovingWindowStep,iConfiguration)
        write(3,*) iConfiguration, dble(iMovingWindow-1)*dt, &
             modZ(iMovingWindowStep,iConfiguration),modN(iMovingWindowStep,iConfiguration), &
             modE(iMovingWindowStep,iConfiguration)
     enddo
  enddo
  close(1)
  close(2)
  close(3)
     

     
  

  
  
  
end program MarsInversion
