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
    use angles
    use mainparameters

    implicit none

    call pinput

    
    !npDSM = iWindowEnd-iWindowStart+1
  
  !print *, np,ntwin
    allocate(taperDSM(iWindowStart:iWindowEnd))
    allocate(taperOBS(0:npData))


  !print *, iMovingWindowStart(1), iMovingWindowStart(2)
  !print *, iMovingWindowEnd(1), iMovingWindowEnd(2)
  !print *, nTimeCombination, npData, npDSM, ntStep


    if(calculMode.eq.2) allocate(ata(1:nmt,1:nmt))
    if(calculMode.eq.3) allocate(ata(1:nmt*nConfiguration*nTimeCombination,1:nmt*nConfiguration*nTimeCombination))
    if(calculMode.eq.4) allocate(ata(1:nmt*nTimeCombination,1:nmt*nTimeCombination))
    !allocate(atainv(1:nmt,1:nmt))
    if(calculMode.eq.2) allocate(atd(1:nmt))
    if(calculMode.eq.3) allocate(atd(1:nmt*nConfiguration*nTimeCombination))
    if(calculMode.eq.4) allocate(atd(1:nmt*nTimeCombination))
    allocate(mtInverted(1:nmt,1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.3) allocate(mtInverted_total(1:nmt*nConfiguration*nTimeCombination))
    if(calculMode.eq.4) allocate(mtInverted_total(1:nmt*nTimeCombination))
    allocate(misfitTaper(1:nmt,1:nTimeCombination,1:nConfiguration))
    allocate(misfitRaw(1:nmt,1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.2) allocate(varZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(varN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(varE(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modE(1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.10) allocate(varZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(varN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(varE(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modE(1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.4) allocate(varZ(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(varN(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(varE(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(modZ(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(modN(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(modE(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(xcorrZ(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(xcorrN(1:1,1:nConfiguration))
    if(calculMode.eq.4) allocate(xcorrE(1:1,1:nConfiguration))

    if(calculMode.eq.2) allocate(varRawZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(varRawN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(varRawE(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modRawZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modRawN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(modRawE(1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.10) allocate(varRawZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(varRawN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(varRawE(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modRawZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modRawN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.10) allocate(modRawE(1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.2) allocate(xcorrRawZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(xcorrRawN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(xcorrRawE(1:nTimeCombination,1:nConfiguration))

    if(calculMode.eq.2) allocate(xcorrZ(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(xcorrN(1:nTimeCombination,1:nConfiguration))
    if(calculMode.eq.2) allocate(xcorrE(1:nTimeCombination,1:nConfiguration))
    

    allocate(conf_depth(1:nConfiguration))
    allocate(conf_lat(1:nConfiguration))
    allocate(conf_lon(1:nConfiguration))
    allocate(conf_gcarc(1:nConfiguration))
    allocate(conf_azimuth(1:nConfiguration))

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


    ! For the observed data we taper the whole signal of a length of npData
    
    ! A priori taper for stabilising bwfilt
    taperOBS=1.d0
    do it=0,npData/20
        xfwin=dsin(0.5d0*pi*dble(it-1)/dble(npData/20-1))
        taperOBS(it)=xfwin*xfwin
    enddo
    do it=npData-npData/20,npData
        xfwin=dsin(0.5d0*pi*dble(it-npData)/dble(-npData/20))
        taperOBS(it)=xfwin*xfwin
    enddo
    
    ! For the observed data we filter the whole signal of a length of npData

    do icomp=1,3
        obsRaw(0:npData,icomp)=obsRaw(0:npData,icomp)*taperOBS(0:npData)
       call bwfilt(obsRaw(0:npData,icomp),obsFilt(0:npData,icomp),dt,npData+1,1,npButterworth,fmin,fmax)
    enddo
    

    ! Raw data and filtered data are written as fort.11-13 for references
    
    call obsFiltWriting

    ! NF will soon remove the following ASCII fils (see above)
    open(21,file="obsZtreated.txt",status='unknown')
    open(22,file="obsNtreated.txt",status='unknown')
    open(23,file="obsEtreated.txt",status='unknown')
    do it=0,npData
       write(21,*) dble(it)*dt,obsRaw(it,1),obsFilt(it,1)
       write(22,*) dble(it)*dt,obsRaw(it,2),obsFilt(it,2)
       write(23,*) dble(it)*dt,obsRaw(it,3),obsFilt(it,3)
    enddo
    close(21)
    close(22)
    close(23)

  
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

    ! we taper the filtered OBS

    ! For the observed data we filter the whole signal of a length of npData

    do icomp=1,3
        obsFiltTapered(0:npData,icomp)=obsFilt(0:npData,icomp)*taperOBS(0:npData)
    enddo
    
    



    ! GreenArray will be the filtered Green's function of 3 x 6 components

    allocate(tmparray(iWindowStart:iWindowEnd,1:3,1:nmt))
    !if(calculMode.eq.3) allocate(tmparrayI(iWindowStart:iWindowEnd,1:3,1:nmt))
    allocate(northTemp(iWindowStart:iWindowEnd))
    allocate(eastTemp(iWindowStart:iWindowEnd))
    allocate(GreenArray(iWindowStart:iWindowEnd,1:3,1:nmt))
    if(calculMode.eq.3) allocate(GreenArrayK(iWindowStart:iWindowEnd,1:3,1:nmt))
    allocate(GreenArrayShifted(0:npData,1:3,1:nmt)) ! Attention this is ok (0:npData) because we shift SYN to OBS
    allocate(GreenArrayShiftedTapered(0:npData,1:3,1:nmt))
    !! NF should think how to do this
    allocate(obsArray(0:npData,1:3),obsRawArray(0:npData,1:3),obsFiltTaperedRotated(0:npData,1:3))
    allocate(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd))

    if(calculMode.eq.2) allocate(modArray(0:npData,1:3))
    if(calculMode.eq.2) allocate(modRawArray(0:npData,1:3))
    if(calculMode.eq.10) allocate(modArray(0:npData,1:3))
    if(calculMode.eq.10) allocate(modRawArray(0:npData,1:3))
    if(calculMode.eq.3) allocate(modArray_total(iWindowStart:iWindowEnd+ntStep*(totalNumberInWindowDimension(1)-1),1:3))
    if(calculMode.eq.4) allocate(modArray(0:npData,1:3))
    allocate(rsgtomega(1:num_rsgtPSV,imin:imax,1:theta_n))
    allocate(rsgtTime(iWindowStart:iWindowEnd,1:num_rsgtPSV))
        !rsgtomega=dcmplx(0.d0)
        !u=0.d0
    allocate(rsgtomegatmp(1:num_rsgtPSV,imin:imax))

   
    

    obsArray=obsFiltTapered
    obsRawArray=obsFilt

    mtInverted=0.d0


    if(calculMode.eq.3) then

        allocate(rsgtomegaK(1:num_rsgtPSV,imin:imax,1:theta_n))
        allocate(rsgtTimeK(iWindowStart:iWindowEnd,1:num_rsgtPSV))

    endif

    

if(calculMode.eq.2) then

    do iConfR=1,nr
        print *, "depth is ", r_(iradiusD(iConfR))
        rsgtomega=dcmplx(0.d0)
        
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)
       

        
        do iConfTheta=1,ntheta
            print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
            rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
            !call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp(1:num_rsgtPSV,imin:imax),rsgtTime(iWindowStart:iWindowEnd,1:num_rsgtPSV),omegai,tlenFull,iWindowStart,iWindowEnd)
            
            call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
                tlenFull,iWindowStart,iWindowEnd)
            
           ! do iloop=1,num_rsgtPSV
            !    call vectorFFT_double(imin,imax,np1,rsgtomegatmp(iloop,imin:imax),rsgtTime(iWindowStart:iWindowEnd,iloop),omegai,tlenFull,iWindowStart,iWindowEnd)
            !enddo

            do iConfPhi=1,nphi
                print *, "source location is ", latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
                iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi
                
                conf_depth(iConfiguration)=r_(iradiusD(iConfR))
                conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
                conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
                conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
                conf_azimuth(iConfiguration)=azimuth(iConfPhi)
                call rsgt2h3time_adhoc(iConfPhi,iConfTheta)
                        
    
                !print *, "iConfR, iConfTheta,iConfPhi,iConfiguration=", iConfR, iConfTheta, iConfPhi, iConfiguration
                    
                ! Here we have to rotate from ZRT to ZNE

                if(ZRTorZNE.eq."ZNE") then
                    obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                    do mtcomp=1,nmt
                        northTemp(iWindowStart:iWindowEnd) = &
                            +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        eastTemp(iWindowStart:iWindowEnd) = &
                            +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                        tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                    enddo
                elseif(ZRTorZNE.eq."ZRT") then
                    obsFiltTaperedRotated(0:npData,2) = &
                        +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                        +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                    obsFiltTaperedRotated(0:npData,3) = &
                        +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                        -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                endif
            
                
                ! Here we first filter Green's function as a whole and taper them
                    
                do mtcomp=1,nmt
                    do icomp=1,3
                        filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                        filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                        
                        call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                            dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                        tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                        GreenArray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                        !do it=iWindowStart,iWindowEnd
                        !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                        !enddo
                   enddo
                enddo
                
                
                
                ! Here we start moving window procedure
                
                do jloop=1,nTimeCombination
                    ! This is for each time moving window set (independent or fixed)
                    iMovingWindowStep=jloop
                    GreenArrayShifted=0.d0
                    GreenArrayShiftedTapered=0.d0
                    !print  *, "now time combination", jloop
                    do iloop=1,ntwinObs
                        ! First shift the GreenArray in order to get to the same timing as OBS
                        GreenArrayShifted(itwinObs(1,iloop):itwinObs(4,iloop),1:3,1:nmt) &
                            =GreenArray(iEachWindowStart(jloop,iloop):iEachWindowEnd(jloop,iloop),1:3,1:nmt)
                    enddo
                    
                    ! Applying the same taper as OBS
                        
                    do mtcomp=1,nmt
                        do icomp=1,3
                            GreenArrayShiftedTapered(0:npData,icomp,mtcomp) &
                                =GreenArrayShifted(0:npData,icomp,mtcomp)*taperOBS(0:npData)
                        enddo
                    enddo


                    ! Construct AtA
                    
                    ata=0.d0
                    do mtcomp=1,nmt
                       do jmtcomp=1,mtcomp
                          ata(mtcomp,jmtcomp) &
                            = sum(GreenArrayShiftedTapered(0:npData,1:3,mtcomp) &
                                *GreenArrayShiftedTapered(0:npData,1:3,jmtcomp))
                       enddo
                    enddo
                    
                    ! AtA is symmetric
                    
                    do mtcomp=1,nmt
                       do jmtcomp=mtcomp,nmt
                          ata(mtcomp,jmtcomp)=ata(jmtcomp,mtcomp)
                       enddo
                    enddo

                    ! Atd construction
                    atd=0.d0
                    do mtcomp=1,nmt
                       atd(mtcomp) &
                        =sum(GreenArrayShiftedTapered(0:npData,1:3,mtcomp) &
                            *obsFiltTaperedRotated(0:npData,1:3))
                    enddo


                    ! MT inversion by CG
                    call invbyCG(nmt,ata,atd,eps,mtInverted(1:nmt,iMovingWindowStep,iConfiguration))

                    ! MT inversion by CG for 5 components (without Mpp)

                    !call invbyCG(5,ata(1:5,1:5),atd(1:5),eps,mtInverted(1:5,iMovingWindowStep,iConfiguration))
                    !mtInverted(6,iMovingWindowStep,iConfiguration)= &
                    !     -(mtInverted(1,iMovingWindowStep,iConfiguration) &
                    !     + mtInverted(4,iMovingWindowStep,iConfiguration))
                    
                    !NF: note that for the shifting window version we don't really have any difference between modRaw and mod but
                    !    I keep it deliberately so that I can use for the further use.

                    ! residual evaluation with/without tapering
                    !modRawArray=0.d0
                    modArray=0.d0
                    do mtcomp=1,nmt
                        modRawArray(0:npData,1:3)=modRawArray(0:npData,1:3) &
                            +GreenArrayShifted(0:npData,1:3,mtcomp)*mtInverted(mtcomp,iMovingWindowStep,iConfiguration)
                        modArray(0:npData,1:3)=modArray(0:npData,1:3) &
                            +GreenArrayShiftedTapered(0:npData,1:3,mtcomp)*mtInverted(mtcomp,iMovingWindowStep,iConfiguration)
                    enddo

                    write(list,'(I7,".",I7)') iConfiguration,iMovingWindowStep
                    do jjj=1,15
                       if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
                    enddo

                    tmpfile=trim(resultDir)//'/'//trim(list)//"modRaw.dat"
                    open(unit=21,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
                    
                    tmpfile=trim(resultDir)//'/'//trim(list)//"mod.dat"
                    open(unit=22,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)


                    tmpfile=trim(resultDir)//'/'//trim(list)//"obs.dat"
                    open(unit=32,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
                    
                    !tmpfile=trim(resultDir)//'/'//trim(list)//"obsRaw.dat"
                    !open(unit=23,file=tmpfile,status='unknown')

                    !tmpfile=trim(resultDir)//'/'//trim(list)//"obs.dat"
                    !open(unit=24,file=tmpfile,status='unknown')

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
                    xcorrZ(iMovingWindowStep,iConfiguration)=0.d0
                    xcorrN(iMovingWindowStep,iConfiguration)=0.d0
                    xcorrE(iMovingWindowStep,iConfiguration)=0.d0
                    xcorrRawZ(iMovingWindowStep,iConfiguration)=0.d0
                    xcorrRawN(iMovingWindowStep,iConfiguration)=0.d0
                    xcorrRawE(iMovingWindowStep,iConfiguration)=0.d0
    
                
                    normaliseModZ=0.d0
                    normaliseModRawN=0.d0
                    normaliseModE=0.d0
                    normaliseModRawZ=0.d0
                    normaliseModN=0.d0
                    normaliseModRawE=0.d0

                        
                    obsArray=obsFiltTaperedRotated
                        
                    do it=0,npData
                        tmpfloat(1)=sngl(dt*dble(it))
                        tmpfloat(2)=sngl(modRawArray(it,1))
                        tmpfloat(3)=sngl(modRawArray(it,2))
                        tmpfloat(4)=sngl(modRawArray(it,3))
                        write(21,rec=it+1) tmpfloat(1:4)
                        tmpfloat(2)=sngl(modArray(it,1))
                        tmpfloat(3)=sngl(modArray(it,2))
                        tmpfloat(4)=sngl(modArray(it,3))
                        write(22,rec=it+1) tmpfloat(1:4)
                        tmpfloat(1)=sngl(dt*dble(it))
                        tmpfloat(2)=sngl(obsArray(it,1))
                        tmpfloat(3)=sngl(obsArray(it,2))
                        tmpfloat(4)=sngl(obsArray(it,3))
                        write(32,rec=it+1) tmpfloat(1:4)
                        !write(23,*) dt*dble(it), obsRawArray(it,1), obsRawArray(it,2), obsRawArray(it,3)
                        !write(24,*) dt*dble(it), obsArray(it,1), obsArray(it,2), obsArray(it,3)

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
                        normaliseModZ=normaliseModZ+modArray(it,1)**2
                        normaliseModRawZ=normaliseModRawZ+modRawArray(it,1)**2
                        normaliseModN=normaliseModN+modArray(it,2)**2
                        normaliseModRawN=normaliseModRawN+modRawArray(it,2)**2
                        normaliseModE=normaliseModE+modArray(it,3)**2
                        normaliseModRawE=normaliseModRawE+modRawArray(it,3)**2
                        xcorrZ(iMovingWindowStep,iConfiguration)= &
                            xcorrZ(iMovingWindowStep,iConfiguration)+obsArray(it,1)*modArray(it,1)
                        xcorrN(iMovingWindowStep,iConfiguration)= &
                            xcorrN(iMovingWindowStep,iConfiguration)+obsArray(it,2)*modArray(it,2)
                        xcorrE(iMovingWindowStep,iConfiguration)= &
                            xcorrE(iMovingWindowStep,iConfiguration)+obsArray(it,3)*modArray(it,3)
                        xcorrRawZ(iMovingWindowStep,iConfiguration)= &
                            xcorrRawZ(iMovingWindowStep,iConfiguration)+obsRawArray(it,1)*modRawArray(it,1)
                        xcorrRawN(iMovingWindowStep,iConfiguration)= &
                            xcorrRawN(iMovingWindowStep,iConfiguration)+obsRawArray(it,2)*modRawArray(it,2)
                        xcorrRawE(iMovingWindowStep,iConfiguration)= &
                            xcorrRawE(iMovingWindowStep,iConfiguration)+obsRawArray(it,3)*modRawArray(it,3)
                    enddo

                    close(21)
                    close(22)
                    close(32)
                    xcorrZ(iMovingWindowStep,iConfiguration) &
                        =xcorrZ(iMovingWindowStep,iConfiguration)/sqrt(normaliseModZ) &
                        /sqrt(varZ(iMovingWindowStep,iConfiguration))
                    xcorrN(iMovingWindowStep,iConfiguration) &
                        =xcorrN(iMovingWindowStep,iConfiguration)/sqrt(normaliseModN) &
                        /sqrt(varN(iMovingWindowStep,iConfiguration))
                    xcorrE(iMovingWindowStep,iConfiguration) &
                        =xcorrE(iMovingWindowStep,iConfiguration)/sqrt(normaliseModE) &
                        /sqrt(varE(iMovingWindowStep,iConfiguration))
                    xcorrRawZ(iMovingWindowStep,iConfiguration) &
                        =xcorrRawZ(iMovingWindowStep,iConfiguration)/sqrt(normaliseModRawZ) &
                        /sqrt(varRawZ(iMovingWindowStep,iConfiguration))
                    xcorrRawN(iMovingWindowStep,iConfiguration) &
                        =xcorrRawN(iMovingWindowStep,iConfiguration)/sqrt(normaliseModRawN) &
                        /sqrt(varRawN(iMovingWindowStep,iConfiguration))
                    xcorrRawE(iMovingWindowStep,iConfiguration) &
                        =xcorrRawE(iMovingWindowStep,iConfiguration)/sqrt(normaliseModRawE) &
                        /sqrt(varRawE(iMovingWindowStep,iConfiguration))


                enddo ! for each time shift

                


            enddo ! iConfPhi loop
        enddo ! iConfTheta loop
    enddo ! iConfR loop


    open(unit=1,file=trim(inversionName)//".inv_result",status='unknown')
    open(unit=2,file=trim(inversionName)//".raw_var",status='unknown')
    open(unit=3,file=trim(inversionName)//".tap_var",status='unknown')
    open(unit=4,file=trim(inversionName)//".conf_info",status='unknown')
    open(unit=5,file=trim(inversionName)//".tap_xcorr",status='unknown')
    open(unit=6,file=trim(inversionName)//".raw_xcorr",status='unknown')
    open(unit=7,file=trim(inversionName)//".shift_info",status='unknown')
    do jloop=1,nTimeCombination
        write(7,*) jloop,fEachShift(jloop,1:ntwinObs)
    enddo
    do iConfiguration=1,nConfiguration
        write(4,*) iConfiguration, conf_depth(iConfiguration), conf_lat(iConfiguration), &
            conf_lon(iConfiguration), conf_gcarc(iConfiguration), conf_azimuth(iConfiguration)
        do jloop=1,nTimeCombination
            ! This is for each time moving window set (independent or fixed)
            iMovingWindowStep=jloop
            write(1,*) iConfiguration, iMovingWindowStep, &
                mtInverted(1:nmt,iMovingWindowStep,iConfiguration)
            write(2,*) iConfiguration, iMovingWindowStep, &
                modRawZ(iMovingWindowStep,iConfiguration),modRawN(iMovingWindowStep,iConfiguration), &
                modRawE(iMovingWindowStep,iConfiguration)
            write(3,*) iConfiguration,iMovingWindowStep,  &
                modZ(iMovingWindowStep,iConfiguration),modN(iMovingWindowStep,iConfiguration), &
                modE(iMovingWindowStep,iConfiguration)
            write(5,*) iConfiguration,iMovingWindowStep,  &
                xcorrZ(iMovingWindowStep,iConfiguration),xcorrN(iMovingWindowStep,iConfiguration), &
                xcorrE(iMovingWindowStep,iConfiguration)
            write(6,*) iConfiguration,iMovingWindowStep,  &
                xcorrRawZ(iMovingWindowStep,iConfiguration),xcorrRawN(iMovingWindowStep,iConfiguration), &
                xcorrRawE(iMovingWindowStep,iConfiguration)
        enddo
    enddo
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)
    close(7)

elseif(calculMode.eq.10) then

do iConfR=1,nr
    print *, "depth is ", r_(iradiusD(iConfR))
    rsgtomega=dcmplx(0.d0)
    
    call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
    call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)
   

    
    do iConfTheta=1,ntheta
        print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
        rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
        !call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp(1:num_rsgtPSV,imin:imax),rsgtTime(iWindowStart:iWindowEnd,1:num_rsgtPSV),omegai,tlenFull,iWindowStart,iWindowEnd)
        
        call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
            tlenFull,iWindowStart,iWindowEnd)
        
       ! do iloop=1,num_rsgtPSV
        !    call vectorFFT_double(imin,imax,np1,rsgtomegatmp(iloop,imin:imax),rsgtTime(iWindowStart:iWindowEnd,iloop),omegai,tlenFull,iWindowStart,iWindowEnd)
        !enddo

        do iConfPhi=1,nphi
            print *, "source location is ", latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
            iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi
            
            conf_depth(iConfiguration)=r_(iradiusD(iConfR))
            conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
            conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
            conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
            conf_azimuth(iConfiguration)=azimuth(iConfPhi)
            call rsgt2h3time_adhoc(iConfPhi,iConfTheta)
                    

            !print *, "iConfR, iConfTheta,iConfPhi,iConfiguration=", iConfR, iConfTheta, iConfPhi, iConfiguration
                
            ! Here we have to rotate from ZRT to ZNE
            

            if(ZRTorZNE.eq."ZNE") then
                obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                do mtcomp=1,nmt
                    northTemp(iWindowStart:iWindowEnd) = &
                        +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                    eastTemp(iWindowStart:iWindowEnd) = &
                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                        -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                    tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                    tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                enddo
            elseif(ZRTorZNE.eq."ZRT") then
                obsFiltTaperedRotated(0:npData,2) = &
                    +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                obsFiltTaperedRotated(0:npData,3) = &
                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                    -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
            endif
            
            ! Here we first filter Green's function as a whole and taper them
                
            do mtcomp=1,nmt
                do icomp=1,3
                    filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                    filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                    
                    call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                        dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                    tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                    GreenArray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                    !do it=iWindowStart,iWindowEnd
                    !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                    !enddo
                    write(list,'(I7,".",I7)') iConfiguration,mtcomp
                    do jjj=1,15
                       if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
                    enddo

                    tmpfile=trim(resultDir)//'/'//trim(list)//".syn.dat"
                    open(unit=21,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
                    do it=iWindowStart,iWindowEnd
                        tmpfloat(1)=sngl(dt*dble(it))
                        tmpfloat(2)=sngl(GreenArray(it,1,mtcomp))
                        tmpfloat(3)=sngl(GreenArray(it,2,mtcomp))
                        tmpfloat(4)=sngl(GreenArray(it,3,mtcomp))
                        write(21,rec=it+1) tmpfloat(1:4)
                    enddo
                    close(21)
                enddo
            enddo
            
            tmpfile=trim(workingDir)//'/'//trim(list)//".synZ.dat"
            open(unit=31,file=tmpfile,status='unknown',form='formatted')
            do it=iWindowStart,iWindowEnd
                write(31,*) dot_product(GreenArray(it,1,1:6),Mij_synthetic(1:6))
            enddo
            close(31)

            if(ZRTorZNE.eq."ZNE") tmpfile=trim(workingDir)//'/'//trim(list)//".synN.dat"
            if(ZRTorZNE.eq."ZRT") tmpfile=trim(workingDir)//'/'//trim(list)//".synR.dat"
            open(unit=31,file=tmpfile,status='unknown',form='formatted')
            do it=iWindowStart,iWindowEnd
                write(31,*) dot_product(GreenArray(it,2,1:6),Mij_synthetic(1:6))
            enddo
            close(31)

            if(ZRTorZNE.eq."ZNE") tmpfile=trim(workingDir)//'/'//trim(list)//".synE.dat"
            if(ZRTorZNE.eq."ZRT") tmpfile=trim(workingDir)//'/'//trim(list)//".synT.dat"
            open(unit=31,file=tmpfile,status='unknown',form='formatted')
            do it=iWindowStart,iWindowEnd
                write(31,*) dot_product(GreenArray(it,3,1:6),Mij_synthetic(1:6))
            enddo
            close(31)

        enddo
    enddo
enddo

open(unit=4,file=trim(inversionName)//".conf_info",status='unknown')
  
   do iConfiguration=1,nConfiguration
       write(4,*) iConfiguration, conf_depth(iConfiguration), conf_lat(iConfiguration), &
           conf_lon(iConfiguration), conf_gcarc(iConfiguration), conf_azimuth(iConfiguration)
   enddo
close(4)
            
  
elseif(calculMode.eq.3) then

    ata=0.d0
    atd=0.d0
    ! big ata and atd construction (maybe we should paralellise this)
    do iConfR=1,nr

        rsgtomega=dcmplx(0.d0)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)

        rsgtomegaI=rsgtomega
        
        do kConfR=1,iConfR

            rsgtomega=dcmplx(0.d0)
            call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtSH,num_rsgtPSV,10)
            call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtPSV,num_rsgtPSV,20)

            rsgtomegaK=rsgtomega ! all the rsgt in freq. for kConfR depth are stored

            do iConfTheta=1,ntheta

                !print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
                rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                
                call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
                    tlenFull,iWindowStart,iWindowEnd) ! rsgtTime is for iConfR and iConfTheta


                do kConfTheta=1,iConfTheta
                    
                    rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomegaK(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                    
                    call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTimeK,omegai, &
                             tlenFull,iWindowStart,iWindowEnd) ! rsgtTimeK is for kConfR and kConfTheta

                    do iConfPhi=1,nphi

                        print *, "source location I is ", r_(iradiusD(iConfR)),latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
                                
                        iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi
                        
                        conf_depth(iConfiguration)=r_(iradiusD(iConfR))
                        conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
                        conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
                        conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
                        conf_azimuth(iConfiguration)=azimuth(iConfPhi)

                        call rsgt2h3time_adhoc(iConfPhi,iConfTheta) ! tmparray is for iConfR, iConfTheta, iConfPhi
                        
                        
                        ! Here we have to rotate from ZRT to ZNE

                        if(ZRTorZNE.eq."ZNE") then
                            obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                            do mtcomp=1,nmt
                                northTemp(iWindowStart:iWindowEnd) = &
                                    +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                    +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                eastTemp(iWindowStart:iWindowEnd) = &
                                    +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                    -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                                tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                            enddo
                        elseif(ZRTorZNE.eq."ZRT") then
                            obsFiltTaperedRotated(0:npData,2) = &
                                +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                                +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                            obsFiltTaperedRotated(0:npData,3) = &
                                +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                                -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                        endif
                    

                

                        
                            
                        ! Here we first filter Green's function as a whole and taper them
                                
                        do mtcomp=1,nmt
                            do icomp=1,3
                                filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                                filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                                    
                                call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                                    dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                                tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                                GreenArray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                                !do it=iWindowStart,iWindowEnd
                                !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                                !enddo
                            enddo
                        enddo
                            


                        do kConfPhi=1,iConfPhi


                            print *, "source location K is ",r_(iradiusD(kConfR)), latgeo(kConfPhi,kConfTheta), longeo(kConfPhi,kConfTheta)
                            kConfiguration=(kConfR-1)*(nphi*ntheta)+(kConfTheta-1)*nphi+kConfPhi
                                    
                            rsgtTime=rsgtTimeK
                            call rsgt2h3time_adhoc(kConfPhi,kConfTheta) ! tmparray is for kConfR, kConfTheta, kConfPhi
                            
                            ! Here we have to rotate from ZRT to ZNE

                           
                            if(ZRTorZNE.eq."ZNE") then
                                obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                                do mtcomp=1,nmt
                                    northTemp(iWindowStart:iWindowEnd) = &
                                        +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                    eastTemp(iWindowStart:iWindowEnd) = &
                                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                        -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                    tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                                    tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                                enddo
                            elseif(ZRTorZNE.eq."ZRT") then
                                obsFiltTaperedRotated(0:npData,2) = &
                                    +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                                obsFiltTaperedRotated(0:npData,3) = &
                                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                                    -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                            endif
                        

                            
                            


                                
                            ! Here we first filter Green's function as a whole and taper them
                                    
                            do mtcomp=1,nmt
                                do icomp=1,3
                                    filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                                    filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                                        
                                    call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                                        dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                                    tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                                    GreenArrayK(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                                    !do it=iWindowStart,iWindowEnd
                                    !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                                    !enddo
                                enddo
                            enddo

                            ! normally all the GreenArray and GreenArrayK are fulfilled
                            do jloop=1,totalNumberInWindowDimension(1) ! we fix to the first timeshift for the K-th green
                                do jmtcomp=1,nmt
                                    do kmtcomp=1,jmtcomp
                                        iBig=(jloop-1)*nConfiguration*nmt+(iConfR-1)*nmt*nphi*ntheta &
                                            +(iConfTheta-1)*nmt*nphi+(iConfPhi-1)*nmt+jmtcomp
                                        kBig=(1-1)+(kConfR-1)*nmt*nphi*ntheta &
                                        +(kConfTheta-1)*nmt*nphi+(kConfPhi-1)*nmt+kmtcomp
                                      
                                        do it=iWindowStart+(jloop-1)*ntStep,iWindowEnd
                                            do icomp=1,3
                                                ata(iBig,kBig)= ata(iBig,kBig)+ &
                                                    GreenArray(it-(jloop-1)*ntStep,icomp,jmtcomp)* &
                                                    GreenArrayK(it,icomp,kmtcomp)
                                            enddo ! icomp
                                        enddo ! time series

                                        if(kConfR*kConfTheta*kConfPhi*kmtcomp.eq.1) then
                                            do it=iWindowStart,iWindowEnd
                                                do icomp=1,3
                                                    atd(iBig)=atd(iBig)+GreenArray(it,icomp,jmtcomp)* &
                                                                obsFiltTaperedRotated(it+(jloop-1)*ntStep,icomp)
                                                enddo
                                            enddo
                                        endif
                                    enddo ! jmtcomp
                                enddo ! mtcomp
                            enddo ! kloop: moving window
                            ! NF should put the other contributions (when I-th green is for the other timeshift)

                            do jloop=2,totalNumberInWindowDimension(1)
                                do kloop=2,jloop
                                    do jmtcomp=1,nmt
                                        do kmtcomp=1,jmtcomp
                                            iBig=(jloop-1)*nConfiguration*nmt+(iConfR-1)*nmt*nphi*ntheta &
                                                +(iConfTheta-1)*nmt*nphi+(iConfPhi-1)*nmt+jmtcomp
                                            kBig=(kloop-1)*nConfiguration*nmt+(kConfR-1)*nmt*nphi*ntheta &
                                                +(kConfTheta-1)*nmt*nphi+(kConfPhi-1)*nmt+kmtcomp
                                            iBigEquivalent=(jloop-kloop)*nConfiguration*nmt+(iConfR-1)*nmt*nphi*ntheta &
                                                +(iConfTheta-1)*nmt*nphi+(iConfPhi-1)*nmt+jmtcomp
                                            kBigEquivalent=(1-1)*nConfiguration*nmt+(kConfR-1)*nmt*nphi*ntheta &
                                                +(kConfTheta-1)*nmt*nphi+(kConfPhi-1)*nmt+kmtcomp
                                            
                                            ata(iBig,kBig)=ata(iBigEquivalent,kBigEquivalent)
                                            ! NF have to verify all above NF
                                        enddo
                                    enddo
                                enddo
                            enddo

                            !!!

                        enddo ! kConfPhi
                    enddo !iConfPhi
                enddo !kConfTheta
            enddo !iConfTheta
        enddo ! kConfR
    enddo ! iConfR


    ! ata is symmetric : fulfil the other half!


    ! ata is symmetric
    do iBig=1,nTimeCombination*nConfiguration*nmt
        do kBig=iBig,nTimeCombination*nConfiguration*nmt
            ata(iBig,kBig)=ata(kBig,iBig)
        enddo
    enddo

    ! MT inversion by CG
    call invbyCG(nTimeCombination*nConfiguration*nmt,ata,atd,eps,mtInverted_total)
    
    ! NF should verify the above equation mtInverted???
    ! then mod waveforms!!
    
    !!! Here is how we calculate the whole time series of modeled waveforms !!!!
    modArray_total=0.d0
    do iConfR=1,nr
        rsgtomega=dcmplx(0.d0)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)

        do iConfTheta=1,ntheta
            rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                       
            call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
                tlenFull,iWindowStart,iWindowEnd)

            do iConfPhi=1,nphi
                conf_depth(iConfiguration)=r_(iradiusD(iConfR))
                conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
                conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
                conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
                conf_azimuth(iConfiguration)=azimuth(iConfPhi)

                call rsgt2h3time_adhoc(iConfPhi,iConfTheta) ! tmparray is for iConfR, iConfTheta, iConfPhi
                
                
                ! Here we have to rotate from ZRT to ZNE

                if(ZRTorZNE.eq."ZNE") then
                    obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                    do mtcomp=1,nmt
                        northTemp(iWindowStart:iWindowEnd) = &
                            +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        eastTemp(iWindowStart:iWindowEnd) = &
                            +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                        tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                    enddo
                elseif(ZRTorZNE.eq."ZRT") then
                    obsFiltTaperedRotated(0:npData,2) = &
                        +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                        +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                    obsFiltTaperedRotated(0:npData,3) = &
                        +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                        -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                endif
            


                
                    
                ! Here we first filter Green's function as a whole and taper them
                        
                do mtcomp=1,nmt
                    do icomp=1,3
                        filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                        filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                            
                        call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                            dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                        tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                        GreenArray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                        !do it=iWindowStart,iWindowEnd
                        !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                        !enddo
                    enddo
                enddo
                    
                
                do jloop=1,totalNumberInWindowDimension(1)
                    do jmtcomp=1,nmt
                        iBig=(jloop-1)*nConfiguration*nmt+(iConfR-1)*nmt*nphi*ntheta &
                                +(iConfTheta-1)*nmt*nphi+(iConfPhi-1)*nmt+jmtcomp
                        do it=iWindowStart,iWindowEnd
                            do icomp=1,3
                                modArray_total(it+(jloop-1)*ntStep,icomp)= &
                                    modArray_total(it+(jloop-1)*ntStep,icomp)+ &
                                    GreenArray(it,icomp,jmtcomp)*mtInverted_total(iBig)
                            enddo ! icomp
                        enddo ! it
                    enddo ! jmtcomp
                enddo ! jloop
            enddo ! iConfPhi
        enddo ! iConfTheta
    enddo ! iConfR

    !write(list,'(I7,".",I7)') iConfiguration,iMovingWindowStep
    !do jjj=1,15
    !    if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
    !enddo

    variance_total=0.d0
    tmpfile=trim(resultDir)//'/mod_waveform_total.dat'
    open(unit=21,file=tmpfile,status='unknown')
    tmpfile=trim(resultDir)//'/obs_waveform_total.dat'
    open(unit=22,file=tmpfile,status='unknown')
    do it=iWindowStart,iWindowEnd+ntStep*(totalNumberInWindowDimension(1)-1)
        write(21,*) dt*dble(it),modArray_total(it,1),modArray_total(it,2),modArray_total(it,3)
        write(22,*) dt*dble(it),obsFiltTapered(it,1),obsFiltTapered(it,2),obsFiltTapered(it,3)
        do icomp=1,3
            variance_total(icomp)=variance_total(icomp)+dt*(modArray_total(it,icomp)-obsFiltTapered(it,icomp))**2
        enddo
    enddo ! it
    close(21)


    mtInverted=0.d0

    open(unit=1,file=trim(inversionName)//".inv_result",status='unknown')
   ! open(unit=2,file=trim(inversionName)//".raw_var",status='unknown')
   ! open(unit=3,file=trim(inversionName)//".tap_var",status='unknown')
    open(unit=4,file=trim(inversionName)//".conf_info",status='unknown')
   ! open(unit=5,file=trim(inversionName)//".tap_xcorr",status='unknown')
   ! open(unit=6,file=trim(inversionName)//".raw_xcorr",status='unknown')
    open(unit=7,file=trim(inversionName)//".shift_info",status='unknown')
    do jloop=1,totalNumberInWindowDimension(1)
        write(7,*) jloop,dt*dble(iMovingWindowStart(1))*dble(ntStep)*dble(jloop)
    enddo

    do iConfR=1,nr
        do iConfTheta=1,ntheta
            do iConfPhi=1,nphi
                do jloop=1,totalNumberInWindowDimension(1)
                    iConfiguration=(jloop-1)*nConfiguration+(iConfR-1)*nphi*ntheta &
                        +(iConfTheta-1)*nphi+(iConfPhi-1)
                    do jmtcomp=1,nmt
                        iBig=(jloop-1)*nConfiguration*nmt+(iConfR-1)*nmt*nphi*ntheta &
                            +(iConfTheta-1)*nmt*nphi+(iConfPhi-1)*nmt+jmtcomp
                        mtInverted(jmtcomp,jloop,iConfiguration)=mtInverted_total(iBig)
                    enddo ! jmtcomp
                enddo ! jloop
            enddo ! iConfPhi
        enddo ! iConfTheta
    enddo ! iConfR

    

    do iConfiguration=1,nConfiguration
        write(4,*) iConfiguration, conf_depth(iConfiguration), conf_lat(iConfiguration), &
            conf_lon(iConfiguration), conf_gcarc(iConfiguration), conf_azimuth(iConfiguration)
        do jloop=1,totalNumberInWindowDimension(1)
            ! This is for each time moving window set (independent or fixed)
            iMovingWindowStep=jloop

            write(1,*) iConfiguration, iMovingWindowStep, &
                mtInverted(1:nmt,iMovingWindowStep,iConfiguration)
           ! write(2,*) iConfiguration, iMovingWindowStep, &
           !     modRawZ(iMovingWindowStep,iConfiguration),modRawN(iMovingWindowStep,iConfiguration), &
           !        modRawE(iMovingWindowStep,iConfiguration)
           !    write(3,*) iConfiguration,iMovingWindowStep,  &
           !        modZ(iMovingWindowStep,iConfiguration),modN(iMovingWindowStep,iConfiguration), &
           !        modE(iMovingWindowStep,iConfiguration)
           !    write(5,*) iConfiguration,iMovingWindowStep,  &
           !        xcorrZ(iMovingWindowStep,iConfiguration),xcorrN(iMovingWindowStep,iConfiguration), &
           !        xcorrE(iMovingWindowStep,iConfiguration)
           !    write(6,*) iConfiguration,iMovingWindowStep,  &
           !        xcorrRawZ(iMovingWindowStep,iConfiguration),xcorrRawN(iMovingWindowStep,iConfiguration), &
           !        xcorrRawE(iMovingWindowStep,iConfiguration)
        enddo
    enddo
    close(1)
    !   close(2)
    !   close(3)
    close(4)
    !   close(5)
    !   close(6)
    close(7)

elseif(calculMode.eq.4) then

mtInverted=0.d0
 ! big ata and atd construction (maybe we should paralellise this)
 do iConfR=1,nr

    rsgtomega=dcmplx(0.d0)
    call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
    call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)
    
    do iConfTheta=1,ntheta

        !print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
        rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
             
        call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
            tlenFull,iWindowStart,iWindowEnd) ! rsgtTime is for iConfR and iConfTheta
        do iConfPhi=1,nphi
            !print *, "source location is ", r_(iradiusD(iConfR)),latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
                            
            iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi
                     
            conf_depth(iConfiguration)=r_(iradiusD(iConfR))
            conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
            conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
            conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
            conf_azimuth(iConfiguration)=azimuth(iConfPhi)

            call rsgt2h3time_adhoc(iConfPhi,iConfTheta) ! tmparray is for iConfR, iConfTheta, iConfPhi
                     
            ! Here we have to rotate from ZRT to ZNE
            if(ZRTorZNE.eq."ZNE") then
                obsFiltTaperedRotated(0:npData,1:3)=obsFiltTapered(0:npData,1:3)
                do mtcomp=1,nmt
                    northTemp(iWindowStart:iWindowEnd) = &
                        +cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                    eastTemp(iWindowStart:iWindowEnd) = &
                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                        -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                    tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                    tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                enddo
            elseif(ZRTorZNE.eq."ZRT") then
                obsFiltTaperedRotated(0:npData,2) = &
                    +cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
                obsFiltTaperedRotated(0:npData,3) = &
                    +sqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,2) &
                    -cqr(iConfPhi,iConfTheta)*obsFiltTapered(0:npData,3)
            endif

                     
                         
            ! Here we first filter Green's function as a whole and taper them
                             
            do mtcomp=1,nmt
                do icomp=1,3
                    filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                    filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                                 
                    call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                        dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                    tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                    GreenArray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                    !do it=iWindowStart,iWindowEnd
                    !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                    !enddo
                enddo
            enddo

            ata=0.d0
            atd=0.d0

            ! normally all the GreenArray and GreenArrayK are fulfilled
            do jloop=1,nTimeCombination
                do jmtcomp=1,nmt
                    iBig=(jloop-1)*nmt+jmtcomp
                    do icomp=1,3
                        do it=iWindowStart,iWindowEnd
                            atd(iBig)=atd(iBig)+GreenArray(it,icomp,jmtcomp)* &
                                obsFiltTaperedRotated(it+(jloop-1)*ntStep,icomp)
                        enddo
                    enddo
                    do kmtcomp=1,jmtcomp
                        kBig=kmtcomp
                        do icomp=1,3
                            do it=iWindowStart+(jloop-1)*ntStep,iWindowEnd
                                ata(iBig,kBig)= ata(iBig,kBig)+ &
                                    GreenArray(it-(jloop-1)*ntStep,icomp,jmtcomp)* &
                                    GreenArray(it,icomp,kmtcomp)
                            enddo ! icomp
                        enddo ! time series
                    enddo ! jmtcomp
                enddo ! jmtcomp
            enddo ! jloop: moving window
                
            ! NF should put the other contributions (when I-th green is for the other timeshift)

            do jloop=2,nTimeCombination
                do kloop=2,jloop
                    do jmtcomp=1,nmt
                        do kmtcomp=1,jmtcomp
                            iBig=(jloop-1)*nmt+jmtcomp
                            kBig=(kloop-1)*nmt+kmtcomp
                            iBigEquivalent=(jloop-kloop)*nmt+jmtcomp
                            kBigEquivalent=kmtcomp
                                         
                            ata(iBig,kBig)=ata(iBigEquivalent,kBigEquivalent)
                            ! NF have to verify all above NF
                        enddo
                    enddo
                enddo
            enddo !jloop for at compilation
            ! ata is symmetric : fulfil the other half!


            ! ata is symmetric
            do iBig=1,nTimeCombination*nmt
                do kBig=iBig,nTimeCombination*nmt
                    ata(iBig,kBig)=ata(kBig,iBig)
                enddo
            enddo



            mtInverted_total=0.d0
            ! MT inversion by CG
            call invbyCG(nTimeCombination*nmt,ata,atd,eps,mtInverted_total)
            
            modArray=0.d0
            do jloop=1,nTimeCombination
                do jmtcomp=1,nmt
                    iBig=(jloop-1)*nmt+jmtcomp

                    mtInverted(jmtcomp,jloop,iConfiguration)=mtInverted_total(iBig)
                    do it=iWindowStart,iWindowEnd
                        modArray(it+(jloop-1)*ntStep,jmtcomp)=modArray(it+(jloop-1)*ntStep,jmtcomp) &
                            +GreenArray(it,icomp,jmtcomp)*mtInverted_total(iBig)
                    enddo
                enddo
            enddo
            write(list,'(I7)') iConfiguration
            do jjj=1,7
                if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
            enddo

            tmpfile=trim(resultDir)//'/'//trim(list)//".mod.dat"
            open(unit=22,file=tmpfile,status='unknown')
                            
            iMovingWindowStep=1

            varZ(iMovingWindowStep,iConfiguration)=0.d0
            varN(iMovingWindowStep,iConfiguration)=0.d0
            varE(iMovingWindowStep,iConfiguration)=0.d0
            modZ(iMovingWindowStep,iConfiguration)=0.d0
            modN(iMovingWindowStep,iConfiguration)=0.d0
            modE(iMovingWindowStep,iConfiguration)=0.d0
                           
            xcorrZ(iMovingWindowStep,iConfiguration)=0.d0
            xcorrN(iMovingWindowStep,iConfiguration)=0.d0
            xcorrE(iMovingWindowStep,iConfiguration)=0.d0
                            
            
            normaliseModZ=0.d0
                           
            normaliseModE=0.d0
        
            normaliseModN=0.d0
                     
                                
            do it=0,npData
                write(22,*) dt*dble(it), modArray(it,1), modArray(it,2), modArray(it,3)
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
               
                normaliseModZ=normaliseModZ+modArray(it,1)**2
                normaliseModN=normaliseModN+modArray(it,2)**2
                normaliseModE=normaliseModE+modArray(it,3)**2
                               
                xcorrZ(iMovingWindowStep,iConfiguration)= &
                    xcorrZ(iMovingWindowStep,iConfiguration)+obsArray(it,1)*modArray(it,1)
                xcorrN(iMovingWindowStep,iConfiguration)= &
                    xcorrN(iMovingWindowStep,iConfiguration)+obsArray(it,2)*modArray(it,2)
                xcorrE(iMovingWindowStep,iConfiguration)= &
                    xcorrE(iMovingWindowStep,iConfiguration)+obsArray(it,3)*modArray(it,3)
            enddo

            close(22)

            xcorrZ(iMovingWindowStep,iConfiguration) &
                =xcorrZ(iMovingWindowStep,iConfiguration)/sqrt(normaliseModZ) &
                /sqrt(varZ(iMovingWindowStep,iConfiguration))
            xcorrN(iMovingWindowStep,iConfiguration) &
                =xcorrN(iMovingWindowStep,iConfiguration)/sqrt(normaliseModN) &
                /sqrt(varN(iMovingWindowStep,iConfiguration))
            xcorrE(iMovingWindowStep,iConfiguration) &
                =xcorrE(iMovingWindowStep,iConfiguration)/sqrt(normaliseModE) &
                /sqrt(varE(iMovingWindowStep,iConfiguration))
        enddo !iConfPhi
    enddo !iConfTheta
 enddo ! iConfR

open(unit=1,file=trim(inversionName)//".inv_result",status='unknown')
!open(unit=2,file=trim(inversionName)//".raw_var",status='unknown')
open(unit=3,file=trim(inversionName)//".tap_var",status='unknown')
open(unit=4,file=trim(inversionName)//".conf_info",status='unknown')
open(unit=5,file=trim(inversionName)//".tap_xcorr",status='unknown')
!open(unit=6,file=trim(inversionName)//".raw_xcorr",status='unknown')
open(unit=7,file=trim(inversionName)//".shift_info",status='unknown')
do jloop=1,nTimeCombination
    write(7,*) jloop,dt*dble(iMovingWindowStart(1))*dble(ntStep)*dble(jloop)
enddo
do iConfiguration=1,nConfiguration
    write(4,*) iConfiguration, conf_depth(iConfiguration), conf_lat(iConfiguration), &
        conf_lon(iConfiguration), conf_gcarc(iConfiguration), conf_azimuth(iConfiguration)
    do jloop=1,nTimeCombination
        ! This is for each time moving window set (independent or fixed)
        iMovingWindowStep=jloop
        write(1,*) iConfiguration, iMovingWindowStep, &
            mtInverted(1:nmt,iMovingWindowStep,iConfiguration)
       ! write(2,*) iConfiguration, iMovingWindowStep, &
        !    modRawZ(iMovingWindowStep,iConfiguration),modRawN(iMovingWindowStep,iConfiguration), &
        !    modRawE(iMovingWindowStep,iConfiguration)
        write(3,*) iConfiguration,iMovingWindowStep,  &
            modZ(iMovingWindowStep,iConfiguration),modN(iMovingWindowStep,iConfiguration), &
            modE(iMovingWindowStep,iConfiguration)
        write(5,*) iConfiguration,iMovingWindowStep,  &
            xcorrZ(iMovingWindowStep,iConfiguration),xcorrN(iMovingWindowStep,iConfiguration), &
            xcorrE(iMovingWindowStep,iConfiguration)
        !write(6,*) iConfiguration,iMovingWindowStep,  &
        !    xcorrRawZ(iMovingWindowStep,iConfiguration),xcorrRawN(iMovingWindowStep,iConfiguration), &
        !    xcorrRawE(iMovingWindowStep,iConfiguration)
    enddo
enddo
close(1)
!close(2)
close(3)
close(4)
close(5)
!close(6)
close(7)
 

 tmpfile=trim(resultDir)//'/mod_waveform_total.dat'
 open(unit=21,file=tmpfile,status='unknown')
 tmpfile=trim(resultDir)//'/obs_waveform_total.dat'
 open(unit=22,file=tmpfile,status='unknown')
 do it=0,npData
     write(21,*) dt*dble(it),modArray(it,1),modArray(it,2),modArray(it,3)
     write(22,*) dt*dble(it),obsFiltTaperedRotated(it,1),obsFiltTapered(it,2),obsFiltTapered(it,3)
     
 enddo ! it
 close(21)
 close(22)


endif
  
  
end program MarsInversion
