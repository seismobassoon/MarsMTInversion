!
!   mainMPI.f90
!   MarsMTInversion
!
!   Created by fuji on 25/03/2020.
!   Copyright 2020 nfuji. All rights reserved.
!
program MarsInversion
    use mpi
    use parameters
    use tmpSGTs
    use angles
    use mainparameters

    implicit none

    ! Here we put only the parameters required for MPI: otherwise we declare all the variables in mainparameters

    integer :: ierr, my_rank, nproc
    integer :: tmpint

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

    if(my_rank.eq.0) then
        print *, "my_rank is", my_rank
        call pinput
    endif

    call sendAllParameters  ! sending all the necessary parameters just after pinput
    call allocatingLocalArrays
    call preprocessing

    if(my_rank.eq.0) then
        call obsFiltWriting
    endif

    
    mtInverted_total=0.d0
    mtInverted_total_previous_iteration=0.d0 ! this is important to calculate synthetics for the previous solution
    
    do lIteration=0,NumberIteration
        
        if(lIteration.ne.0) then
            call MPI_BCAST(mtInverted_total(1:nmt*nTimeCombination*nConfiguration),nmt*nTimeCombination*nConfiguration, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            mtInverted_total_previous_iteration=mtInverted_total  ! Since we prefer the Gauss-Siedel to Jacobi,
                                                                  !  we store this previous solution for variance calculation but
                                                                  !  otherwise mtInverted_total is updated for every iConfiguration
        endif

        modArray_total=0.d0
        modArray_local=0.d0
        do iConfR=1,nr

            rsgtomega=dcmplx(0.d0)
            call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
            call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)
            rsgtomegaI=rsgtomega
        

            do iConfTheta=1,ntheta
            
                !print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
                rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomegaI(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
        
                call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
                    tlenFull,iWindowStart,iWindowEnd) ! rsgtTime is for iConfR and iConfTheta

                do iConfPhi=1,nphi

                
                   
                    
                    iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi

                    if(mod(iConfiguration,nproc).ne.my_rank) cycle
                    !print *, "my_rank:", my_rank, "; source location I is "
                    !print *, r_(iradiusD(iConfR)),latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
            
                    ! update the mtInverted_total for each iConfiguration
                    if(lIteration.ne.0) then
                        call MPI_BCAST(mtInverted_total(1:nmt*nTimeCombination*nConfiguration),nmt*nTimeCombination*nConfiguration, &
                            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                    endif


                    conf_depth(iConfiguration)=r_(iradiusD(iConfR))
                    conf_lat(iConfiguration)=latgeo(iConfPhi,iConfTheta)
                    conf_lon(iConfiguration)=longeo(iConfPhi,iConfTheta)
                    conf_gcarc(iConfiguration)=thetaD(ithetaD(iConfTheta))
                    conf_azimuth(iConfiguration)=azimuth(iConfPhi)

                    call rsgt2h3time_adhoc(iConfPhi,iConfTheta) ! tmparray is for iConfR, iConfTheta, iConfPhi
            

                    ! Here we have to rotate from ZRT to ZNE

                    do mtcomp=1,nmt
                        northTemp(iWindowStart:iWindowEnd) = &
                            -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        eastTemp(iWindowStart:iWindowEnd) = &
                            -sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                            -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                        tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                        tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                    enddo

            
                
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
            
                    ! Here is the diagonal part of Jacobi method

                    ata=0.d0
                    atd=0.d0

                    ! normally all the GreenArray and GreenArrayK are fulfilled
                    do jloop=1,nTimeCombination
                        do jmtcomp=1,nmt
                            iBig=(jloop-1)*nmt+jmtcomp
                            do icomp=1,3
                                do it=iWindowStart,iWindowEnd
                                    atd(iBig)=atd(iBig)+GreenArray(it,icomp,jmtcomp)* &
                                        obsFiltTapered(it+(jloop-1)*ntStep,icomp)
                                    modArray_local(it+(jloop-1)*ntStep,icomp)= &
                                        modArray_local(it+(jloop-1)*ntStep,icomp)+ &
                                        GreenArray(it,icomp,jmtcomp)*mtInverted_total_previous_iteration(iBig)
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



                    ! Here we are dealing with non-diagonal parts of Hessian matrix to conduct Jacobi method
                    ! The first iteration does not need to pass this step
                    
                    do kConfR=1,nr !iConfR
                        if(lIteration.eq.0) cycle
                        if(abs(r_(iradiusD(iConfR))-r_(iradiusD(kConfR)))>toleranceDistance) cycle
                        rsgtomega=dcmplx(0.d0)
                        call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtSH,num_rsgtPSV,10)
                        call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtPSV,num_rsgtPSV,20)
                        rsgtomegaK=rsgtomega ! all the rsgt in freq. for kConfR depth are stored
                    
                        do kConfTheta=1,ntheta !kConfTheta

                            ! we don't take into account the cross-talks between points
                            distanceKm=sqrt(r_(iradiusD(iConfR))**2+r_(iradiusD(kConfR))**2+ &
                                2.d0*r_(iradiusD(iConfR))*r_(iradiusD(kConfR))* &
                                dcos(degree2radian*(gcarc(iConfTheta)-gcarc(kConfTheta))) )
                            if(distanceKm>toleranceDistance) cycle
                                      
                            rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomegaK(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                                  
                            call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTimeK,omegai, &
                                tlenFull,iWindowStart,iWindowEnd) ! rsgtTimeK is for kConfR and kConfTheta




                            do kConfPhi=1,nphi !kConfPhi
                                kConfiguration=(kConfR-1)*(nphi*ntheta)+(kConfTheta-1)*nphi+kConfPhi

                                ata_nondiagonal=0.d0

                                if(iConfiguration.eq.kConfiguration) cycle  ! We are just dealing with non-diagonal parts
                                ! we don't take into account the cross-talks between points
                                distanceKm=sqrt(r_(iradiusD(iConfR))**2+r_(iradiusD(kConfR))**2+ &
                                    2.d0*r_(iradiusD(iConfR))*r_(iradiusD(kConfR))* &
                                    (dsin(gcarc(iConfTheta)*degree2radian)*dsin(gcarc(kConfTheta)*degree2radian) * &
                                    dcos(degree2radian*(azimuth(iConfPhi)-azimuth(kConfPhi))) + &
                                    dcos(gcarc(iConfTheta)*degree2radian)*dcos(gcarc(kConfTheta)*degree2radian)))
                                if(distanceKm>toleranceDistance) cycle
                                print *, "source location K is ",r_(iradiusD(kConfR)), latgeo(kConfPhi,kConfTheta), &
                                    longeo(kConfPhi,kConfTheta)
                                rsgtTime=rsgtTimeK
                                call rsgt2h3time_adhoc(kConfPhi,kConfTheta) ! tmparray is for kConfR, kConfTheta, kConfPhi
                                          
                                ! Here we have to rotate from ZRT to ZNE

                                do mtcomp=1,nmt
                                    northTemp(iWindowStart:iWindowEnd) = &
                                        -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                        +sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                    eastTemp(iWindowStart:iWindowEnd) = &
                                        -sqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,2,mtcomp) &
                                        -cqr(iConfPhi,iConfTheta)*tmparray(iWindowStart:iWindowEnd,3,mtcomp)
                                    tmparray(iWindowStart:iWindowEnd,2,mtcomp)=northTemp(iWindowStart:iWindowEnd)
                                    tmparray(iWindowStart:iWindowEnd,3,mtcomp)=eastTemp(iWindowStart:iWindowEnd)
                                enddo
                      
                                ! Here we first filter Green's function as a whole and taper them
                                                    
                                do mtcomp=1,nmt
                                    do icomp=1,3
                                        filtbefore(iWindowStart:iWindowEnd)=tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)
                                        filtbefore(iWindowStart:iWindowEnd)=filtbefore(iWindowStart:iWindowEnd)*taperDSM(iWindowStart:iWindowEnd)
                                                          
                                        call bwfilt(filtbefore(iWindowStart:iWindowEnd),filtafter(iWindowStart:iWindowEnd), &
                                            dt,iWindowEnd-iWindowStart+1,1,npButterworth,fmin,fmax)
                                        tmparray(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                                        GreenArrayK(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)
                                    enddo
                                enddo
                    
                                ! here is the A_{ij} m_j for Jacobi since ata_{ij} is not symmetric ... we have to pay attentiion
                                
                                ! normally all the GreenArray and GreenArrayK are fulfilled
                                do jloop=1,nTimeCombination
                                    do jmtcomp=1,nmt
                                        iBig=(jloop-1)*nmt+jmtcomp
                                            
                                        do kmtcomp=1,jmtcomp
                                            kBig=kmtcomp
                                            do icomp=1,3
                                                do it=iWindowStart+(jloop-1)*ntStep,iWindowEnd
                                                    ata_nondiagonal(iBig,kBig)= ata_nondiagonal(iBig,kBig)+ &
                                                        GreenArray(it-(jloop-1)*ntStep,icomp,jmtcomp)* &
                                                        GreenArrayK(it,icomp,kmtcomp)
                                                enddo ! icomp
                                            enddo ! time series
                                        enddo ! kmtcomp

                                        do kmtcomp=1,jmtcomp
                                            kBig=kmtcomp
                                            do icomp=1,3
                                                do it=iWindowStart+(jloop-1)*ntStep,iWindowEnd
                                                    ata_nondiagonal(kBig,iBig)= ata_nondiagonal(kBig,iBig)+ &
                                                        GreenArrayK(it-(jloop-1)*ntStep,icomp,kmtcomp)* &
                                                        GreenArray(it,icomp,jmtcomp)
                                                enddo ! icomp
                                            enddo ! time series
                                        enddo ! kmtcomp
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
                                                                
                                                ata_nondiagonal(iBig,kBig)=ata_nondiagonal(iBigEquivalent,kBigEquivalent)
                                                ata_nondiagonal(kBig,iBig)=ata_nondiagonal(kBigEquivalent,iBigEquivalent)
                                                ! NF have to verify all above NF
                                            enddo
                                        enddo
                                    enddo
                                enddo !jloop for at compilation

                                ! gradient direction modification using Jacobian method
                                atd(1:nmt*nTimeCombination) &
                                    =atd(1:nmt*nTimeCombination)+ &
                                    matmul(ata_nondiagonal, &
                                    mtInverted_total(nmt*nTimeCombination*(kConfiguration-1)+1:nmt*nTimeCombination*kConfiguration))

                            enddo  !kConfPhi
                        enddo !kConfTheta
                    enddo !kConfR
    


                
                    mtInverted_local=0.d0
                    ! MT inversion by CG
                    call invbyCG(nTimeCombination*nmt,ata,atd,eps,mtInverted_local)
                    
                    
                    ! Our strategy here is rather Gauss-Seidel globally and locally (inside iConfiguration)
                    !     rather Jacobi method
                    call MPI_GATHER(mtInverted_local(1:nmt*nTimeCombination),nmt*nTimeCombination, &
                        MPI_DOUBLE_PRECISION,&
                        mtInverted_total(nmt*nTimeCombination*(iConfiguration-1)+1:nmt*nTimeCombination*iConfiguration), &
                        nmt*nTimeCombination,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                    !mtInverted_total(nmt*nTimeCombination*(iConfiguration-1)+1:nmt*nTimeCombination*iConfiguration)&
                    !    =mtInverted_local(1:nmt*nTimeCombination)
                    
                    
                    
                enddo ! iConfPhi
            enddo ! iConfTheta
        enddo ! iConfR
        ! we should wait for all the iConf to finished
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! summing up the modified waveforms
        call MPI_REDUCE(modArray_local,modArray_total, &
            ((iWindowEnd-iWindowStart+1)+ntStep*(nTimeCombination-1))*3, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        
        if(my_rank.eq.0) then
            !! write data and calculate xcorr and variance
            write(list,'(I7)') lIteration
            do jjj=1,7
                if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
            enddo
            tmpfile=trim(resultDir)//'/'//trim(list)//"."//trim(modelname)//".mod"
            open(unit=22,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
            if(lIteration.eq.0) then
                tmpfile=trim(resultDir)//'/data.obs'
                open(unit=23,file=tmpfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*4)
            endif
            ampsq_total=0.d0
            variance_total=0.d0
            xcorr_total=0.d0
            do it=iWindowStart,iWindowEnd+ntStep*(nTimeCombination-1)
                tmpfloat(1)=sngl(dt*dble(it))
                tmpfloat(2)=sngl(modArray_total(it,1))
                tmpfloat(3)=sngl(modArray_total(it,2))
                tmpfloat(4)=sngl(modArray_total(it,3))
                write(22,rec=it+1-iWindowStart) tmpfloat(1:4)
                if(lIteration.eq.0) then
                    tmpfloat(1)=sngl(dt*dble(it))
                    tmpfloat(2)=sngl(obsFiltTapered(it,1))
                    tmpfloat(3)=sngl(obsFiltTapered(it,2))
                    tmpfloat(4)=sngl(obsFiltTapered(it,3))
                    write(23,rec=it+1-iWindowStart) tmpfloat(1:4)
                endif
                do icomp=1,3
                    variance_total(icomp)=variance_total(icomp)+dt*(modArray_total(it,icomp)-obsFiltTapered(it,icomp))**2
                    ampsq_total(icomp,0)=ampsq_total(icomp,0)+obsFiltTapered(it,icomp)**2
                    ampsq_total(icomp,1)=ampsq_total(icomp,1)+modArray_total(it,icomp)**2
                    xcorr_total(icomp)=xcorr_total(icomp)+obsFiltTapered(it,icomp)*modArray_total(it,icomp)
                enddo
            enddo
            if(lIteration.ne.0) then
                do icomp=1,3
                    xcorr_total(icomp)=xcorr_total(icomp)/sqrt(ampsq_total(icomp,0))/sqrt(ampsq_total(icomp,1))
                enddo
            endif
            close(22)
            if(lIteration.eq.0) close(23)
            if(lIteration.eq.0) then
                open(unit=10,file=trim(inversionName)//".var",status='unknown',form='formatted')
                open(unit=11,file=trim(inversionName)//".xcorr",status='unknown',form='formatted')
            else
                open(unit=10,file=trim(inversionName)//".var",status = 'old',access='append', form = 'formatted')
                open(unit=11,file=trim(inversionName)//".xcorr",status = 'old',access='append', form = 'formatted')
            endif
            write(10,*) lIteration,variance_total(1:3)
            write(11,*) lIteration,xcorr_total(1:3)
            close(10)
            close(11)
            ! so we write the variance and xcorr here

            ! write inversion result

            allocate(mtInverted_total_single(1:nmt*nTimeCombination*nConfiguration))
            mtInverted_total_single=sngl(mtInverted_total)
            write(list,'(I7)') lIteration
            do jjj=1,7
                if(list(jjj:jjj).eq.' ') list(jjj:jjj)='0'
            enddo
            tmpfile=trim(resultDir)//'/'//trim(list)//"."//trim(modelname)//".whole_inv"
            open(unit=21,file=tmpfile,status='unknown',form='unformatted',access='direct', &
                recl=kind(0e0)*nmt*nTimeCombination*nConfiguration)
            write(21,rec=1) mtInverted_total_single
            close(21)
            deallocate(mtInverted_total_single)
        endif
        
    enddo ! lIteration
                
    if(my_rank.eq.0) then

        open(unit=4,file=trim(inversionName)//".conf_info",status='unknown')
        !open(unit=5,file=trim(inversionName)//".tap_xcorr",status='unknown')
        !open(unit=6,file=trim(inversionName)//".raw_xcorr",status='unknown')
        open(unit=7,file=trim(inversionName)//".shift_info",status='unknown')
        do jloop=1,nTimeCombination
            write(7,*) jloop,dt*dble(iMovingWindowStart(1))*dble(ntStep)*dble(jloop)
        enddo
        do iConfiguration=1,nConfiguration
            write(4,*) iConfiguration, conf_depth(iConfiguration), conf_lat(iConfiguration), &
                conf_lon(iConfiguration), conf_gcarc(iConfiguration), conf_azimuth(iConfiguration)
        enddo
        close(4)
        close(7)

    endif

    call MPI_FINALIZE(ierr)
    stop
end program MarsInversion
