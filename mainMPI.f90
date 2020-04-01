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

    mtInverted=0.d0
    mtInverted_total=0.d0
    
    do lIteration=0,NumberIteration

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

                
                    print *, "source location I is ", r_(iradiusD(iConfR)),latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
                    
                    iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi

                    if(mod(iConfiguration,nproc).ne.my_rank) cycle
            
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
                    
                    



                    ! MT inversion by CG
                    call invbyCG(nTimeCombination*nConfiguration*nmt,ata,atd,eps,mtInverted_local)
                    






            do kConfR=1,nr !iConfR

                if(abs(r_(iradiusD(iConfR))-r_(iradiusD(kConfR)))>toleranceDistance) cycle
                rsgtomega=dcmplx(0.d0)
                call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtSH,num_rsgtPSV,10)
                call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtPSV,num_rsgtPSV,20)
                rsgtomegaK=rsgtomega ! all the rsgt in freq. for kConfR depth are stored

            
            


                do kConfTheta=1,ntheta !iConfTheta

                    ! we don't take into account the cross-talks between points
                    distanceKm=sqrt(r_(iradiusD(iConfR))**2+r_(iradiusD(kConfR))**2+ &
                        2.d0*r_(iradiusD(iConfR))*r_(iradiusD(kConfR))* &
                        dcos(degree2radian*(gcarc(iConfTheta)-gcarc(kConfTheta))) )
                    if(distanceKm>toleranceDistance) cycle
                    
                    rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomegaK(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                
                    call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTimeK,omegai, &
                        tlenFull,iWindowStart,iWindowEnd) ! rsgtTimeK is for kConfR and kConfTheta




                        do kConfPhi=1,nphi !iConfPhi

                            ! we don't take into account the cross-talks between points
                            distanceKm=sqrt(r_(iradiusD(iConfR))**2+r_(iradiusD(kConfR))**2+ &
                                2.d0*r_(iradiusD(iConfR))*r_(iradiusD(kConfR))* &
                                (dsin(gcarc(iConfTheta)*degree2radian)*dsin(gcarc(kConfTheta)*degree2radian) * &
                                dcos(degree2radian*(azimuth(iConfPhi)-azimuth(kConfPhi))) + &
                                dcos(gcarc(iConfTheta)*degree2radian)*dcos(gcarc(kConfTheta)*degree2radian)))
                            if(distanceKm>toleranceDistance) cycle

                           


                            print *, "source location K is ",r_(iradiusD(kConfR)), latgeo(kConfPhi,kConfTheta), longeo(kConfPhi,kConfTheta)
                            kConfiguration=(kConfR-1)*(nphi*ntheta)+(kConfTheta-1)*nphi+kConfPhi
                            
                            if((lIteration.eq.0).and.(iConfiguration.ne.kConfiguration)) cycle
                            
                            
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
                                    GreenArrayK(iWindowStart:iWindowEnd,icomp,mtcomp)=filtafter(iWindowStart:iWindowEnd)! *taperDSM(1:npDSM)


                                    !do it=iWindowStart,iWindowEnd
                                    !    write(15,*) GreenArray(it,1,mtcomp), GreenArray(it,2,mtcomp),GreenArray(it,3,mtcomp)
                                    !enddo
                                enddo
                            enddo

                            
                            if(iConfiguration.eq.kConfiguration) then
                                
                                

                            else
                                    
                                ! here is the A_{ij} m_j for Jacobi since ata_{ij} is not symmetric ... we have to pay attentiion
                                ! it's not completed here
                                ! NF
                                ! normally all the GreenArray and GreenArrayK are fulfilled
                                do jloop=1,nTimeCombination
                                    do jmtcomp=1,nmt
                                        iBig=(jloop-1)*nmt+jmtcomp
                                        do icomp=1,3
                                            do it=iWindowStart,iWindowEnd
                                                atd(iBig)=atd(iBig)+GreenArray(it,icomp,jmtcomp)* &
                                                    obsFiltTapered(it+(jloop-1)*ntStep,icomp)
                                            enddo
                                        enddo
                                        do kmtcomp=1,jmtcomp
                                            kBig=kmtcomp
                                            do icomp=1,3
                                                do it=iWindowStart+(jloop-1)*ntStep,iWindowEnd
                                                    ata(iBig,kBig)= ata(iBig,kBig)+ &
                                                        GreenArray(it-(jloop-1)*ntStep,icomp,jmtcomp)* &
                                                        GreenArrayK(it,icomp,kmtcomp)
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


                                   ! ata is symmetric for iConfiguration = kConfiguration
                                do iBig=1,nTimeCombination*nmt
                                    do kBig=iBig,nTimeCombination*nmt
                                           ata(iBig,kBig)=ata(kBig,iBig)
                                    enddo
                                enddo
                                

                                
                            endif




enddo ! lIteration
                


    

    call MPI_FINALIZE(ierr)
    stop
end program MarsInversion
