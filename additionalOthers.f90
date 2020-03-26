!
!   additionalOthers.f90
!   MarsMTInversion
!
!   Created by fuji on 26/03/2020.
!   Copyright 2020 nfuji. All rights reserved.
!
subroutine sendAllParameters
    use mpi
    use parameters
    use angles
    implicit none
    integer :: ierr,my_rank
    integer :: iloop

    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_BCAST(calculMode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    nmt=6
    call MPI_BCAST(SGTinfo,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(parentDir,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(eventName,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(stationName,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(stla,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(stla,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rlat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rlon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(gcarcmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(gcarcmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dgcarc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntheta,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(my_rank.ne.0) then
        allocate(gcarc(1:ntheta))
        allocate(ithetaD(1:ntheta))
        do iloop=1,ntheta
            gcarc(iloop) = gcarcmin + dble(iloop-1)*dgcarc
        enddo
    endif
    call MPI_BCAST(azimuthmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(azimuthmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dazimuth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nphi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(my_rank.ne.0) then
        allocate(azimuth(1:nphi))
        do iloop=1,nphi
            azimuth(iloop) = azimuthmin + dble(iloop-1)*dazimuth
        enddo
    endif
    call MPI_BCAST(radiusmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(radiusmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dradius,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(my_rank.ne.0) then
        allocate(radius(1:nr))
        allocate(iradiusD(1:nr))
        do iloop=1,nr
            radius(iloop) = radiusmin + dble(iloop-1)*dradius
        enddo
    endif
    call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(samplingHz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tlenData,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(npData,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(movingWindowStep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntStep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(npButterworth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(fmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(start,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(end,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntwin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(my_rank.ne.0) then
        allocate(twin(1:4,1:ntwin))
        allocate(itwin(1:4,1:ntwin))
    endif
    call MPI_BCAST(twin,4*ntwin,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(itwin,4*ntwin,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(my_rank.ne.0) then
        allocate(obsRaw(0:npData,1:3))
        allocate(obsFilt(0:npData,1:3))
        allocate(obsFiltTapered(0:npData,1:3))
    endif
    call MPI_BCAST(obsRaw,3*(npData+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ntwinObs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    NmovingWindowDimension=1

    if(my_rank.ne.0) then
        allocate(twinObs(1:4,1:ntwinObs))
        allocate(itwinObs(1:4,1:ntwinObs))
    endif
    call MPI_BCAST(twinObs,4*ntwinObs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(itwinObs,4*ntwinObs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
    call MPI_BCAST(DSMconfFile,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(PoutputDir,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(psvmodel,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(modelname,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tlenFull,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rmin_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rmax_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rdelta_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(r0min,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(r0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(r0delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(thetamin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(thetamax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(thetadelta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(imin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(imax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(rsgtswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tsgtswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(synnswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ratc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ratl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(omegai,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(maxlmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(np0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(np1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(lsmooth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dtn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iWindowStart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iWindowEnd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
    if(my_rank.ne.0) then
        call makingIndependentWindow
    endif

    INFO_TSGT = trim(parentDir)//"/INFO_TSGT.TXT"
    INFO_RSGT = trim(parentDir)//"/INFO_RSGT.TXT"
    rsampletxt = trim(parentDir)//"/rsample.txt"
    modelcard = trim(parentDir)//"/"//trim(modelname)//".card"

    call MPI_BCAST(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(theta_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(my_rank.ne.0) then
        allocate(r_(1:r_n))
        allocate(thetaD(1:theta_n))
        
    endif

    call MPI_BCAST(r_,r_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(thetaD,theta_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        


        ! check r(:) and r_(:) because we don't interpolate for depth
        iradiusD=0
        do iloop=1,nr
            icheck=0
            do jloop=1,r_n
                if(dabs(r_(jloop)-radius(iloop))< 1.d-5*radius(iloop)) then
                    icheck=1
                    iradiusD(iloop)=jloop
                endif
            enddo
            if(icheck.eq.0) then
                print *, radius(iloop), "is not in the catalogue, sorry"
                stop
            endif
        enddo
        
        ! check gcarc(:) and thetaD(:) because we don't interpolate for depth
        ithetaD=0
        do iloop=1,ntheta
            icheck=0
            do jloop=1,theta_n
                if(dabs(thetaD(jloop)-gcarc(iloop))<1.d-5*gcarc(iloop)) then
                    icheck=1
                    ithetaD(iloop)=jloop
                endif
            enddo
            if(icheck.eq.0) then
                print *, gcarc(iloop), "is not in the catalogue, sorry"
                stop
            endif
        enddo
     
        allocate(latgeo(1:nphi,1:ntheta),longeo(1:nphi,1:ntheta))
        allocate(crq(1:nphi,1:ntheta),crq2(1:nphi,1:ntheta))
        allocate(srq(1:nphi,1:ntheta),srq2(1:nphi,1:ntheta))
        allocate(cqr(1:nphi,1:ntheta),sqr(1:nphi,1:ntheta))
        latgeo=0.d0
        longeo=0.d0
    
        do iloop=1,ntheta
            do jloop=1,nphi
                call geoCoordinates(stla*degree2radian,stlo*degree2radian, &
                    latgeo(jloop,iloop),longeo(jloop,iloop),azimuth(jloop)*degree2radian, &
                    fdummy,gcarc(iloop)*degree2radian) !!! geoCooridnates uses radian
            enddo
        enddo

        latgeo=latgeo*radian2degree
        longeo=longeo*radian2degree

        print *, "InSight: ", stla, stlo

        do iloop=1,ntheta
            do jloop=1,nphi
                !print *, gcarc(iloop),azimuth(jloop), latgeo(jloop,iloop),longeo(jloop,iloop)
                !azimuth(jloop)=-azimuth(jloop)
                call azimth(0,latgeo(jloop,iloop),longeo(jloop,iloop),stla,stlo,fdummy2,phirq,fdummy)
               
                !print *, " inverse: ",fdummy2,phirq,fdummy
                
                !print *, phirq,fdummy
                phirq=pi-phirq*degree2radian
                fdummy=fdummy*degree2radian
           
                if(phirq.lt.0.d0) phirq=phirq+2.d0*pi
                if(phirq.gt.(2.d0*pi)) phirq=phirq-2.d0*pi
                
                crq(jloop,iloop)=dcos(phirq)
                srq(jloop,iloop)=dsin(phirq)
                crq2(jloop,iloop)=dcos(2.d0*phirq)
                srq2(jloop,iloop)=dsin(2.d0*phirq)
                cqr(jloop,iloop)=dcos(fdummy)
                sqr(jloop,iloop)=dsin(fdummy)
            enddo
        enddo

      
      


       



        ! allocate SGTs, synthetics in frequency
        allocate(omega(imin:imax))
        do iloop = imin, imax
            omega(iloop) = 2.d0*pi*dble(iloop)/tlenFull
        enddo
     
        !nConfiguration=r_n*theta_n*phi_n

        nConfiguration=nr*ntheta*nphi
        print *, "source locations to be tested (in r, theta, phi, total): "
        print *, "      ", nr,ntheta,nphi,nConfiguration

    
    
end subroutine sendAllParameters
