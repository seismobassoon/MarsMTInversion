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
    

    do iConfR=1,nr

        rsgtomega=dcmplx(0.d0)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtSH,num_rsgtPSV,10)
        call rdsgtomega(r_(iradiusD(iConfR)),num_rsgtPSV,num_rsgtPSV,20)
        rsgtomegaI=rsgtomega

        do kConfR=1,iConfR

            if(abs(r_(iradiusD(iConfR))-r_(iradiusD(kConfR)))>toleranceDistance) cycle
            rsgtomega=dcmplx(0.d0)
            call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtSH,num_rsgtPSV,10)
            call rdsgtomega(r_(iradiusD(kConfR)),num_rsgtPSV,num_rsgtPSV,20)
            rsgtomegaK=rsgtomega ! all the rsgt in freq. for kConfR depth are stored

            


        do iConfTheta=1,ntheta

            !print *,"distance is", thetaD(ithetaD(iConfTheta)),theta_n, ntheta,ithetaD(iConfTheta)
            rsgtomegatmp(1:num_rsgtPSV,imin:imax)=rsgtomega(1:num_rsgtPSV,imin:imax,ithetaD(iConfTheta))
                
            call tensorFFT_double(num_rsgtPSV,imin,imax,np1,rsgtomegatmp,rsgtTime,omegai, &
                tlenFull,iWindowStart,iWindowEnd) ! rsgtTime is for iConfR and iConfTheta
            do iConfPhi=1,nphi
                !print *, "source location is ", r_(iradiusD(iConfR)),latgeo(iConfPhi,iConfTheta), longeo(iConfPhi,iConfTheta)
                               
                iConfiguration=(iConfR-1)*(nphi*ntheta)+(iConfTheta-1)*nphi+iConfPhi
                if(mod(iConfiguration,nproc).ne.my_rank) cycle
                
                


    

    call MPI_FINALIZE(ierr)
    stop
end program MarsInversion
