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

    call sendAllParameters

    if(calculMode.ne.5) then
        print *, "MPI version accepts only heavyMonitor option."
        stop
    endif
    
    tmpint=1+my_rank
    call MPI_BCAST(tmpint,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    print *, "my_rank and tmpint", my_rank, tmpint






    call MPI_FINALIZE(ierr)
    stop
end program MarsInversion
