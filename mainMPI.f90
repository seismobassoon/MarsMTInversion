!
!   mainMPI.f90
!   MarsMTInversion
!
!   Created by fuji on 25/03/2020.
!   Copyright 2020 nfuji. All rights reserved.
!
program MarsInversion
    use mpi

    implicit none
    integer :: ierr, my_rank, nproc
    integer :: tmpint

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

    if(my_rank.eq.0) then
        print *, "my_rank is", my_rank
    endif
    tmpint=1+my_rank
    call MPI_BCAST(tmpint,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    print *, "my_rank and tmpint", my_rank, tmpint
    
    stop
end program MarsInversion