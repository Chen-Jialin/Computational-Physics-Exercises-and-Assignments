program simple
    ! required MPI include file
    use mpi
    implicit none

    integer :: numtasks, rank, ierr

    ! initialize MPI
    call MPI_INIT(ierr)

    ! get number of tasks
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    ! get my rank
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    write(*,*) 'My rank is', rank

    ! done with MPI
    call MPI_FINALIZE(ierr)

end program simple
