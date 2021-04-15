program Parallel_Sum
    use MPI
    implicit none

    ! variable for parallel environment
    integer :: rank, ntasks, ierr
    integer :: master = 0
    integer, allocatable :: status(:)

    ! local vars
    integer :: i, divider, Sum_local, Sum_global
    
    ! prepare MPI environment
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
        write(*,*) 'MPI initialization failed'
        stop
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    divider = 100 / ntasks
    Sum_local = 0
    do i = rank * divider + 1, (rank + 1) * divider
        Sum_local = Sum_local + i
    end do

    call MPI_REDUCE(Sum_local, Sum_global, 1, MPI_INTEGER, MPI_SUM, master, MPI_COMM_WORLD, ierr)
    if (rank == master) then
        write(*,*) 'Total sum is', Sum_global
    end if

    ! clean up MPI
    deallocate(status)
    call MPI_FINALIZE(ierr)
end program Parallel_Sum

