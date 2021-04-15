program main
    use mpi
    implicit none
    integer :: ntasks, id, rc
    integer, allocatable :: status(:)
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    real(8) :: x

    ! initialize the MPI environment
    call MPI_INIT(rc)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)
    allocate(status(MPI_STATUS_SIZE))

    ! initialize the random number for different processes
    if (id == 0) then
        call SYSTEM_CLOCK(clock)
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = clock + 37 * i
        end do
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
        do i = 1, ntasks - 1
            call RANDOM_NUMBER(x)
            clock = clock + Int(x * 1000000)
            call MPI_SEND(clock, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, rc)
        end do
    else
        call MPI_RECV(clock, 1, MPI_INTEGER, 0, id, MPI_COMM_WORLD, status, rc)
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = clock + 37 * i
        end do
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
    end if

    call RANDOM_NUMBER(x)
    write(*,*) x

    ! done with MPI
    call MPI_FINALIZE(rc)
end program main
