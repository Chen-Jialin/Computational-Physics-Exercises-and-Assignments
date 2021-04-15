program main
    use mpi
    implicit none
    integer :: ntasks, rank, ierr
    integer, allocatable :: status(:)
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    integer(4) :: n_in = 0, n_tot = 100000000
    real(8), parameter :: l = 0.d0, r = 1.d0, b = -1.d0, t = 2.d0
    real(8) :: x, y, y_real
    real(8) :: integral_local, integral

    ! initialize the MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    if (rank == 0) then
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
            call MPI_SEND(clock, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, ierr)
        end do
    else
        call MPI_RECV(clock, 1, MPI_INTEGER, 0, rank, MPI_COMM_WORLD, status, ierr)
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = clock + 37 * i
        end do
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
    end if

    do i = 1, n_tot
        call RANDOM_NUMBER(x)
        x = x * (r - l) + l
        call RANDOM_NUMBER(y)
        y = y * (t - b) + b
        call func(x, y_real)
        if (y < y_real) then
            n_in = n_in + 1
        end if
    end do
    integral_local = (t - b) * (r - l) * dble(n_in) / dble(n_tot) + b * (r - l)

    call MPI_REDUCE(integral_local, integral, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        integral = integral / dble(ntasks)
        write(*,'(f10.5)') integral
    end if

    call MPI_FINALIZE(ierr)
end program main

subroutine func(x, y)
    ! the function to be integrated
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: y

    y = x**2
end subroutine func
