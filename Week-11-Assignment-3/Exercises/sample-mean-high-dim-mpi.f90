program main
    use mpi
    implicit none
    integer :: ntasks, rank, ierr
    integer, allocatable :: status(:)
    integer :: clock, n
    integer, allocatable :: seed(:)
    integer(4) :: i, n_tot = 100000000
    real(8), parameter :: l = 0.d0, r = 1.d0, b = 0.d0, t = 1.d0
    real(8) :: x, y, z, H
    real(8) :: integral_local = 0.d0, integral

    ! initialize the MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    ! initialize the random number for different processes
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
        call func(x, y, z)
        call inArea(x, y, H)
        integral_local = integral_local + z * H
    end do
    integral_local = (r - l) * (t - b) / dble(n_tot) * integral_local

    call MPI_REDUCE(integral_local, integral, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        integral = integral / dble(ntasks)
        write(*,'(f10.5)') integral
    end if

    call MPI_FINALIZE(ierr)
end program main

subroutine func(x, y, z)
    ! the function to be integrated
    implicit none
    real(8), intent(in) :: x, y
    real(8), intent(out) :: z

    z = x * y
end subroutine func

subroutine inArea(x, y, H)
    ! judge whether the dot is in the integral area
    real(8), intent(in) :: x, y
    real(8), intent(out) :: H

    H = 1.d0
end subroutine inArea
