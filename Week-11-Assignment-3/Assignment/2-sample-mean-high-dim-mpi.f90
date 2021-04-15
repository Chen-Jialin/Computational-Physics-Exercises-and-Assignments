program main
    use mpi
    implicit none
    integer :: ntasks, rank, ierr
    integer, allocatable :: status(:)
    integer :: clock, n
    integer, allocatable :: seed(:)
    integer(4) :: i, n_tot = 125000000
    real(8), parameter :: pi = acos(-1.d0)
    real(8), parameter :: phi_l = 0.d0, phi_u = 2.d0 * pi, rho_l = 0.d0, rho_u = 3.d0, z_l = -3.d0, z_u = 3.d0
    real(8) :: phi, rho, z, f, H
    real(8) :: integral_local = 0.d0, integral
    real :: start, finish, time
    real(8) :: sigma

    ! initialize the MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    call CPU_TIME(start)

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
            call RANDOM_NUMBER(z)
            clock = clock + Int(z * 1000000)
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
        call RANDOM_NUMBER(phi)
        phi = phi * (phi_u - phi_l) + phi_l
        call RANDOM_NUMBER(rho)
        rho = rho * (rho_u - rho_l) + rho_l
        call RANDOM_NUMBER(z)
        z = z * (z_u - z_l) + z_l
        call func(phi, rho, z, f)
        call inArea(phi, rho, z, H)
        integral_local = integral_local + f * H
    end do
    integral_local = (phi_u - phi_l) * (rho_u - rho_l) * (z_u - z_l) / dble(n_tot) * integral_local

    call MPI_REDUCE(integral_local, integral, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        integral = integral / dble(ntasks)
        write(*,'(f10.5)') integral
    end if

    call CPU_TIME(finish)
    call MPI_REDUCE(finish - start, time, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        write(*,'("Time consumed: ",f10.5)') time
    end if
    call MPI_REDUCE(((integral_local - dble(648) * pi / dble(5))**2)**2, sigma, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        write(*,'("Error: ",f10.5)') sqrt(sigma / dble(ntasks))
    end if

    call MPI_FINALIZE(ierr)
end program main

subroutine func(phi, rho, z, f)
    ! the function to be integrated
    implicit none
    real(8), intent(in) :: phi, rho, z
    real(8), intent(out) :: f

    f = rho**3
end subroutine func

subroutine inArea(phi, rho, z, H)
    ! judge whether the dot is in the integral area
    real(8), intent(in) :: phi, rho, z
    real(8), intent(out) :: H

    if ((z > -sqrt(9.d0 - rho**2)) .and. (z < sqrt(9.d0 - rho**2))) then
        H = 1.d0
    else
        H = 0.d0
    end if
end subroutine inArea
