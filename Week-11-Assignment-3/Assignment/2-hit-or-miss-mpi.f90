program main
    use mpi
    implicit none
    integer :: ntasks, rank, ierr
    integer, allocatable :: status(:)
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    real(8), parameter :: pi = acos(-1.d0)
    integer(4) :: n_in_local = 0, n_in, n_tot = 125000000
    real(8), parameter :: phi_l = 0.d0, phi_u = 2 * pi, rho_l = 0.d0, rho_u = 3.d0, z_l = -3.d0, z_u = 3.d0,&
        f_l = 0.d0, f_u = 27.d0
    integer :: H
    real(8) :: phi, rho, z, f, f_real
    real(8) :: integral_local, integral
    real :: start, finish, time
    real(8) :: sigma

    ! initialize the MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    call CPU_TIME(start)

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
        call RANDOM_NUMBER(f)
        f = f * (f_u - f_l) + f_l
        call func(phi, rho, z, f_real)
        call inArea(phi, rho, z, H)
        if (f < f_real) then
            n_in = n_in + H
        end if
    end do
    integral_local = (phi_u - phi_l) * (rho_u - rho_l) * (z_u - z_l) * (f_u - f_l) * dble(n_in) / dble(n_tot)

    call MPI_REDUCE(integral_local, integral, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        integral = integral / dble(ntasks)
        write(*,'(f10.5)') integral
    end if

    call CPU_TIME(finish)
    call MPI_REDUCE(finish - start, time, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        write(*,'("Time consumed: "f10.5)') time
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
    integer, intent(out) :: H

    if ((z > -sqrt(9.d0 - rho**2)) .and. (z < sqrt(9.d0 - rho**2))) then
        H = 1
    else
        H = 0
    end if
end subroutine inArea
