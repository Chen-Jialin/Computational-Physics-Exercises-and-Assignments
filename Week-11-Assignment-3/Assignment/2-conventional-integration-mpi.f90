program main
    use mpi
    implicit none
    integer :: ntasks, rank, ierr
    integer, allocatable :: status(:)
    integer, parameter :: n = 500
    integer :: i, j, k
    real(8), parameter :: pi = acos(-1.d0)
    real(8), parameter :: phi_l = 0.d0, phi_u = 2 * pi, rho_l = 0.d0, rho_u = 3.d0, z_l = -3.d0, z_u = 3.d0
    real(8), parameter :: d_phi = (phi_u - phi_l) / dble(n), d_rho = (rho_u - rho_l) / dble(n), d_z = (z_u - z_l) / dble(n)
    real(8) :: phi, rho, z, f, H
    real(8) :: integral_local = 0.d0, integral
    real :: start, finish, time


    ! initialize the MPI environment
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    allocate(status(MPI_STATUS_SIZE))

    call CPU_TIME(start)

    phi = phi_l - d_phi / dble(2)
    do i = 1, n
        phi = phi + d_phi
        rho = rho_l - d_rho + d_rho / dble(ntasks) / dble(2) + d_rho / dble(ntasks) * dble(rank)
        do j = 1, n
            rho = rho + d_rho
            z = z_l - d_z / dble(2)
            do k = 1, n
                z = z + d_z
                call func(phi, rho, z, f)
                call inArea(phi, rho, z, H)
                integral_local = integral_local + f * H
            end do
        end do
    end do
   
    call MPI_REDUCE(integral_local, integral, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
        integral = (phi_u - phi_l) * (rho_u - rho_l) * (z_u - z_l) / dble(ntasks) / dble(n**3) * integral
        write(*,'(f10.5)') integral
    end if

    call CPU_TIME(finish)
    call MPI_REDUCE(finish - start, time, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (rank == 0) then
        write(*,'(f10.5)') time
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
