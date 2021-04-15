program main
    ! solve the Schrodinger equation of a particle in 1-dimension infinite square potential well with Numerov's method

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: eps = 1.d-5, dE = 1.d-2

    ! local vars
    real(dp) :: L = 1.d0
    integer :: Nx = 100, i, n = 1, nmax = 100
    real(dp) :: E = 7.d0    ! trial energy, lower than the first energy level
    real(dp) :: dx, ySquareSum
    real(dp), allocatable :: x(:), Vpot(:), y(:), g(:), f(:)
    logical :: yNSignChange, yNSign
    integer :: looptime
    real(dp) :: Eprev1, Eprev2, yprev1, yprev2

    ! executable
    ! write(*, '(a10,a20,a20)') 'Iteration', 'Energy', 'Boundary value'

    ! memory dynamical allocation
    allocate(x(-Nx:Nx))
    allocate(Vpot(-Nx:Nx))
    allocate(y(-Nx:Nx))
    allocate(g(-Nx:Nx))
    allocate(f(-Nx:Nx))

    ! discretization and define the potential Vpot
    dx = L / dble(Nx)
    do i = -Nx, Nx
        x(i) = dx * i
        if (abs(x(i)) <= 0.5) then
            Vpot(i) = 10.d0
        else
            Vpot(i) = 0.d0
        end if
    end do

    open(2, file = '3-energy-Numerov.txt', status = 'unknown')
    do while (n <= nmax)
        ! write(*,*) n
        looptime = 1
        yNSignChange = .false.
        do while (.true.)
            ! set boundary condition: y(-Nx) and y(-Nx+1)
            y(-Nx) = 0.d0
            y(-Nx + 1) = 1.d-4

            do i = -Nx, Nx
                g(i) = 2.d0 * (E - Vpot(i))
                f(i) = 1.d0 + g(i) / 12.d0 * dx**2
            end do

            ! Numerov's method
            do i = -Nx + 1, Nx - 1
                y(i + 1) = ((12.d0 - 10.d0 * f(i)) * y(i) - f(i - 1) * y(i - 1)) / f(i + 1)
            end do

            ! renormalize y(x)
            call Simpson(2 * Nx + 1, y(:), dx, ySquareSum)
            y = y / sqrt(ySquareSum)

            ! output to screen
            ! write(*,'(i10,2f20.8)') looptime, E, y(Nx)

            ! change E
            if (abs(y(Nx)) < eps) then    ! if precision requirement is satisfied
                ! output to file
                ! open(1, file = '3-wavefunction.txt', status = 'unknown')
                ! do i = -Nx, Nx
                    ! write(1, '(2f20.8)') x(i), y(i)
                ! end do
                ! close(1)
                write(2, '(i4,f20.8)') n, E
                exit
            else
                if (.not. yNSignChange) then    ! if the sign of yN has not changed
                    E = E + dE
                    if (looptime == 1) then    ! if this is the first loop
                        yNSign = (y(Nx) >= 0)
                        yprev1 = y(Nx)
                    else
                        if ((y(Nx) >=0) .neqv. (yNSign)) then    ! if the sign of yN changed this time
                            yNSignChange = .true.
                            E = E - dE
                            Eprev1 = E
                            Eprev2 = E - dE
                            E = E - dE / 2.d0
                        end if
                        yprev2 = yprev1
                        yprev1 = y(Nx)
                    end if
                else    ! once the sign of yN has changed
                    if (abs(yprev1) <= abs(yprev2)) then
                        Eprev2 = E
                        ! Eprev1 = Eprev1
                        E = (Eprev1 + Eprev2) / 2.d0
                        yprev2 = y(Nx)
                        ! yprev1 = yprev1
                    else
                        ! Eprev2 = Eprev2
                        Eprev1 = E
                        E = (Eprev1 + Eprev2) / 2.d0
                        ! yprev2 = yprev2
                        yprev1 = y(Nx)
                    end if
                end if
            end if
            looptime = looptime + 1
        end do
        ! write(*,*)
        E = E + 2 * dE
        n = n + 1
    end do
    close(2)
end program main

subroutine Simpson(N, y, dx, ySquareSum)
    ! calculate the integral of the square of the wavefunction from -Xmax to Xmax with compound simpson formula

    implicit none
    integer :: N
    real(8), intent(in) :: y(N), dx
    real(8), intent(out) :: ySquareSum

    ! local vars
    integer :: i
    real(8) :: Work(N)

    if (mod(N,2) == 0) then
    ! Simpson does not work for integral of array with even number of elements
        write(*, *) 'Array with even elements, Simpson does not know how to work!'
    end if

    do i = 1,N
        Work(i) = y(i)**2
    end do

    ySquareSum = Work(1) + Work(N)
    do i = 2, N - 1, 2
        ySquareSum = ySquareSum + 4.d0 * Work(i)
    end do
    do i = 3, N - 2, 2
        ySquareSum = ySquareSum + 2.d0 * Work(i)
    end do
    ySquareSum = ySquareSum * dx / 3.d0

end subroutine Simpson
