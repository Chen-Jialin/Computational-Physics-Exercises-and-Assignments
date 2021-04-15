program main
    ! solve the Schrodinger equation of a 1-dimension harmonic oscillator with Numerov's method

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: eps = 1.d-5

    ! local vars
    real(dp) :: Xmax
    integer :: Nx, i
    real(dp) :: E    ! trial energy
    real(dp) :: dx, ySquareSum, dE = 1.d-2
    real(dp), allocatable :: x(:), Vpot(:), y(:), g(:), f(:)
    logical :: yNSignChange = .false., yNSignPrev
    integer :: looptime = 1

    ! executable
    Xmax = -1.d0
    do while (Xmax < 0.d0)
        write(*,*) 'Please indicate the range of x, [-Xmax, Xmax]; Xmax ='
        read(*,*) Xmax
    end do

    Nx = -1
    do while (Nx <= 0)
        write(*,*) 'Please specify the # of slices between [0,Xmax], Nx ='
        read(*,*) Nx
    end do

    write(*,*) "Please give a trial energy E0, the first allowed state with E>E0 will be calculated! E0 ="
    read(*,*) E

    write(*, '(a10,a20,a20)') 'Iteration', 'Energy', 'Boundary value'

    ! memory dynamical allocation
    allocate(x(-Nx:Nx))
    allocate(Vpot(-Nx:Nx))
    allocate(y(-Nx:Nx))
    allocate(g(-Nx:Nx))
    allocate(f(-Nx:Nx))

    dx = Xmax / dble(Nx)
    Vpot = 0.d0
    do i = -Nx, Nx
        x(i) = dx * i
        Vpot(i) = .5d0 * x(i)**2
    end do

    do while (.true.)
        ! set boundary condition at y(-Nx) and y(-Nx+1)
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
        write(*,'(i10,2f20.8)') looptime, E, y(Nx)

        ! change E
        if (abs(y(Nx)) >= eps) then    ! if precision requirement not satisfied
            if (.not. yNSignChange) then    ! if the sign of yN not changed
                if (looptime == 1) then    ! if first loop
                    yNSignPrev = (y(Nx) >= 0)
                elseif ((y(Nx) >= 0) .neqv. yNSignPrev) then    ! if not first loop and the sign of yN changed
                    yNSignChange = .true.
                    dE = - dE / 2.d0
                    yNSignPrev = (y(Nx) >= 0)
                end if
            else    ! once the sign of yN changed
                if ((y(Nx) >= 0) .neqv. yNSignPrev) then
                    yNSignPrev = (y(Nx) >= 0)
                    dE = -dE / 2.d0
                else
                    dE = dE / 2.d0
                end if
            end if
            E = E + dE
        else    ! if precision requirement satisfied
            ! output to file
            open(1, file = 'wavefunction1.txt', status = 'unknown')
            do i = -Nx, Nx
                write(1, '(2f20.8)') x(i), y(i)
            end do
            close(1)
            open(2, file = 'energy1.txt', status = 'unknown')
            write(2, *) E
            close(2)
            exit
        end if

        looptime = looptime + 1
    end do
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
