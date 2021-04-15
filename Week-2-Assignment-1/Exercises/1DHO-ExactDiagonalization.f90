program main
    ! solve the Schrodinger equation of a harmonic oscillator with exact diagonalization

    implicit none
    integer, parameter :: dp = selected_real_kind(8)

    ! local vars
    integer :: i, j, Nx, LWORK, INFO
    real(dp) :: Xmax, dx
    real(dp), allocatable :: x(:), Vpot(:), Ham(:,:), W(:), WORK(:), y(:)

    ! executable
    Xmax = -1.d0
    do while (Xmax < 0.d0)
        write(*,*) 'Please indicate the range of x, [-Xmax,Xmax]; Xmax ='
        read(*,*) Xmax
    end do

    Nx = -1
    do while (Nx <= 0)
        write(*,*) 'Please specify the # of slices between [0,Xmax], Nx ='
        read(*,*) Nx
    end do

    ! memory dynamical allocation
    allocate(x(-Nx:Nx))
    allocate(Vpot(-Nx:Nx))

    dx = Xmax / dble(Nx)
    do i = -Nx, Nx
        x(i) = dx * dble(i)
        Vpot(i) = .5d0 * x(i)**2
    end do

    allocate(Ham(2 * Nx - 1, 2 * Nx - 1))
    Ham = 0.d0
    do i = 1, 2 * Nx - 1
        Ham(i,i) = Vpot(i - Nx) + 1.d0 / dx**2
        if (i > 1) Ham(i,i - 1) = -.5d0 / dx**2
        if (i < 2 * Nx - 1) Ham(i,i + 1) = -.5d0 / dx**2
    end do

    allocate(W(2 * Nx - 1))
    LWORK = 10 * (Nx - 1)
    allocate(WORK(LWORK))
    call DSYEV('V', 'U', 2 * Nx - 1, Ham, 2 * Nx - 1, W, WORK, LWORK, INFO)

    if (INFO == 0) then
        open(unit = 1, file = 'eigenvalue2.txt', status = 'unknown')
        open(unit = 2, file = 'wavefunction2.txt', status = 'unknown')

        ! output results
        do i = 1, 2 * Nx + 1
            write(1, '(i4,f20.12)') i, W(i)
            write(2, '(f12.6)', advance = 'no') x(i - Nx -1)
            do j = 1, 2 * Nx + 1
                write(2, '(f12.6)', advance = 'no') Ham(i,j)
            end do
            write(2, *)
        end do
        close(1)
        close(2)
    else
        write(*, *) 'Diagonalization went wrong! aborting ...'
        stop
    end if

    deallocate(Ham)
    deallocate(W, WORK)
    deallocate(x, Vpot)

end program main
