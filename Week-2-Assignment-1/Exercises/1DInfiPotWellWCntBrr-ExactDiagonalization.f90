program main
    ! solve the Schrodinger equation of a particle in an infinite potential well with central barrier with Exact Diagonalization

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: pi = acos(-1.d0)

    ! local vars
    integer :: i, j, Nx, LWORK, INFO
    real(dp) :: Xmax, dx, V0
    real(dp), allocatable :: x(:), Vpot(:), Ham(:,:), basis(:,:), W(:), WORK(:)

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
        if (abs(x(i)) <= .5d0) then
            Vpot(i) = 10.0d0
        else
            Vpot(i) = 0.d0
        end if
    end do

    allocate(basis(2 * Nx + 1, 2 * Nx +1))
    basis = 0.d0
    do i = 1, 2 * Nx + 1    ! order of basis function
        if (mod(i, 2) /= 0) then
            ! i is odd
            do j = 1, 2 * Nx + 1    ! coordinate index
                basis(i,j) = sqrt(1.d0 / Xmax) * cos(dble(i) * pi * x(j - Nx - 1) / 2.d0 / Xmax)
            end do
        else
            ! i is even
            do j = 1, 2 * Nx + 1    ! coordinate index
                basis(i,j) = sqrt(1.d0 / Xmax) * sin(dble(i) * pi * x(j - Nx - 1) / 2.d0 / Xmax)
            end do
        end if
    end do

    allocate(Ham(2 * Nx + 1, 2 * Nx + 1))
    Ham = 0.d0
    do i = 1, 2 * Nx + 1
        do j = i, 2 * Nx + 1
            call Simpson(2 * Nx + 1, basis(i, 1:2 * Nx + 1), basis(j, 1:2 * Nx + 1), Vpot(-Nx:Nx), dx, V0)
            Ham(i, j) = V0
            if (i == j) then
                Ham(i, j) = Ham(i, i) + (dble(i) * pi)**2 / 8.d0 / Xmax**2
            end if
            Ham(j, i) = Ham(i, j)
        end do
    end do

    allocate(W(2 * Nx + 1))
    LWORK = 10 * Nx
    allocate(WORK(LWORK))
    call DSYEV('V', 'U', 2 * Nx + 1, Ham, 2 * Nx + 1, W, WORK, LWORK, INFO)

    if (INFO == 0) then
        open(unit = 1, file = 'eigenvalue.dat', status = 'unknown')
        open(unit = 2, file = 'wavefunction.dat', status = 'unknown')

        ! multiply the coefficient (which is stored in Ham now) with the basis function to get the wave function
        Ham = matmul(transpose(basis), Ham)

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
    deallocate(basis, W, WORK)
    deallocate(x, Vpot)

end program main

subroutine Simpson(N, u, v, Vpot, dx, V0)
    ! calculate the integral of the product of u, v, and Vpot from -Xmax to Xmax with compound simpson formula

    implicit none
    integer :: N
    real(8), intent(in) :: u(N), v(N), Vpot(N), dx
    real(8), intent(out) :: V0

    ! local vars
    integer :: i
    real(8) :: Work(N)

    if (mod(N, 2) == 0) then
    ! Simpson does not work for integral of array with even number of elements
        write(*, *) 'Array with even elements, Simpson does not know how to work!'
    end if

    do i = 1, N
        Work(i) = Vpot(i) * u(i) * v(i)
    end do
     
    V0 = Work(1) + Work(N)
    do i = 2, N - 1, 2
        V0 = V0 + 4.d0 * Work(i)
    end do
    do i = 3, N - 2, 2
        V0 = V0 + 2.d0 * Work(i)
    end do
    V0 = V0 * dx / 3.d0

end subroutine Simpson
