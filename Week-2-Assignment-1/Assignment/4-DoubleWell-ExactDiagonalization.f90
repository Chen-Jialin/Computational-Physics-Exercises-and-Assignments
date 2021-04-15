program main
    ! solve the Schrodinger equation of a particle in an infinite potential well with central barrier with Exact Diagonalization

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: pi = acos(-1.d0)

    ! local vars
    integer :: i, j, Nx, LWORK, INFO
    real(dp) :: Xmax, dx, V0
    real(dp), allocatable :: x(:), Vpot(:), Ham(:,:), basis(:,:), W(:), WORK(:)

    integer :: LWORK1, INFO1
    real(dp), allocatable :: x1(:), Vpot1(:), Ham1(:,:), W1(:), WORK1(:)

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
    allocate(x1(-Nx - 1:Nx + 1))
    allocate(Vpot1(-Nx - 1:Nx + 1))

    dx = Xmax / dble(Nx)
    do i = -Nx, Nx
        x(i) = dx * dble(i)
        Vpot(i) = .5d0 * ((x(i)**2 - 1.d0)**2 / 4.d0 - x(i)**2)
    end do

    ! calculate the wavefunction of 1D HO
    do i = -Nx - 1, Nx + 1
        x1(i) = dx * dble(i)
        Vpot1(i) = .5d0 * x1(i)**2
    end do

    allocate(Ham1(2 * Nx + 1, 2 * Nx + 1))
    Ham1 = 0.d0
    do i = 1, 2 * Nx + 1
        Ham1(i,i) = Vpot1(i - Nx) + 1.d0 / dx**2
        if (i > i) Ham1(i,i - 1) = -.5d0 / dx**2
        if (i < 2 * Nx + 1) Ham1(i,i + 1) = -.5d0 / dx**2
    end do

    allocate(W1(2 * Nx + 1))
    LWORK1 = 10 * (Nx + 1)
    allocate(WORK1(LWORK1))
    call DSYEV('V','U',2 * Nx + 1, Ham1, 2 * Nx + 1, W1, WORK1, LWORK1, INFO1)

    allocate(basis(2 * Nx + 1, 2 * Nx +1))
    if (INFO1 == 0) then
        do i = 1, 2 * Nx + 1
            do j = 1, 2 * Nx + 1
                basis(i,j) = Ham1(j,i)
            end do
        end do
    else
        write(*, *) 'Diagonalization went wrong! aborting ...'
        stop
    end if

    allocate(Ham(2 * Nx + 1, 2 * Nx + 1))
    Ham = 0.d0
    do i = 1, 2 * Nx + 1
        do j = i, 2 * Nx + 1
            call Simpson(2 * Nx + 1, basis(i, 1:2 * Nx + 1), basis(j, 1:2 * Nx + 1), Vpot(-Nx:Nx) - .5d0, dx, V0)
            Ham(i, j) = V0
            if (i == j) then
                Ham(i, j) = Ham(i, j) + (i - .5d0)
            end if
            Ham(j, i) = Ham(i, j)
        end do
    end do

    allocate(W(2 * Nx + 1))
    LWORK = 10 * Nx
    allocate(WORK(LWORK))
    call DSYEV('V', 'U', 2 * Nx + 1, Ham, 2 * Nx + 1, W, WORK, LWORK, INFO)

    if (INFO == 0) then
        open(unit = 1, file = '4-energy-ExactDiagonalization.txt', status = 'unknown')
        open(unit = 2, file = '4-wavefunction-ExactDiagonalization.txt', status = 'unknown')

        ! multiply the coefficient (which is stored in Ham now) with the basis function to get the wave function
        Ham = matmul(transpose(basis), Ham)

        ! output results
        do i = 1, 2 * Nx + 1
            write(1, '(i4,f20.12)') i, W(i)
            write(2, '(f12.8)', advance = 'no') x(i - Nx -1)
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
