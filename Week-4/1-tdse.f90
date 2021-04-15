program main
    ! solve the time-dependent schrodinger equation

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: pi = acos(-1.d0)

    ! local vars
    integer :: Jmax = 200, Nmax = 100, j, n
    real(dp) :: Xmax = 10.d0, dx, Tmax = 2.d0, dt, psiSquareSum
    complex(dp) :: alpha
    real(dp), allocatable :: x(:)
    complex(dp), allocatable :: psi(:,:), M1(:,:), M2(:,:), M(:,:)

    integer :: LWORK, INFO
    complex(dp), allocatable :: IPIV(:), WORK(:)

    integer :: i
    complex(dp), allocatable :: a(:), b(:), c(:), d(:), w(:)

    ! excutable
    ! discretization
    dx = Xmax / dble(Jmax)
    dt = Tmax / dble(Nmax)
    allocate(x(0:Jmax))
    allocate(psi(0:Jmax,0:Nmax))
    do j = 0,Jmax
        x(j) = dble(j) * dx
        psi(j,0) = (1.d0 / pi)**(1.d0 / 4.d0) * exp(cmplx(0.d0,1.d0) * 5.d0 * x(j) - (x(j) - 5.d0)**2 / 2)
    end do
    psi(0,0) = 0.d0
    psi(Jmax,0) = 0.d0
    call Simpson(Jmax + 1,psi(:,0),dx,psiSquareSum)
    psi(:,0) = psi(:,0) / sqrt(psiSquareSum)

    ! define M_1, M_2
    alpha = cmplx(0.d0,1.d0) * dt / 2.d0 / dx**2
    allocate(M1(0:Jmax,0:Jmax))
    M1 = 0.d0
    allocate(M2(0:Jmax,0:Jmax))
    M1 = 0.d0
    do j = 0,Jmax
        M1(j,j) = 1.d0 + cmplx(0.d0,1.d0) * dt / 2.d0 * (2.d0 / dx**2)    ! + V(j) in last ()
        M2(j,j) = 1.d0 - cmplx(0.d0,1.d0) * dt / 2.d0 * (2.d0 / dx**2)    ! + V(j) in last ()
        if (j < Jmax) then
            M1(j,j + 1) = -alpha
            M2(j,j + 1) = alpha
        end if
        if (j > 1) then
            M1(j,j - 1) = -alpha
            M2(j,j - 1) = alpha
        end if
    end do

    ! calculate M_1^{-1}M_2
    allocate(M(0:Jmax,0:Jmax))
    ! with LU factorization
    ! allocate(IPIV(0:Jmax))
    ! call ZGETRF(Jmax + 1,Jmax + 1,M1,Jmax + 1,IPIV,INFO)
    ! if (INFO /= 0) then
    !     write(*,*) 'ZGETRE went wrong! (INFO = ', INFO, ') Aborting ...'
    !     stop
    ! end if
    ! LWORK = 10 * (Jmax + 1)
    ! allocate(WORK(LWORK))
    ! call ZGETRI(Jmax + 1,M1,Jmax + 1,IPIV,WORK,LWORK,INFO)
    ! if (INFO /= 0) then
    !     write(*,*) 'ZGETRI went wrong! (INFO = ', INFO, ') Aborting ...'
    !     stop
    ! end if
    ! M = matmul(M1,M2)

    ! with chasing method
    allocate(a(1:Jmax))
    allocate(b(0:Jmax))
    allocate(c(0:Jmax - 1))
    allocate(d(0:Jmax))
    allocate(w(1:Jmax))
    do j = 0,Jmax
        b(j) = M1(j,j)
        if (j > 0) a(j) = M1(j,j - 1)
        if (j < Jmax) c(j) = M1(j,j + 1)
    end do
    do j = 1,Jmax
        w(j) = a(j) / b(j - 1)
        b(j) = b(j) - w(j) * c(j - 1)
    end do
    do i = 0,Jmax
        d = 0.d0
        d(i) = 1.d0
        do j = 1,Jmax
            d(j) = d(j) - w(j) * d(j - 1)
        end do
        M(Jmax,i) = d(Jmax) / b(Jmax)
        do j = Jmax - 1,0,-1
            M(j,i) = (d(j) - c(j) * M(j + 1,i)) / b(j)
        end do
    end do
    M = matmul(M,M2)

    ! time iteration
    do n = 1,Nmax
        psi(:,n) = matmul(M,psi(:,n - 1))
        call Simpson(Jmax + 1,psi(:,n),dx,psiSquareSum)
        psi(:,n) = psi(:,n) / sqrt(psiSquareSum)
    end do

    ! output to file
    open(unit = 1, file = '1-pdf.txt', status = 'unknown')
    do j = 0,Jmax
        write(1,'(f24.16)', advance = 'no') x(j)
        do n = 0,Nmax
            write(1,'(f24.16)', advance = 'no') (real(psi(j,n))**2 + aimag(psi(j,n))**2)
        end do
        write(1,*)
    end do
    close(1)

    ! output to screen (animate)
    CALL system('gnuplot "1-plt.gnu"')
end program main

subroutine Simpson(N,f,dx,fSquareSum)
    ! calculate the integral of f^2 using Simpson method
    
    implicit none
    integer :: N
    complex(8), intent(in) :: f(N)
    real(8) :: dx
    real(8), intent(out) :: fSquareSum

    ! local vars
    integer :: i
    complex(8) :: f2(N)

    if (mod(N,2) == 0) then
        ! Simpson does not work for integral of array with even number of elements
        write(*,*) 'Array with even elements, Simpson does not know how to work!'
    end if

    f2 = real(f)**2 + aimag(f)**2
    fSquareSum = f(1) + f(N)
    fSquareSum = fSquareSum + sum(4.d0 * f2(2:N - 1:2))
    fSquareSum = fSquareSum + sum(2.d0 * f2(3:N - 2:2))
    fSquareSum = fSquareSum * dx / 3.d0
end subroutine Simpson
