program main
    implicit none
    real(8), parameter :: pi = acos(-1.d0)
    integer :: clock, n, i
    integer, allocatable :: seed(:)

    integer :: n_sample = 10000000, n_accept = 0
    real(8), parameter :: Xmax = 5.d0
    real(8) :: x = 0.d0, r1, r2, A, Integral = 0.d0

    integer, allocatable :: bin(:)
    integer :: Nbin = 1000
    real(8) :: binwidth

    real(8) :: mu, theta

    !real(8) ::

    ! seed initialization
    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(put = seed)

    ! bin initialization
    binwidth = 2.d0 * Xmax / dble(Nbin)
    allocate(bin(Nbin))
    bin = 0

    ! Metropolis
    do i = 1, n_sample
        call RANDOM_NUMBER(r1)
        r1 = (r1 - .5d0) * (2.d0 * Xmax)
        A = exp(-(r1**2 - x**2) / 2.d0)
        call RANDOM_NUMBER(r2)
        if (r2 <= A) then
            x = r1
            bin(floor((x + Xmax) / binwidth)) = bin(floor((x + Xmax) / binwidth)) + 1
            n_accept = n_accept + 1
        end if
        integral = integral + x**2
    end do
    integral = sqrt(2.d0 * pi) * integral / dble(n_sample)
    write(*,'("Metropolis : ",f10.5," Acceptance rate : ",f10.5)') integral, dble(n_accept) / dble(n_sample)
    open(unit = 1, file = 'samples.txt', status = 'unknown')
    do i = 1, Nbin
        write(1,'(f10.5,i8)') -Xmax + (i - .5d0) * binwidth, bin(i)
    end do
    close(1)

    ! Importance sampling - inverse transform
    integral = 0.d0
    do i = 1, n_sample
        call RANDOM_NUMBER(mu)
        call RANDOM_NUMBER(theta)
        theta = 2.d0 * pi * theta
        x = sqrt(-2.d0 * log(mu)) * cos(theta)
        integral = integral + x**2
    end do
    integral = sqrt(2 * pi) * integral / n_sample
    write(*,'("Importance sampling (inverse transform) : ",f10.5)') integral

    ! Importance sampling - acceptance & rejection
    integral = 0.d0
    n_accept = 0
    do i = 1,n_sample
        call RANDOM_NUMBER(x)
        x = (x - .5d0) * (2 * Xmax)
        call WEIGHT(x,r1)
        call random_NUMBER(r2)
        if (r2 <= r1) then
            integral = integral + x**2
            n_accept = n_accept + 1
        end if
    end do
    integral = sqrt(2.d0 * pi) * integral / dble(n_accept)
    write(*,'("Important sampling (acceptance & rejection) : ",f10.5," Acceptance rate : ",f10.5)')&
        integral, dble(n_accept) / dble(n_sample)
end program main

subroutine WEIGHT(x,p)
    implicit none
    real(8), parameter :: pi = acos(-1.d0)
    real(8), intent(in) :: x
    real(8), intent(out) :: p

    p = exp(-x**2 / 2.d0) / sqrt(2.d0 * pi)
end subroutine WEIGHT
