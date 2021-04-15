program main
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    integer :: n_in = 0, n_tot = 100000000
    real(8), parameter :: l = 0.d0, r = 1.d0, b = -1.d0, t = 2.d0
    real(8) :: x, y, y_real
    real(8) :: integral

    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(PUT = seed)
    deallocate(seed)

    do i = 1, n_tot
        call RANDOM_NUMBER(x)
        x = x * (r - l) + l
        call RANDOM_NUMBER(y)
        y = y * (t - b) + b
        call func(x, y_real)
        if (y < y_real) then
            n_in = n_in + 1
        end if
    end do

    integral = (t - b) * (r - l) * dble(n_in) / dble(n_tot) + b * (r - l)
    write(*,'(f10.5)') integral
end program main

subroutine func(x, y)
    ! the function to be integrated
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: y

    y = x**2
end subroutine func
