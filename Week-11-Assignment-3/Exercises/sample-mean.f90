program main
    implicit none
    integer :: clock, n
    integer, allocatable :: seed(:)
    integer(4) :: i, n_tot = 100000000
    real(8), parameter :: l = 0.d0, r = 1.d0
    real(8) :: x, y
    real(8) :: integral = 0.d0

    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(put = seed)
    do i = 1, n_tot
        call RANDOM_NUMBER(x)
        x = x * (r - l) + l
        call func(x, y)
        integral = integral + y
    end do
    integral = (r - l) / dble(n_tot) * integral
    write(*,'(f10.5)') integral
end program main

subroutine func(x, y)
    ! the function to be integrated
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: y

    y = x**2
end subroutine func
