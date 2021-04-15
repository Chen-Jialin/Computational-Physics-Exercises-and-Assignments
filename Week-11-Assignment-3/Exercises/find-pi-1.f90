program main
    implicit none
    integer :: clock, i, n
    integer, allocatable :: seed(:)
    integer(4) :: n_in, n_tot = 100000000
    real(8) :: x, y
    real(8) :: pi

    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(put = seed)
    do i = 1, n_tot
        call RANDOM_NUMBER(x)
        x = 2 * x - 1
        call RANDOM_NUMBER(y)
        y = 2 * y - 1
        if (x**2 + y**2 <= 1.d0) then
            n_in = n_in + 1
        end if
    end do
    pi = 4 * dble(n_in) / dble(n_tot)
    write(*,'(f10.8)') pi
end program main
