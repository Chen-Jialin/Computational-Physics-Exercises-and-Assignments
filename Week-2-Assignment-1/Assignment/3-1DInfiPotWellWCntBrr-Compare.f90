program main
    ! compare the energy calculated from matrix diagonalization method and Numerov's method

    implicit none
    integer, parameter :: dp = selected_real_kind(8)

    ! local vars
    integer :: i, n(100)
    real(dp) :: EMD(100), EN(100)

    open(unit = 1, file = '3-energy-ExactDiagonalization.txt', status = 'unknown')
    do i = 1, 100
        read(1, '(i4,f20.8)') n(i), EMD(i)
    end do
    close(1)

    open(unit = 2, file = '3-energy-Numerov.txt', status = 'unknown')
    do i = 1, 100
        read(2, '(i4,f20.8)') n(i), EN(i)
    end do
    close(2)

    open(unit = 3, file = '3-energy-Compare.txt', status = 'unknown')
    do i = 1, 10
        write(3, '(i4,4f20.8)') n(i), EMD(i), EN(i), EMD(i) - EN(i), (EMD(i) - EN(i)) / EN(i)
    end do
    write(3, *)
    do i = 91, 100
        write(3, '(i4,4f20.8)') n(i), EMD(i), EN(i), EMD(i) - EN(i), (EMD(i) - EN(i)) / EN(i)
    end do
    close(3)
end program main
