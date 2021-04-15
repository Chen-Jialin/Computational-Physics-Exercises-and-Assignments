program main
    ! compare the calculate and exact solution of the wavefunction of a particle in in 1-dimension infinite square potential well

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    
    ! local vars
    integer :: i
    real(dp) :: x(201), yCalculated(201), yExact(201), err = .0d0

    ! get the data of calculated and exact solution
    open(unit = 1, file = '1-wavefunction.txt', status = 'unknown')
    do i = 1,201
        read(1, '(2f20.8)') x(i), yCalculated(i)
    end do
    close(1)
    open(unit = 2, file = '1-wavefunction-exact.txt', status = 'unknown')
    do i = 1,201
        read(2, '(2f20.8)') x(i), yExact(i)
    end do
    close(2)

    do i = 1,201
        err = err + abs(yCalculated(i) - yExact(i))
    end do
    err = err / 201.d0
    write(*, *) err
end program main
