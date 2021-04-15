program main
    ! calculate the exact solution of the Schrodinger equation of a particle in 2-dimension infinite square potential well

    implicit none
    integer, parameter :: dp = selected_real_kind(8)
    real(dp), parameter :: pi = acos(-1.d0)

    ! local vars
    real(dp) :: L = 1.d0
    integer :: Nx
    integer :: n
    real(dp) :: dx, x, y

    ! excutable
    Nx = -1
    do while (Nx <= 0)
        write(*,*) 'Please specify the # number of slices between [0,L], Nx ='
        read(*,*) Nx
    end do

    write(*,*) 'The # of the energy level, n ='
    read(*,*) n

    x = 0.d0
    dx = L / Nx
    open(unit = 1,file = '1-wavefunction-exact.txt',status = 'unknown')
    do while (x < L + dx)
        y = sqrt(2.d0 / L) * sin(n * pi * x / L)
        write(1, '(2f20.8)') x, y 
        x = x + dx
    end do
    close(1)

end program main
