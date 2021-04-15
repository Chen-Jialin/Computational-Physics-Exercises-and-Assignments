program main
    implicit none
    real(8), parameter :: pi = acos(-1.d0), kB = 1.d0
    integer :: clock, n, i
    integer, allocatable :: seed(:)

    real(8) :: T = 1.d0, dT = .1d0, T_final = 5.d0, beta
    integer, parameter :: L_x = 32, L_y = 32
    integer :: lattice(0:L_x - 1, 0:L_y - 1)
    integer, parameter :: n_warmup = 100000, n_evol = 100000
    real(8) :: J = 1.d0, B = 0.d0
    real(8) :: M

    integer :: x, y

    ! seed initialization
    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(put = seed)

    open(unit = 1, file = 'M-T.txt', status = 'unknown')

    do while (T < T_final)

    ! spin initialization
    beta = 1.d0 / kB / T
    lattice = 1
    ! warm up
    do i = 1, n_warmup
        call EVOLUTION(lattice, L_x, L_y, beta, B, J)
    end do
    open(unit = 2, file = 'lattice.txt', status = 'unknown')
    do x = 0, L_x - 1
        do y = 0, L_y - 1
            write(2,'(3i4)') x, y, lattice(x,y)
        end do
    end do
    close(2)

    ! do while (T < T_final)
        beta = 1.d0 / kB / T
        M = 0.d0

        ! evolution
        do i = 1, n_evol
            call EVOLUTION(lattice, L_x, L_y, beta, B, J)
            M = M + dble(sum(lattice))
        end do
        M = M / dble(n_evol)
        write(*,'(2f15.5)') T, abs(M)
        write(1,'(2f15.5)') T, abs(M)
        T = T + dT
    end do
    close(1)
end program main

subroutine EVOLUTION(lattice, L_x, L_y, beta, B, J)
    implicit none
    integer, intent(in) :: L_x, L_y
    integer, intent(inout) :: lattice(0:L_x - 1, 0:L_y - 1)
    real(8), intent(in) :: beta, B, J
    real(8) :: r
    integer :: x, y
    real(8) :: dE

    call RANDOM_NUMBER(r)
    x = floor(r * dble(L_x))
    call RANDOM_NUMBER(r)
    y = floor(r * dble(L_y))

    dE = 0.d0
    dE = dE + 2.d0 * J * lattice(x,y) * (lattice(modulo(x - 1, L_x), y)&
        + lattice(modulo(x + 1, L_x), y)&
        + lattice(x, modulo(y - 1, L_y))&
        + lattice(x, modulo(y + 1, L_y)))&
        + 2.d0 * B * dble(lattice(x,y))

    call RANDOM_NUMBER(r)
    if (r < exp(-beta * dE)) then
        lattice(x,y) = -lattice(x,y)
    end if
end subroutine EVOLUTION
