program main
    use mpi
    implicit none
    real(8), parameter :: pi = acos(-1.d0), kB = 1.d0
    integer :: ntasks, id, rc
    integer, allocatable :: status(:)
    integer :: i, n, clock
    integer, allocatable :: seed(:)
    real(8) :: r

    real(8) :: T = 1.d0, dT = .01d0, T_final = 5.d0, beta
    integer, parameter :: L_x = 32, L_y = 32
    integer :: lattice(0:L_x - 1,0:L_y - 1) = 1
    integer :: x, y
    integer, parameter :: n_warmup = 10000, n_evol = 100000
    real(8) :: J = 1.d0, B = 0.d0
    real(8) :: M, M_sqr, E_tmp, E, E_sqr, M_ave, M_sqr_ave, chi, E_ave, E_sqr_ave, C

    ! initialize MPI environment
    call MPI_INIT(rc)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)
    allocate(status(MPI_STATUS_SIZE))

    ! initialize seeds for different processes
    if (id == 0) then
        call SYSTEM_CLOCK(clock)
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = clock + 37 * i
        end do
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
        do i = 1, ntasks - 1
            call RANDOM_NUMBER(r)
            clock = clock + Int(r * 1000000)
            call MPI_SEND(clock, 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, rc)
        end do
    else
        call MPI_RECV(clock, 1, MPI_INTEGER, 0, id, MPI_COMM_WORLD, status, rc)
        call RANDOM_SEED(size = n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = clock + 37 * i
        end do
        call RANDOM_SEED(PUT = seed)
        deallocate(seed)
    end if

    if (id == 0) then
        open(unit = 1, file = 'data.txt', status = 'unknown')
        write(*,'(4a15)') 'T', 'M / M_max', 'chi', 'C'
    end if

    do while (T < T_final)
    ! warm up
    beta = 1.d0 / kB / T
    do i = 1, n_warmup
        call EVOLUTION(lattice, L_x, L_y, beta, B, J)
    end do

    if (id == 0) then
        open(unit = 2, file = 'lattice-mpi.txt', status = 'unknown')
        do x = 0, L_x - 1
            do y = 0, L_y - 1
                write(2,'(3i4)') x, y, lattice(x,y)
            end do
        end do
        close(2)
    end if

    ! do while (T < T_final)
        ! beta = 1.d0 / kB / T
        M = 0.d0
        M_sqr = 0.d0
        E = 0.d0
        E_sqr = 0.d0

        ! evolution
        do i = 1, n_evol
            call EVOLUTION(lattice, L_x, L_y, beta, B, J)
            M = M + abs(sum(lattice))
            M_sqr = M_sqr + sum(lattice)**2
            E_tmp = 0.d0
            do x = 0, L_x - 2
                do y = 0, L_y - 2
                    E_tmp = E_tmp + lattice(x, y) * (lattice(x + 1, y) + lattice(x, y + 1))
                end do
            end do
            do x = 0, L_x - 2
                E_tmp = E_tmp + lattice(x, L_y - 1) * (lattice(x + 1, L_y - 1) + lattice(x, 0))
            end do
            do y = 0, L_y - 2
                E_tmp = E_tmp + lattice(L_x - 1, y) * (lattice(L_x - 1, y + 1) + lattice(0, y))
            end do
            E_tmp = E_tmp + lattice(L_x - 1, L_y - 1) * (lattice(0, L_y - 1) + Lattice(L_x - 1, 0))
            E_tmp = - E_tmp * J - B * sum(lattice)
            E = E + E_tmp
            E_sqr = E_sqr + E_tmp**2
        end do
        call MPI_REDUCE(M, M_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, rc)
        call MPI_REDUCE(M_sqr, M_sqr_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, rc)
        call MPI_REDUCE(E, E_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, rc)
        call MPI_REDUCE(E_sqr, E_sqr_ave, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, rc)
        if (id == 0) then
            M_ave = M_ave / dble(ntasks * n_evol)
            M_sqr_ave = M_sqr_ave / dble(ntasks * n_evol)
            chi = beta * (M_sqr_ave - M_ave**2)
            E_ave = E_ave / dble(ntasks * n_evol)
            E_sqr_ave = E_sqr_ave / dble(ntasks * n_evol)
            C = kB * beta**2 / dble((L_x) * (L_y)) * (E_sqr_ave - E_ave**2)
            write(*,'(4f15.5)') T, M_ave / dble(L_x * L_y), chi, C
            write(1,'(4f15.5)') T, M_ave / dble(L_x * L_y), chi, C
        end if
        T = T + dT
    end do

    if (id == 0) then
        close(1)
    end if

    ! done with MPI
    call MPI_FINALIZE(rc)
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
    if (r <= exp(-beta * dE)) then
        lattice(x,y) = -lattice(x,y)
    end if
end subroutine EVOLUTION
