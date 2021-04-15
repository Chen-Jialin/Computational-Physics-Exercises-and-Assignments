program main
    ! 2D model of crystal growth based on diffusion-limited aggregation
    implicit none
    integer :: clock, n, i
    integer, allocatable :: seed(:)
    integer, parameter :: map_size_x = 1024, map_size_y = 1024, n_particle = 1000
    integer :: map(0:map_size_x - 1, 0:map_size_y - 1)
    integer :: x, y

    ! seed initialization
    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    do i = 1, n
        seed(i) = clock + 37 * i
    end do
    call RANDOM_SEED(put = seed)

    ! map initialization
    do x = 0, map_size_x - 1
        do y = 0, map_size_y - 1
            map(x, y) = 0
        end do
    end do
    map(map_size_x / 2, map_size_y / 2) = 1
    open(unit = 1, file = 'crystal.txt', status = 'unknown')
    write(1, '(2i5)') map_size_x / 2, map_size_y / 2
    close(1)

    ! particle aggregation
    do i = 1, n_particle
        do while (.true.)
            call INITIAL_POSITION(x, y, map_size_x, map_size_y)
            if (map(x, y) == 0) exit
        end do
        do while (.true.)
            call RANDOM_WALK(x, y, map_size_x, map_size_y)
            if ((map(modulo(x - 1, map_size_x), y) == 1) &
                .or. (map(modulo(x + 1, map_size_x), y) == 1) &
                .or. (map(x, modulo(y - 1, map_size_y)) == 1) &
                .or. (map(x, modulo(y + 1, map_size_y))) == 1) then
                map(x, y) = 1
                exit
            end if
        end do
        open(unit = 1, file = 'crystal.txt', status = 'old', position = 'append')
        write(1, '(2i5)') x, y
        close(1)
    end do
    close(1)
end program main

subroutine INITIAL_POSITION(x, y, map_size_x, map_size_y)
    real :: r
    integer, intent(inout) :: x, y
    integer, intent(in) :: map_size_x, map_size_y

    call RANDOM_NUMBER(r)
    if (r < .25d0) then
        x = 0
        y = floor(r * 4.d0 * dble(map_size_y))
    elseif (r < .5d0) then
        x = map_size_x - 1
        y = floor((r - .25d0) * 4.d0 * dble(map_size_y))
    elseif (r < .75d0) then
        x = floor((r - .5d0) * 4.d0 * dble(map_size_x))
        y = 0
    else
        x = floor((r - .75d0) * 4.d0 * dble(map_size_x))
        y = map_size_y - 1
    end if
end subroutine INITIAL_POSITION

subroutine RANDOM_WALK(x, y, map_size_x, map_size_y)
    real :: r
    integer, intent(inout) :: x, y
    integer, intent(in) :: map_size_x, map_size_y

    call RANDOM_NUMBER(r)
    if (r < .25d0) then
        x = modulo(x - 1, map_size_x)
    elseif (r < .5d0) then
        x = modulo(x + 1, map_size_x)
    elseif (r < .75d0) then
        y = modulo(y - 1, map_size_y)
    else
        y = modulo(y + 1, map_size_y)
    end if
end subroutine RANDOM_WALK
