program main
    use MPI
    implicit none

    ! variables for parallel environment
    integer :: rank, ntasks, ierr
    integer :: master = 0
    integer, allocatable :: status(:)

    ! local vars
    integer, parameter :: row = 4, col = 4, np = 4
    integer(2) :: i, r, c
    real(8) :: x(1, row), A(row, col), A_local(row, 1)
    real(8) :: prod_local(1,1), prod_global(1, col)

    ! prepare MPI environment
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
        write(*,*) 'MPI initialization failed!'
        stop
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (ntasks /= np) then
        write(*,*) 'Must specify', np, 'processors. Terminating.'
        stop
    end if

    ! set & distribute data
    if (rank == master) then
        ! master node
        ! set data
        x = reshape((/ 1.d0, 2.d0, 3.d0, 4.d0 /), (/ 1, row /))
        A = reshape((/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0, 11.d0, 12.d0, 13.d0, 14.d0, 15.d0, 16.d0 /),&
            (/ row, col /))

        ! write(*,*) 'Matrix A:'
        ! do r = 1,row
        !     do c = 1,col
        !         write(*,'("A(",I1,",",I1,")=")') r,c
        !         read(*,*) A(r,c)
        !     end do
        ! end do

        ! write(*,*) 'Vector x:'
        ! do r = 1,row
        !     write(*,'("x(",I1,",1)=")') r
        !     read(*,*) x(r,1)
        ! end do

        ! output data
        ! write(*,*) 'Vector x'
        ! do c = 1, row
        !     write(*,'(f8.2)') x(1,c)
        ! end do

        ! write(*,*) 'Matrix A:'
        ! do c = 1, col
        !     do r = 1, row
        !         write(*,'(f8.2)',advance = 'no') A(r,c)
        !     end do
        !     write(*,*)
        ! end do

        ! write(*,*) 'The result by MPI is'
        ! do c = 1, row
        !     write(*,'(f8.2)', advance = 'no') prod_global(1,c)
        ! end do
    end if

    call MPI_BCAST(x, np, MPI_REAL8, master, MPI_COMM_WORLD, ierr)
    call MPI_SCATTER(A, np, MPI_REAL8, A_local, np, MPI_REAL8, master, MPI_COMM_WORLD, ierr)

    ! calculate local result
    prod_local = matmul(x, A_local)

    ! gather result
    call MPI_GATHER(prod_local, 1, MPI_REAL8, prod_global, 1, MPI_REAL8, master, MPI_COMM_WORLD, ierr)

    ! output result
    if (rank == master) then
        write(*,*) 'The result by MPI is'
        do c = 1, col
            write(*, '(f8.2)', advance = 'no') prod_global(1, c)
        end do
        write(*,*)
        prod_global = matmul(x, A)
        write(*,*) 'The correct result should be'
        do c = 1,col
            write(*, '(f8.2)', advance = 'no') prod_global(1, c)
        end do
        write(*,*)
    end if

    ! clean up MPI
    ! deallocate(status)
    call MPI_FINALIZE(ierr)
end program main

