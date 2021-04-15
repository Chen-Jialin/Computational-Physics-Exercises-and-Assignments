program double_quantum_well

   implicit none

   integer,  parameter :: dp = selected_real_kind(8)
   real(dp), parameter :: Zero = 0.d0, Half = 0.5d0, One = 1.d0, Two = 2.d0
   real(dp), parameter :: pi = acos(-One)

   integer,  parameter :: N = 200         ! discretization of x
   integer,  parameter :: N_basis = 500   ! number of basis functions
   real(dp), parameter :: Range_of_X = 6.d0, x0 = 1.5d0
   real(dp) :: V_pot(-N:N), X(-N:N), dx, E, Y(-N:N)

   integer :: i, Idx

   ! ... executable ...
   dx = Range_of_X/N

   ! generate the potential data
   do i = -N, N
      X(i) = dx*i
      V_pot(i) = Half*( (x(i)**2 - x0**2)**2/Two/Two/x0**2 - X(i)**2 )
   end do  

   open(unit=1, file='Vpot.dat', status = 'unknown')
   do i = -N, N
      write(1, '(2f12.6)') X(i), V_pot(i)
   end do
   close(1)

   ! calculate with numeral approach
   E = -1.5d0
   call Numerov(N, X, V_pot, E, Y)
   Y = Zero
   ! 

   do Idx = 1, 2
      call Matrix_diagonalization(Idx, N, N_basis, X, V_pot)
   end do   

contains
   !----------------------------------------------------------------------
   subroutine Numerov(Nx, X, Vpot, E, Y)
      !
      ! Purpose
      ! =======
      !   Solving one-dimensional Schrodinger equation with Numerov method.
      !
      integer,  intent(in)  :: Nx
      real(dp), intent(in)  :: X(-Nx:Nx), Vpot(-Nx:Nx)
      real(dp), intent(out) :: Y(-Nx:Nx)
      real(dp), intent(inout) :: E

      ! ... local vars ...
      integer  :: ite, i
      real(dp) :: f(-Nx:Nx), g(-Nx:Nx), dx, fac, Y_old, Emin, Emax, dE, temp(-Nx:Nx)
      logical  :: bisection = .false.

      ! ... executable ...
      write(*, *) 'Numerov calculation starts ...'
      dx = X(2) - X(1)
      Y = Zero

      ! set boundary condition at Y(-Nx) and Y(-Nx+1) 

      Emin = -200.d0
      Emax =  200.d0

   iteration: do ite = 1, 10000000
 
         Y = Zero
         Y(-Nx) = 0.d0
         Y(-Nx+1) = 0.00001d0
         dE   = 0.0001d0

         do i = -Nx, Nx
            g(i) = 2.d0*(E - Vpot(i))
            f(i) = 1.d0 + g(i)/12.d0*dx*dx
         end do

         ! Nomerov's method
         do i = -Nx+1, Nx-1
            Y(i+1) = ((12.d0 - 10.d0*f(i))*Y(i) - f(i-1)*Y(i-1))/f(i+1)
         end do

         ! normalize the wave function Y(x)
         do i = -Nx, Nx
            Temp(i) = Y(i)**2
         end do
         call Simpson(Temp(-Nx:Nx), 2*Nx+1, dx, fac)
         Y = Y/sqrt(fac)


         if (ite == 1) Y_old = Y(Nx)
      
         ! monitor the boundary value of Y(Nx)
         if (ite == 1) write(*, *) 'Iteration, energy and boundary value'
         write(*, '(i10, f20.12, 4x, f20.12)') ite, E, Y(Nx)
 
         if (.NOT. bisection) then
            if (Y_old*Y(Nx) >= 0.d0) then
               Emin = E
               E = E + dE
               Y_old = Y(Nx)
            else
               Emax = E
               Emin = E-de
               bisection = .true.
               E = 0.5d0*(Emin+Emax)
            end if
         else
            if (Y_old*Y(Nx) >= 0 ) then
               Emin = E
               Y_old = Y(Nx)
            else
               Emax = E
            end if
            E = 0.5d0*(Emax + Emin)
            if (abs(Y(Nx)) < 1.d-8) exit
         end if

      end do iteration

      open(unit = 1, file='wave_Numerov.dat', status='unknown')
      do i = -Nx, Nx
         write(1, '(2f12.6)') X(i), Y(i)
      end do
      close(1)

   end subroutine Numerov

   !----------------------------------------------------------------------
   subroutine Matrix_diagonalization(Idx, Nx, N_basis, X, Vpot)
      !
      ! Purpose
      ! =======
      !   By using basis funtion to express the Hamiltonian in matrix form and 
      ! calculate the ground state energy by diagonalizing it.
      !
      integer,  intent(in)  :: Idx, Nx, N_basis
      real(dp), intent(in)  :: X(-Nx:Nx), Vpot(-Nx:Nx)

      ! local vars 
      integer   :: i, j, k, INFO, LWORK
      real(dp)  :: basis(2*Nx+1, N_basis), Ham(N_basis, N_basis), WORK(10*N_basis), W(N_basis) 
      real(dp)  :: Xmax, V0, Temp(2*Nx+1)
      character(len=30) :: FLE  


      Xmax = X(Nx)
      dx   = X(2) - X(1)

      ! step[1]: generating the basis functions
      basis = 0.d0
      call basis_generation(Idx, Nx, N_basis, basis)

      ! step[2]: generating the Hamiltonian matrix
      Ham = 0.d0
      do i = 1, N_basis
         do j = i, N_basis

            select case (Idx)
            case(1)
               ! eigenvectors of inifite quantum well as basis
               do k = -Nx, Nx
                  Temp(k+Nx+1) = Vpot(k)*basis(k+Nx+1, i)*basis(k+Nx+1, j)  
               end do 
               call Simpson(Temp, 2*Nx+1, dx, V0)
               Ham(i, j) = V0

               if (i == j) then
                  Ham(i, i) = Ham(i, i) + (dble(i)*pi)**2/8.d0/Xmax/Xmax
               end if
            case(2)
               ! eigenvectors of quantum harmonic oscillator as basis
               do k = -Nx, Nx
                  Temp(k+Nx+1) = (Vpot(k)-Half*X(k)**2)*basis(k+Nx+1, i)*basis(k+Nx+1, j)  
               end do 
               call Simpson(Temp, 2*Nx+1, dx, V0)
               Ham(i, j) = V0

               if (j == i)   Ham(i, i) = Ham(i, i) + 0.5d0 + dble(j-1)    
            end select
            Ham(j, i) = Ham(i, j)
         end do
      end do

      LWORK = 10*N_basis
      call DSYEV( 'V', 'U', N_basis, Ham, N_basis, W, WORK, LWORK, INFO)
    
      if (INFO == 0) then
         write(FLE, '(a,i1,a)') 'eigenvalue-', Idx, '.dat'
         open(unit=1, file=FLE, status='unknown')
    
         ! multiply the coefficient (which is stored in Ham now) with the basis function
         ! to get the wave functions. 
         basis = matmul(basis, Ham)  ! use basis to store Psi(x) = Sum_n c_n * Phi_n(x)
         
         ! renormalize the wave function
         do i = 1, N_basis
            call Simpson(basis(1:2*Nx+1, i)**2, 2*Nx+1, dx, V0)
            basis(1:2*Nx+1, i) = basis(1:2*Nx+1, i)/sqrt(V0)
         end do

         ! output results
         do i = 1, N_basis
            write(1, '(i4, f20.12)') i, W(i)
         end do
         close(1)

         write(FLE, '(a,i1,a)') 'wavefunction-', Idx, '.dat'
         open(unit=2, file=FLE, status='unknown')
         do i = 1, 2*Nx+1
            write(2, '(f12.6)', advance='no') X(i-Nx-1)
            do j = 1, N_basis
               write(2, '(f12.6)', advance='no') basis(i, j)
            end do
            write(2, *)
         end do
         close(2)
      else
         write(*, *) 'Diagonalization went wrong! aborting ...'
         stop   
      end if

   end subroutine Matrix_diagonalization

   !--------------------------------------------------------------------
   subroutine basis_generation(Idx, Nx, N_basis, basis)

      integer,  intent(in)  :: Idx
      integer,  intent(in)  :: Nx, N_basis
      real(dp), intent(out) :: basis(2*Nx+1, N_basis) 
      !
      ! Purposes
      ! =========
      !   This subroutine generates the basis functions for the Hamiltonian matrix.
      ! Here, two different eigenvectors can be selected by sepecifying the switcher
      ! Idx = 1  --- eigenvectors of the infinite quantum Well
      ! Idx = 2  --- eigenvectors of the quantum harmonic oscillator
      !

      ! ... local vars ...
      integer  :: i, j 

      open(unit = 1, file = 'basis_function.dat', status='unknown')
 
      select case(Idx)
      case(1)
         write(*, *) 'case 1: eigenvector of infinite quantum well as basis functions'
         do i = 1, 2*Nx+1   ! coordinate index
            do j = 1, N_basis  ! order of basis function
               if (mod(j, 2) /= 0) then
                  ! j is odd
                  basis(i, j) = sqrt(1.d0/X(Nx))*cos(dble(j)*pi*X(i-Nx-1)/2.d0/X(Nx))
               else
                  ! j is even
                  basis(i, j) = sqrt(1.d0/X(Nx))*sin(dble(j)*pi*X(i-Nx-1)/2.d0/X(Nx))
               end if
            end do
         end do
      case(2)
         write(*, *) 'case 2: eigenvectors of quantum harmonic oscillator as basis functions'
         basis = 0.d0
         do i = 1, 2*Nx + 1   ! order of basis function
            basis(i, 1) = (1.d0/pi)**0.25*exp(-X(i-Nx-1)**2/2.d0)
            basis(i, 2) = (1.d0/pi)**0.25*exp(-X(i-Nx-1)**2/2.d0)*X(i-Nx-1)*sqrt(Two)
         end do         

         do i = 1, 2*Nx + 1
            write(1, '(3f12.6)', advance='no') X(i-Nx-1), basis(i, 1), basis(i, 2)
            do j = 3, N_basis ! coordinate index
               basis(i, j) = (Two*X(i-Nx-1)*basis(i, j-1) - sqrt(Two*(j-2))*basis(i,j-2))/sqrt(Two*(j-1)) 
               write(1, '(f12.6)', advance='no') basis(i, j)
            end do
            write(1, *)
         end do   
      end select 

      close(1)

   end subroutine basis_generation

   !--------------------------------------------------------------------
   subroutine Simpson(Func, nL, h, solution)

      implicit none
      integer, intent(in)   :: nL
      real(dp), intent(in)  :: Func(nL), h
      real(dp), intent(out) :: solution

      integer :: i

     if (Mod(nL, 2) == 0) then
         print*, "nL must be odd!"
         stop
      endif

      solution = (Func(1) + Func(nL))
      do i = 2, nL-1, 2
         solution = solution + 4.d0 * Func(i)
      enddo
      do i = 3, nL-2, 2
         solution = solution + 2.d0 * Func(i)
      enddo
      solution = solution * h/3.d0

      end subroutine Simpson

end program double_quantum_well