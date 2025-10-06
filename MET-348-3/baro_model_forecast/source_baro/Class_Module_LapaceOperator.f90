!  $Author: pkubota 						$
!  $Date: 2008/09/23 17:51:54 					$
!  $Revision: 1.9 						$
!  $Revisions are currently made by the class's students.	$
!  $Update Date: 01/11/2023 10:19 AM				$
!  
!  Implementações: 
!  	1) Colocar INIT e FINALIZE - OK 
! 	2) Colocar as equações dos campos no "Class_Module_Dynamics" - OK
! 	7) Amortecimento das condições de contorno - OK (funcionando parcialmente)
MODULE Class_Module_LapaceOperator
  PRIVATE
  ! Selecting Kinds
  INTEGER,PUBLIC, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER,PUBLIC, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER,PUBLIC, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  type :: mg_level
    integer  :: nx, ny
    real(r8), allocatable :: dx(:), dy(:)
    real(r8), allocatable :: phi(:,:), rhs(:,:), res(:,:), corr(:,:)
  contains
    procedure :: allocate_level   => mg_allocate_level
    procedure :: deallocate_level => mg_deallocate_level
    procedure :: set_grid         => mg_set_grid
    procedure :: zero_fields      => mg_zero_fields
    procedure :: compute_residual_norm => mg_compute_residual_norm
    procedure :: get_phi          => mg_get_phi
  end type mg_level


  PUBLIC :: init_rhs_and_boundary
  PUBLIC :: sor_solve
  PUBLIC :: multigrid_solver
  PUBLIC :: report_error

CONTAINS
    subroutine init_rhs_and_boundary(phi, f, nx, ny, xmin, ymin, dx, dy)
    implicit none
    integer, intent(in) :: nx, ny
    real(r8), intent(in) :: xmin, ymin, dx(nx), dy(ny)
    real(r8), intent(out) :: phi(nx,ny), f(nx,ny)
    integer :: i,j
    real(r8) :: x,y,pi
    pi=function_pi()
    do j = 1, ny
      y = ymin + (j-1) * dy(j)
      do i = 1, nx
        x = xmin + (i-1) * dx(i)
        ! Manufactured solution: phi_exact = sin(pi x) sin(pi y)
        phi(i,j) = sin(pi*x) * sin(pi*y)
        f(i,j) = -2.0_r8 * pi*pi * sin(pi*x) * sin(pi*y)
      end do
    end do

    ! Apply Dirichlet boundary conditions: keep phi on boundary as exact
    ! For solver initial guess we want interior initial guess maybe zero

    do j = 1, ny
      phi(1,j)  = sin(pi*(xmin +  0    *dy(j))) * sin(pi*(ymin + (j-1)*dy(j)))
      phi(nx,j) = sin(pi*(xmin + (nx-1)*dy(j))) * sin(pi*(ymin + (j-1)*dy(j)))
    end do

    do i = 1, nx
      phi(i, 1) = sin(pi*(xmin + (i-1)*dx(i))) * sin(pi*(ymin +   0   *dx(i)))
      phi(i,ny) = sin(pi*(xmin + (i-1)*dx(i))) * sin(pi*(ymin + (ny-1)*dx(i)))
    end do

    ! Interior initial guess set to zero (or keep exact? set to 0 to test convergence)
    do j = 2, ny-1
      do i = 2, nx-1
        phi(i,j) = 0.0_r8
      end do
    end do
  end subroutine init_rhs_and_boundary


  subroutine sor_solve(f, phi, nx, ny, dx, dy, tol, maxiter, omega)
    implicit none
    integer, intent(in) :: nx, ny, maxiter
    real(r8), intent(in) :: dx(nx), dy(ny), tol, omega
    real(r8), intent(in) :: f(nx,ny)
    real(r8), intent(inout) :: phi(nx,ny)
    integer :: iter, i, j
    real(r8) :: dx2, dy2, denom, diff, error

    do iter = 1, maxiter
      error = 0.0_r8
      ! Red-black or simple Gauss-Seidel ordering; here simple
      do j = 2, ny-1
        do i = 2, nx-1
          dx2   = dx(i)*dx(i); dy2 = dy(j)*dy(j)
          denom = 2.0_r8*(1.0_r8/dx2 + 1.0_r8/dy2)
          ! compute Gauss-Seidel update then apply SOR
          diff = (((phi(i+1,j) + phi(i-1,j))/dx2) + &
                  ((phi(i,j+1) + phi(i,j-1))/dy2) - f(i,j)) / denom
          diff = diff - phi(i,j)  ! new - old (Gauss-Seidel new would be diff + phi(i,j))
          phi(i,j) = phi(i,j) + omega * diff
          error = error + abs(omega * diff)
print*,error,diff
        end do
      end do
      !if (mod(iter,100) == 0) then
      !  print '(A,I6,A,ES12.5)', 'SOR iter=', iter, ' residual=', error
      !end if
      if (error < tol) then
        print '(A,I6)', 'SOR converged in iterations=', iter
        exit
      end if
    end do
    if (iter >= maxiter) print *, 'SOR: reached maxiter=', maxiter

  end subroutine sor_solve


  subroutine multigrid_solver(f, phi, nx, ny, dx, dy, tol)
    implicit none
    integer, intent(in) :: nx, ny
    real(r8), intent(in) :: dx(nx), dy(ny), tol
    real(r8), intent(in) :: f(nx,ny)
    real(r8), intent(inout) :: phi(nx,ny)

    ! Build multigrid hierarchy: level 1 = fine
    real(r8) :: res_norm
    integer :: cycles, maxcycles

    integer :: maxlevels, lvl
    integer, allocatable :: nxl(:), nyl(:)
    real(r8), allocatable :: dxl(:,:), dyl(:,:)
    type(mg_level), allocatable :: levels(:)

    maxlevels = floor(log(real(min(nx,ny)-1,r8)) / log(2.0_r8))

    if (maxlevels < 1) maxlevels = 1

    allocate(nxl(maxlevels), nyl(maxlevels), dxl(maxlevels,nx), dyl(maxlevels,ny))
    allocate(levels(maxlevels))

    ! initialize fine level sizes
    nxl(1) = nx; nyl(1) = ny; dxl(1,1:nx) = dx(1:nx); dyl(1,1:ny) = dy(1:ny)

    do lvl = 2, maxlevels
      nxl(lvl) = (nxl(lvl-1) - 1) / 2 + 1
      nyl(lvl) = (nyl(lvl-1) - 1) / 2 + 1
      dxl(lvl,1:ny) = dxl(lvl-1,1:ny) * 2.0_r8
      dyl(lvl,1:ny) = dyl(lvl-1,1:ny) * 2.0_r8
    end do
    ! allocate level arrays
    do lvl = 1, maxlevels
      call levels(lvl)%allocate_level(nxl(lvl), nyl(lvl),dxl(lvl,1:nx), dyl(lvl,1:ny))
    end do

    ! fill fine level f and phi
    call levels(1)%set_grid(f, phi)
    ! initialize coarse rhs/phi to zero
    do lvl = 2, maxlevels
      call levels(lvl)%zero_fields()
    end do

    ! Multigrid cycles until residual norm < tol
    maxcycles = 1000

    do cycles = 1, maxcycles

      call v_cycle(levels, 1, maxlevels)

      call levels(1)%compute_residual_norm(res_norm)

      !if (mod(cycles,1) == 0) print '(A,I4,A,ES12.5)', 'V-cycle=', cycles, ' residual norm=', res_norm
      if (res_norm < tol) then
        print '(A,I4)', 'Multigrid converged in cycles=', cycles
        exit
      end if
    end do
    if (cycles >= maxcycles) print *, 'Multigrid: reached max cycles=', maxcycles

    ! copy back fine phi
    call levels(1)%get_phi(phi)

    ! deallocate
    do lvl = 1, maxlevels
      call levels(lvl)%deallocate_level()
    end do
    deallocate(nxl,nyl,dxl,dyl,levels)
  end subroutine multigrid_solver

  subroutine mg_allocate_level(this, nx_in, ny_in,dx_in,dy_in)
    class(mg_level), intent(inout) :: this
    integer, intent(in)  :: nx_in, ny_in
    real(r8), intent(in) :: dx_in(nx_in),dy_in(ny_in)
    this%nx = nx_in; this%ny = ny_in
    allocate(this%dx(nx_in))
    this%dx = dx_in;
    allocate(this%dy(ny_in))
    this%dy = dy_in
    allocate(this%phi(nx_in, ny_in))
    allocate(this%rhs(nx_in, ny_in))
    allocate(this%res(nx_in, ny_in))
    allocate(this%corr(nx_in, ny_in))
  end subroutine mg_allocate_level

  subroutine mg_deallocate_level(this)
    class(mg_level), intent(inout) :: this
    if (allocated(this%dx)) deallocate(this%dx)
    if (allocated(this%dy)) deallocate(this%dy)
    if (allocated(this%phi)) deallocate(this%phi)
    if (allocated(this%rhs)) deallocate(this%rhs)
    if (allocated(this%res)) deallocate(this%res)
    if (allocated(this%corr)) deallocate(this%corr)
  end subroutine mg_deallocate_level

  subroutine mg_set_grid(this, rhs_in, phi_in)
    class(mg_level), intent(inout) :: this
    real(r8), intent(in) :: rhs_in(:,:), phi_in(:,:)
    ! assume shapes match
    this%rhs = rhs_in
    this%phi = phi_in
  end subroutine mg_set_grid

  subroutine mg_get_phi(this, out_phi)
    class(mg_level), intent(in) :: this
    real(r8), intent(out) :: out_phi(this%nx, this%ny)
    ! copy internal phi to external array
    out_phi = this%phi
  end subroutine mg_get_phi

  subroutine mg_zero_fields(this)
    class(mg_level), intent(inout) :: this
    this%phi = 0.0_r8
    this%rhs = 0.0_r8
    this%res = 0.0_r8
    this%corr = 0.0_r8
  end subroutine mg_zero_fields


  subroutine mg_compute_residual_norm(this, norm)
    class(mg_level), intent(inout) :: this
    real(r8), intent(out) :: norm
    integer :: i,j
    real(r8) :: dx2, dy2, rloc

    ! residual r = rhs - A phi
    norm = 0.0_r8
    do j = 2, this%ny-1
      do i = 2, this%nx-1
        dx2 = (this%dx(i))**2
        dy2 = (this%dy(j))**2

        rloc = this%rhs(i,j) - ((this%phi(i+1,j) - 2.0_r8*this%phi(i,j) + this%phi(i-1,j))/dx2 + &
                                (this%phi(i,j+1) - 2.0_r8*this%phi(i,j) + this%phi(i,j-1))/dy2)
        norm = norm + abs(rloc)
      end do
    end do
  end subroutine mg_compute_residual_norm

  recursive subroutine v_cycle(levels, ilevel, maxlevels)
    type(mg_level), intent(inout) :: levels(:)
    integer, intent(in) :: ilevel, maxlevels
    ! Pre-smoothing
    call gauss_seidel(levels(ilevel)%phi, levels(ilevel)%rhs, levels(ilevel)%nx, levels(ilevel)%ny, &
                                                              levels(ilevel)%dx, levels(ilevel)%dy, 2)

    ! compute residual
    call compute_residual(levels(ilevel)%phi, levels(ilevel)%rhs, &
                          levels(ilevel)%res, levels(ilevel)%nx, &
                          levels(ilevel)%ny, levels(ilevel)%dx, &
                          levels(ilevel)%dy)


    if (ilevel == maxlevels) then
      ! On coarsest level solve directly with a few GS iterations
      levels(ilevel)%phi = 0.0_r8

      call gauss_seidel(levels(ilevel)%phi, levels(ilevel)%rhs, levels(ilevel)%nx, levels(ilevel)%ny, &
                                                                levels(ilevel)%dx, levels(ilevel)%dy, 50)
    else
      ! Restrict residual to coarser rhs

      call restrict_full_weighting(levels(ilevel)%res, levels(ilevel+1)%rhs, levels(ilevel)%nx, levels(ilevel)%ny)
      ! initialize coarse phi to zero
      levels(ilevel+1)%phi = 0.0_r8
      ! recursive V-cycle on coarse grid
      call v_cycle(levels, ilevel+1, maxlevels)
      ! Prolongate coarse correction and add
      call prolongate_and_add(levels(ilevel+1)%phi, levels(ilevel)%phi, levels(ilevel+1)%nx, levels(ilevel+1)%ny)
      ! Post-smoothing
      call gauss_seidel(levels(ilevel)%phi, levels(ilevel)%rhs, levels(ilevel)%nx, levels(ilevel)%ny, &
                                                                levels(ilevel)%dx, levels(ilevel)%dy, 2)
    end if
  end subroutine v_cycle


  subroutine gauss_seidel(phi, rhs, nx, ny,dx, dy, niter)
    implicit none
    integer , intent(in) :: nx, ny, niter
    real(r8), intent(inout) :: phi(nx,ny)
    real(r8), intent(in) :: rhs(nx,ny)
    real(r8), intent(in) :: dx(nx), dy(ny)
    integer :: iter, i, j
    real(r8) :: dx2, dy2, denom
    do iter = 1, niter
      do j = 2, ny-1
        do i = 2, nx-1
           dx2 = dx(i)*dx(i); dy2 = dy(j)*dy(j)
           denom = 2.0_r8*(1.0_r8/dx2 + 1.0_r8/dy2)
           phi(i,j) = ((phi(i+1,j) + phi(i-1,j))/dx2 + (phi(i,j+1) + phi(i,j-1))/dy2 - rhs(i,j)) / denom
        end do
      end do
    end do
  end subroutine gauss_seidel


  subroutine compute_residual(phi, rhs, res, nx, ny,dx,dy)
    implicit none
    integer, intent(in) :: nx, ny
    real(r8), intent(in) :: dx(nx), dy(ny)
    real(r8), intent(in) :: phi(nx,ny), rhs(nx,ny)
    real(r8), intent(out) :: res(nx,ny)
    integer :: i,j
    real(r8) :: dx2, dy2
    res = 0.0_r8
    do j = 2, ny-1
      do i = 2, nx-1
         dx2 = dx(i)*dx(i); dy2 = dy(j)*dy(j)
         res(i,j) = rhs(i,j) - (((phi(i+1,j) - 2.0_r8*phi(i,j) + phi(i-1,j))/dx2) + &
                                ((phi(i,j+1) - 2.0_r8*phi(i,j) + phi(i,j-1))/dy2)   )
      end do
    end do
  end subroutine compute_residual

  subroutine restrict_full_weighting(fine, coarse, nx_f, ny_f)
    implicit none
    real(r8), intent(in) :: fine(nx_f, ny_f)
    real(r8), intent(out) :: coarse((nx_f-1)/2+1, (ny_f-1)/2+1)
    integer, intent(in) :: nx_f, ny_f
    integer :: i_c, j_c, i_f, j_f
    integer :: nx_c, ny_c
    nx_c = (nx_f-1)/2 + 1
    ny_c = (ny_f-1)/2 + 1
    coarse = 0.0_r8

    do j_c = 2, ny_c-1
      j_f = 2*(j_c-1)
      do i_c = 2, nx_c-1
        i_f = 2*(i_c-1)
        coarse(i_c,j_c) = ( fine(i_f-1,j_f-1) + 2.0_r8*fine(i_f,j_f-1) + fine(i_f+1,j_f-1) + &
                             2.0_r8*fine(i_f-1,j_f) + 4.0_r8*fine(i_f,j_f) + 2.0_r8*fine(i_f+1,j_f) + &
                             fine(i_f-1,j_f+1) + 2.0_r8*fine(i_f,j_f+1) + fine(i_f+1,j_f+1) ) / 16.0_r8
      end do
    end do
    ! boundaries: simple injection (copy)
    do j_c = 1, ny_c
      coarse(1,j_c) = fine(1, 2*(j_c-1)+1)
      coarse(nx_c,j_c) = fine(nx_f, 2*(j_c-1)+1)
    end do
    do i_c = 1, nx_c
      coarse(i_c,1) = fine(2*(i_c-1)+1, 1)
      coarse(i_c,ny_c) = fine(2*(i_c-1)+1, ny_f)
    end do
  end subroutine restrict_full_weighting


  subroutine prolongate_and_add(coarse_phi, fine_phi, nx_c, ny_c)
   implicit none
   real(r8), intent(in) :: coarse_phi(nx_c, ny_c)
   real(r8), intent(inout) :: fine_phi( (nx_c-1)*2+1, (ny_c-1)*2+1 )
   integer, intent(in) :: nx_c, ny_c
   integer :: i_c, j_c, i_f, j_f, nx_f, ny_f
   integer :: ic, jc
   real(r8) :: val


   nx_f = (nx_c-1)*2 + 1
   ny_f = (ny_c-1)*2 + 1


   ! inject coarse grid points (odd indices on fine grid)
   do j_c = 1, ny_c
      j_f = 2*(j_c-1) + 1
      do i_c = 1, nx_c
         i_f = 2*(i_c-1) + 1
         fine_phi(i_f, j_f) = fine_phi(i_f, j_f) + coarse_phi(i_c, j_c)
      end do
   end do


   ! interpolate for fine points where i is even, j is odd (horizontal edges)
   do j_c = 1, ny_c
      j_f = 2*(j_c-1) + 1
      do i_f = 2, nx_f-1, 2
         ic = nint(real(i_f)/2.0)
         fine_phi(i_f, j_f) = fine_phi(i_f, j_f) + 0.5_r8 * (coarse_phi(ic, j_c) + coarse_phi(ic+1, j_c))
       end do
   end do


   ! interpolate for fine points where i is odd, j is even (vertical edges)
   do i_c = 1, nx_c
      i_f = 2*(i_c-1) + 1
      do j_f = 2, ny_f-1, 2
         jc = nint(real(j_f)/2.0)
         fine_phi(i_f, j_f) = fine_phi(i_f, j_f) + 0.5_r8 * (coarse_phi(i_c, jc) + coarse_phi(i_c, jc+1))
      end do
   end do


   ! interpolate for fine points where i and j are even (cell centers)
   do j_f = 2, ny_f-1, 2
      do i_f = 2, nx_f-1, 2
         ic = nint(real(i_f)/2.0)
         jc = nint(real(j_f)/2.0)
         val = 0.25_r8 * (coarse_phi(ic, jc  ) + coarse_phi(ic+1, jc  ) + &
                          coarse_phi(ic, jc+1) + coarse_phi(ic+1, jc+1))
         fine_phi(i_f, j_f) = fine_phi(i_f, j_f) + val
      end do
   end do
  end subroutine prolongate_and_add
  
  
  subroutine report_error(phi, nx, ny, xmin, ymin, dx, dy)
    implicit none
    real(r8), intent(in) :: phi(nx,ny)
    integer, intent(in) :: nx, ny
    real(r8), intent(in) :: xmin, ymin, dx(nx), dy(ny)
    integer :: i,j
    real(r8) :: x,y,err,errL2,pi
    pi=function_pi()
    errL2 = 0.0_r8
    do j = 1, ny
      y = ymin + (j-1)*dy(j)
      do i = 1, nx
        x = xmin + (i-1)*dx(i)
        err = abs(phi(i,j) - sin(pi*x)*sin(pi*y))
        errL2 = errL2 + err*err
      end do
    end do
    errL2 = sqrt(errL2 / real((nx)*(ny), r8))
    print '(A,ES12.5)', 'L2 error norm =', errL2
  end subroutine report_error  

  real(r8) function function_pi()
    function_pi = 4.0_r8*atan(1.0_r8)
  end function function_pi


END MODULE Class_Module_LapaceOperator


!PROGRAM Main
!USE Class_Module_LapaceOperator, only : r8,&
!                                        init_rhs_and_boundary,&
!                                        sor_solve,&
!                                        multigrid_solver,&
!                                        report_error
! ! Parameters
! integer, parameter :: nx = 129         ! number of grid points in x (must be 2^k + 1)
! integer, parameter :: ny = 129
! real(r8), parameter :: xmin = 0.0_r8, xmax = 1.0_r8
! real(r8), parameter :: ymin = 0.0_r8, ymax = 1.0_r8
! integer, parameter :: maxiter_sor = 20000
! real(r8), parameter :: tol = 1.0e-8_r8
! real(r8), parameter :: omega = 1.7_r8   ! SOR relaxation

! ! Arrays
! real :: start, finish
! real(r8),allocatable :: dx(:), dy(:)
! real(r8),allocatable :: t0, t1
! character(len=1) :: choice
! real(r8), allocatable :: phi(:,:), f(:,:)
! allocate(phi(nx,ny), f(nx,ny))
! allocate(dx(nx), dy(ny))

! dx(1:nx) = (xmax - xmin) / real(nx-1, r8)
! dy(1:ny) = (ymax - ymin) / real(ny-1, r8)

! call init_rhs_and_boundary(phi, f, nx, ny, xmin, ymin, dx, dy)
! print *, 'Poisson solver (2D) - options: (s) SOR, (m) Multigrid'
! print *, 'Grid: ', nx, ' x ', ny, ' dx=', dx(1), ' dy=', dy(1)
! print *, 'Enter choice (s/m):'
! !read(*,'(A)') choice
! choice = 'm'
! choice = adjustl(choice(1:1))
! !  print '("Time = ",f6.3," seconds.")',finish-start
! if (choice == 's' .or. choice == 'S') then
!   print *, 'Running SOR ...'
!   call cpu_time(start)
!   t0 = start
!   call sor_solve(f, phi, nx, ny, dx, dy, tol, maxiter_sor, omega)
!   call cpu_time(finish)
!   t1 = finish
!   print *, 'SOR finished. Time (s)=', t1-t0
!   call report_error(phi, nx, ny, xmin, ymin, dx, dy)
! else if (choice == 'm' .or. choice == 'M') then
!   print *, 'Running Multigrid V-cycles ...'
!   call cpu_time(start)
!   t0 = start
!   call multigrid_solver(f, phi, nx, ny, dx, dy, tol)
!   call cpu_time(finish)
!   t1 = finish
!   print *, 'Multigrid finished. Time (s)=', t1-t0
!   call report_error(phi, nx, ny, xmin, ymin, dx, dy)
! else
!   print *, 'Invalid choice. Exiting.'
! end if
!END PROGRAM Main
