module drying_problem

  use pdeone_mod, only: wp
  implicit none
  private

  ! public :: dopri5

  public :: drying_diff
  public :: drying_rhs
  public :: drying_bndry

contains

  !> Diffusivity function
  subroutine drying_diff(t,x,u,dval,npde)
      integer, intent(in) :: npde
      real(wp), intent(in) :: t, x
      real(wp), intent(in) :: u(npde)
      real(wp), intent(out) :: dval(npde,npde)
      print *, "Hello from diff"
      dval(1,1) = 0.1_wp
  end subroutine

  !> PDE function
  subroutine drying_rhs(t,x,u,ux,duxx,fval,npde)
      integer, intent(in) :: npde
      real(wp), intent(in) :: t
      real(wp), intent(in) :: x
      real(wp), intent(in) :: u(npde),ux(npde),duxx(npde,npde)
      real(wp), intent(out) :: fval(npde)
      print *, "Hello from RHS"
      fval(1) = duxx(1,1)
  end subroutine

  !> Boundary function
  subroutine drying_bndry(left,t,x,u,alpha,beta,gamm,npde)
      logical, intent(in) :: left
      integer, intent(in) :: npde
      real(wp), intent(in) :: t, x
      real(wp), intent(in) :: u(npde)
      real(wp), intent(out) :: alpha(npde), beta(npde), gamm(npde)


      if (left) then
      print *, "Hello from LEFT bndry"
        alpha(1) = 0.0_wp
        beta(1) = 1.0_wp
        gamm(1) = 1.0_wp
      else
      print *, "Hello from RIGHT bndry"
        alpha(1) = 1.0_wp
        beta(1) = 0.0_wp
        gamm(1) = 0.12_wp
      end if
  end subroutine

end module

program main_drying_example

  use iso_fortran_env, only: output_unit
  use pdeone_mod, only: wp, pdeone_type
  use drying_problem, only: drying_rhs, drying_diff, drying_bndry
  use dopri5_ode, only: dopri5
  implicit none

  integer, parameter :: npde = 1
  integer, parameter :: icord = 0 ! SLAB

  integer, parameter :: n = 10
  real(wp) :: L0
  real(wp), target :: x(n), y(n)

  type(pdeone_type) :: my_pde

  real(wp) :: dx, t, tend
  real(wp) :: rtol(1), atol(1)
  integer :: itol

  integer, parameter :: iout = 1

  real(wp) :: rdum(1)
  integer :: idum(1)

  integer :: lwork
  real(wp), allocatable :: work(:)

  integer :: liwork
  integer, allocatable :: iwork(:)

  integer :: i, idid

  ! external :: dopri5

  L0 = 0.001_wp
  dx = L0/real(n,wp)
  x = [((i-1)*dx,i=1,n)]

  ! The parabolic partial differential equation
  my_pde = pdeone_type(npde,x,icord, &
    drying_diff, &
    drying_rhs, &
    drying_bndry)

  ! Initial moisture concentration
  y = 0.3_wp

  ! Time interval
  t = 0.0_wp
  tend = 10.0_wp

  ! Tolerances
  rtol(1) = 0.0001_wp
  atol(1) = 0.00001_wp
  itol = 0

  ! Integer work array
  liwork = 21
  allocate(iwork(liwork))
  iwork(:) = 0

  iwork(3) = output_unit

  ! Real work arrays
  lwork = 8*n + 5*iwork(5) + 21
  allocate(work(lwork))
  work(:) = 0

  call dopri5(n,fcn,t,y,tend,rtol,atol,itol,solout,iout, &
        work,lwork,iwork,liwork,rdum,idum,idid)

contains

  !> dopri5 compatible interface to PDEONE
  subroutine fcn(n,x,y,f,rpar,ipar)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    real(wp), intent(inout) :: y(n)
    real(wp), intent(inout) :: f(n)
    real(wp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    print *, "In RHS FCN"
    call my_pde%pdeone(x, y, f)
    print *, "Exiting RHS FCN"
  end subroutine

  subroutine solout(nr,told,t,y,n,con,icomp,nd,rpar,ipar,irtrn)
    ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    ! IMPLICIT INTEGER (I-N)
    ! DIMENSION Y(N),CON(5*ND),ICOMP(ND)
    ! DIMENSION RPAR(*),IPAR(*)
    integer, intent(in) :: nr
    real(wp), intent(in) :: told, t
    real(wp), intent(inout) :: y(n)
    integer, intent(in) :: n
    real(wp), intent(inout) :: con(5*nd)
    integer, intent(inout) :: icomp(nd)
    integer, intent(in) :: nd
    real(wp), intent(inout) :: rpar(*)
    integer, intent(inout) :: ipar(*)
    integer, intent(inout) :: irtrn

    character(len=20) :: stepstr
    integer :: output_unit, i

    print *, "In SOLOUT"

    write(stepstr,'(I0.8)') nr

    open(newunit=output_unit,file="drying"//trim(stepstr)//".txt")

    do i = 1, n
      write(output_unit,*) x(i), y(i)
    end do

    close(output_unit)

    print *, "Exiting SOLOUT"

  end subroutine

end program