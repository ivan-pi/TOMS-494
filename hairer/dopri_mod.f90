MODULE DOPRI_MOD

  IMPLICIT NONE
  PUBLIC

  DOUBLE PRECISION :: XOLD
  DOUBLE PRECISION :: HOUT

  integer, parameter :: dp = kind(1.d0)

  abstract interface
    subroutine fcn_proc(n,x,y,f,rpar,ipar)
        import dp
      integer, intent(in) :: n
      real(dp), intent(in) :: x
      real(dp), intent(inout) :: y(n)
      real(dp), intent(inout) :: f(n)
      real(dp), intent(inout) :: rpar(*)
      integer, intent(inout) :: ipar(*)
    end subroutine
    subroutine solout_proc(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
        import dp
      integer, intent(in) :: nr
      real(dp), intent(in) :: xold, x
      real(dp), intent(inout) :: y(n)
      integer, intent(in) :: n
      real(dp), intent(inout) :: con(5*nd)
      integer, intent(inout) :: icomp(nd)
      integer, intent(in) :: nd
      real(dp), intent(inout) :: rpar(*)
      integer, intent(inout) :: ipar(*)
      integer, intent(inout) :: irtrn
    end subroutine
  end interface

  ! interface
  !   subroutine dopri5(n,fcn,x,y,xend,rtol,atol,itol,solout,iout, &
  !       work,lwork,iwork,liwork,rpar,ipar,idid)
  !       import dp
  !     integer, intent(in) :: n
  !       !! dimension of the system
  !     procedure(fcn_proc) :: fcn
  !     real(dp), intent(inout) :: x
  !       !! initial x-value
  !     real(dp), intent(inout) :: y(n)
  !       !! initial y-value
  !     real(dp), intent(in) :: xend
  !     real(dp), intent(in) :: rtol(*), atol(*)
  !     integer, intent(in) :: itol
  !     procedure(solout_proc) :: solout
  !     integer, intent(in) :: iout
  !     real(dp), intent(inout) :: work(lwork)
  !     integer, intent(in) :: lwork
  !     integer, intent(inout) :: iwork(liwork)
  !     integer, intent(in) :: liwork
  !     real(dp), intent(inout) :: rpar(*)
  !     integer, intent(inout) :: ipar(*)
  !     integer, intent(out) :: idid
  !   end subroutine
  ! end interface

CONTAINS

  SUBROUTINE CDOPRI(C2,C3,C4,C5,E1,E3,E4,E5,E6,E7, &
          A21,A31,A32,A41,A42,A43,A51,A52,A53,A54, &
          A61,A62,A63,A64,A65,A71,A73,A74,A75,A76, &
          D1,D3,D4,D5,D6,D7)
! ----------------------------------------------------------
!     RUNGE-KUTTA COEFFICIENTS OF DORMAND AND PRINCE (1980)
! ----------------------------------------------------------
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    
    C2 = 0.2_dp
    C3 = 0.3_dp
    C4 = 0.8_dp
    C5 = 8._dp/9._dp

    A21 = 0.2_dp
    A31 = 3._dp/40._dp
    A32 = 9._dp/40._dp
    A41 = 44._dp/45._dp
    A42 = -56._dp/15._dp
    A43 = 32._dp/9._dp
    A51 = 19372._dp/6561._dp
    A52 = -25360._dp/2187._dp
    A53 = 64448._dp/6561._dp
    A54 = -212._dp/729._dp
    A61 = 9017._dp/3168._dp
    A62 = -355._dp/33._dp
    A63 = 46732._dp/5247._dp
    A64 = 49._dp/176._dp
    A65 = -5103._dp/18656._dp
    A71 = 35._dp/384._dp
    A73 = 500._dp/1113._dp
    A74 = 125._dp/192._dp
    A75 = -2187._dp/6784._dp
    A76 = 11._dp/84._dp
    E1 = 71._dp/57600._dp
    E3 = -71._dp/16695._dp
    E4 = 71._dp/1920._dp
    E5 = -17253._dp/339200._dp
    E6 = 22._dp/525._dp
    E7 = -1._dp/40._dp  
! ---- DENSE OUTPUT OF SHAMPINE (1986)
    D1 = -12715105075._dp/11282082432._dp
    D3 = 87487479700._dp/32700410799._dp
    D4 = -10690763975._dp/1880347072._dp
    D5 = 701980252875._dp/199316789632._dp
    D6 = -1453857185._dp/822651844._dp
    D7 = 69997945._dp/29380423._dp
  END SUBROUTINE

END MODULE