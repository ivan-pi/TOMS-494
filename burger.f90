module mod_pdeone

    implicit none
    public

    integer, parameter :: npde = 1 
    integer, parameter :: npts = 201
    integer, parameter :: icord = 0
    real :: x(npts)

end module

module mod_user
    
    implicit none
    private

    public :: mypde,jac ! used in lsode
    public :: f,d,bndry ! used by pdeone
    public :: trusol    ! true solution

contains

    subroutine mypde(neq,t,y,ydot)
        use mod_pdeone, only: npde,npts
        integer, intent(in) :: neq
        real, intent(in) :: t
        real, intent(in) :: y(neq)
        real, intent(in) :: ydot(neq)
        external :: pdeone
        call pdeone(t,y,ydot,npde,npts)
    end subroutine

    subroutine jac(neq,t,y,ml,mu,pd,nrowpd)
        integer neq, ml, mu, nrowpd
        real t, y(*), pd(nrowpd,*)
        ! dummy
    end subroutine

!
!  F ROUTINE FOR BURGERS EQUATION
!
    subroutine f(t,x,u,ux,duxx,fval,npde)
        integer, intent(in) :: npde
        real, intent(in) :: t
        real, intent(in) :: x
        real, intent(in) :: u(npde),ux(npde),duxx(npde,npde)
        real, intent(out) :: fval(npde)
        fval(1) = duxx(1,1) - u(1) * ux(1)
    end subroutine

!
!  D ROUTINE FOR BURGERS EQUATION
!
    subroutine d(t,x,u,dval,npde)
        integer, intent(in) :: npde
        real, intent(in) :: t, x
        real, intent(in) :: u(npde)
        real, intent(out) :: dval(npde,npde)
        dval(1,1) = .003
    end subroutine

!
!  BNDRY ROUTINE FOR BURGERS EQUATION
!
    subroutine bndry(t,x,u,alpha,beta,gamma_,npde)
        use mod_pdeone, only: xm => x
        integer, intent(in) :: npde
        real, intent(in) :: t, x
        real, intent(in) :: u(npde)
        real, intent(out) :: alpha(npde), beta(npde),gamma_(npde)
        alpha(1) = 1.
        beta(1) = 0.
        gamma_(1) = trusol(t,x)
    end subroutine

!
!  TRUSOL IS A FUNCTION ROUTINE WHICH COMPUTES THE KNOWN TRUE SOLUTION
!
    elemental real function trusol(t,x)
        real, intent(in) :: t, x
        real :: anu,a,b,c
        anu = 3.e-3
        a = -(x-0.5 + 4.95*t)/(20.*anu)
        b = -(x-0.5 + 0.75*t)/(4.*anu)
        c = -(x-0.375)/(2.*anu)
        a = exp(a)
        b = exp(b)
        c = exp(c)
        trusol = ( 0.1*a + 0.5*b + c) / ( a + b + c)
    end function

end module

!>  *Burgers' Equation*
!   We consider the numerical solution of the equation
!   $$ \partial u/\partialt = -u \frac{\partial u}{\partial x} + .003\frac{\partial^2 u}{\partial x^2}$$
!   on the $x$ interval $[0,1]$ and for $0 < t \leq 1.1$. We assume initial
!   and boundary conditions of the form
!   $$ u(0,x)=\phi(0,x), \quad u(t,0) = \phi(t,0), \quad u(t,1) = \phi(t,1) $$
!   where
!   $$ \phi(t,x) = (.1e^{-A} + .5e^{-B} + e^{-C})/(e^{-A} + e^{-B} + e^{-C}) $$
!   and $A = \frac{50}{3}(x -.5 + 4.95*t)$, $B = \frac{250}{3}(x - .5 + 0.75t)$, $C=\frac{500}{3}(x-.375)$.
!
program burger
    
    use mod_pdeone, only: npde, npts, icord, x
    use mod_user, only: mypde,jac,trusol
    implicit none

    real, allocatable :: u(:)
    integer :: i, unit

    integer :: lrw, liw, mf, ml, mu
    integer :: neq, itol, itask, istate, iopt
    integer, allocatable ::iwork(:)

    real :: t,tout,rtol,atol
    real, allocatable :: rwork(:)

    open(newunit=unit,file="pdeout.txt")

    !
    ! Define problem parameters
    !
    ! mf = 10     ! Method flag 
    ! mf = 22
    mf = 25

    neq = npde*npts ! Number of first order ODEs
    allocate(u(neq))
    
    select case(mf)
    case(10)
        ! Nonstiff (Adams) method, no Jacobian used
        lrw = 20 + 16*neq 
        liw = 20
    case(22)
        ! Stiff method, internally generated full Jacobian
        lrw = 22 + 9*neq + neq**2
        liw = 20 + neq
    case(25)
        ! Stiff method, internally generated banded Jacobian
        ml = 2*npde-1; mu = 2*npde-1;
        lrw = 22 + 10*neq + (2*ml + mu)*neq
        liw = 20 + neq
    end select
    allocate(rwork(lrw))
    allocate(iwork(liw))
    if (mf == 25) iwork(1:2) = [ml,mu]

    itol = 1        ! Scalar absolute eror control
    rtol = 1.e-5    ! Allowed relative error
    atol = 0.0      ! Pure relative error control
    itask = 1       ! Normal computation
    istate = 1      ! First call to the problem
    iopt = 0        ! No optional inputs used

    !
    ! Define mesh and initial data
    !
    t = 0.0
    tout = 0.0
    x = [(real(i-1)*0.005,i=1,npts)]
    u = trusol(t,x)

    do while (t < 1.1)
        write(unit,25) t
        write(unit,30) (u(i),i=1,npts)
        tout = tout + 0.1

        !
        ! Call ODE Integrator
        !
        call slsode(mypde,neq,u,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
        if (istate /= 2) then
            write(unit,35) istate
        end if
        print *, t, iwork(11), iwork(12), iwork(13)
    end do

    close(unit)
    
25  format("t =  ",e10.3)
30  format(5e20.6)
35  format("istate = ",i3)
    
end program


