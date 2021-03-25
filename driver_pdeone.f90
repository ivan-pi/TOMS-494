module burger_example

    use pdeone_mod, only: wp
    implicit none
    private

    public :: f, d, bndry, trusol

contains
!
!  F ROUTINE FOR BURGERS EQUATION
!
    subroutine f(t,x,u,ux,duxx,fval,npde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: t
        real(wp), intent(in) :: x
        real(wp), intent(in) :: u(npde),ux(npde),duxx(npde,npde)
        real(wp), intent(out) :: fval(npde)
        fval(1) = duxx(1,1) - u(1) * ux(1)
    end subroutine

!
!  D ROUTINE FOR BURGERS EQUATION
!
    subroutine d(t,x,u,dval,npde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: t, x
        real(wp), intent(in) :: u(npde)
        real(wp), intent(out) :: dval(npde,npde)
        dval(1,1) = .003_wp
    end subroutine

!
!  BNDRY ROUTINE FOR BURGERS EQUATION
!
    subroutine bndry(t,x,u,alpha,beta,gamma_,npde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: t, x
        real(wp), intent(in) :: u(npde)
        real(wp), intent(out) :: alpha(npde), beta(npde),gamma_(npde)
        alpha(1) = 1._wp
        beta(1) = 0._wp
        gamma_(1) = trusol(t,x)
    end subroutine

!
!  TRUSOL IS A FUNCTION ROUTINE WHICH COMPUTES THE KNOWN TRUE SOLUTION
!
    elemental real(wp) function trusol(t,x)
        real(wp), intent(in) :: t, x
        real(wp) :: anu,a,b,c
        anu = 3.e-3_wp
        a = -(x-0.5_wp + 4.95_wp*t)/(20._wp*anu)
        b = -(x-0.5_wp + 0.75_wp*t)/(4._wp*anu)
        c = -(x-0.375_wp)/(2._wp*anu)
        a = exp(a)
        b = exp(b)
        c = exp(c)
        trusol = ( 0.1_wp*a + 0.5_wp*b + c) / ( a + b + c)
    end function

end module

module matlab_example

    use pdeone_mod, only: wp
    implicit none
    private
    real(wp), parameter :: pi = 4._wp*atan(1.0_wp)

contains

    subroutine f(t,x,u,ux,duxx,fval,npde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: t
        real(wp), intent(in) :: x
        real(wp), intent(in) :: u(npde),ux(npde),duxx(npde,npde)
        real(wp), intent(out) :: fval(npde)
        fval(1) = duxx(1,1)/pi**2
    end subroutine

    subroutine d(t,x,u,dval,npde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: t, x
        real(wp), intent(in) :: u(npde)
        real(wp), intent(out) :: dval(npde,npde)
        dval(1,1) = 1._wp
    end subroutine

    ! subroutine bndry(t,x,u,alpha,beta,gamma_,npde)
    !     integer, intent(in) :: npde
    !     real(wp), intent(in) :: t, x
    !     real(wp), intent(in) :: u(npde)
    !     real(wp), intent(out) :: alpha(npde), beta(npde),gamma_(npde)
    !     alpha(1) = 1._wp
    !     beta(1) = 0._wp
    !     gamma_(1) = trusol(t,x)
    ! end subroutine

end module

module nonlinear_diffusion

    use pdeone_mod, only: wp, pde_t
    implicit none
    private

    type, extends(pde_t), public :: diffeq_t
    contains
        procedure :: f
        procedure :: d
        procedure :: bndry
    end type

contains

    subroutine f(this,t,x,u,ux,duxx,fval)
        class(diffeq_t), intent(in) :: this
        real(wp), intent(in) :: t
        real(wp), intent(in) :: x
        real(wp), intent(in) :: u(this%npde),ux(this%npde),duxx(this%npde,this%npde)
        real(wp), intent(out) :: fval(this%npde)
        fval(1) = duxx(1,1) - u(1)**2
    end subroutine

!
!  D ROUTINE FOR BURGERS EQUATION
!
    subroutine d(this,t,x,u,dval)
        class(diffeq_t), intent(in) :: this
        real(wp), intent(in) :: t, x
        real(wp), intent(in) :: u(this%npde)
        real(wp), intent(out) :: dval(this%npde,this%npde)
        dval(1,1) = u(1)
    end subroutine

!
!  BNDRY ROUTINE FOR BURGERS EQUATION
!
    subroutine bndry(this,t,x,u,alpha,beta,gamm)
        class(diffeq_t), intent(in) :: this
        real(wp), intent(in) :: t, x
        real(wp), intent(in) :: u(this%npde)
        real(wp), intent(out) :: alpha(this%npde), beta(this%npde),gamm(this%npde)
        
        if (x == this%xm(1)) then
            alpha(1) = 1._wp
            beta(1) = 0._wp
            gamm(1) = 50._wp
        else
            alpha(1) = 0.0_wp
            beta(1) = 1._wp
            gamm = 1._wp - sin(u(1))
        end if
    end subroutine

end module

program main

    use pdeone_mod, only: wp, pdeone_type, CARTESIAN
    use burger_example, only: d, f, bndry, trusol
    use nonlinear_diffusion, only: diffeq_t
    implicit none

    type(pdeone_type) :: pde
    type(diffeq_t) :: pde_g

    integer, parameter :: icord = CARTESIAN
    integer, parameter :: npts = 51
    integer, parameter :: npde = 1
    real(wp), allocatable :: x(:)

    real(wp), allocatable :: u(:)
    integer :: i, unit

    integer :: lrw, liw, mf, ml, mu
    integer :: neq, itol, itask, istate, iopt, nnz, lwm
    integer, allocatable ::iwork(:)

    real(wp) :: t,tout,rtol,atol
    real(wp), allocatable :: rwork(:)

    open(newunit=unit,file="pdeout.txt")

    !
    ! Define problem parameters
    !
    ! mf = 10     ! Method flag 
    ! mf = 22
    mf = 25
    ! mf = 222 ! SLSODES

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
        ml = 2*npde-1; mu = 2*npde-1; ! upper and lower bandwidth
        lrw = 22 + 10*neq + (2*ml + mu)*neq
        liw = 20 + neq
    ! case(222) ! SLSODES
    !     nnz = 4*neq ! estimate of non-zero elements in jacobian
    !     lwm = 2*nnz + 2*neq + (nnz + 10*neq)
    !     lrw = 20 + 9*neq + lwm
    !     liw = neq
    end select
    allocate(rwork(lrw))
    allocate(iwork(liw))
    if (mf == 25) iwork(1:2) = [ml,mu]

    itol = 1        ! Scalar absolute eror control
    rtol = 1.e-5_wp    ! Allowed relative error
    atol = 0.0_wp      ! Pure relative error control
    itask = 1       ! Normal computation
    istate = 1      ! First call to the problem
    iopt = 0        ! No optional inputs used

    !
    ! Define mesh and initial data
    !
    t = 0.0_wp
    tout = 0.0_wp
    allocate(x(npts))
    x = [(real(i-1)*0.005_wp,i=1,npts)]
    u = trusol(t,x)


    !
    ! Define PDE
    !
    pde = pdeone_type(1,x,icord,d,f,bndry) ! makes a copy of the grid

    do while (t < 1.1_wp)
        ! write(unit,25) t
        ! write(unit,30) (u(i),i=1,npts)
        call output(unit,x,u)
        tout = tout + 0.1_wp

        !
        ! Call ODE Integrator
        !
        call dlsode(mypde,neq,u,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
        if (istate /= 2) then
            write(unit,35) istate
        end if
        print *, t, iwork(11), iwork(12), iwork(13)
    end do

    close(unit)

    !
    ! EXAMPLE G
    !

    open(newunit=unit,file="example_g.out")

    itol = 1        ! Scalar absolute eror control
    rtol = 1.e-4_wp    ! Allowed relative error
    atol = 0.0_wp      ! Pure relative error control
    itask = 1       ! Normal computation
    istate = 1      ! First call to the problem
    iopt = 0        ! No optional inputs used

    t = 0.0_wp
    tout = 0.0_wp

    call example_g(npts,x,u)
    pde_g%icord = CARTESIAN
    pde_g%npde = npde
    pde_g%npts = npts
    allocate(pde_g%xm,source=x)


    do while(t < 0.11_wp)

        ! call output(unit,x,u)
        call output(unit,x,u)
        tout = tout + 1.e-5_wp
        
        !
        ! Call ODE Integrator
        !
        call dlsode(mypde_g,neq,u,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,mf)
        if (istate /= 2) then
            write(unit,35) istate
        end if
        print *, t, iwork(11), iwork(12), iwork(13)
    end do

    close(unit)
    
25  format("t =  ",e10.3)
30  format(5e20.6)
35  format("istate = ",i3)

contains

    subroutine output(unit,x,u)
        integer, intent(in) :: unit
        real(wp), intent(in) :: x(:)
        real(wp), intent(in) :: u(:)

        do i = 1, min(size(x),size(u))
            write(unit,*), x(i), u(i)
        end do
        write(unit,*)
    end subroutine

    ! required by ODEPACK
    subroutine mypde(neq,t,y,ydot)
        integer, intent(in) :: neq
        real(wp), intent(in) :: t
        real(wp) :: y(neq)
        real(wp) :: ydot(neq)
        ! external :: pdeone
        ! call pdeone(t,y,ydot,npde,npts)
        call pde%pdeone(t,y,ydot)
    end subroutine

    subroutine jac(neq,t,y,ml,mu,pd,nrowpd)
        integer neq, ml, mu, nrowpd
        real t, y(*), pd(nrowpd,*)
        ! dummy
    end subroutine


    subroutine example_g(n,x,u)
        integer, intent(in) :: n
        real(wp), intent(out) :: x(n), u(n)
        real(wp) :: dx

        dx = 1._wp/(n-1)
        x = [((i-1)*dx,i=1,n)]
        u = 100._wp
    end subroutine

    ! required by ODEPACK
    subroutine mypde_g(neq,t,y,ydot)
        integer, intent(in) :: neq
        real(wp), intent(in) :: t
        real(wp) :: y(neq)
        real(wp) :: ydot(neq)
        ! external :: pdeone
        ! call pdeone(t,y,ydot,npde,npts)
        call pde_g%pdeone(t,y,ydot)
    end subroutine

end program