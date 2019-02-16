module pdeone_mod

    implicit none
    private

    integer, public, parameter :: wp = kind(1.0d0)

    integer, public, parameter :: CARTESIAN = 0
    integer, public, parameter :: SLAB = CARTESIAN
    integer, public, parameter :: CYLINDRICAL = 1
    integer, public, parameter :: SPHERICAL = 2

    abstract interface
        subroutine d_interface(t,x,u,dval,npde)
            import wp
            integer, intent(in) :: npde
            real(wp), intent(in) :: t, x
            real(wp), intent(in) :: u(npde)
            real(wp), intent(out) :: dval(npde,npde)
        end subroutine
        subroutine f_interface(t,x,u,ux,duxx,fval,npde)
            import wp
            integer, intent(in) :: npde
            real(wp), intent(in) :: t
            real(wp), intent(in) :: x
            real(wp), intent(in) :: u(npde),ux(npde),duxx(npde,npde)
            real(wp), intent(out) :: fval(npde)
        end subroutine
        subroutine bndry_interface(t,x,u,alpha,beta,gamm,npde)
            import wp
            integer, intent(in) :: npde
            real(wp), intent(in) :: t, x
            real(wp), intent(in) :: u(npde)
            real(wp), intent(out) :: alpha(npde), beta(npde), gamm(npde)
        end subroutine
    end interface

    type, public :: pdeone_type
        private
        integer :: npde, npts
        real(wp), allocatable :: x(:)
        integer :: icord
        procedure(d_interface), pointer, nopass :: d => null()
        procedure(f_interface), pointer, nopass :: f => null()
        procedure(bndry_interface), pointer, nopass :: bndry => null()
    contains
        procedure :: pdeone
    end type

    type, abstract, public :: pde_t
        integer :: npde
        integer :: npts
        real(wp), pointer :: xm(:)
        integer :: icord
    contains
        procedure(oo_d), deferred :: d
        procedure(oo_f), deferred :: f
        procedure(oo_bndry), deferred :: bndry
        procedure :: pdeone => pdeone_oo
    end type

    abstract interface
        subroutine oo_d(this,t,x,u,dval)
            import wp, pde_t
            class(pde_t), intent(in) :: this
            real(wp), intent(in) :: t, x
            real(wp), intent(in) :: u(this%npde)
            real(wp), intent(out) :: dval(this%npde,this%npde)
        end subroutine
        subroutine oo_f(this,t,x,u,ux,duxx,fval)
            import wp, pde_t
            class(pde_t), intent(in) :: this
            real(wp), intent(in) :: t, x
            real(wp), intent(in) :: u(this%npde),ux(this%npde),duxx(this%npde,this%npde)
            real(wp), intent(out) :: fval(this%npde)
        end subroutine
        subroutine oo_bndry(this,t,x,u,alpha,beta,gamm)
            import wp, pde_t
            class(pde_t), intent(in) :: this
            real(wp), intent(in) :: t, x
            real(wp), intent(in) :: u(this%npde)
            real(wp), intent(out) :: alpha(this%npde), beta(this%npde), gamm(this%npde)
        end subroutine
    end interface    


    interface pdeone_type
        module procedure :: new_pdeone
    end interface

contains

    function new_pdeone(npde,x,icord,d,f,bndry) result(pde)
        integer, intent(in) :: npde
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: icord
        procedure(d_interface) :: d
        procedure(f_interface) :: f
        procedure(bndry_interface) :: bndry
        type(pdeone_type) :: pde

        pde%npde = npde
        pde%npts = size(x)
        allocate(pde%x,source=x)
        pde%icord = icord

        pde%d => d
        pde%f => f
        pde%bndry => bndry
    end function

    subroutine pdeone(this, t, u, udot)
        class(pdeone_type), intent(in) :: this
        real(wp), intent(in) :: t
        real(wp), intent(inout) :: u(this%npde,this%npts)
        real(wp), intent(out) :: udot(this%npde,this%npts)
!  pdeone is an interface subroutine which uses centered dif-
!  ference approximations to convert one-dimensional systems
!  of partial differential equations into a system of ordinary
!  differential equations, udot = f(t,x,u).  this routine is
!  intended to be used with a robust ode integrator.
!  input..
!   npde = number of partial differential equations.
!   npts = number of spatial grid points.
!      t = current value of time.
!      u = an npde by npts array containing the computed
!          solution at time t.
!  output..
!   udot = an npde by npts array containing the right hand
!          side of the resulting system of ode*s, f(t,x,u),
!          obtained by discretizing the given pde*s.
!  the user must insert a dimension statement of the follow-
!  ing form..
        real(wp) :: dval(this%npde,this%npde,2), ux(this%npde), uavg(this%npde)
        real(wp) :: alpha(this%npde), beta(this%npde), gamma_(this%npde)
        integer :: k, l, i, icord1, itest, ibck, ifwd, ilim
        real(wp) :: dxi, dxil, dxir, dxic, xavgl, xavgr
! the symbols ** denote the actual numerical value of npde
! for the problem being solved.
! common block mesh contains the user specified spatial
! grid points.
! common block coord contains 0,1, or 2 depending on whether
! the problem is in cartesian, cylindrical, or spherical
! coordinates, respectively.
!      common /mesh/ x(1)
!      common /coord/ icord

    associate(icord => this%icord,x => this%x, npts => this%npts, npde => this%npde)

        icord1 = icord + 1
    ! update the left boundary values
        call this%bndry(t, x(1), u, alpha, beta, gamma_, npde)
        itest = 0
        dxi = 1._wp/(x(2) - x(1))
        do k = 1, npde
            if (beta(k) /= 0.0_wp) cycle
            u(k,1) = gamma_(k)/alpha(k)
            itest = itest + 1
        end do
        if (itest.eq.npde) go to 50
        if (itest.eq.0) go to 20
        call this%bndry(t, x(1), u, alpha, beta, gamma_, npde)
    ! evaluate d coefficients at the left boundary
   20   call this%d(t, x(1), u, dval, npde)
    ! form approximations to du/dx at the left boundary
        do 40 k=1,npde
            if (beta(k).ne.0.0) go to 30
            ux(k) = dxi*(u(k,2)-u(k,1))
            go to 40
   30       ux(k) = (gamma_(k)-alpha(k)*u(k,1))/beta(k)
   40   continue
    ! evaluate u-average in the first interval
   50   do k = 1, npde
            uavg(k) = 0.5_wp*(u(k,2) + u(k,1))
        end do
    ! evaluate d coefficients in the first interval
        xavgr = 0.5_wp*(x(2) + x(1))
        call this%d(t,xavgr,uavg,dval(1,1,2),npde)
        dxil = 1.0_wp
        dxir = dxi
        if (icord.eq.0) go to 70
        dxil = x(1)**icord
        dxir = dxir*xavgr**icord
    ! evaluate duxx at the left boundary
   70   if (itest.eq.npde) go to 100
        dxic = real(icord1,wp)/(xavgr**icord1-x(1)**icord1)
        do l = 1, npde
            do k = 1, npde
                dval(k,l,1) = dxic*(dval(k,l,2)*(u(l,2)-u(l,1))*dxir-dval(k,l,1)*ux(l)*dxil)
            end do
        end do
    ! evaluate righthand side of pde*s at the left boundary
        call this%f(t,x(1),u,ux,dval,udot,npde)
    ! set udot = 0 for known left boundary values
  100   do k = 1, npde
            if (beta(k) == 0.0_wp) udot(k,1) = 0.0_wp
        end do
    ! update the right boundary values
        call this%bndry(t,x(npts),u(1,npts),alpha,beta,gamma_,npde)
        itest = 0
        do k = 1, npde
            if (beta(k) /= 0.0_wp) cycle
            u(k,npts) = gamma_(k)/alpha(k)
            itest = itest + 1
        end do
        ibck = 1
        ifwd = 2
        ilim = npts - 1
    ! main loop to form ordinary differential equations at the
    ! grid points
        do 160 i=2,ilim
            k = ibck
            ibck = ifwd
            ifwd = k
            xavgl = xavgr
            xavgr = .5*(x(i+1)+x(i))
            dxi = 1./(x(i+1)-x(i-1))
            dxil = dxir
            dxir = 1./(x(i+1)-x(i))
            if (icord.ne.0) dxir = dxir*xavgr**icord
            dxic = real(icord1,wp)/(xavgr**icord1-xavgl**icord1)
    ! evaluate du/dx and u-average at the i-th grid point and
    ! interval respectively
            do k = 1, npde
                ux(k) = dxi*(u(k,i+1) - u(k,i-1))
                uavg(k) = 0.5_wp*(u(k,i+1) + u(k,i))
            end do
    ! evaluate d coefficients in the i-th interval
            call this%d(t,xavgr,uavg,dval(1,1,ifwd),npde)
    ! evaluate duxx at the i-th grid point
             do 150 l=1,npde
                do 140 k=1,npde
                    dval(k,l,ibck) = (dval(k,l,ifwd)*(u(l,i+1)-u(l,i))*dxir &
                        - dval(k,l,ibck)*(u(l,i)-u(l,i-1))*dxil)*dxic
  140           continue
  150       continue
! evaluate righthand side of pde*s at the i-th grid point
            call this%f(t,x(i),u(1,i),ux,dval(1,1,ibck),udot(1,i),npde)
  160   continue
! finish updating the right boundary if necessary
        if (itest.eq.npde) go to 220
        if (itest.eq.0) go to 170
        call this%bndry(t,x(npts),u(1,npts),alpha,beta,gamma_,npde)
    ! evaluate d coefficients at the right boundary
  170   call this%d(t, x(npts), u(1,npts), dval(1,1,ibck), npde)
    ! form approximations to du/dx at the right boundary
        dxi = 1./(x(npts)-x(ilim))
        do 190 k=1,npde
            if (beta(k).ne.0.0) go to 180
            ux(k) = dxi*(u(k,npts)-u(k,ilim))
            go to 190
  180       ux(k) = (gamma_(k)-alpha(k)*u(k,npts))/beta(k)
  190   continue
        dxil = dxir
        dxir = 1.0_wp
        if (icord /= 0) dxir = x(npts)**icord
        dxic = real(icord1,wp)/(x(npts)**icord1 - xavgr**icord1)
    ! evaluate duxx at the right boundary
        do l = 1, npde
            do k = 1, npde
                dval(k,l,ibck) = dxic*(dval(k,l,ibck)*ux(l)*dxir-dval(k,l,ifwd)*(u(l,npts)-u(l,ilim))*dxil)
            end do
        end do
    ! evaluate righthand side of pde*s at the right boundary
        call this%f(t,x(npts),u(1,npts),ux,dval(1,1,ibck),udot(1,npts),npde)
    ! set udot = 0 for known right boundary values
  220   do k=1,npde
            if (beta(k) == 0.0_wp) udot(k,npts) = 0.0_wp
        end do
    end associate
    end subroutine

subroutine pdeone_oo(this, t, u, udot)
        class(pde_t), intent(in) :: this
        real(wp), intent(in) :: t
        real(wp), intent(inout) :: u(this%npde,this%npts)
        real(wp), intent(out) :: udot(this%npde,this%npts)
!  pdeone is an interface subroutine which uses centered dif-
!  ference approximations to convert one-dimensional systems
!  of partial differential equations into a system of ordinary
!  differential equations, udot = f(t,x,u).  this routine is
!  intended to be used with a robust ode integrator.
!  input..
!   npde = number of partial differential equations.
!   npts = number of spatial grid points.
!      t = current value of time.
!      u = an npde by npts array containing the computed
!          solution at time t.
!  output..
!   udot = an npde by npts array containing the right hand
!          side of the resulting system of ode*s, f(t,x,u),
!          obtained by discretizing the given pde*s.
!  the user must insert a dimension statement of the follow-
!  ing form..
        real(wp) :: dval(this%npde,this%npde,2), ux(this%npde), uavg(this%npde)
        real(wp) :: alpha(this%npde), beta(this%npde), gamma_(this%npde)
        integer :: k, l, i, icord1, itest, ibck, ifwd, ilim
        real(wp) :: dxi, dxil, dxir, dxic, xavgl, xavgr
! the symbols ** denote the actual numerical value of npde
! for the problem being solved.
! common block mesh contains the user specified spatial
! grid points.
! common block coord contains 0,1, or 2 depending on whether
! the problem is in cartesian, cylindrical, or spherical
! coordinates, respectively.
!      common /mesh/ x(1)
!      common /coord/ icord

    associate(icord => this%icord,x => this%xm, npts => this%npts, npde => this%npde)

        icord1 = icord + 1
    ! update the left boundary values
        call this%bndry(t, x(1), u, alpha, beta, gamma_)
        itest = 0
        dxi = 1._wp/(x(2) - x(1))
        do k = 1, npde
            if (beta(k) /= 0.0_wp) cycle
            u(k,1) = gamma_(k)/alpha(k)
            itest = itest + 1
        end do
        if (itest.eq.npde) go to 50
        if (itest.eq.0) go to 20
        call this%bndry(t, x(1), u, alpha, beta, gamma_)
    ! evaluate d coefficients at the left boundary
   20   call this%d(t, x(1), u, dval)
    ! form approximations to du/dx at the left boundary
        do 40 k=1,npde
            if (beta(k).ne.0.0) go to 30
            ux(k) = dxi*(u(k,2)-u(k,1))
            go to 40
   30       ux(k) = (gamma_(k)-alpha(k)*u(k,1))/beta(k)
   40   continue
    ! evaluate u-average in the first interval
   50   do k = 1, npde
            uavg(k) = 0.5_wp*(u(k,2) + u(k,1))
        end do
    ! evaluate d coefficients in the first interval
        xavgr = 0.5_wp*(x(2) + x(1))
        call this%d(t,xavgr,uavg,dval(1,1,2))
        dxil = 1.0_wp
        dxir = dxi
        if (icord.eq.0) go to 70
        dxil = x(1)**icord
        dxir = dxir*xavgr**icord
    ! evaluate duxx at the left boundary
   70   if (itest.eq.npde) go to 100
        dxic = real(icord1,wp)/(xavgr**icord1-x(1)**icord1)
        do l = 1, npde
            do k = 1, npde
                dval(k,l,1) = dxic*(dval(k,l,2)*(u(l,2)-u(l,1))*dxir-dval(k,l,1)*ux(l)*dxil)
            end do
        end do
    ! evaluate righthand side of pde*s at the left boundary
        call this%f(t,x(1),u,ux,dval,udot)
    ! set udot = 0 for known left boundary values
  100   do k = 1, npde
            if (beta(k) == 0.0_wp) udot(k,1) = 0.0_wp
        end do
    ! update the right boundary values
        call this%bndry(t,x(npts),u(1,npts),alpha,beta,gamma_)
        itest = 0
        do k = 1, npde
            if (beta(k) /= 0.0_wp) cycle
            u(k,npts) = gamma_(k)/alpha(k)
            itest = itest + 1
        end do
        ibck = 1
        ifwd = 2
        ilim = npts - 1
    ! main loop to form ordinary differential equations at the
    ! grid points
        do 160 i=2,ilim
            k = ibck
            ibck = ifwd
            ifwd = k
            xavgl = xavgr
            xavgr = .5*(x(i+1)+x(i))
            dxi = 1./(x(i+1)-x(i-1))
            dxil = dxir
            dxir = 1./(x(i+1)-x(i))
            if (icord.ne.0) dxir = dxir*xavgr**icord
            dxic = real(icord1,wp)/(xavgr**icord1-xavgl**icord1)
    ! evaluate du/dx and u-average at the i-th grid point and
    ! interval respectively
            do k = 1, npde
                ux(k) = dxi*(u(k,i+1) - u(k,i-1))
                uavg(k) = 0.5_wp*(u(k,i+1) + u(k,i))
            end do
    ! evaluate d coefficients in the i-th interval
            call this%d(t,xavgr,uavg,dval(1,1,ifwd))
    ! evaluate duxx at the i-th grid point
             do 150 l=1,npde
                do 140 k=1,npde
                    dval(k,l,ibck) = (dval(k,l,ifwd)*(u(l,i+1)-u(l,i))*dxir &
                        - dval(k,l,ibck)*(u(l,i)-u(l,i-1))*dxil)*dxic
  140           continue
  150       continue
! evaluate righthand side of pde*s at the i-th grid point
            call this%f(t,x(i),u(1,i),ux,dval(1,1,ibck),udot(1,i))
  160   continue
! finish updating the right boundary if necessary
        if (itest.eq.npde) go to 220
        if (itest.eq.0) go to 170
        call this%bndry(t,x(npts),u(1,npts),alpha,beta,gamma_)
    ! evaluate d coefficients at the right boundary
  170   call this%d(t, x(npts), u(1,npts), dval(1,1,ibck))
    ! form approximations to du/dx at the right boundary
        dxi = 1./(x(npts)-x(ilim))
        do 190 k=1,npde
            if (beta(k).ne.0.0) go to 180
            ux(k) = dxi*(u(k,npts)-u(k,ilim))
            go to 190
  180       ux(k) = (gamma_(k)-alpha(k)*u(k,npts))/beta(k)
  190   continue
        dxil = dxir
        dxir = 1.0_wp
        if (icord /= 0) dxir = x(npts)**icord
        dxic = real(icord1,wp)/(x(npts)**icord1 - xavgr**icord1)
    ! evaluate duxx at the right boundary
        do l = 1, npde
            do k = 1, npde
                dval(k,l,ibck) = dxic*(dval(k,l,ibck)*ux(l)*dxir-dval(k,l,ifwd)*(u(l,npts)-u(l,ilim))*dxil)
            end do
        end do
    ! evaluate righthand side of pde*s at the right boundary
        call this%f(t,x(npts),u(1,npts),ux,dval(1,1,ibck),udot(1,npts))
    ! set udot = 0 for known right boundary values
  220   do k=1,npde
            if (beta(k) == 0.0_wp) udot(k,npts) = 0.0_wp
        end do
    end associate
    end subroutine

end module

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