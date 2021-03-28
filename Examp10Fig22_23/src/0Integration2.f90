    Module AdaptedIntegration
    use constants
    use RandUtils 
    !use BLcoordinate 
    Implicit none
    real(mcp), parameter :: EPS = 1.D-20
    real(mcp), parameter :: alpha = dsqrt( 2.D0 / 3.D0 ) 
    real(mcp), parameter :: beta = 1.D0 / dsqrt( 5.D0 )
    real(mcp), parameter :: x1 = 0.942882415695480D0
    real(mcp), parameter :: x2 = 0.641853342345781D0
    real(mcp), parameter :: x3 = 0.236383199662150D0
    real(mcp), parameter :: x(1: 12) = (/ 0.D0, -x1, -alpha, -x2, -beta, -x3, &
                                           0.D0, x3, beta, x2, alpha, x1 /)

    type, public :: AdaptInteg
        real(mcp) :: TOL, toler
        logical :: terminate = .true., out_of_tolerance = .false.
    contains
        procedure, public :: Integrate
        procedure, public :: adaptlob
     
    end type AdaptInteg
 
    contains

!*************************************************************************************
    real(mcp) function Integrate( this, func, a, b )
!*************************************************************************************
    Implicit none
    class(AdaptInteg) :: this
    real(mcp), intent(inout) :: a, b
    real(mcp), external :: func
    real(mcp) :: m, h, i1, i2, fa, fb, is, erri1, erri2, r, y(0: 12), toler
    integer :: i
!*************************************************************************************

    if( this%TOL < 10.0 * EPS )this%TOL = 10.D0 * EPS
    m = 0.5D0 * ( a + b )
    h = 0.5D0 * ( b - a )
    fa = func( a )
    fb = func( b )
    y(0) = fa
    y(12) = fb
    do i = 1, 11
        y( i ) = func( m + x(i) * h )
    enddo 
    i2 = ( h / 6.D0 ) * ( y(0) + y(12) + 5.D0 * ( y(4) + y(8) ) )   !4-point Gauss-Lobatto formula.
    i1 = ( h / 1470.D0 ) * ( 77.D0 * ( y(0) + y(12) ) + 432.D0 * ( y(2) + y(10) ) + &
    625.D0 * ( y(4) + y(8) ) + 672.D0 * y(6) )    !7-point Kronrod extension.

    is = h * ( 0.0158271919734802D0 * ( y(0) + y(12) ) + 0.0942738402188500D0 * &
    ( y(1) + y(11) ) + 0.155071987336585D0 * ( y(2) + y(10) ) + &
    0.188821573960182D0 * ( y(3) + y(9) ) + 0.199773405226859D0 * &
    ( y(4) + y(8) ) + 0.224926465333340D0 * ( y(5) + y(7) ) +  &
    0.242611071901408D0 * y(6) )   !13-point Kronrod extension.

    erri1 = dabs( i1 - is )
    erri2 = dabs( i2 - is )
    if( erri2 /= 0.D0 )then
        r = erri1 / erri2
    else
        r = 1.D0
    endif 
    if( r > 0.D0 .and. r < 1.D0 )then
        this%toler = this%TOL / r
    else 
        this%toler = this%TOL
    endif
    if( is == zero ) is = b - a
    is = dabs( is )
    
    Integrate = this%adaptlob( func, a, b, fa, fb, is )
    end function Integrate

!*************************************************************************************
    recursive real(mcp) function adaptlob( this, func, a, b, fa, fb, is)
!*************************************************************************************
    Implicit none
    class(AdaptInteg) :: this
    real(mcp), intent(in) :: a, b, fa, fb, is
    real(mcp), external :: func
    real(mcp) :: m, h, mll, ml, mr, mrr, fmll, fml, fm, fmrr, fmr, i1, i2
!*************************************************************************************

    m = 0.5D0 * ( a + b )
    h = 0.5D0 * ( b - a )
    mll = m - alpha * h
    ml = m - beta * h
    mr = m + beta * h
    mrr = m + alpha * h
    fmll = func( mll )
    fml = func( ml )
    fm = func( m )
    fmr = func( mr )
    fmrr = func( mrr )
    i2 = h / 6.D0 * ( fa + fb + 5.0D0 * ( fml + fmr ) )   !4-point Gauss-Lobatto formula.
    i1 = h / 1470.D0 * ( 77.D0 * ( fa + fb ) + 432.D0 * ( fmll + fmrr ) &
                      + 625.D0 * ( fml + fmr ) + 672.D0 * fm )  ! 7-point Kronrod extension.

    if( dabs( i1 - i2 ) <= this%toler * is .and. mll <= a .and. b <= mrr )then
        if( (mll <= a .and. b <= mrr ) .and. this%terminate )then
            this%out_of_tolerance = .true.
            this%terminate = .false.
        else
            adaptlob = i1
            return
        endif
    else
        adaptlob = this%adaptlob( func, a, mll, fa, fmll, is ) + &
                   this%adaptlob( func, mll, ml, fmll, fml, is) + &
                   this%adaptlob( func, ml, m, fml, fm, is) + &
                   this%adaptlob( func, m, mr, fm, fmr, is) + &
                   this%adaptlob( func, mr, mrr, fmr, fmrr, is) + &
                   this%adaptlob( func, mrr, b, fmrr, fb, is)
        return 
    endif
    end function adaptlob
!*************************************************************************************

    End Module AdaptedIntegration












 
