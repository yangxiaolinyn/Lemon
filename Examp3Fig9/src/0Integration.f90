      module Integrations
      use constants 
      implicit none
      real(mcp), parameter :: c1 = 20.D0
      real(mcp), parameter :: c2 = 2.D0**(11.D0/6.D0)
      real(mcp), parameter :: c3 = 2.D0**(11.D0/12.D0) * 8.D0
      real(mcp), parameter :: c4 = 2.D0**(11.D0/12.D0)
      integer(kind=8), parameter :: n_all = 2**11 - 2 
      real(mcp) :: x5(0: 4), w5(0: 4)
      real(mcp) :: x7(0: 6), w7(0: 6), x10(0: 9), w10(0: 9), x20(0: 19), w20(0: 19)
      real(mcp) :: x50(0: 49), w50(0: 49), x100(0: 99), w100(0: 99), x200(0: 199), w200(0: 199)
      real(mcp) :: x2(0: 1), w2(0: 1), x4(0: 3), w4(0: 3), x15(0: 14), w15(0: 14)
      real(mcp) :: x8(0: 7), w8(0: 7), x16(0: 15), w16(0: 15)
      real(mcp) :: x32(0: 31), w32(0: 31), x64(0: 63), w64(0: 63)
      real(mcp) :: x128(0: 127), w128(0: 127), x256(0: 255), w256(0: 255)
      real(mcp) :: x512(0: 511), w512(0: 511), x1024(0: 1023), w1024(0: 1023)
      real(mcp) :: x500(0: 499), w500(0: 499), x1000(0: 999), w1000(0: 999)
      real(mcp) :: xall(0: n_all), wall(0: n_all)
      real(mcp) :: x0la1000(0: 361), w0la1000(0: 361)
      real(mcp) :: x1la100(0: 99), w1la100(0: 99)
      real(mcp) :: x2la100(0: 99), w2la100(0: 99)
      real(mcp) :: x10la100(0: 99), w10la100(0: 99)

      contains
!*******************************************************************************************************
        real(mcp) function Integrals(a1, b1) 
!*******************************************************************************************************
        implicit none 
        real(mcp) :: ant
        real(mcp), intent(in) :: a1, b1
        integer :: i, j, N=300000 
        real(mcp) :: eps, d, temp, delta, a, b
        !real(mcp), external :: ff
!*******************************************************************************************************
        
        eps = 1.D-5
        d = 1.D-5    
     
        ant = zero
        !a = one
        delta = (b1 - a1) / 1000.D0 !1.D-2
        a = a1
        b = one + delta
        do   
            temp = fpts_gama(a, b, eps, d, f_temp)
            !write(*,*)ant,temp,b
            ant = ant + temp
            if ( temp <= 1.D-15 ) exit
            a = b
            b = b + delta
            if( b > b1 )exit
        enddo 
        Integrals = ant
        return
        end function Integrals

!******************************************************************************************************* 
        function f_temp(y) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: y
        real(mcp) :: eps, d, theta
        integer, parameter :: n = 2 
 
        ant = ( y**2 + c2 + c4 * two * y ) * y**3 * dexp( - y )
        !write(*,*)ant
        return
        end function f_temp  

!*******************************************************************************************************
        function fpts_gama(a, b, eps, d, f_temp) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: a, b, eps, d
        real(mcp) :: h, f0, f1, t0
        real(mcp) :: integal
        real(mcp), external :: f_temp

        h = b - a
        integal = zero 
        f0 = f_temp( a ) 
        f1 = f_temp( b )
        t0 = h*(f0 + f1)/two 
            
        call ppp1(a, b, h, f0, f1, t0, eps, d, integal, f_temp)
        ant = integal
        return
        end function fpts_gama

!*******************************************************************************************************
        recursive subroutine ppp1(x0, x1, h, f0, f1, t0, eps, d, integal, f_temp)
!*******************************************************************************************************
        implicit none
        real(mcp), intent(in) :: x0, x1, h, f0, f1, t0, eps, d
        real(mcp), intent(inout) :: integal
        real(mcp) :: x, f, t1, t2, p, g, eps1
        real(mcp), external :: f_temp

        x = x0 + h/two 
        f = f_temp(x)
        t1 = h*(f0 + f)/four 
        t2 = h*(f + f1)/four
        p = dabs(t0 - (t1 + t2))  
        if ( (p < eps) .or. (h/two < d) ) then
            integal = integal + t1 + t2 
            return
        else 
            g = h/two
            eps1 = eps/1.4D0
            call ppp1(x0, x, g, f0, f, t1, eps1, d, integal, f_temp)
            call ppp1(x, x1, g, f, f1, t2, eps1, d, integal, f_temp)
        endif
        end subroutine ppp1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration2( a, b, n )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      real(mcp), intent(in) :: a, b
      real(mcp) :: Del_y, ints, y, Del_ints, Del_ints_0, Del_ints_n
      integer, intent(in) :: n
      integer :: i 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          ints = zero  
          i = 0 
          Del_y = ( b - a ) / (n + 1) 
          y = a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_ints = ( y + c2 * y**(one/three) + c4 * two * &
                      y**(two/three) ) * dexp( - y**(one/three) ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              ints = ints + Del_ints    
              If ( i == 0 ) then 
                  Del_ints_0 = Del_ints
              Endif  
              If ( i > n ) then
                  Del_ints_n = Del_ints
                  exit
              Endif  
              i = i + 1  
              y = Del_y + y
              !write(*,*)'sd', y, Del_y, ints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          Integration2 = ( ints - Del_ints_0 / two - Del_ints_n / two )*Del_y
      end function Integration2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function f_Legendre_Gauss( x )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      real(mcp), intent(in) :: x
      
      f_Legendre_Gauss = ( x + c2 * x**(one/three) + c4 * two * &
                      x**(two/three) ) * dexp( - x**(one/three) ) / x
      end function f_Legendre_Gauss

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration_Legendre_Gauss( a, b, eps )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      real(mcp), intent(in) :: a, b, eps
      real(mcp) :: h, s, p, ep, aa, bb, w, x, g
      real(mcp), parameter :: t(0: 4)=(/- 0.9061798459D0, - 0.5384693101D0, 0.D0, &
                 0.5384693101D0, 0.9061798459D0/), &
                 l(0: 4) = (/ 0.2369268851D0, 0.4786286705D0, 0.5688888889D0, &
                 0.4786286705D0, 0.2369268851D0 /)
      integer :: i, j, m 
  
      m = 1
      h = b - a
      s = dabs(0.001D0 * h)
      p = 1.D35
      ep = eps + one
      Do
          if( ep < eps .and. dabs(h) < s )exit 
          g = zero
          do i = 1, m
              aa = a + ( i - one ) * h
              bb = a + i * h
              w = zero
              do j = 0, 4
                  x = ( ( bb - aa ) * t(j) + ( bb + aa ) ) / two
                  w = w + f_Legendre_Gauss( x ) * l(j)
              enddo
              g = g + w
          enddo
          g = g * h / two
          ep = dabs( g - p ) / ( one + dabs( g ) )
          p = g
          m = m + 1
          h = ( b - a ) / m
      ENDDO
      Integration_Legendre_Gauss = g
      end function Integration_Legendre_Gauss

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration_Legendre_Gauss_xw( a, b, eps, t, w, n )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      integer, intent(in) :: n
      real(mcp), intent(in) :: a, b, eps, t(0: n-1), w(0: n-1)
      real(mcp) :: h, s, p, ep, aa, bb, w1, x, g 
      integer :: i, j, m 
  
      m = 1
      h = ( b - a ) / m
      !s = dabs(0.001D0 * h)
      !p = 1.D35
      !ep = eps + one
      !Do
          !if( ep < eps .and. m >=4 )exit !.and. dabs(h) < s )exit 
          g = zero
          do i = 1, m
              aa = a + ( i - one ) * h
              bb = a + i * h
              w1 = zero
              do j = 0, n - 1
                  x = ( ( bb - aa ) * t(j) + ( bb + aa ) ) / two
                  w1 = w1 + f_Legendre_Gauss( x ) * w(j)
              enddo
              g = g + w1
          enddo
          g = g * h / two
          !ep = dabs( g - p ) / ( one + dabs( g ) )
          !write(*,*)'wt = ', ep, m, g
          !p = g
          !m = m + 1
          !h = ( b - a ) / m
      !ENDDO
      Integration_Legendre_Gauss_xw = g
      end function Integration_Legendre_Gauss_xw

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration_Chebyshev( a, b, eps )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      real(mcp), intent(in) :: a, b, eps
      real(mcp) :: h, d, p, ep, aa, bb, w, x, g
      real(mcp), parameter :: t(0: 4)=(/- 0.8324975D0, - 0.3745414D0, 0.D0, &
                 0.3745414D0, 0.8324975D0/)
      integer :: i, j, m 
  
      m = 1
      h = b - a
      d = dabs(0.001D0 * h)
      p = 1.D35
      ep = eps + one
      Do
          g = zero
          if( ep < eps .and. dabs(h) < d )exit 
          do i = 1, m
              aa = a + ( i - one ) * h
              bb = a + i * h
              w = zero
              do j = 0, 4
                  x = ( ( bb - aa ) * t(j) + ( bb + aa ) ) / two
                  w = w + f_Legendre_Gauss( x )
              enddo
              g = g + w
          enddo
          g = g * h / two
          ep = dabs( g - p ) / ( one + dabs( g ) )
          p = g
          m = m + 1
          h = ( b - a ) / m
      ENDDO
      Integration_Chebyshev = g
      end function Integration_Chebyshev

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      subroutine gauleg_x_w( x1, x2, x, w, n )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      real(mcp), intent(in) :: x1, x2
      integer, intent(in) :: n
      real(mcp), intent(inout) :: x(0:n-1), w(0:n-1)
      real(mcp), parameter :: eps = 1.D-14
      real(mcp) :: z1, z, xm, xl, pp, p3, p2, p1
      integer :: m, i, j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      m = (n + 1) / 2
      xm = 0.5D0 * ( x2 + x1 )
      xl = 0.5D0 * ( x2 - x1 )
      do i = 0, m-1
          z = dcos( pi * (i + 0.75D0) / (n + 0.5D0) )
          do
              p1 = one
              p2 = zero
              do j = 0, n-1
                  p3 = p2
                  p2 = p1
                  p1 = ( ( two * j + one ) * z * p2 - j * p3 ) / ( j + 1 )
              enddo
              pp = n * ( z * p1 - p2 ) / ( z**2 - one )
              z1 = z
              z = z1 - p1 / pp
              if( dabs(z-z1) < eps )exit
          enddo
          x(i) = xm - xl * z
          x(n-1-i) = xm + xl * z
          w(i) = two * xl / ( ( one - z**2 ) * pp**2 )
          w(n-1-i) = w(i)
      enddo
      end subroutine gauleg_x_w

!*******************************************************************
      subroutine Set_xi_wi_all()
!*******************************************************************
      implicit none
      real(mcp), allocatable :: x_temp(:), w_temp(:)
      integer :: i, m, i_end, i_start

      call gauleg_x_w( -one, one, x10,   w10,   10 )
      call gauleg_x_w( -one, one, x5,    w5,    5 )
      call gauleg_x_w( -one, one, x7,    w7,    7 )
      call gauleg_x_w( -one, one, x20,   w20,   20 )
      call gauleg_x_w( -one, one, x50,   w50,   50 )
      call gauleg_x_w( -one, one, x100,  w100,  100 )
      call gauleg_x_w( -one, one, x200,  w200,  200 )
      call gauleg_x_w( -one, one, x500,  w500,  500 )
      call gauleg_x_w( -one, one, x1000, w1000, 1000 )
 
      call gauleg_x_w( -one, one, x2, w2, 2 )
      call gauleg_x_w( -one, one, x4, w4, 4 )
      call gauleg_x_w( -one, one, x8, w8, 8 )
      call gauleg_x_w( -one, one, x16, w16, 16 )
      call gauleg_x_w( -one, one, x32, w32, 32 )
      call gauleg_x_w( -one, one, x64, w64, 64 )
      call gauleg_x_w( -one, one, x128, w128, 128 )
      call gauleg_x_w( -one, one, x256, w256, 256 )
      call gauleg_x_w( -one, one, x512, w512, 512 )
      call gauleg_x_w( -one, one, x1024, w1024, 1024 ) 
      call gaulag( x10la100, w10la100, 100, ten )
      call gaulag( x2la100, w2la100, 100, two )
      call gaulag( x1la100, w1la100, 100, one )
      call gaulag( x0la1000, w0la1000, 362, zero ) 
      do i = 1, 10
          m = 2**i
          allocate(x_temp(0: m-1))
          allocate(w_temp(0: m-1))
          call gauleg_x_w( -one, one, x_temp, w_temp, m )
          i_end = 2 * ( m - 1 ) - 1
          i_start = i_end - m + 1
          xall(i_start: i_end) = x_temp
          wall(i_start: i_end) = w_temp
          deallocate(x_temp)
          deallocate(w_temp)
      enddo
      end subroutine Set_xi_wi_all
!*******************************************************************************************************

!*******************************************************************************************************
      subroutine gaulag( x, w, n, alf )
!*******************************************************************************************************
      implicit none
      integer, intent(in) :: n
      real(mcp), intent(out) :: x(0: n - 1), w(0: n - 1)
      real(mcp), intent(in) :: alf
      real(mcp), parameter :: EPS = 1.0D-14
      real(mcp) :: ai, p1, p2, p3, pp, z, z1
      integer, parameter :: MAXIT = 10
      integer :: i, its, j
!*******************************************************************************************************

      do i = 0, n - 1
          if(i == 0)then    !Initial guess for the smallest root.
              z = (1.D0 + alf ) * (3.D0 + 0.92D0 * alf ) / ( 1.D0 + 2.4D0 * n + 1.8D0 * alf )
          else if(i == 1)then    !Initial guess for the second root.
              z = z + ( 15.D0 + 6.25D0 * alf ) / ( 1.D0+0.9D0 * alf + 2.5D0 * n )
          else         !Initial guess for the other roots.
              ai = i - 1
              z = z + ( ( 1.D0 + 2.55D0 * ai ) / ( 1.9D0 * ai ) + 1.26D0 * ai * alf / &
                      ( 1.D0 + 3.5D0 * ai ) ) * ( z - x(i-2) ) / ( 1.D0 + 0.3D0 * alf )
          endif
          do its = 0, MAXIT - 1
              p1 = 1.D0
              p2 = 0.D0
              do j = 0, n-1     ! Loop up the recurrence relation to get the
                  p3 = p2       ! Laguerre polynomial evaluated at z.
                  p2 = p1
                  p1 = ( ( 2 * j + 1 + alf - z ) * p2 - ( j + alf ) * p3 ) / ( j + 1 )
              enddo
              pp = ( n * p1 - ( n + alf ) * p2 ) / z 
              z1 = z 
              z = z1 - p1 / pp         !Newton’s formula.
              if( dabs( z - z1 ) <= EPS )exit
          enddo
          if(its >= MAXIT)write(*, '(a)')"too many iterations in gaulag"
          x(i) = z           !Store the root and the weight.
          w(i) = - dexp( gammln( alf + n ) - gammln( dfloat(n) ) ) / ( pp * n * p2 )
          !write(*, *)'ss=qq', z, w(i), gammln(dfloat(n)), pp, p2, n
      enddo
      end subroutine gaulag
!*******************************************************************************************************

!*******************************************************************************************************
      real(mcp) function gammln( xx )
!*******************************************************************************************************
!     Returns the value lnŒ .xx/ for xx > 0.
      implicit none
      real(mcp), intent(in) :: xx
      Integer :: j 
      real(mcp) :: x, tmp, y, ser
      real(mcp), parameter :: cof(0: 13) = (/ 57.1562356658629235D0, -59.5979603554754912D0, &
      14.1360979747417471D0, -0.491913816097620199D0, .339946499848118887D-4, &
      .465236289270485756D-4, -.983744753048795646D-4, .158088703224912494D-3, & 
      -.210264441724104883D-3, .217439618115212643D-3, -.164318106536763890D-3, &
      .844182239838527433D-4, -.261908384015814087D-4, .368991826595316234D-5 /)
      if( xx <= 0 )write(*, '(a)')"bad arg in gammln"
      y = xx
      x = xx
      tmp = x + 5.24218750000000000D0            !Rational 671/128.D0
      tmp = ( x + 0.5D0 ) * log( tmp ) - tmp
      ser = 0.999999999999997092D0
      do j=0, 13
          y = y + 1.D0
          ser = ser + cof(j) / y
      enddo
      gammln = tmp + dlog( 2.5066282746310005D0 * ser / x )
      return
      end function gammln

      end module Integrations




