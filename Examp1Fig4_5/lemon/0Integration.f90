      module Integrations
      use constants 
      implicit none

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
        real(mcp), parameter :: c1 = 20.D0
        real(mcp), parameter :: c2 = 2.D0**(11.D0/6.D0)
        real(mcp), parameter :: c3 = 2.D0**(11.D0/12.D0) * 8.D0
        real(mcp), parameter :: c4 = 2.D0**(11.D0/12.D0)
 
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
      real(mcp), parameter :: c1 = 20.D0
      real(mcp), parameter :: c2 = 2.D0**(11.D0/6.D0)
      real(mcp), parameter :: c3 = 2.D0**(11.D0/12.D0) * 8.D0
      real(mcp), parameter :: c4 = 2.D0**(11.D0/12.D0)

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

!*******************************************************************************************************

       end module Integrations




