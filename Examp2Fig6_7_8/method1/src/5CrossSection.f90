      module CrossSection
      use constants 
      implicit none

      contains
!*******************************************************************************************************
        function sigma_a(T_e, E) result(ant)
!*******************************************************************************************************
        implicit none 
        real(mcp) :: ant
        integer :: i, j, N=300000
        real(mcp), intent(in) :: T_e, E
        real(mcp) :: epsi, mu, delta_mu, gama, gama1, gama2, delta_gama,&
               v_speed, theta, sum_integal_mu, sum_integal_gama, eps, d, temp, a, b,&
               delta
        !real(mcp), external :: ff
!*******************************************************************************************************
        
        eps = 1.D-5
        d = 1.D-5
        theta = T_e/mec2
        delta_gama = (10.D5-one)/dfloat(N)
        delta_mu = two/dfloat(N)
        sum_integal_gama = 0.D0
        !do i=0,N
       !     gama = one + i*delta_gama 
        !    sum_integal_mu = fpts(-one, one, eps, d, ff, gama, T_e, E)
            !write(*,*)i,  sum_integal_mu
  
        !    sum_integal_gama = sum_integal_gama + sum_integal_mu*N_e(gama, theta)
        !enddo
        ant = zero
        a = one
        delta = 1.D-2
        b = one + delta
        do   
            temp = fpts_gama(a, b, eps, d, f_temp, T_e, E)
            !write(*,*)ant,temp,b
            ant = ant + temp
            if ( temp <= 1.D-15 ) exit
            a = b
            b = b + delta
        enddo 
        ant = ant*(three*sigma_T/four)!/barn
        return
        end function sigma_a

!******************************************************************************************************* 
        function ff(mu, gama, E) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: mu, gama, E
        real(mcp) :: epsi, v_speed, temp
  
        v_speed = dsqrt(one - one/gama/gama) 
        epsi = (two*E/mec2)*gama*(one - mu*v_speed)
        if (epsi > 1.D-4) then 
            temp = ( (one - four/epsi - eight/epsi/epsi)*dlog(one + epsi)&
                  + half + eight/epsi - half/(one + epsi)**two )/epsi!*(three*sigma_T/four)
        else
            temp = (one - epsi) * four / three
        endif
        ant = (one - mu*v_speed)*temp
        return
        end function ff

!******************************************************************************************************* 
        function ff1(xx, gama, T_e, E) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: xx, gama, T_e, E 
  
        ant = one/(one + xx**2*25.D0)
        return
        end function ff1

!******************************************************************************************************* 
        function f_temp(gama, T_e, E) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: gama, T_e, E 
        real(mcp) :: eps, d, theta
        integer, parameter :: n = 2

        theta = T_e/mec2
        eps = 1.D-7
        d = 1.D-7
        ant = fpts_mu(-one, one, eps, d, ff, gama, E)/(two*theta*bessk(n,one/theta))&
                      *dsqrt(one-one/gama/gama)*gama**2*dexp(-gama/theta) 
        !write(*,*)ant
        return
        end function f_temp

!*******************************************************************************************************
        function fpts_mu(a, b, eps, d, ff, gama, E) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: a, b, eps, d, gama, E
        real(mcp) :: h, f0, f1, t0
        real(mcp) :: integal
        real(mcp), external :: ff

        h = b - a
        integal = zero
        f0 = ff(a, gama, E)
        f1 = ff(b, gama, E)
        t0 = h*(f0 + f1)/two 
            
        call ppp(a, b, h, f0, f1, t0, eps, d, integal, ff, gama, E)
        ant = integal
        return
        end function fpts_mu


!*******************************************************************************************************
        recursive subroutine ppp(x0, x1, h, f0, f1, t0, eps, d, integal, ff, gama, E)
!*******************************************************************************************************
        implicit none
        real(mcp), intent(in) :: x0, x1, h, f0, f1, t0, eps, d, gama, E
        real(mcp), intent(inout) :: integal
        real(mcp) :: x, f, t1, t2, p, g, eps1
        real(mcp), external :: ff

        x = x0 + h/two
        f = ff(x, gama, E)
        t1 = h*(f0 + f)/four 
        t2 = h*(f + f1)/four
        p = dabs(t0 - (t1 + t2))  
        if ( (p < eps) .or. (h/two < d) ) then
            integal = integal + t1 + t2
            return
        else 
            g = h/two
            eps1 = eps/1.4D0
            call ppp(x0, x, g, f0, f, t1, eps1, d, integal, ff, gama, E)
            call ppp(x, x1, g, f, f1, t2, eps1, d, integal, ff, gama, E)
        endif
        end subroutine ppp

!*******************************************************************************************************
        function fpts_gama(a, b, eps, d, f_temp, T_e, E) result(ant)
!*******************************************************************************************************
        implicit none
        real(mcp) :: ant
        real(mcp), intent(in) :: a, b, eps, d, T_e, E
        real(mcp) :: h, f0, f1, t0
        real(mcp) :: integal
        real(mcp), external :: f_temp

        h = b - a
        integal = zero 
        f0 = f_temp(a, T_e, E) 
        f1 = f_temp(b, T_e, E)
        t0 = h*(f0 + f1)/two 
            
        call ppp1(a, b, h, f0, f1, t0, eps, d, integal, f_temp, T_e, E)
        ant = integal
        return
        end function fpts_gama

!*******************************************************************************************************
        recursive subroutine ppp1(x0, x1, h, f0, f1, t0, eps, d, integal, f_temp, T_e, E)
!*******************************************************************************************************
        implicit none
        real(mcp), intent(in) :: x0, x1, h, f0, f1, t0, eps, d, T_e, E
        real(mcp), intent(inout) :: integal
        real(mcp) :: x, f, t1, t2, p, g, eps1
        real(mcp), external :: f_temp

        x = x0 + h/two 
        f = f_temp(x, T_e, E)
        t1 = h*(f0 + f)/four 
        t2 = h*(f + f1)/four
        p = dabs(t0 - (t1 + t2))  
        if ( (p < eps) .or. (h/two < d) ) then
            integal = integal + t1 + t2 
            return
        else 
            g = h/two
            eps1 = eps/1.4D0
            call ppp1(x0, x, g, f0, f, t1, eps1, d, integal, f_temp, T_e, E)
            call ppp1(x, x1, g, f, f1, t2, eps1, d, integal, f_temp, T_e, E)
        endif
        end subroutine ppp1
 

!*******************************************************************************************************
        DOUBLE PRECISION FUNCTION bessi0(x)
!*******************************************************************************************************
        double precision, intent(in) :: x
        !Returns the modified Bessel function I0 (x) for any real x.
        double precision :: ax
        double precision :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,&
                q8,q9,y   ! Accumulate polynomials in double precision.
        SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
        DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,&
                0.2659732d0,0.360768d-1,0.45813d-2/
        DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
                0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,&
                0.2635537d-1,-0.1647633d-1,0.392377d-2/
        if (abs(x).lt.3.75) then
            y=(x/3.75)**2
            bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
        else
            ax=abs(x)
            y=3.75/ax
            bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4 &
                +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        endif
        return
        END FUNCTION bessi0
!*******************************************************************************************************
 
!*******************************************************************************************************
        double precision FUNCTION bessk0(x)
!*******************************************************************************************************
        double precision, intent(in) :: x
        !USES bessi0
        !Returns the modified Bessel function K0 (x) for positive real x. 
        double precision p1,p2,p3,p4,p5,p6,p7,q1,&
               q2,q3,q4,q5,q6,q7,y    ! Accumulate polynomials in double precision.
        SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
        DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,&
        0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
        DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,&
        -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
        if (x.le.2.0) then   !    Polynomial fit.
            y=x*x/4.0
            bessk0=(-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+&
                       y*(p4+y*(p5+y*(p6+y*p7))))))
        else
            y=(2.0/x)
            bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+&
                       y*(q4+y*(q5+y*(q6+y*q7))))))
        endif
        return
        END FUNCTION bessk0

!*******************************************************************************************************
        double precision FUNCTION bessi1(x)
!*******************************************************************************************************
        double precision, intent(in) :: x
        !Returns the modified Bessel function I1 (x) for any real x.
        double precision ax
        double precision p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,&
                q8,q9,y  !Accumulate polynomials in double precision.
        SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
        DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,&
                0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
        DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,&
                -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,&
                -0.2895312d-1,0.1787654d-1,-0.420059d-2/
        if (abs(x) < 3.75) then  ! Polynomial fit.
            y=(x/3.75)**2
            bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
        else
            ax=abs(x)
            y=3.75/ax
            bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+&
                  y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
            if(x.lt.0.)bessi1=-bessi1
        endif
        return
        END FUNCTION bessi1

!*******************************************************************************************************
        double precision FUNCTION bessk1(x)
!*******************************************************************************************************
        double precision, intent(in) :: x
        !USES bessi1
        !Returns the modified Bessel function K1 (x) for positive real x. 
        double precision p1,p2,p3,p4,p5,p6,p7,q1,&
                q2,q3,q4,q5,q6,q7,y ! Accumulate polynomials in double precision.
        SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
        DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,&
                -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
        DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,&
                0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
        if (x <= 2.D0) then !Polynomial fit.
            y=x*x/4.0
            bessk1=(log(x/2.0)*bessi1(x))+(1.0/x)*(p1+y*(p2+&
                    y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
        else
            y=2.0/x
            bessk1=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+&
                 y*(q4+y*(q5+y*(q6+y*q7))))))
        endif
        return
        END FUNCTION bessk1

 
!*******************************************************************************************************
        double precision FUNCTION bessk(n,x)
!*******************************************************************************************************
        integer, intent(in) :: n
        double precision, intent(in) :: x
        !USES bessk0,bessk1
        !Returns the modified Bessel function Kn (x) for positive x and n ≥ 2.
        integer j
        double precision bk,bkm,bkp,tox!
        if (n < 2) pause 'bad argument n in bessk'
        tox=2.0/x
        bkm=bessk0(x)
        bk=bessk1(x)
        do  j=1,n-1
            bkp=bkm+j*tox*bk
            bkm=bk
            bk=bkp
        enddo 
        bessk=bk
        return
        END FUNCTION bessk
 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      end module CrossSection
!********************************************************************************
