      module BCS_simulations
      use ModuleForCompScatEsti 
      implicit none 
      integer, parameter :: N_BCS = 100

      type, public, extends(Photon_ForEstimation) :: BCS_photons 
      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Func_Sigma1_power
          procedure, public :: Func_Sigma1 
          procedure, public :: Func_m_gamma
          procedure, public :: Func_Sigma2_power
          procedure, public :: Func_Sigma2
          procedure, public :: Func_m_gamma2
          procedure, public :: BCS_analytical_formula_Powerlaw  =>  &
                               BCS_analytical_formula_Powerlaw_Sub
          procedure, public :: BCS_analytical_formula_HotElectron  =>  &
                               BCS_analytical_formula_HotElectron_Sub
          procedure, public :: BCS_IQUV2xytheta   =>   BCS_IQUV2xytheta_Sub
          procedure, public :: BCS_xytheta2IQUV   =>   BCS_xytheta2IQUV_Sub
          procedure, public :: BCS_Get_c1_c2   =>   BCS_Get_c1_c2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end type BCS_photons
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      private :: BCS_analytical_formula_Powerlaw_Sub 
      private :: BCS_analytical_formula_HotElectron_Sub
      private :: BCS_IQUV2xytheta_Sub
      private :: BCS_xytheta2IQUV_Sub
      private :: BCS_Get_c1_c2_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
      contains       
!*******************************************************************************************************
      subroutine BCS_analytical_formula_Powerlaw_Sub( this, S_in, mu_obs, filename )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: S_in(1: 4), mu_obs
      character*80, intent(inout) :: filename
      real(mcp) :: tau, T_bb, T_elec, Emin, temp_v1, freq_s, J_I, J_Q, J_U, J_V, J_I1, J_I2
      real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS), IQUV(1: 4), xythe(1: 4), cx, cy, &
                   IQUV1(1: 4), xythe1(1: 4), coef_1(1:2, 1: 3), coef_2(1:2, 1: 3), sigma1, sigma2
      real(mcp) :: cos_p1ap2, sin_p1ap2, cos_p1mp2, sin_p1mp2, alp, xt, x1, x2
      integer :: i


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IQUV(1) = S_in(1)
        IQUV(2) = S_in(2)
        IQUV(3) = S_in(3)
        IQUV(4) = S_in(4)
        call this%BCS_IQUV2xytheta( IQUV, xythe )

        IQUV1(1) = one
        IQUV1(2) = zero
        IQUV1(3) = zero
        IQUV1(4) = one
        call this%BCS_IQUV2xytheta( IQUV1, xythe1 )

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        coef_1(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        coef_2(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2

        write(unit=6, fmt=*)coef_1, coef_2, xythe, xythe1
 
        call this%BCS_Get_c1_c2( IQUV, coef_1, coef_2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
        temp_v1 = dsqrt( two * ( one - mu_obs ) )
        open(unit=13, file = filename, status="replace")  

        alp = this%alp 
        do i = 0, N_BCS
            freq_s = 10.D0**( this%vy1 + (this%vy2 - this%vy1) / N_BCS * i )
            Emin = dsqrt( freq_s ) / temp_v1

            !sigma1 = this%Func_Sigma1( Emin, Theta_e, x0la1000, w0la1000, 362 )
            !sigma2 = this%Func_Sigma2( Emin, Theta_e, x0la1000, w0la1000, 362 )
            !sigma1 = this%Func_Sigma1_power( Emin, Theta_e, x1000, w1000, 1000 )
            !sigma2 = this%Func_Sigma2_power( Emin, Theta_e, x1000, w1000, 1000 ) 
            !sigma1 = Emin ** (-6.D0)/8.D0 - Emin ** (-6.D0) / 4.D0 + Emin**()
            x1 = Emin / this%gama2
            if(this%gama1 <= Emin)then
                x2 = one
            else
                x2 = Emin / this%gama1
            endif
            !write(unit=6, fmt=*)'lls=',Emin, this%gama1, this%gama2, x1,x2
!~~~~~ Here the expressions of sigma1 and sigma2 are different from the Equation (144) of our paper,
!~~~~~ Since in our paper, gama2 is infinity, then x1 = 0.D0, the below expressions of 
!~~~~~ sigma1 and sigma2 will equal to the Equation (144) of our paper. Here we choose gama2 is a 
!~~~~~ finite number, we get the following expressions for sigma1 and sigma2.
            sigma1 = this%N_coef / Emin**(alp + two) * ( ( x2**(alp+5.D0) - x1**(alp+5.D0) ) /& 
                 (alp + 5.D0) - ( x2**(alp+1.D0) - x1**(alp+1.D0) ) /(alp+one) + &
                two*( x2**(alp+3.D0) - x1**(alp+3.D0) )/(alp+3.D0) )

            sigma2 = this%N_coef / Emin**(alp + two) * ( ( x2**(alp+5.D0) - x1**(alp+5.D0) ) /& 
                 (alp + 5.D0) + ( x2**(alp+1.D0) - x1**(alp+1.D0) ) /(alp+one) - &
                two*( x2**(alp+3.D0) - x1**(alp+3.D0) )/(alp+3.D0) )
 
            J_I = ( sigma1 + 3.D0 * sigma2 ) * Emin * freq_s

            J_Q = ( coef_1(1, 1) * sigma1 + coef_2(1, 1) * sigma2 - &
                   (coef_1(2, 1) * sigma1 + coef_2(2, 1) * sigma2) ) * Emin * freq_s
  
            J_U = ( coef_1(1, 2) * sigma1 + coef_2(1, 2) * sigma2 - &
                   (coef_1(2, 2) * sigma1 + coef_2(2, 2) * sigma2) ) * Emin * freq_s

            J_V = ( coef_1(1, 3) * sigma1 + coef_2(1, 3) * sigma2 - &
                   (coef_1(2, 3) * sigma1 + coef_2(2, 3) * sigma2) ) * Emin * freq_s
 
            write(unit=13, fmt=100)J_I , J_Q / J_I, J_U / J_I, J_V / J_I 
            !write(unit=6, fmt=*)'mms= ', coef_1(1, 3), coef_1(2, 3), coef_2(1, 3), coef_2(2, 3)
        enddo
        100 FORMAT(' ', '    ', ES20.10, '   ', ES20.10, '   ', ES20.10, '   ', ES20.10)
        200 FORMAT(' ', '    ', ES20.10) 
        write(unit=6, fmt=*)'lls=', two * ( one - mu_obs ), dlog10( two * ( one - mu_obs ) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      end subroutine BCS_analytical_formula_Powerlaw_Sub


!*******************************************************************************************************
      subroutine BCS_analytical_formula_HotElectron_Sub( this, Theta_e, S_in, filename )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: Theta_e, S_in(1: 4)
      character*80, intent(inout) :: filename
      real(mcp) :: tau, T_bb, T_elec, gam, temp_v1, freq_s, J_I, J_Q, J_U, J_V, J_I1, J_I2
      real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS), IQUV(1: 4), xythe(1: 4), cx, cy, &
                   IQUV1(1: 4), xythe1(1: 4), coef_1(1:2, 1: 3), coef_2(1:2, 1: 3), sigma1, sigma2
      real(mcp) :: cos_p1ap2, sin_p1ap2, cos_p1mp2, sin_p1mp2
      integer :: i


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IQUV(1) = S_in(1)
        IQUV(2) = S_in(2)
        IQUV(3) = S_in(3)
        IQUV(4) = S_in(4)
        call this%BCS_IQUV2xytheta( IQUV, xythe )

        IQUV1(1) = one
        IQUV1(2) = zero
        IQUV1(3) = zero
        IQUV1(4) = one
        call this%BCS_IQUV2xytheta( IQUV1, xythe1 )

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        coef_1(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        coef_2(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2

        write(unit=6, fmt=*)coef_1, coef_2, xythe, xythe1
 
        call this%BCS_Get_c1_c2( IQUV, coef_1, coef_2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call Set_xi_wi_all()
        !Theta_e = 100.D0
        temp_v1 = dsqrt( two * ( one - dcos( pi * 85.D0 / 180.D0 ) ) )
        open(unit=13, file = filename, status="replace") 
        !open(unit=14, file = './spectrum/bcs2.txt', status="replace") 

        do i = 0, N_BCS
            freq_s = 10.D0**( 1.D0 + (7.D0 - 1.D0) / N_BCS * i )
            gam = dsqrt( freq_s ) / temp_v1

            sigma1 = this%Func_Sigma1( gam, Theta_e, x0la1000, w0la1000, 362 )
            sigma2 = this%Func_Sigma2( gam, Theta_e, x0la1000, w0la1000, 362 )
            !write(unit=6, fmt=*)'lls=',gam, sigma1, sigma2
 
            J_I = ( sigma1 + 3.D0 * sigma2 ) * gam * freq_s

            J_Q = ( coef_1(1, 1) * sigma1 + coef_2(1, 1) * sigma2 - &
                   (coef_1(2, 1) * sigma1 + coef_2(2, 1) * sigma2) ) * gam * freq_s
  
            J_U = ( coef_1(1, 2) * sigma1 + coef_2(1, 2) * sigma2 - &
                   (coef_1(2, 2) * sigma1 + coef_2(2, 2) * sigma2) ) * gam * freq_s

            J_V = ( coef_1(1, 3) * sigma1 + coef_2(1, 3) * sigma2 - &
                   (coef_1(2, 3) * sigma1 + coef_2(2, 3) * sigma2) ) * gam * freq_s
 
            write(unit=13, fmt=100)J_I, J_Q / J_I, J_U / J_I, J_V / J_I 
  
        enddo
        100 FORMAT(' ', '    ', ES20.10, '   ', ES20.10, '   ', ES20.10, '   ', ES20.10)
        200 FORMAT(' ', '    ', ES20.10) 
        write(unit=6, fmt=*)'lls=', two * ( one - dcos( pi * 85.D0 / 180.D0 ) ), &
                              dlog10( two * ( one - dcos( pi * 85.D0 / 180.D0 ) ) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ! stop

      end subroutine BCS_analytical_formula_HotElectron_Sub


!*******************************************************************************************************
      subroutine BCS_IQUV2xytheta_Sub( this, IQUV, xythe )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: IQUV(1: 4)
      real(mcp), intent(out) :: xythe(1: 4)
      real(mcp) :: tau, T_bb, T_elec, gam, temp_v1, freq_s, J_I, J_I1, J_I2
      real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS)
      integer :: i
   
      xythe(1) = dsqrt( ( IQUV(1) + IQUV(2) ) / two )
      xythe(2) = dsqrt( ( IQUV(1) - IQUV(2) ) / two )
      if( IQUV(3) /= zero .and. IQUV(4) /= zero )then
          xythe(3) = IQUV(3) / dsqrt( IQUV(3)**2 + IQUV(4)**2 )
          xythe(4) = IQUV(4) / dsqrt( IQUV(3)**2 + IQUV(4)**2 )
      else
          xythe(3) = one
          xythe(4) = one
      endif

      end subroutine BCS_IQUV2xytheta_Sub


!*******************************************************************************************************
      subroutine BCS_xytheta2IQUV_Sub( this, IQUV, xythe )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(out) :: IQUV(1: 4)
      real(mcp), intent(in) :: xythe(1: 4)
      real(mcp) :: tau, T_bb, T_elec, gam, temp_v1, freq_s, J_I, J_I1, J_I2
      real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS)
      integer :: i

      IQUV(1) = xythe(1)**2 + xythe(2)**2 
      IQUV(2) = xythe(1)**2 - xythe(2)**2 
      IQUV(3) = two * xythe(1) * xythe(2) * xythe(3)
      IQUV(4) = two * xythe(1) * xythe(2) * xythe(4) 

      end subroutine BCS_xytheta2IQUV_Sub



!*******************************************************************************************************
      subroutine BCS_Get_c1_c2_Sub( this, IQUV, c1, c2 )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: IQUV(1: 4) 
      real(mcp), intent(out) :: c1(1: 2, 1: 3), c2(1: 2, 1: 3)
      real(mcp) :: xythe(1: 4), xythe1(1: 4)
      real(mcp) :: cos_p1ap2, sin_p1ap2, cos_p1mp2, sin_p1mp2 
      integer :: i

        call this%BCS_IQUV2xytheta( IQUV, xythe ) 
  
        xythe1(1) = one
        xythe1(2) = zero
        xythe1(3) = one
        xythe1(4) = zero

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        c1(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        c2(1, 1) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2

        xythe1(1) = zero
        xythe1(2) = one
        xythe1(3) = one
        xythe1(4) = zero

        c1(2, 1) = ( xythe(2)*xythe1(2) )**2
        c2(2, 1) = ( xythe(2)*xythe1(2) )**2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        xythe1(1) = one / dsqrt(two)
        xythe1(2) = one / dsqrt(two)
        xythe1(3) = one
        xythe1(4) = zero

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        c1(1, 2) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        c2(1, 2) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2


        xythe1(1) = one / dsqrt(two)
        xythe1(2) = one / dsqrt(two)
        xythe1(3) = - one
        xythe1(4) = zero

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        c1(2, 2) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        c2(2, 2) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xythe1(1) = one / dsqrt(two)
        xythe1(2) = one / dsqrt(two)
        xythe1(3) = zero
        xythe1(4) = one

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        c1(1, 3) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        c2(1, 3) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2

        xythe1(1) = one / dsqrt(two)
        xythe1(2) = one / dsqrt(two)
        xythe1(3) = zero
        xythe1(4) = - one

        cos_p1ap2 = xythe(3) * xythe1(3) - xythe(4) * xythe1(4)
        sin_p1ap2 = xythe(4) * xythe1(3) + xythe(3) * xythe1(4)
        cos_p1mp2 = xythe(3) * xythe1(3) + xythe(4) * xythe1(4)
        sin_p1mp2 = xythe(4) * xythe1(3) - xythe(3) * xythe1(4)
        
        c1(2, 3) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1mp2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1mp2 )**2

        c2(2, 3) = ( xythe(1)*xythe1(1) - xythe(2)*xythe1(2)*cos_p1ap2 )**2 + &
                 ( xythe(2)*xythe1(2)*sin_p1ap2 )**2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      end subroutine BCS_Get_c1_c2_Sub



!*******************************************************************************************************
      real(mcp) function Func_Sigma1( this, gam, Theta_e, t0la, w0la, nla )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      integer, intent(in) :: nla
      real(mcp), intent(in) :: gam, Theta_e, t0la(0: nla-1), w0la(0: nla-1)
      integer :: i, j
      real(mcp) :: w0, Coef_A, gama0, Ints 

      w0 = zero  
      Coef_A = gam * dexp( - gam / Theta_e ) / bessk( 2, one / Theta_e )

      do j = 0, nla - 1
          gama0 = t0la(j) * Theta_e + gam
          if(gama0 < one)gama0=one
          Ints = this%Func_m_gamma( gama0, gam ) 
          w0 = w0 + Ints * w0la(j) 
          !write(unit=6, fmt=*)'lls=',gama0, gam
      enddo 
      Func_Sigma1 = w0 * Coef_A
         ! write(unit=6, fmt=*)'lls=1', Func_Sigma1 
   

      end function Func_Sigma1



!*******************************************************************************************************
      real(mcp) function Func_m_gamma( this, ga, gam )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: ga, gam
      integer :: i
      real(mcp) :: x
 
      !if(ga < one)ga=one
      x = gam / ga
  
      Func_m_gamma = dsqrt( one - one / ga**2 ) / ( ga**2 ) * ( x**2 - one / x / x + two )
      !!Notice that we have used a Gaussian integral method to calculate a integral with form as:
      !! f(x) * exp( - x ), and x in (1, infinity), so, in the final expression, exp(-x) will not
      !! appear. this may confuse me latter.

      end function Func_m_gamma


!*******************************************************************************************************
      real(mcp) function Func_Sigma1_power( this, gam, Theta_e, t0la, w0la, nla )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      integer, intent(in) :: nla
      real(mcp), intent(in) :: gam, Theta_e, t0la(0: nla-1), w0la(0: nla-1)
      integer :: i, j
      real(mcp) :: w0, Coef_A, x_int, Ints, gama

      w0 = zero  
      !Coef_A = gam * dexp( - gam / Theta_e ) / bessk( 2, one / Theta_e )

      do j = 0, nla - 1
          x_int = ( t0la(j) + one ) / two
          !if(gama0 < one)gama0=one
          gama = gam / x_int
          Ints = gama ** (-3.D0 - two) * ( x_int**2 - one / x_int / x_int + two )
          w0 = w0 + Ints * w0la(j) 
          !write(unit=6, fmt=*)'lls=',gama0, gam
      enddo 
      Func_Sigma1_power = w0 ! * Coef_A
         ! write(unit=6, fmt=*)'lls=1', Func_Sigma1 
   

      end function Func_Sigma1_power

 
!*******************************************************************************************************
      real(mcp) function Func_Sigma2( this, gam, Theta_e, t0la, w0la, nla )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      integer, intent(in) :: nla
      real(mcp), intent(in) :: gam, Theta_e, t0la(0: nla-1), w0la(0: nla-1)
      integer :: i, j
      real(mcp) :: w0, Coef_A, gama0, Ints 

      w0 = zero  
      Coef_A = gam * dexp( - gam / Theta_e ) / bessk( 2, one / Theta_e )

      do j = 0, nla - 1
          gama0 = t0la(j) * Theta_e + gam
          if(gama0 < one)gama0=one
          Ints = this%Func_m_gamma2( gama0, gam ) 
          w0 = w0 + Ints * w0la(j) 
      enddo 
  
      Func_Sigma2 = w0 * Coef_A

      end function Func_Sigma2



!*******************************************************************************************************
      real(mcp) function Func_m_gamma2( this, ga, gam )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: ga, gam
      integer :: i
      real(mcp) :: x

      !if(ga < one)ga=one
      x = gam / ga
  
      Func_m_gamma2 = dsqrt( one - one / ga**2 ) / ( ga**2 ) * ( one / x - x )**2 

      end function Func_m_gamma2



 
!*******************************************************************************************************
      real(mcp) function Func_Sigma2_power( this, gam, Theta_e, t0la, w0la, nla )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      integer, intent(in) :: nla
      real(mcp), intent(in) :: gam, Theta_e, t0la(0: nla-1), w0la(0: nla-1)
      integer :: i, j
      real(mcp) :: w0, Coef_A, gama0, Ints, x_int, gama

      w0 = zero  
      !Coef_A = gam * dexp( - gam / Theta_e ) / bessk( 2, one / Theta_e )

      do j = 0, nla - 1
          x_int = ( t0la(j) + one ) / two 
          gama = gam / x_int
          Ints = gama ** (-3.D0 - two) * ( one / x_int - x_int )**2 
          w0 = w0 + Ints * w0la(j)   
      enddo 
  
      Func_Sigma2_power = w0 !* Coef_A

      end function Func_Sigma2_power
 

!*******************************************************************************************************
 
      end module BCS_simulations




 
