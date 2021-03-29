      module BCS_simulations
      use ModuleForCompScatEsti 
      implicit none 
      integer, parameter :: N_BCS = 100

      type, public, extends(Photon_ForEstimation) :: BCS_photons
      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !procedure, public :: Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi   =>   &
          !                     Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub 
          !procedure, public :: Set_Initial_Values_For_Photon_Parameters    =>    &
          !                     Set_Initial_Values_For_Photon_Parameters_Sub
          procedure, public :: Func_Sigma1_power
          procedure, public :: Func_Sigma1 
          procedure, public :: Func_m_gamma
          procedure, public :: Func_Sigma2_power
          procedure, public :: Func_Sigma2
          procedure, public :: Func_m_gamma2
          procedure, public :: BCS_analytical_formula  =>  BCS_analytical_formula_Sub
          procedure, public :: BCS_IQUV2xytheta   =>   BCS_IQUV2xytheta_Sub
          procedure, public :: BCS_xytheta2IQUV   =>   BCS_xytheta2IQUV_Sub
          procedure, public :: BCS_Get_c1_c2   =>   BCS_Get_c1_c2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Generate_A_Photon   =>   Generate_A_Photon_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p    =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering    =>    &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE   =>    FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Photon_Electron_Scattering   =>   Photon_Electron_Scattering_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type BCS_photons
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      !private :: Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !private :: Set_Initial_Values_For_Photon_Parameters_Sub
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
      private :: Set_InI_Conditions_For_Next_Scattering_Sub
      private :: FIRST_SCATTERING_OF_PHOT_ELCE_Sub
      private :: Photon_Electron_Scattering_Sub
      private :: BCS_analytical_formula_Sub
      private :: BCS_IQUV2xytheta_Sub
      private :: BCS_xytheta2IQUV_Sub
      private :: BCS_Get_c1_c2_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
      contains       
!*******************************************************************************************************
      subroutine BCS_analytical_formula_Sub( this, Theta_e, S_in, filename )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this 
      real(mcp), intent(in) :: Theta_e, S_in(1: 4) 
      character*80, intent(inout) :: filename
      real(mcp) :: tau, T_bb, T_elec, gam, temp_v1, freq_s, J_I, J_Q, J_U, J_V, J_I1, J_I2
      real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS), IQUV(1: 4), xythe(1: 4), cx, cy, &
                   IQUV1(1: 4), xythe1(1: 4), coef_1(1:2, 1: 3), coef_2(1:2, 1: 3), sigma1, sigma2
      real(mcp) :: cos_p1ap2, sin_p1ap2, cos_p1mp2, sin_p1mp2, alp, xt, N_coef, x1, x2
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

        alp = this%alp
        N_coef = ( this%alp - one ) / ( this%gama1**(one -this%alp) - this%gama2**(one -this%alp) )
        do i = 0, N_BCS
            freq_s = 10.D0**( this%vy1 + (this%vy2 - this%vy1) / N_BCS * i )
            gam = dsqrt( freq_s ) / temp_v1

            !sigma1 = this%Func_Sigma1( gam, Theta_e, x0la1000, w0la1000, 362 )
            !sigma2 = this%Func_Sigma2( gam, Theta_e, x0la1000, w0la1000, 362 )
            !sigma1 = this%Func_Sigma1_power( gam, Theta_e, x1000, w1000, 1000 )
            !sigma2 = this%Func_Sigma2_power( gam, Theta_e, x1000, w1000, 1000 )
            !sigma1 = one / gam ** (alp + two) * ( one / (alp + 5.D0) - one /(alp+one) + two/(alp+3.D0) )
            !sigma2 = one / gam ** (alp + two) * ( one / (alp + 5.D0) + one /(alp+one) - two/(alp+3.D0) )
            !sigma1 = gam ** (-6.D0)/8.D0 - gam ** (-6.D0) / 4.D0 + gam**()
            x1 = gam / this%gama2
            if(this%gama1 <= gam)then
                x2 = one
            else
                x2 = gam / this%gama1
            endif
            !write(unit=6, fmt=*)'lls=',gam, this%gama1, this%gama2, x1,x2
            sigma1 = N_coef / gam**(alp + two) * ( ( x2**(alp+5.D0) - x1**(alp+5.D0) ) / (alp + 5.D0) - & 
                ( x2**(alp+1.D0) - x1**(alp+1.D0) ) /(alp+one) + &
                two*( x2**(alp+3.D0) - x1**(alp+3.D0) )/(alp+3.D0) )

            sigma2 = N_coef / gam**(alp + two) * ( ( x2**(alp+5.D0) - x1**(alp+5.D0) ) / (alp + 5.D0) + & 
                ( x2**(alp+1.D0) - x1**(alp+1.D0) ) /(alp+one) - &
                two*( x2**(alp+3.D0) - x1**(alp+3.D0) )/(alp+3.D0) )
 
            J_I = ( sigma1 + 3.D0 * sigma2 ) * gam * freq_s

            J_Q = ( coef_1(1, 1) * sigma1 + coef_2(1, 1) * sigma2 - &
                   (coef_1(2, 1) * sigma1 + coef_2(2, 1) * sigma2) ) * gam * freq_s
  
            J_U = ( coef_1(1, 2) * sigma1 + coef_2(1, 2) * sigma2 - &
                   (coef_1(2, 2) * sigma1 + coef_2(2, 2) * sigma2) ) * gam * freq_s

            J_V = ( coef_1(1, 3) * sigma1 + coef_2(1, 3) * sigma2 - &
                   (coef_1(2, 3) * sigma1 + coef_2(2, 3) * sigma2) ) * gam * freq_s
 
            write(unit=13, fmt=100)J_I , J_Q / J_I, J_U / J_I, J_V / J_I 
            !write(unit=6, fmt=*)'mms= ', coef_1(1, 3), coef_1(2, 3), coef_2(1, 3), coef_2(2, 3)
        enddo
        100 FORMAT(' ', '    ', ES20.10, '   ', ES20.10, '   ', ES20.10, '   ', ES20.10)
        200 FORMAT(' ', '    ', ES20.10) 
        write(unit=6, fmt=*)'lls=', two * ( one - dcos( pi * 85.D0 / 180.D0 ) ), &
                              dlog10( two * ( one - dcos( pi * 85.D0 / 180.D0 ) ) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       ! stop

      end subroutine BCS_analytical_formula_Sub


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
      real(mcp), intent(in) :: IQUV(1: 4)!, IQUV1(1: 4)
      real(mcp), intent(out) :: c1(1: 2, 1: 3), c2(1: 2, 1: 3)
      real(mcp) :: xythe(1: 4), xythe1(1: 4)
      real(mcp) :: cos_p1ap2, sin_p1ap2, cos_p1mp2, sin_p1mp2
      !real(mcp) :: tau, T_bb, T_elec, gam, temp_v1, freq_s, J_I, J_I1, J_I2
      !real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS)
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
      subroutine Generate_A_Photon_Sub( this, Emitter )
!*******************************************************************************************************
      implicit none
      class(BCS_photons) :: this
      TYPE(Photon_Emitter), intent(inout) :: Emitter  
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL Emitter%get_Phot4k_CtrCF_CovCF_BoundReflec() 
      !Emitter%Vector_of_Momentum(1:3) = Emitter%Phot4k_CtrCF(2:4) / Emitter%Phot4k_CtrCF(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%z_tau = this%z_max
      !this%medium_case = Emitter%medium_case
      this%Vector_of_Momentum_ini = Emitter%Vector_of_Momentum
      !this%Vector_of_position_ini = Emitter%Vector_of_position
      this%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%E_ini = DABS( Emitter%Phot4k_CtrCF(1) ) 
      !write(*, *)'f1 = ', this%E_ini
      this%Sigma_a_E_ini = this%sigma_fn( this%E_ini )! * this%n_e1 
      this%ne_times_Sigma_a = this%Sigma_a_E_ini * this%n_e1 
      !write(*, *)'f1333 = ', this%Sigma_KN_E_ini, this%E_ini
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%w_ini = Emitter%w_ini_em 
      this%w_ini0 = Emitter%w_ini_em
      this%Q_sp = zero
      this%U_sp = zero
      this%V_sp = zero
      this%delta_pd = zero
      this%Psi_I = one
      this%Psi_Q = zero
      this%Psi_U = zero
      this%Psi_V = zero
      this%Vector_Stokes4_CF(1) = one
      this%Vector_Stokes4_CF(2) = zero
      this%Vector_Stokes4_CF(3) = zero
      this%Vector_Stokes4_CF(4) = zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      !call this%get_J_emissivity_for_estimation()
      call this%get_J_emissivity_for_estimation_Phiarr()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Generate_A_Photon_Sub

!************************************************************************************
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub( this ) 
!************************************************************************************
      IMPLICIT NONE
      class(BCS_photons) :: this  
   
      this%z_tau = this%Get_scatter_distance_BoundReflec( )    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%w_ini = this%w_ini * this%NormalA  ! * dexp( - this%CROS_absorption )
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub



!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(BCS_photons) :: this
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      real(mcp) :: Sigma_COH, Sigma_INCOH, P_Rayl, P_Comp, P_fluor, xi_1, vLv
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      sphot%delta_pd = this%delta_pd
      sphot%f4_CF    = this%f4_CF 
      sphot%Vector_Stokes4_CF = this%Vector_Stokes4_CF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Photon_Tetrad_In_CF()
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scat_Elec( T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Elec_Tetrad_In_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Phot4k_In_Elec_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%StokesPara_Rotation_Matrix(sphot%Elec_Phot_phi, sphot%Vector_Stokes4_CF)
      sphot%Vector_Stokes4_ECF = sphot%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Making_An_Estimation_One_MC_Component_1( sphot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL sphot%Compton_Scattering_With_Zero_QU_StokesVec()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(BCS_photons) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Phot4k_CtrCF_ini = sphot%Scattered_Phot4k_CF 
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / dabs( this%Phot4k_CtrCF_ini(1) )
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mmsf=', this%Phot4k_CtrCF_ini

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !this%Q_sp    = sphot%Q_sp_scat
      !this%U_sp    = sphot%U_sp_scat
      !this%delta_pd = sphot%delta_pd_scat
      this%Vector_Stokes4_CF = sphot%Vector_Stokes4_ECF_scat
      this%f4_CF    = sphot%f4_scat_CF 
      !write(*, fmt="(' ', A5, 1ES18.7)")'ss1=', &
      !    Vector3D_Inner_Product( this%Phot4k_CtrCF_ini(2: 4), this%f4_CF(2: 4) ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%E_ini = DABS( this%Phot4k_CtrCF_ini(1) ) 
      !write(*, fmt = "(' ', A10, ES20.6, I10)")'f1 = ', this%E_ini, this%scatter_times
      this%Sigma_a_E_ini = this%sigma_fn( this%E_ini )
      this%ne_times_Sigma_a = this%Sigma_a_E_ini * this%n_e1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 

!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(BCS_photons) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      real(mcp) :: Sigma_COH, Sigma_INCOH, P_Rayl, P_Comp, P_fluor, xi_1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      sphot%delta_pd = this%delta_pd
      sphot%f4_CF    = this%f4_CF 
      sphot%Vector_Stokes4_CF = this%Vector_Stokes4_CF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !write(*, fmt="(' ', A5, 1ES18.7)")'ss2=', &
      !    Vector3D_Inner_Product( sphot%Phot4k_CtrCF(2: 4), sphot%f4_CF(2: 4) ) 
      !write(*, fmt = "(' ', A40)")'******************************************************'
      CALL sphot%Set_Photon_f3_Tetrad_In_CF()
      !write(*, fmt = "(' ', A40)")'******************************************************'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Elec_Tetrad_In_CF()
      CALL sphot%StokesPara_Rotation_Matrix(sphot%Elec_Phot_phi, sphot%Vector_Stokes4_CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Phot4k_In_Elec_CF()
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Set_f4_In_Elec_CF()
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Set_Phot_f4_Tetrad_In_Elec_CF()    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      sphot%Vector_Stokes4_ECF = sphot%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !this%Phot4k_In_Elec_CF = sphot%Phot4k_In_Elec_CF
      !this%f4_In_Elec_CF = sphot%f4_In_Elec_CF
      !this%Phot_f4_AxisX = sphot%Phot_f4_AxisX
      !this%Elec_Phot_mu = sphot%Elec_Phot_mu
      !this%Elec_Phot_sin = sphot%Elec_Phot_sin
      !this%Elec_gama = sphot%Elec_gama
      !this%Elec_V = sphot%Elec_V
      !this%Matrix_Of_Tetrad_Of_ElecAxis = sphot%Matrix_Of_Tetrad_Of_ElecAxis
      !this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis
      !this%Matrix_Of_Tetrad_Of_Phot_f4_Axis = sphot%Matrix_Of_Tetrad_Of_Phot_f4_Axis
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Vector_Stokes4_ECF = sphot%Vector_Stokes4_ECF
      this%Phot4k_In_Elec_CF = sphot%Phot4k_In_Elec_CF
      this%Elec_Phot_mu = sphot%Elec_Phot_mu
      this%Elec_Phot_sin = sphot%Elec_Phot_sin
      this%Elec_Phot_phi = sphot%Elec_Phot_phi
      this%Elec_gama = sphot%Elec_gama
      this%Elec_V = sphot%Elec_V

      this%Elec_Phot_sin_In_Elec_CF = sphot%Elec_Phot_sin_In_Elec_CF 
      this%Elec_Phot_mu_In_Elec_CF = sphot%Elec_Phot_mu_In_Elec_CF

      this%Matrix_Of_Tetrad_Of_ElecAxis = sphot%Matrix_Of_Tetrad_Of_ElecAxis
      this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis  
      this%Matrix_Of_Tetrad1_Of_photAxis = sphot%Matrix_Of_Tetrad1_Of_photAxis
      this%Matrix_ECF_2_ECF1 = sphot%Matrix_ECF_2_ECF1
      this%Matrix_ECF1_2_ECF = sphot%Matrix_ECF1_2_ECF
      !CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation() 
      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Compton_Scattering_With_Zero_QU()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withpol_medium2(1)  
      !CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withpol_medium2(2) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !if( sphot%delta_pd /= zero )then 
      !write(*, fmt = "(' ', A40)")'rrrrr******************************************************'
          CALL sphot%Compton_Scattering_With_Polar_StokesVec()
      !write(*, fmt = "(' ', A40)")'rrrrr******************************************************'
      !else   
      !    CALL sphot%Compton_Scattering_With_Zero_QU_StokesVec()
      !endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module BCS_simulations




 
