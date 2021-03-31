      module Photons_FlatSP
      use ScatDistance_FlatSP  
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type, public, extends(Photon_With_ScatDistance_FlatSP) :: Photon_FlatSP
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
          real(mcp) :: v_L_v_i(1: vL_sc_up)  
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: nu_obs  
          real(mcp) :: PolarArrayd(0: Num_PolDeg)
          real(mcp) :: PolarArrayI(0: Num_PolDeg)
          real(mcp) :: PolarArrayQ(0: Num_PolDeg)
          real(mcp) :: PolarArrayU(0: Num_PolDeg)
          real(mcp) :: PolarArrayV(0: Num_PolDeg) 
          real(mcp) :: R_Matrix(1: 4, 1: 4) 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: d_theta, d_phi, d_tau
          logical :: first_time_recording

      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Calc_Phot_Informations_At_Observor_FLST_IQ  =>   &
                                Calc_Phot_Informations_At_Observor_FLST_IQ_Sub 
          procedure, public :: Matrix_Multiplied
          procedure, public :: Matrix_Multiplication
          procedure, public :: Set_initial_parameter_values 
          procedure, public :: Initiate_The_Scatter_Sequence
          procedure, public :: Initiate_The_Scatter_Sequence2
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p 
          procedure, public :: Determine_Next_Scattering_Site
          procedure, public :: Photon_Electron_Scattering
          procedure, public :: Get_Analytical_Solutions_Of_RTE_QUV
          procedure, public :: Get_Analytical_Solutions_Of_RTE_IQ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_FlatSP
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Calc_Phot_Informations_At_Observor_FLST_IQ_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Matrix_Multiplied( this, I, Q, U, V, F_I, F_Q, F_U, F_V )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
      real(mcp), intent(in) :: I, Q, U, V
      real(mcp), intent(out) :: F_I, F_Q, F_U, F_V 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      F_I = ( this%alpI - this%alp ) * I + this%alpQ * Q + this%alpU * U + this%alpV * V
      F_Q = this%alpQ * I + ( this%alpI - this%alp ) * Q + this%rhoV * U - this%rhoU * V
      F_U = this%alpU * I - this%rhoV * Q + ( this%alpI - this%alp ) * U + this%rhoQ * V
      F_V = this%alpV * I + this%rhoU * Q - this%rhoQ * U + ( this%alpI - this%alp ) * V
  
      RETURN
      END SUBROUTINE Matrix_Multiplied


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      SUBROUTINE Matrix_Multiplication( this, IQUV, FIQUV )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
      real(mcp), intent(in) :: IQUV(1: 4)
      real(mcp), intent(out) :: FIQUV(1: 4)
      integer :: i, j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      FIQUV = zero
      do i = 1, 4
          do j = 1, 4
              FIQUV(i) = FIQUV(i) + this%R_Matrix(i, j) * IQUV(j)
          enddo
      enddo
  
      RETURN
      END SUBROUTINE Matrix_Multiplication



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values( this, tau, &
              jIQUV, alpIQUV, rhoQUV, IQUV0, alpha )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      class(Photon_FlatSP) :: this  
      REAL(mcp), INTENT(IN) :: tau, jIQUV(1: 4), alpIQUV(1: 4), &
                               rhoQUV(1: 3), IQUV0(1: 4), alpha 
      REAL(mcp) :: E_low, E_up
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set Initial conditions for the Photon                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          this%effect_number = 0 
          this%j_I = jIQUV(1)
          this%j_Q = jIQUV(2)
          this%j_U = jIQUV(3)
          this%j_V = jIQUV(4)

          this%alpI = alpIQUV(1)
          this%alpQ = alpIQUV(2)
          this%alpU = alpIQUV(3)
          this%alpV = alpIQUV(4)

          this%rhoQ = rhoQUV(1)
          this%rhoU = rhoQUV(2)
          this%rhoV = rhoQUV(3)

          this%I0 = IQUV0(1)
          this%Q0 = IQUV0(2)
          this%U0 = IQUV0(3)
          this%V0 = IQUV0(4)
  
          this%alp = alpha 
 
          this%rho  = dsqrt(this%rhoQ**2 + this%rhoU**2 + this%rhoV**2)
 
          !write(*, *)'sss = ', this%j_Q = -10.71710  
      !this%alpI = zero
      !this%rho  = dsqrt(this%rhoQ**2 + this%rhoU**2 + this%rhoV**2)

      !this%j_Q = ( -this%rho**2 / this%rhoQ * 0.95D0 - this%j_V * this%rhoV ) / this%rhoQ

      this%alpp = this%alp + this%alpQ
      this%alpm = this%alp - this%alpQ
      this%s_max  = tau

      this%R_Matrix(1, 1) = this%alpI - this%alp 
      this%R_Matrix(1, 2) = this%alpQ
      this%R_Matrix(1, 3) = this%alpU
      this%R_Matrix(1, 4) = this%alpV

      this%R_Matrix(2, 1) = this%alpQ
      this%R_Matrix(2, 2) = this%alpI - this%alp
      this%R_Matrix(2, 3) = this%alpV
      this%R_Matrix(2, 4) = this%alpU

      this%R_Matrix(3, 1) = this%alpU
      this%R_Matrix(3, 2) = this%alpV
      this%R_Matrix(3, 3) = this%alpI - this%alp
      this%R_Matrix(3, 4) = this%alpQ

      this%R_Matrix(4, 1) = this%alpV
      this%R_Matrix(4, 2) = this%alpU
      this%R_Matrix(4, 3) = this%alpQ
      this%R_Matrix(4, 4) = this%alpI - this%alp
 
      RETURN
      END SUBROUTINE Set_initial_parameter_values
 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Initiate_The_Scatter_Sequence2( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon_FlatSP) :: this   
      real(mcp) :: temp_value
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
      this%s_var = this%s_max * ranmar() 
      temp_value = dexp( - 50.D0 * this%s_var * this%alp ) * this%alp * 49.D0

      this%psi_I = this%j_I + this%I0 * temp_value 
      this%psi_Q = this%j_Q + this%Q0 * temp_value 
      this%psi_U = this%j_U + this%U0 * temp_value 
      this%psi_V = this%j_V + this%V0 * temp_value 

      !this%psi_I = this%j_I
      !this%psi_Q = this%j_Q
      !this%psi_U = this%j_U
      !this%psi_V = this%j_V
      !write(*,*)'5551=',  this%x_ini, this%y_ini, this%z_ini 
      Call this%Calc_Phot_Informations_At_Observor_FLST_IQ( ) 
  
      RETURN
      END SUBROUTINE Initiate_The_Scatter_Sequence2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Initiate_The_Scatter_Sequence( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon_FlatSP) :: this  
      real(mcp) :: temp_value, F_I, F_Q, F_U, F_V
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
      this%s_var = this%s_max * ranmar()

      temp_value = dexp( - this%s_var * this%alp )  
 

      CALL this%Matrix_Multiplied( this%I0, this%Q0, this%U0, &
                                  this%V0, F_I, F_Q, F_U, F_V )

      this%psi_I = this%j_I - F_I * temp_value 
      this%psi_Q = this%j_Q - F_Q * temp_value 
      this%psi_U = this%j_U - F_U * temp_value 
      this%psi_V = this%j_V - F_V * temp_value 
      !write(*,*)'ss1=',  F_I, F_Q, F_U, F_V  
  
      Call this%Calc_Phot_Informations_At_Observor_FLST_IQ( ) 
  
      RETURN
      END SUBROUTINE Initiate_The_Scatter_Sequence


!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p( this ) 
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
      integer cases
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      this%s_var = this%Get_scatter_distance_IQ( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      this%psi_I = this%psi_I * this%NormalA
      this%psi_Q = this%psi_Q * this%NormalA
      this%psi_U = this%psi_U * this%NormalA
      this%psi_V = this%psi_V * this%NormalA 
      !write(unit = *, fmt = *)'s5 = ', this%s_var, this%NormalA 
      !stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p
 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this   
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      this%s_var = this%Get_scatter_distance_IQ( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       this%psi_I = this%psi_I * this%NormalA
       this%psi_Q = this%psi_Q * this%NormalA
       this%psi_U = this%psi_U * this%NormalA
       this%psi_V = this%psi_V * this%NormalA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering( this )
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this   
      real(mcp) :: F_I, F_Q, F_U, F_V
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   

      CALL this%Matrix_Multiplied( this%psi_I, this%psi_Q, this%psi_U, &
                                   this%psi_V, F_I, F_Q, F_U, F_V )

      this%psi_I = - F_I / this%alp
      this%psi_Q = - F_Q / this%alp
      this%psi_U = - F_U / this%alp 
      this%psi_V = - F_V / this%alp

      Call this%Calc_Phot_Informations_At_Observor_FLST_IQ( ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Get_Analytical_Solutions_Of_RTE_QUV( this, tau, n, &
                                         AnalyticResult_Fig27_right )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
      real(mcp), intent(in) :: tau
      integer, intent(in) :: n
      character*80, intent(in) :: AnalyticResult_Fig27_right
      real(mcp) :: ds, s_i, I0, Q0, U0, V0, rhoS, rhoJ
      real(mcp) :: alpQ, alpU, alpV, jI, jQ, jU, jV, alpp, alpm
      real(mcp) :: I_s(0: n), Q_s(0: n), U_s(0: n), V_s(0: n)
      real(mcp) :: sg
      integer :: i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      ds = tau / (n * one)

      I_s = zero
      Q_s = zero
      U_s = zero
      V_s = zero
      s_i = zero 

      rhoS = this%Q0  * this%rhoQ + this%U0  * this%rhoU + this%V0  * this%rhoV
      rhoJ = this%j_Q * this%rhoQ + this%j_U * this%rhoU + this%j_V * this%rhoV
      do i = 0, n
            Q_s(i) = this%rhoQ / this%rho**2 * s_i * rhoJ + dcos( this%rho*s_i ) * this%Q0 + &
           ( this%j_Q * this%rho**2 - this%rhoQ * rhoJ ) * dsin(this%rho * s_i) / this%rho**3 + &
         ( this%rhoU * this%j_V - this%rhoV * this%j_U ) / this%rho**2 * ( one - dcos( this%rho*s_i ) ) & 
            + two * this%rhoQ * rhoS * dsin(this%rho * s_i/two)**2 / this%rho**2 + &
            (this%rhoU * this%V0 - this%rhoV * this%U0) * dsin(this%rho * s_i) / this%rho &
              - this%Q0 * dexp(-s_i * this%alp)
        

            U_s(i) = this%rhoU / this%rho**2 * s_i * rhoJ + dcos( this%rho*s_i ) * this%U0 + &
               ( this%j_U * this%rho**2 - this%rhoU * rhoJ ) * dsin(this%rho * s_i) / this%rho**3 + &
         ( this%rhoV * this%j_Q - this%rhoQ * this%j_V ) / this%rho**2 * ( one - dcos( this%rho*s_i ) ) & 
            + two * this%rhoU * rhoS * dsin(this%rho * s_i/two)**2 / this%rho**2 + &
            (this%rhoV * this%Q0 - this%rhoQ * this%V0) * dsin(this%rho * s_i) / this%rho &
              - this%U0 * dexp(-s_i * this%alp)


            V_s(i) = this%rhoV / this%rho**2 * s_i * rhoJ + dcos( this%rho*s_i ) * this%V0 + &
               ( this%j_V * this%rho**2 - this%rhoV * rhoJ ) * dsin(this%rho * s_i) / this%rho**3 + &
         ( this%rhoQ * this%j_U - this%rhoU * this%j_Q ) / this%rho**2 * ( one - dcos( this%rho*s_i ) ) & 
            + two * this%rhoV * rhoS * dsin(this%rho * s_i/two)**2 / this%rho**2 + &
            (this%rhoQ * this%U0 - this%rhoU * this%Q0) * dsin(this%rho * s_i) / this%rho &
              - this%V0 * dexp(-s_i * this%alp)

          s_i = s_i + ds  
          !write(*, *)'ssdd = ', s_i, Q_s(i), U_s(i), V_s(i)
      enddo 

      open(unit=10, file = AnalyticResult_Fig27_right, status="replace")      
      do i = 0, 1000
          write(unit = 10, fmt = 100)Q_s(i), U_s(i), V_s(i)
      enddo
      100  FORMAT(' ', 1F22.12, '   ', 1F22.12, '   ', 1F22.12) 
      close(unit=10)  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Get_Analytical_Solutions_Of_RTE_QUV
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Get_Analytical_Solutions_Of_RTE_IQ( this, tau, n, &
                                         AnalyticResult_Fig27_left )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
      real(mcp), intent(in) :: tau
      integer, intent(in) :: n
      character*80, intent(in) :: AnalyticResult_Fig27_left 
      real(mcp) :: ds, s_i, I0, Q0
      real(mcp) :: alpI, alpQ, jI, jQ, alpp, alpm
      real(mcp) :: I_s(0: n), Q_s(0: n)
      real(mcp) :: sg
      integer :: i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      alpI = this%alpI
      alpQ = this%alpQ
      I0 = this%I0
      Q0 = this%Q0
      jI = this%j_I
      jQ = this%j_Q
      alpp = alpI + alpQ
      alpm = alpI - alpQ
 
      ds = tau / (n * one)

      I_s = zero
      Q_s = zero
      s_i = zero
      do i = 0, n
          I_s(i) = ( (jI + jQ) / alpp * ( one - dexp( - alpp * s_i ) ) + &
                     (jI - jQ) / alpm * ( one - dexp( - alpm * s_i ) ) + &
                     (I0 + Q0) * dexp( - alpp * s_i ) + &
                     (I0 - Q0) * dexp( - alpm * s_i ) ) / two
                   
          Q_s(i) = ( (jI + jQ) / alpp * ( one - dexp( - alpp * s_i ) ) - &
                     (jI - jQ) / alpm * ( one - dexp( - alpm * s_i ) ) + &
                     (I0 + Q0) * dexp( - alpp * s_i ) - &
                     (I0 - Q0) * dexp( - alpm * s_i ) ) / two
          s_i = s_i + ds   
      enddo 

      open(unit=10, file = AnalyticResult_Fig27_left, status="replace")      
      do i = 0, 1000
          write(unit = 10, fmt = 100)I_s(i), Q_s(i)
      enddo
      100  FORMAT(' ', 1F22.12, '   ', 1F22.12) 
      close(unit=10)  
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Get_Analytical_Solutions_Of_RTE_IQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_FLST_IQ_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv, p_out1
      integer :: i_obs 
  
      do i_obs = 0, N_s
          p_out1 = this%alp * ( this%s_max / N_s * i_obs - this%s_var )
          if( p_out1 < zero )cycle
          vLv = dexp( - p_out1 )
          this%I_i( i_obs ) = this%I_i( i_obs ) + vLv * this%psi_I
          this%Q_i( i_obs ) = this%Q_i( i_obs ) + vLv * this%psi_Q 
          this%U_i( i_obs ) = this%U_i( i_obs ) + vLv * this%psi_U 
          this%V_i( i_obs ) = this%V_i( i_obs ) + vLv * this%psi_V 
      enddo
      this%effect_number = this%effect_number + 1
  
      return
      end subroutine Calc_Phot_Informations_At_Observor_FLST_IQ_Sub
 
 
   
!*******************************************************************************************************
 
      end module Photons_FlatSP




 
