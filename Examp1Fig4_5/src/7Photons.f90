      module Photons 
      use PhotonEmitter
      !use ScatDistance
      implicit none 
      integer, parameter :: N_sigma = 1000 

      type, public, extends(Photon_Emitter) :: Photon
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(0:500) = zero
          real(mcp) :: v_L_v_i_ET(1:500, 1:1000)
          real(mcp) :: frequency_v
          real(mcp) :: Optical_Depth_scatter   
          real(mcp) :: logE_low
          real(mcp) :: logE_up
          real(mcp) :: dindexE  
          real(mcp) :: NormalA 
          integer(kind=8) :: effect_number
          integer(kind=8) :: scatter_times

      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
       procedure, public :: Calc_Phot_Informations_At_Observor_2zones => &
                            Calc_Phot_Informations_At_Observor_2zones_Sub 
       procedure, public :: Set_initial_parameter_values =>  &
                            Set_initial_parameter_values_Sub
       procedure, public :: Emitter_A_Photon   =>   Emitter_A_Photon_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Calc_Phot_Informations_At_Observor_2zones_Sub
      private :: Set_initial_parameter_values_Sub
      private :: Emitter_A_Photon_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon) :: this
      !REAL(mcp), INTENT(IN) :: T_e
      !TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter
      !TYPE(Photon), INTENT(INOUT) :: Phot 
      REAL(mcp) :: E_low, E_up
  

      this%ln_nu1 = 8.D0
      this%ln_nu2 = 15.D0
      this%T_e = 100.D0 * mec2
      !this%Theta_e = 100.D0 * mec2
      this%R_out = one

      E_low = 1.D-5
      E_up = 1.D1 
      this%n_e = 1.D20 
      this%logE_low = DLOG10(E_low)
      this%logE_up = DLOG10(E_up)
      this%effect_number = 0

      RETURN
      END SUBROUTINE Set_initial_parameter_values_Sub


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photon) :: this
      !TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r_ini = ranmar()**(one/three) 
      this%mucos_ini = one - two * ranmar()
      this%musin_ini = dsqrt( one - this%mucos**2 )
      this%phi_ini = twopi * ranmar()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%x_ini = this%r * this%musin * dcos( this%phi )
      this%y_ini = this%r * this%musin * dsin( this%phi )
      this%z_ini = this%r * this%mucos
      this%Vector_of_position_ini(1) = this%x_ini
      this%Vector_of_position_ini(2) = this%y_ini
      this%Vector_of_position_ini(3) = this%z_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%j_nu_theta_emissity_sampling( )  
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF(2:4) / this%Phot4k_CtrCF(1)  

      RETURN
      END SUBROUTINE Emitter_A_Photon_Sub
 
 
 
 
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      !integer, intent(in) :: cases   
      real(mcp) :: index_i1, del_x1, r_times_p, p_out1
      integer :: k, N_low1
  
      del_x1 = (this%ln_nu2 - this%ln_nu1) / 500.D0
      !N_low1 = - floor( this%ln_nu1 / del_x1 )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
      index_i1 = dlog10( this%frequency_v / 1.D8 ) 
      k = floor( index_i1 / del_x1 ) + 1 ! + N_low1 + 1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - r_times_p + dsqrt( r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * dexp( - this%Optical_Depth_absorption * p_out1 ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
     
      return
      end subroutine Calc_Phot_Informations_At_Observor_2zones_Sub
!*******************************************************************************************************
!******************************************************************************************************* 
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Photons





