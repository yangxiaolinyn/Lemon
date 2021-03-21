      module PhotonModule
      use PhotonsEstimation
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photons_Esti) :: Photons
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
      contains   
          procedure, public :: Set_initial_parameter_values => Set_initial_parameter_values_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p   =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub   
          procedure, public :: Photon_Electron_Scattering    =>  &
                               Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photons

      private :: Set_initial_parameter_values_Sub 
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub   
      private :: Photon_Electron_Scattering_Sub
     

      contains   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values_Sub( this, tau, mu0, phi0 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      class(Photons) :: this   
      !TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter 
      REAL(mcp), INTENT(INOUT) :: tau, mu0, phi0
      integer :: i 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set Initial conditions for the Photons                    !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%tau_max = tau  
      this%cos_Theta0 = mu0
      this%sin_Theta0 = dsqrt( one - mu0**2 ) 

      this%phi_ini0 = phi0
      this%cosphi_ini0 = dcos( phi0 )
      this%sinphi_ini0 = dsin( phi0 )
      this%cos2phi_ini0 = dcos( two*phi0 )
      this%sin2phi_ini0 = dsin( two*phi0 )


      this%Vector_of_Momentum0(1) = this%sin_Theta0 * this%cosphi_ini0
      this%Vector_of_Momentum0(2) = this%sin_Theta0 * this%sinphi_ini0
      this%Vector_of_Momentum0(3) = this%cos_Theta0
 
      CALL Set_psi_phi_chi_zera_array() 
      this%effect_number = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !this%delta_pds = zero
      this%v_L_v_i = zero
      this%d_theta = one / Num_PolDeg
      this%d_phi =  twopi / Num_Phi 

      this%phi_estimat(1) = zero
      this%phi_estimat(2) = pi / two 
      this%phi_estimat(3) = pi
      do i = 1, Num_PolDeg
          this%mu_estimates( i ) = this%d_theta * i
      enddo
      this%mu_estimates( 0 ) = 1.D-3
 
      this%Psi_I0 = one
      this%Psi_Q0 = one / 4.D0
      this%Psi_U0 = one / 4.D0
      this%Psi_V0 = one
   
      RETURN
      END SUBROUTINE Set_initial_parameter_values_Sub

 

!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub( this ) 
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this  
      integer cases
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      this%z_tau = this%Get_scatter_distance_tau( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !CALL this%Calc_Phot_Informations_At_Observer_Diffuse_Reflec() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      this%Psi_I = this%Psi_I * this%NormalA  
      this%Psi_Q = this%Psi_Q * this%NormalA  
      this%Psi_U = this%Psi_U * this%NormalA  
      this%Psi_V = this%Psi_V * this%NormalA  
      !write(*,*)'ss2=', this%NormalA, this%Vector_of_Momentum_ini(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      !this%Vector_of_Momentum_p = this%Vector_of_Momentum_ini
      !this%Vector_of_position_p = this%Vector_of_position_ini
      !write(unit = *, fmt = *)'************************************************************'

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub



!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this )
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
!~~~~~~~~~Implement estimations for observational quantities before scattering~~~~~~~~~~~~~~  
      CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(1)
      CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(2) 
      CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(3) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~Scatterings~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Tompson_Scat_With_Polarized_Diffuse_Reflection() 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
 
      end module PhotonModule



