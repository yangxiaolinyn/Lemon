      module PhotonModule
      use PhotonsEstimation
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photons_Esti) :: Photons
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 


      contains   
          procedure, public :: Set_initial_parameter_values => Set_initial_parameter_values_Sub
          procedure, public :: Emitter_A_Photon  =>   Emitter_A_Photon_Sub
          procedure, public :: Transmit_Data_And_Parameters_From_Emitter2Photon   =>   &
                               Transmit_Data_And_Parameters_From_Emitter2Photon_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p   =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE   =>   &
                               FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering  =>   &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: Determine_Next_Scattering_Site  =>  &
                               Determine_Next_Scattering_Site_Sub
          procedure, public :: Photon_Electron_Scattering    =>  &
                               Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photons

      private :: Set_initial_parameter_values_Sub
     

      contains   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values_Sub( this, Emitter, tau )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      class(Photons) :: this   
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter 
      REAL(mcp), INTENT(INOUT) :: tau
      integer :: i 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!| Set Initial conditions for the Emitter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      Emitter%tau_max = tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set Initial conditions for the Photon                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%tau_max = tau   
 
      CALL Set_psi_phi_chi_zera_array() 
      this%effect_number = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%delta_pds = zero
      this%v_L_v_i = zero
      this%d_theta = one / Num_PolDeg
      this%d_phi =  twopi / Num_Phi 

      this%phi_estimat(1) = zero
      this%phi_estimat(2) = pi / two
      !this%phi_estimat(3) = pi / four
      this%phi_estimat(3) = pi
      do i = 1, Num_PolDeg
          this%mu_estimates( i ) = this%d_theta * i
      enddo
      this%mu_estimates( 0 ) = 1.D-3
   
      RETURN
      END SUBROUTINE Set_initial_parameter_values_Sub

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon_Sub( this, Emitter )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photons) :: this   
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL Emitter%get_Phot4k_CtrCF_CovCF_Reflection()  

      RETURN
      END SUBROUTINE Emitter_A_Photon_Sub

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon_Sub( this, Emitter )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photons) :: this 
      TYPE(Photon_Emitter), INTENT(IN) :: Emitter 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      this%r_ini      = Emitter%r
      !this%theta_ini  = Emitter%theta
      !this%mucos_ini  = Emitter%mucos
      !this%musin_ini  = Emitter%musin
      !this%phi_ini    = Emitter%phi 
      this%x_ini  = Emitter%x
      this%y_ini  = Emitter%y
      this%z_ini  = Emitter%z
      this%z_tau  = Emitter%z_tau
      this%Vector_of_Momentum_ini = Emitter%Vector_of_Momentum
      this%Vector_of_position_ini = Emitter%Vector_of_position
      !this%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF 
      !this%Phot4k_CovCF_ini = Emitter%Phot4k_CovCF 
      this%cosphi_ini = Emitter%cosphi_ini
      this%sinphi_ini = Emitter%sinphi_ini
      this%cos2phi_ini = Emitter%cos2phi_ini
      this%sin2phi_ini = Emitter%sin2phi_ini
      this%phi_ini = Emitter%phi_ini
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !this%E_ini = DABS( Emitter%Phot4k_CovCF(1) )  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%w_ini = Emitter%w_ini_em 
      this%w_ini0 = Emitter%w_ini_em 
      this%Psi_I = one
      this%Psi_Q = one / 4.D0
      this%Psi_U = one / 4.D0
      this%Psi_V = one
  
      RETURN
      END SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon_Sub

!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub( this ) 
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this  
      integer cases
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      this%z_tau = this%Get_scatter_distance_tau( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Calc_Phot_Informations_At_Observer_Diffuse_Reflec() 
      this%Psi_I = this%Psi_I * this%NormalA  
      this%Psi_Q = this%Psi_Q * this%NormalA  
      this%Psi_U = this%Psi_U * this%NormalA  
      this%Psi_V = this%Psi_V * this%NormalA  
      !write(*,*)'ss2=', this%z_tau, this%p_scattering, this%Vector_of_Momentum_ini(3) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      this%Vector_of_Momentum_p = this%Vector_of_Momentum_ini
      !this%Vector_of_position_p = this%Vector_of_position_ini
      !write(unit = *, fmt = *)'************************************************************'

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub

!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this  
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Vector_of_Momentum_ini = this%Vector_of_Momentum_p
      !sphot%Vector_of_position_ini = this%Vector_of_position_p
      sphot%Psi_I = this%Psi_I
      sphot%Psi_Q = this%Psi_Q
      sphot%Psi_U = this%Psi_U
      sphot%Psi_V = this%Psi_V
      sphot%cosphi_ini = this%cosphi_ini
      sphot%sinphi_ini = this%sinphi_ini
      sphot%cos2phi_ini = this%cos2phi_ini
      sphot%sin2phi_ini = this%sin2phi_ini
      sphot%phi_ini = this%phi_ini  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If( .true. )then
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(1)
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(2)
          !CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(3)
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(3)
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection()   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this  
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Vector_of_Momentum_ini(3) = sphot%Vector_of_Momentum_ini(3)
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mms=', this%Phot4k_CtrCF_ini

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Psi_I = sphot%Psi_I
      this%Psi_Q = sphot%Psi_Q
      this%Psi_U = sphot%Psi_U
      this%Psi_V = sphot%Psi_V
      this%cosphi_ini = sphot%cosphi_ini
      this%sinphi_ini = sphot%sinphi_ini
      this%cos2phi_ini = sphot%cos2phi_ini
      this%sin2phi_ini = sphot%sin2phi_ini
      this%phi_ini = sphot%phi_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !this%E_ini = DABS( this%Phot4k_CovCF_ini(1) ) 
      !write(*,*)'5555==', this%E_ini, this%Phot4k_CovCF_ini(1) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site_Sub( this, sphot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      class(Photons) :: this  
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      this%z_tau = this%Get_scatter_distance_tau( ) 
      !write(*,*)'7777==', this%E_ini, this%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Calc_Phot_Informations_At_Observer_Diffuse_Reflec() 
      !this%w_ini = this%w_ini * this%NormalA   
      this%Psi_I = this%Psi_I * this%NormalA
      this%Psi_Q = this%Psi_Q * this%NormalA  
      this%Psi_U = this%Psi_U * this%NormalA  
      this%Psi_V = this%Psi_V * this%NormalA   
      !write(*,*)'ss3=',this%w_ini, this%NormalA, this%r_one_hvlmec2_one_cosE    
      this%Vector_of_Momentum_p = this%Vector_of_Momentum_ini  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photons) :: this  
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Vector_of_Momentum_ini = this%Vector_of_Momentum_p
      !sphot%Vector_of_position_ini = this%Vector_of_position_p
      sphot%Psi_I = this%Psi_I
      sphot%Psi_Q = this%Psi_Q
      sphot%Psi_U = this%Psi_U
      sphot%Psi_V = this%Psi_V
      sphot%cosphi_ini = this%cosphi_ini
      sphot%sinphi_ini = this%sinphi_ini
      sphot%cos2phi_ini = this%cos2phi_ini
      sphot%sin2phi_ini = this%sin2phi_ini
      sphot%phi_ini = this%phi_ini 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If( .true. )then
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(1)
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(2)
          !CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(3)
          CALL this%Calc_Phot_Inform_At_Observer_with_mu_phi_Given2(3)
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection() 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 
      end module PhotonModule



