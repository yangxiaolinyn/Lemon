      module Photons_FlatSP
      use ModuleForCompScatEsti 
      implicit none 

      type, public, extends(Photon_ForEstimation) :: Photon_FlatSP 
      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Set_Initial_Values_For_Photon_Parameters    =>    &
                               Set_Initial_Values_For_Photon_Parameters_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Generate_A_Photon   =>   Generate_A_Photon_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p    =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering    =>    &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE   =>    FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Photon_Electron_Scattering   =>   Photon_Electron_Scattering_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_FlatSP
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Set_Initial_Values_For_Photon_Parameters_Sub
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
      private :: Set_InI_Conditions_For_Next_Scattering_Sub
      private :: FIRST_SCATTERING_OF_PHOT_ELCE_Sub
      private :: Photon_Electron_Scattering_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
      contains     
  

!*******************************************************************************************************
      subroutine Set_Initial_Values_For_Photon_Parameters_Sub( this, T_elec, &
                            CrossSec_filename, S_in, alp, &
                            gama1, gama2, vy1, vy2, E_ini )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp), intent(in) :: T_elec, S_in(1: 4), alp, gama1, gama2, vy1, vy2, E_ini
      character*80, intent(inout) :: CrossSec_filename
      integer :: i
      real(mcp) :: dy

!~~~~~~~Emitte a photon and take its parameters as initial values~~~~~~~~~~~~~~~ 
      this%E_ini = E_ini

      this%alp = alp
      this%gama1 = gama1
      this%gama2 = gama2   
      this%T_e = T_elec 
      this%log10_Tbb = dlog10( E_ini )  
       
      this%n_e1 = 1.D20 
      this%ne_times_SigmaT = this%n_e1 * Sigma_T 
      this%effect_number = 0 
      this%CrossSectFileName = CrossSec_filename  
 
      CALL this%Set_Cross_Section_3Te()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%y1 = vy1
      this%y2 = vy2
      this%dy = ( this%y2 - this%y1 ) / vL_sc_up 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%mu_estimates(1) = dcos(85.D0 * dtor) 
      this%smu_estimates(1) = dsqrt(one - this%mu_estimates(1)**2) 

      dy = twopi / Num_phi
      do i = 0, Num_phi - 1
          this%phi_estimates(i) = pi / two !dy * i
          this%sin_phi_esti(i) = one  !dsin(this%phi_estimates(i))
          this%cos_phi_esti(i) = zero !dcos(this%phi_estimates(i)) 
      enddo 
      this%P_mu_normal = one / twopi 
      this%Vector_Stokes4_CF(1) = S_in(1)
      this%Vector_Stokes4_CF(2) = S_in(2)
      this%Vector_Stokes4_CF(3) = S_in(3)
      this%Vector_Stokes4_CF(4) = S_in(4)
       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end subroutine Set_Initial_Values_For_Photon_Parameters_Sub


!*******************************************************************************************************
      subroutine Generate_A_Photon_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this   
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%get_Phot4k_CtrCF_CovCF_BoundReflec()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      this%Q_sp = zero
      this%U_sp = zero
      this%V_sp = zero
      this%delta_pd = zero
      this%Psi_I = one
      this%Psi_Q = zero
      this%Psi_U = zero
      this%Psi_V = zero 
      
!~~~ Set the initial polarization vector: f = (1, 0, 0), i.e., along the x-axis.
      this%f4_CF(1) = zero
      this%f4_CF(2) = one
      this%f4_CF(3) = zero
      this%f4_CF(4) = zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Generate_A_Photon_Sub

!************************************************************************************
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub( this ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this  
   
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
      class(Photon_FlatSP) :: this
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      real(mcp) :: Sigma_COH, Sigma_INCOH, P_Rayl, P_Comp, P_fluor, xi_1, vLv
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      sphot%alp = this%alp
      sphot%gama1 = this%gama1
      sphot%gama2 = this%gama2 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_ini
      !sphot%delta_pd = this%delta_pd
      sphot%f4_CF    = this%f4_CF 
      sphot%Vector_Stokes4_CF = this%Vector_Stokes4_CF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !write(*, fmt="(' ', A5, 1ES18.7)")'ss2=', &
      !    Vector3D_Inner_Product( sphot%Phot4k_CtrCF(2: 4), sphot%f4_CF(2: 4) )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Photon_f3_Tetrad_In_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
      !call sphot%Get_gama_mu_phi_Of_Scat_Elec( this%T_e ) 
      CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_Power( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Elec_Tetrad_In_CF()
      CALL sphot%StokesPara_Rotation_Matrix(sphot%Elec_Phot_phi, sphot%Vector_Stokes4_CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      sphot%Vector_Stokes4_ECF = sphot%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Compton_Scattering_With_Polar_StokesVec( ) 
 

      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
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
      class(Photon_FlatSP) :: this 
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
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
      CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_Power( this%T_e )
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
 
      end module Photons_FlatSP




 
