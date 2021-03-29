      module Photons_FlatSP
      use ModuleForCompScatEsti 
      USE MPI
      implicit none 

      type, public, extends(Photon_ForEstimation) :: Photon_FlatSP 
      contains 
!*******************************************************************************************************  
          procedure, public :: Set_Initial_Values_For_Photon_Parameters    =>    &
                               Set_Initial_Values_For_Photon_Parameters_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Generate_A_Photon   =>   Generate_A_Photon_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p    =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering    =>    &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE   =>    &
                               FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Photon_Electron_Scattering   =>   &
                               Photon_Electron_Scattering_Sub 
          procedure, public :: Photon_Electron_Scattering2   =>   &
                               Photon_Electron_Scattering2_Sub 
          procedure, public :: Photon_Electron_Scattering3   =>   &
                               Photon_Electron_Scattering3_Sub 
          procedure, public :: Set_InI_Conditions_For_Next_Scattering2  =>  &
                               Set_InI_Conditions_For_Next_Scattering2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_FlatSP
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Set_Initial_Values_For_Photon_Parameters_Sub
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
      private :: Set_InI_Conditions_For_Next_Scattering_Sub
      private :: FIRST_SCATTERING_OF_PHOT_ELCE_Sub
      private :: Photon_Electron_Scattering_Sub 
      private :: Photon_Electron_Scattering2_Sub
      private :: Photon_Electron_Scattering3_Sub 
      private :: Set_InI_Conditions_For_Next_Scattering2_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
      contains
!************************************************************************************* 
      subroutine Set_Initial_Values_For_Photon_Parameters_Sub( this, T_elec, &
                 T_bb, tau, E1_scat, E2_scat, y_obs1, y_obs2, mu_esti, sin_esti, &
                 Num_mu_esti, CrossSec_filename )
!************************************************************************************* 
      implicit none
      class(Photon_FlatSP) :: this  
      real(mcp), intent(in) :: T_elec, T_bb, tau, E1_scat, E2_scat, y_obs1, &
                       y_obs2, mu_esti(1: Num_mu_esti), sin_esti(1: Num_mu_esti)
      character*80, intent(inout) :: CrossSec_filename
      integer, intent(in) :: Num_mu_esti
      integer :: i, ierr
      real(mcp) :: dy

!~~~~~~~~~~~~Emitte a photon and take its parameters as initial values~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!    Set Initial conditions for the Emitter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%tau_max = tau
      this%T_s = T_bb 
      this%T_e = T_elec  
      CALL this%Set_Emin_Emax()   
  
      this%logE_low = DLOG10( E1_scat )
      this%logE_up = DLOG10( E2_scat )
      this%log10_Tbb = dlog10( mec2 )   
      
      this%n_e1 = 1.D20 
      this%ne_times_SigmaT = this%n_e1 * Sigma_T
      this%z_max = tau / this%ne_times_SigmaT
      !this%z_max1 = tau * 100000.00001D0 / this%ne_times_SigmaT
      this%effect_number = 0 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%CrossSectFileName = CrossSec_filename 

      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma)
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if( this%myid == this%num_np-1 )then
          CALL this%Set_Cross_Section()
      endif
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
      CALL MPI_BCAST( this%sigmaaTeE_FST, N_sigma + 1, MPI_DOUBLE_PRECISION, &
                        this%num_np-1, MPI_COMM_WORLD, ierr )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL Set_psi_phi_chi_zera_array() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%y1 = y_obs1
      this%y2 = y_obs2
      this%dy = ( this%y2 - this%y1 ) / vL_sc_up 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%num_mu_esti = Num_mu_esti
      this%mu_estimates(1: Num_mu_esti) = mu_esti
      this%smu_estimates(1: Num_mu_esti) = sin_esti
 
      this%smu_estimates = dsqrt(one - this%mu_estimates **2)
 
      dy = twopi / Num_phi
      do i = 0, Num_phi-1 
          this%phi_estimates(i) = zero 
          this%sin_phi_esti(i) = dsin(this%phi_estimates(i))
          this%cos_phi_esti(i) = dcos(this%phi_estimates(i)) 
      enddo 
      this%P_mu_normal = one / twopi 
      this%P_nu_normal = one / this%dy_em
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end subroutine Set_Initial_Values_For_Photon_Parameters_Sub


!**************************************************************************** 
      subroutine Generate_A_Photon_Sub( this )
!**************************************************************************** 
      implicit none
      class(Photon_FlatSP) :: this  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%get_Phot4k_CtrCF_CovCF_BoundReflec()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      this%Sigma_a_E_ini = this%sigma_fn( this%E_ini ) 
      this%ne_times_Sigma_a = this%Sigma_a_E_ini * this%n_e1 
      !write(*, *)'f1333 = ', this%Sigma_KN_E_ini, this%E_ini
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
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
      this%Vector_Stokes4_CF(4) = one
!~~~~~~~~~~~~~Implement initial estimation~~~~~~~~~~~~~~~~~~~~~~~~~ 
      call this%get_J_emissivity_for_estimation_Phiarr()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Generate_A_Photon_Sub


!************************************************************************************
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub( this ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this  
   
      this%z_tau = this%Get_scatter_distance_BoundReflec( )    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%w_ini = this%w_ini * this%NormalA  
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%StokesPara_Rotation_Matrix(sphot%Elec_Phot_phi, sphot%Vector_Stokes4_CF)
      !sphot%Vector_Stokes4_ECF = sphot%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Phot4k_In_Elec_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Making_An_Estimation_One_MC_Component_1( sphot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL sphot%Compton_Scattering_With_Zero_QU_StokesVec()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  


!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, Scattered_Phot4k_CF, &
                       Vector_Stokes4_ECF_scat, f4_scat_CF )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
      real(mcp), intent(in) :: Scattered_Phot4k_CF(1: 4), f4_scat_CF(1: 4), &
                               Vector_Stokes4_ECF_scat(1: 4)
      !TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot 
      !type( Photon_ForEstimation ), INTENT(INOUT) :: sphot


      this%Phot4k_CtrCF_ini = Scattered_Phot4k_CF 
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / dabs( this%Phot4k_CtrCF_ini(1) )
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mmsf=', this%Phot4k_CtrCF_ini

      this%Vector_Stokes4_CF = Vector_Stokes4_ECF_scat
      this%f4_CF    = f4_scat_CF 
      !write(*, fmt="(' ', A5, 1ES18.7)")'ss1=', &
      !    Vector3D_Inner_Product( this%Phot4k_CtrCF_ini(2: 4), this%f4_CF(2: 4) ) 

      this%E_ini = DABS( this%Phot4k_CtrCF_ini(1) ) 
      !write(*, fmt = "(' ', A10, ES20.6, I10)")'f1 = ', this%E_ini, this%scatter_times
      this%Sigma_a_E_ini = this%sigma_fn( this%E_ini )
      this%ne_times_Sigma_a = this%Sigma_a_E_ini * this%n_e1

      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering2_Sub( this )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Phot4k_CtrCF_ini = this%Scattered_Phot4k_CF 
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / &
                                       dabs( this%Phot4k_CtrCF_ini(1) )
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mmsf=', this%Phot4k_CtrCF_ini

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Vector_Stokes4_CF = this%Vector_Stokes4_ECF_scat
      this%f4_CF = this%f4_scat_CF 
      !write(*, fmt="(' ', A5, 1ES18.7)")'ss1=', &
      !    Vector3D_Inner_Product( this%Phot4k_CtrCF_ini(2: 4), this%f4_CF(2: 4) ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%E_ini = DABS( this%Phot4k_CtrCF_ini(1) ) 
      !write(*, fmt = "(' ', A10, ES20.6, I10)")'f1 = ', this%E_ini 
      this%Sigma_a_E_ini = this%sigma_fn( this%E_ini )
      this%ne_times_Sigma_a = this%Sigma_a_E_ini * this%n_e1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 

!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering2_Sub( this )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this  
      real(mcp) :: Sigma_COH, Sigma_INCOH, P_Rayl, P_Comp, P_fluor, xi_1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      this%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p    
      CALL this%Set_Photon_f3_Tetrad_In_CF()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Get_gama_mu_phi_Of_Scat_Elec( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Set_Elec_Tetrad_In_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%StokesPara_Rotation_Matrix( this%Elec_Phot_phi, this%Vector_Stokes4_CF )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%Vector_Stokes4_ECF = this%Vector_Stokes4_CF 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      CALL this%Compton_Scattering_With_Polar_StokesVec()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Photon_Electron_Scattering2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      !type( Photon_ForEstimation ), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      !sphot%delta_pd = this%delta_pd
      sphot%f4_CF    = this%f4_CF 
      sphot%Vector_Stokes4_CF = this%Vector_Stokes4_CF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Photon_f3_Tetrad_In_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
      CALL sphot%Get_gama_mu_phi_Of_Scat_Elec( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Elec_Tetrad_In_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%StokesPara_Rotation_Matrix(sphot%Elec_Phot_phi, sphot%Vector_Stokes4_CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Phot4k_In_Elec_CF()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      sphot%Vector_Stokes4_ECF = sphot%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Compton_Scattering_With_Polar_StokesVec()  
      !write(*, *)'f1 = ', sphot%Vector_Stokes4_ECF_scat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 




!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering3_Sub( this, &
                   Phot4k_CtrCF_At_p, f4_CF, Vector_Stokes4_CF )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
      real(mcp), intent(in) :: Phot4k_CtrCF_At_p(1: 4), f4_CF(1: 4), Vector_Stokes4_CF(1: 4)
      !TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      !type(Photon_FlatSP) :: sphot
      !real(mcp) :: Sigma_COH, Sigma_INCOH, P_Rayl, P_Comp, P_fluor, xi_1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%Phot4k_CtrCF = Phot4k_CtrCF_At_p
      !sphot%delta_pd = this%delta_pd
      this%f4_CF    = f4_CF 
      this%Vector_Stokes4_CF = Vector_Stokes4_CF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Set_Photon_f3_Tetrad_In_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
      CALL this%Get_gama_mu_phi_Of_Scat_Elec( this%T_e )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Set_Elec_Tetrad_In_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%StokesPara_Rotation_Matrix(this%Elec_Phot_phi, this%Vector_Stokes4_CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%Set_Phot4k_In_Elec_CF()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Set_Phot_Tetrad1_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%Vector_Stokes4_ECF = this%Vector_Stokes4_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !this%Vector_Stokes4_ECF = sphot%Vector_Stokes4_ECF
      !this%Phot4k_In_Elec_CF = sphot%Phot4k_In_Elec_CF
      !this%Elec_Phot_mu = sphot%Elec_Phot_mu
      !this%Elec_Phot_sin = sphot%Elec_Phot_sin
      !this%Elec_Phot_phi = sphot%Elec_Phot_phi
      !this%Elec_gama = sphot%Elec_gama
      !this%Elec_V = sphot%Elec_V

      !this%Elec_Phot_sin_In_Elec_CF = sphot%Elec_Phot_sin_In_Elec_CF 
      !this%Elec_Phot_mu_In_Elec_CF = sphot%Elec_Phot_mu_In_Elec_CF

      !this%Matrix_Of_Tetrad_Of_ElecAxis = sphot%Matrix_Of_Tetrad_Of_ElecAxis
      !this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis  
      !this%Matrix_Of_Tetrad1_Of_photAxis = sphot%Matrix_Of_Tetrad1_Of_photAxis
      !this%Matrix_ECF_2_ECF1 = sphot%Matrix_ECF_2_ECF1
      !this%Matrix_ECF1_2_ECF = sphot%Matrix_ECF1_2_ECF 
      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Compton_Scattering_With_Polar_StokesVec()  
      !write(*, *)'f1 = ', this%Vector_Stokes4_ECF_scat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      END SUBROUTINE Photon_Electron_Scattering3_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!*************************************************************************** 
 
      end module Photons_FlatSP




 
