      module Photons_FlatSP
      use EstimationsModule
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type, public, extends(PhotonEstimate) :: Photon_FlatSP
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   

      contains 
!**************************************************************************************** 
          procedure, public :: Set_Initial_Values_For_Photon_Parameters    =>    &
                               Set_Initial_Values_For_Photon_Parameters_Sub 
          procedure, public :: Generate_A_Photon   =>   Generate_A_Photon_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p   =>   &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering    =>    &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE   =>   &
                               FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Photon_Electron_Scattering   =>   &
                               Photon_Electron_Scattering_Sub  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
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
      subroutine Set_Initial_Values_For_Photon_Parameters_Sub( this, &
                     T_elec, T_bb, tau, E1_scat, E2_scat, CrossSec_filename )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp), intent(in) :: T_elec, T_bb, tau, E1_scat, E2_scat
      character*80, intent(inout) :: CrossSec_filename
      integer :: i
      real(mcp) :: dy

!~~~~~~~~~Emitte a photon and take its parameters as initial values~~~~~~~~~~~~  
      this%tau_max = tau  
      this%T_e = T_elec 
      this%T_s = T_bb
      !CALL Emitter%Set_Emin_Emax() 
      !this%nu_low = Emitter%nu_low
      !this%nu_up = Emitter%nu_up
      !this%ln_nu1 = Emitter%ln_nu1
      !this%ln_nu2 = Emitter%ln_nu2
      this%logE_low = DLOG10( E1_scat )
      this%logE_up = DLOG10( E2_scat )
 
      this%n_e = 1.D16
      this%tau_max = tau / sigma_T / this%n_e
      this%effect_number = 0
      
      this%CrossSectFileName = CrossSec_filename

      write(*, *)'fsd', this%CrossSectFileName
      CALL this%Set_Cross_Section_3Te()
      CALL Set_psi_phi_chi_zera_array() 
     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%mu_esti(1) = 0.05D0
      this%mu_esti(2) = 0.1D0
      this%mu_esti(3) = 0.15D0
      this%mu_esti(4) = 0.25D0
      this%mu_esti(5) = 0.50D0
      this%mu_esti(6) = 0.75D0 
      this%mu_esti(7) = 0.95D0
      this%mu_esti(8) = 1.D0
      this%smu_esti(1: 8) = dsqrt(one - this%mu_esti(1: 8)**2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%y1 = 0.D0
      this%y2 = 0.6D0
      this%dy = ( this%y2 - this%y1 ) / vL_sc_up
      do i = 0, vL_sc_up
          this%E_array_esti(i) =  this%y1 + this%dy * i 
      enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      dy = twopi / N_esti_phi
      do i = 0, N_esti_phi
          this%phi_esti(i) = dy * i
          this%sin_phi_esti(i) = dsin(this%phi_esti(i))
          this%cos_phi_esti(i) = dcos(this%phi_esti(i))
      enddo   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end subroutine Set_Initial_Values_For_Photon_Parameters_Sub


!*******************************************************************************************************
      subroutine Generate_A_Photon_Sub( this, Emitter )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this
      TYPE(Photon_Emitter_BB), intent(inout) :: Emitter  
      
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL Emitter%get_Phot4k_CtrCF_CovCF()   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%z_tau  = Emitter%z_tau
      this%Vector_of_Momentum_ini = Emitter%Vector_of_Momentum 
      this%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF 
      this%Phot4k_CovCF_ini = Emitter%Phot4k_CovCF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%E_ini = DABS( Emitter%Phot4k_CovCF(1) )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%w_ini = Emitter%w_ini_em 
      this%w_ini0 = Emitter%w_ini_em 
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
      !write(*, *)'fs=',this%w_ini, this%NormalA
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini
      this%Phot4k_CovCF_At_p = this%Phot4k_CovCF_ini  

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_sub



!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      integer :: mu_i
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Set_Photon_Tetrad_In_CF() 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis
      !write(*, *)'ss1=', this%scatter_times
      CALL Inverse_Matrix_Of_Matrix_3X3_Sub(this%Matrix_Of_Tetrad_Of_PhotAxis, this%Matrix_CF_2_Elec)
      do  mu_i = 1, 8
          CALL this%Get_K_P1_P2_Scat_Kernel_for_Estimation( mu_i )
      enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Static_Compton_Scattering_WithOut_Polarizations() 
      !write(*, *)'PXfddddd = ', Vector4D_Inner_Product_Mski( sphot%Scattered_Phot4k_CF, &
      !!                        sphot%f4_scat_CF ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Phot4k_CtrCF_ini = sphot%Scattered_Phot4k_CF
      this%Phot4k_CovCF_ini = sphot%Scattered_Phot4k_CovCF
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / dabs( this%Phot4k_CtrCF_ini(1) )
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mmsf=', this%Phot4k_CtrCF_ini
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      this%E_ini = DABS( this%Phot4k_CovCF_ini(1) )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !phot%vector_of_position_ini = phot%vector_of_position_p
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 

!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot
      integer :: mu_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p
      sphot%delta_pd = this%delta_pd
      !sphot%f4_CF    = this%f4_CF
      !sphot%f4_CovCF = this%f4_CovCF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Set_Photon_Tetrad_In_CF()
      !CALL sphot%Set_Phot_f4_Tetrad_In_Phot_CF()
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Get_gama_mu_phi_Of_Scat_Elec(T_e)
      !CALL sphot%Get_gama_mu_phi_Of_Scatter_Electron_HXM( this%T_e )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Set_Elec_Tetrad_In_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL sphot%Set_Phot4k_In_Elec_CF()
      !CALL sphot%Set_f4_In_Elec_CF()
      !CALL sphot%Set_Phot_f4_Tetrad_In_Elec_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !this%Phot4k_In_Elec_CF = sphot%Phot4k_In_Elec_CF
      !this%f4_In_Elec_CF = sphot%f4_In_Elec_CF
      !this%Phot_f4_AxisX = sphot%Phot_f4_AxisX
      !this%Elec_Phot_mu = sphot%Elec_Phot_mu
      !this%Elec_Phot_sin = sphot%Elec_Phot_sin
      !this%Elec_gama = sphot%Elec_gama
      !this%Elec_V = sphot%Elec_V
      !this%Matrix_Of_Tetrad_Of_ElecAxis = sphot%Matrix_Of_Tetrad_Of_ElecAxis
      this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis
      CALL Inverse_Matrix_Of_Matrix_3X3_Sub(this%Matrix_Of_Tetrad_Of_PhotAxis, this%Matrix_CF_2_Elec)
      do  mu_i = 1, 8
          CALL this%Get_K_P1_P2_Scat_Kernel_for_Estimation( mu_i )
      enddo
      !this%Matrix_Of_Tetrad_Of_Phot_f4_Axis = sphot%Matrix_Of_Tetrad_Of_Phot_f4_Axis
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol(1)  
      !CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withpol(2)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !write(*,*)'555sss3=',  sphot%f4_In_Elec_CF, sphot%Phot4k_In_Elec_CF
      !stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Static_Compton_Scattering_WithOut_Polarizations() 
      !if( sphot%delta_pd /= zero )then 
      !    CALL sphot%Compton_Scattering_With_Polarization() 
      !else  
      !    CALL sphot%Compton_Scattering_With_Zero_QU() 
      !endif
      !write(*,*)'555sss3=',  sphot%delta_pd
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Photons_FlatSP




 
