      module Photons_FlatSP
      use ScatDistance_FlatSP
      !use PhotonEmitterSyn
      !use PhotonEmitter
      !use ScatPhotSequences
      implicit none 
      integer, parameter :: Num_mu = 100

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photon_With_ScatDistance_FlatSP) :: Photon_FlatSP 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: nu_obs 
          real(mcp) :: mu_estis(0: Num_PolDeg) = zero 
          real(mcp) :: PolarArrayd(0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayI(0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayQ(0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayU(0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayV(0: Num_PolDeg) = zero
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: d_theta, d_phi, d_tau 

      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Set_initial_parameter_values   =>   Set_initial_parameter_values_Sub
          procedure, public :: Emitter_A_Photon  =>   Emitter_A_Photon_Sub 
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p  =>  &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub
          procedure, public :: FIRST_SCATTERING_OF_PHOT_ELCE  =>   FIRST_SCATTERING_OF_PHOT_ELCE_Sub
          procedure, public :: Set_InI_Conditions_For_Next_Scattering  =>   &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: Determine_Next_Scattering_Site  =>   Determine_Next_Scattering_Site_Sub
          procedure, public :: Photon_Electron_Scattering   =>   Photon_Electron_Scattering_Sub
          procedure, public :: Calc_Phot_Informations_At_Observor_FLST_IQ  =>   &
                               Calc_Phot_Informations_At_Observor_FLST_IQ_Sub 
          procedure, public :: Calc_Esti_Informations_For_mu_Given   =>   &
                               Calc_Esti_Informations_For_mu_Given_Sub
          procedure, public :: Calc_Esti_Informations_For_mu_Given_J   =>   &
                               Calc_Esti_Informations_For_mu_Given_J_Sub
          
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_FlatSP
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Calc_Phot_Informations_At_Observor_FLST_IQ_Sub
      private :: Calc_Esti_Informations_For_mu_Given_Sub
      private :: Calc_Esti_Informations_For_mu_Given_J_Sub
      private :: Photon_Electron_Scattering_Sub
      private :: Determine_Next_Scattering_Site_Sub
      private :: Set_InI_Conditions_For_Next_Scattering_Sub
      private :: FIRST_SCATTERING_OF_PHOT_ELCE_Sub
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub 
      private :: Emitter_A_Photon_Sub
      private :: Set_initial_parameter_values_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values_Sub( this, tau )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      class(Photon_FlatSP) :: this  
      REAL(mcp), INTENT(INOUT) :: tau
      REAL(mcp) :: dy
      integer :: i
      
      this%effect_number = 0

      dy = one / Num_PolDeg
      do i = 1, Num_PolDeg
          this%mu_estis(i) = i * dy
      enddo
      this%mu_estis(0) = dy * 0.5D0

      this%tau_max = tau 
      this%PolarArrayI = zero
      this%PolarArrayQ = zero   
      this%PolarArrayU = zero
   
      RETURN
      END SUBROUTINE Set_initial_parameter_values_Sub


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon_FlatSP) :: this  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      CALL this%get_Phot4k_CtrCF_CovCF_IQ_only()  
 
      this%I_IQ = this%w_ini_em 
      this%Q_IQ = zero   
      CALL this%Calc_Esti_Informations_For_mu_Given_J()

      RETURN
      END SUBROUTINE Emitter_A_Photon_Sub
 

!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub( this ) 
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      this%z_tau = this%Get_scatter_distance_IQ( )    
      this%I_IQ = this%I_IQ * this%NormalA
      this%Q_IQ = this%Q_IQ * this%NormalA 
      !write(*,*)'ss2=', this%z_tau, this%p_scattering, this%Vector_of_Momentum_ini(3) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini  

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub

!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      sphot%I_IQ     = this%I_IQ
      sphot%Q_IQ     = this%Q_IQ
      sphot%Vector_of_Momentum_ini(3) = this%Vector_of_Momentum_ini(3) 
      sphot%z_tau = this%z_tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Calc_Esti_Informations_For_mu_Given()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL sphot%Tompson_Scattering_WithOut_Polarization_IQ()   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Phot4k_CtrCF_ini = sphot%Scattered_Phot4k_CF
      !this%Phot4k_CovCF_ini = sphot%Scattered_Phot4k_CovCF
      this%Vector_of_Momentum_ini(3) = this%Phot4k_CtrCF_ini(4) 
      if( isnan( this%Vector_of_Momentum_ini(3) ) )write(*, *)'mms=', this%Phot4k_CtrCF_ini

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%I_IQ     = sphot%I_IQ
      this%Q_IQ     = sphot%Q_IQ 
      !this%Sin_Theta_Scat = sphot%Sin_Theta_Scat
      !this%Cos_Theta_Scat = sphot%Cos_Theta_Scat
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
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !write(*,*)'6666==', this%E_ini, this%Phot4k_CovCF_ini(1) 
      !this%p_scattering = this%Get_scatter_distance_IQ( )  
      this%z_tau = this%Get_scatter_distance_IQ( )  
      !write(*,*)'7777==', this%E_ini, this%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !Call this%Calc_Phot_Informations_At_Observor_FLST_IQ( )
      !write(*,*)'ss3=',this%w_ini, this%NormalA, this%r_one_hvlmec2_one_cosE 
      this%I_IQ = this%I_IQ * this%NormalA
      this%Q_IQ = this%Q_IQ * this%NormalA
      !write(*,*)'ss4=',this%w_ini  
      !this%z_tau = this%z_tau - this%p_scattering * this%Vector_of_Momentum_ini(3)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !write(unit = *, fmt = *)'************************************************************'
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini
      !this%Phot4k_CovCF_At_p = this%Phot4k_CovCF_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, sphot )
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this 
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p
      sphot%I_IQ     = this%I_IQ
      sphot%Q_IQ     = this%Q_IQ
      sphot%Vector_of_Momentum_ini(3) = this%Vector_of_Momentum_ini(3) 
      sphot%z_tau = this%z_tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Calc_Esti_Informations_For_mu_Given()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL sphot%Tompson_Scattering_WithOut_Polarization_IQ()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*******************************************************************************************************
      subroutine Calc_Esti_Informations_For_mu_Given_J_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv, p_out1, mu_ini, mup2, mu2, mu, mup, Psi_I, Psi_Q
      integer :: i_obs , i_mu
  
      !mup = this%Vector_of_Momentum_ini(3)
      !mup2 = mup**2

          do i_mu = 0, Num_PolDeg
              mu = this%mu_estis(i_mu)
              mu2 = mu**2
              p_out1 = this%z_tau / mu
              vLv = dexp( - p_out1 ) / mu 

              Psi_I = one !mu !( this%I_IQ * ( 3.D0 - mup2 - mu2 + 3.D0 * mup2 * mu2 - mup * mu ) + &
                         ! this%Q_IQ * ( one - mup2 ) * ( one - 3.D0 * mu2 )  ) * three / 16.D0
              Psi_Q = zero !( this%I_IQ * ( one - mu2 ) * ( one - three * mup2 ) + this%Q_IQ * three * &
                           !( one - mu2 ) * ( one - mup2 ) ) * three / 16.D0

              this%PolarArrayI( i_mu ) = this%PolarArrayI( i_mu ) + vLv * Psi_I
  
              !this%PolarArrayQ( i_mu ) = this%PolarArrayQ( i_mu ) + vLv * Psi_Q
          enddo
 
      !Endif 
      return
      end subroutine Calc_Esti_Informations_For_mu_Given_J_Sub


!*******************************************************************************************************
      subroutine Calc_Esti_Informations_For_mu_Given_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv, p_out1, mu_ini, mup2, mu2, mu, mup, Psi_I, Psi_Q
      integer :: i_obs , i_mu
  
      mup = this%Vector_of_Momentum_ini(3)
      mup2 = mup**2

          do i_mu = 0, Num_PolDeg 
              mu = this%mu_estis(i_mu)
              mu2 = mu**2
              p_out1 =  this%z_tau / mu
              vLv = dexp( - p_out1 ) / mu 

              Psi_I = ( this%I_IQ * ( 3.D0 - mup2 - mu2 + 3.D0 * mup2 * mu2 ) + &
                        this%Q_IQ * ( one - mup2 ) * ( one - 3.D0 * mu2 )  )! * three / 16.D0
              Psi_Q = ( this%I_IQ * ( one - mu2 ) * ( one - three * mup2 ) + this%Q_IQ * three * &
                      ( one - mu2 ) * ( one - mup2 ) )! * three / 16.D0

              !write(*, *)'ffs==', vLv, Psi_I, Psi_Q, vLv * Psi_I, vLv * Psi_Q
              this%PolarArrayI( i_mu ) = this%PolarArrayI( i_mu ) + vLv * Psi_I
              !if(i_mu==2) write(*, *)'ffs2==', this%PolarArrayI( 2 ), vLv * Psi_I, this%I_IQ
              !if(i_mu==90) write(*, *)'ffs9==', this%PolarArrayI( 90 ), vLv * Psi_I, this%I_IQ
  
              this%PolarArrayQ( i_mu ) = this%PolarArrayQ( i_mu ) + vLv * Psi_Q
          enddo
      !stop
      !Endif 
      return
      end subroutine Calc_Esti_Informations_For_mu_Given_Sub



!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_FLST_IQ_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv
      integer :: i_obs 

      this%d_theta = one / Num_PolDeg

      if( this%InterSection_Cases == -1 .or. this%InterSection_Cases == -2 )then
          return
      else
   
          i_obs = floor( dabs( this%Vector_of_Momentum_ini(3) ) / this%d_theta )  
          !vLv = this%w_ini * dexp( - this%Optical_Depth_scatter ) / dabs(this%Vector_of_Momentum_ini(3))
          vLv = dexp( - this%Optical_Depth_scatter ) / dabs(this%Vector_of_Momentum_ini(3))

          this%PolarArrayI( i_obs ) = this%PolarArrayI( i_obs ) + vLv * this%I_IQ
  
          this%PolarArrayQ( i_obs ) = this%PolarArrayQ( i_obs ) + vLv * this%Q_IQ
          !write(*, *)'mmsff==', vLv * this%Q_IQ, vLv, dabs(this%Vector_of_Momentum_ini(3))
          !write(*, *)'mms44==', this%PolarArrayI( i_obs ), this%PolarArrayQ( i_obs ), vLv * this%I_IQ, & 
           !     vLv, this%I_IQ, this%Q_IQ  
          this%effect_number = this%effect_number + 1
 
      Endif 
      return
      end subroutine Calc_Phot_Informations_At_Observor_FLST_IQ_Sub 
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Photons_FlatSP




 
