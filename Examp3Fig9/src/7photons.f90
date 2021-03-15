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
          real(mcp) :: mu_estis_sq(0: Num_PolDeg) = zero 
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
          this%mu_estis_sq(i) = this%mu_estis(i)**2
      enddo
      this%mu_estis(0) = dy * 0.5D0
      this%mu_estis_sq(0) = this%mu_estis(0)**2

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
      !this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini  
      this%mu_zp_p = this%mu_zp_ini

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub 
 
  
!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this )
!************************************************************************************
      IMPLICIT NONE 
      class(Photon_FlatSP) :: this  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p
      !sphot%I_IQ     = this%I_IQ
      !sphot%Q_IQ     = this%Q_IQ
      !sphot%Vector_of_Momentum_ini(3) = this%Vector_of_Momentum_ini(3) 
      !sphot%z_tau = this%z_tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Calc_Esti_Informations_For_mu_Given()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL this%Tompson_Scattering_WithOut_Polarization_IQ()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*******************************************************************************************************
      subroutine Calc_Esti_Informations_For_mu_Given_J_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv, p_out1, mu_ini, mu, Psi_I !, Psi_Q
      integer :: i_obs , i_mu
 
          do i_mu = 0, Num_PolDeg
              mu = this%mu_estis(i_mu) 
              p_out1 = this%z_tau / mu
              vLv = dexp( - p_out1 ) / mu 

              Psi_I = one 
              !Psi_Q = zero 

              this%PolarArrayI( i_mu ) = this%PolarArrayI( i_mu ) + vLv * Psi_I 
              !this%PolarArrayQ( i_mu ) = this%PolarArrayQ( i_mu ) + vLv * Psi_Q
          enddo 
          return
      end subroutine Calc_Esti_Informations_For_mu_Given_J_Sub


!*******************************************************************************************************
      subroutine Calc_Esti_Informations_For_mu_Given_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_FlatSP) :: this 
      real(mcp) ::  vLv, p_out1, mu_ini, mup2, mu2, mu, mup, Psi_I, Psi_Q
      integer :: i_obs , i_mu
  
      !mup = this%Vector_of_Momentum_ini(3)
      mup = this%mu_zp_p
      mup2 = mup**2

          do i_mu = 0, Num_PolDeg 
              mu = this%mu_estis(i_mu)
              !mu2 = mu**2
              mu2 = this%mu_estis_sq(i_mu)
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




 
