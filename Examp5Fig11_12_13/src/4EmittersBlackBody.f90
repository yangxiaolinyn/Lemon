      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter
          real(mcp) :: R_in
          real(mcp) :: R_out
          real(mcp) :: w_ini_em
          real(mcp) :: T_s
          real(mcp) :: max_nu
          real(mcp) :: E_low1
          real(mcp) :: E_up1
          real(mcp) :: E_max
          real(mcp) :: f_max 
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: dnu 
          real(mcp) :: sin_Theta0, cos_Theta0
          !real(mcp) :: sin_phi0, cos_phi0
          real(mcp) :: Vector_of_Momentum0(1: 3)
          real(mcp) :: Psi_I0
          real(mcp) :: Psi_Q0
          real(mcp) :: Psi_U0
          real(mcp) :: Psi_V0
          real(mcp) :: phi_ini0 
          real(mcp) :: cosphi_ini0
          real(mcp) :: sinphi_ini0
          real(mcp) :: cos2phi_ini0
          real(mcp) :: sin2phi_ini0
          real(mcp) :: w_ini0


      contains   
          procedure, public :: Emitter_A_Photon  =>   Emitter_A_Photon_Sub 
          !procedure, public :: max_plankexp    =>  max_plankexp_fn
          !procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn 

      end type Photon_Emitter
 
      private :: Emitter_A_Photon_Sub
      !private :: max_plankexp_fn
      !private :: scatter_photon_blackbody_fn 

      contains    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      SUBROUTINE Emitter_A_Photon_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon_Emitter) :: this    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
      this%z_tau = zero 
      this%w_ini_em = one
      this%w_ini0 = one 

      this%r_ini = zero
      !this%theta_ini  = Emitter%theta
      !this%mucos_ini  = Emitter%mucos
      !this%musin_ini  = Emitter%musin
      !this%phi_ini    = Emitter%phi 
      !this%x_ini  = Emitter%x
      !this%y_ini  = Emitter%y
      !this%z_ini  = Emitter%z
      !this%z_tau  = Emitter%z_tau
      this%Vector_of_Momentum_ini = this%Vector_of_Momentum0
      this%Vector_of_position_ini = zero
      !this%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF 
      !this%Phot4k_CovCF_ini = Emitter%Phot4k_CovCF 
      this%cosphi_ini = this%cosphi_ini0
      this%sinphi_ini = this%sinphi_ini0
      this%cos2phi_ini = this%cos2phi_ini0
      this%sin2phi_ini = this%sin2phi_ini0
      this%phi_ini = this%phi_ini0
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !this%E_ini = DABS( Emitter%Phot4k_CovCF(1) )  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Psi_I = this%Psi_I0
      this%Psi_Q = this%Psi_Q0
      this%Psi_U = this%Psi_U0
      this%Psi_V = this%Psi_V0

      RETURN
      END SUBROUTINE Emitter_A_Photon_Sub


      end module PhotonEmitter





