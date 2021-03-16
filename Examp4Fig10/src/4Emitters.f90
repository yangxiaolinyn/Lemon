      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter
          real(mcp) :: w_ini_em
          real(mcp) :: mu_zp_ini
          real(mcp) :: mu_zp_p
          real(mcp) :: mu_zp_sc
      contains  
          procedure, public :: get_Phot4k_CtrCF_CovCF_IQ_only  =>  get_Phot4k_CtrCF_CovCF_IQ_only_Sub 
          !procedure, public :: max_plankexp    =>  max_plankexp_fn
          !procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn 

      end type Photon_Emitter
 
      private :: get_Phot4k_CtrCF_CovCF_IQ_only_Sub 
      !private :: max_plankexp_fn
      !private :: scatter_photon_blackbody_fn 

      contains   
!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: v,  index_nu, phi
      real(mcp) :: sin_theta, cos_theta, C_tau
    
      do
          this%mu_zp_ini = one - two * ranmar() !dabs( dcos( pi / two * ranmar() ) )
          if( this%mu_zp_ini /= zero )exit
      enddo 
 
      this%z_tau = this%tau_max / two
      this%w_ini_em = one 
      this%w_ini = this%w_ini_em 

      this%I_IQ = this%w_ini_em 
      this%Q_IQ = zero   
      
      end subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub

      end module PhotonEmitter





