      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter 
          real(mcp) :: w_ini_em
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
   
      this%Phot4k_CtrCF_ini(1) = one
 
      cos_theta = ranmar() !one - two * ranmar() !dcos( pi /two * ranmar() ) 
      sin_theta = dsqrt( one - cos_theta**2 )
      phi = twopi * ranmar()
  
      this%Vector_of_Momentum_ini(1) = sin_theta * DCOS(phi)
      this%Vector_of_Momentum_ini(2) = sin_theta * DSIN(phi)
      this%Vector_of_Momentum_ini(3) = cos_theta

      this%Phot4k_CtrCF_ini(2: 4) = this%Phot4k_CtrCF_ini(1) * this%Vector_of_Momentum_ini(1: 3)
 
      this%z_tau = this%tau_max 
      this%w_ini_em = one 
      this%w_ini = this%w_ini_em 
      
      end subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub

      end module PhotonEmitter





