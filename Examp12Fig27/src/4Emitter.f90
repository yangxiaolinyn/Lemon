      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter  
          real(mcp) :: w_ini_em   

      contains  
          procedure, public :: get_Phot4k_CtrCF_CovCF_IQ_only  =>  &
                               get_Phot4k_CtrCF_CovCF_IQ_only_Sub  

      end type Photon_Emitter 
 
      private :: get_Phot4k_CtrCF_CovCF_IQ_only_Sub  

      contains 
  

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: v,  index_nu, phi
      real(mcp) :: sin_theta, cos_theta, C_tau
 
      this%s_var = this%s_max * ranmar()
      
      end subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub

      end module PhotonEmitter





