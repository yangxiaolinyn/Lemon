      module PhotonEmitterBB
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter_BB
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

      contains 
          procedure, public :: get_Phot4k_CtrCF_CovCF   =>  &
                               get_Phot4k_CtrCF_CovCF_Sub 
      end type Photon_Emitter_BB
   
      private :: get_Phot4k_CtrCF_CovCF_Sub

      contains 
!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: maxn, theta, phi, w_ini, v, fn, nu, index_nu
      real(mcp) :: p0(1: 3), p2(1: 3), cos_BigTheta, r1, sin_theta, cos_theta, &
                   cos_Theta1, sin_Theta1, chi_phi, mu0, sin_phi, cos_phi, Chi
      integer :: i
  
      v = mec2
      this%Phot4k_CtrCF(1) = v 
      cos_theta = - ranmar()
      sin_theta = dsqrt( one - cos_theta**2 )
      phi =  twopi * ranmar()
      sin_phi = dsin( phi )
      cos_phi = dcos( phi )

      this%Vector_of_Momentum(1) = sin_theta * cos_phi
      this%Vector_of_Momentum(2) = sin_theta * sin_phi
      this%Vector_of_Momentum(3) = cos_theta 
      this%Phot4k_CtrCF(2: 4) = this%Phot4k_CtrCF(1) * this%Vector_of_Momentum(1:3)

      this%Phot4k_CovCF = this%Phot4k_CtrCF
      this%Phot4k_CovCF(1) = - this%Phot4k_CovCF(1)

      this%w_ini_em = one 

      this%z_tau = zero 
  
      end subroutine get_Phot4k_CtrCF_CovCF_Sub

      end module PhotonEmitterBB





