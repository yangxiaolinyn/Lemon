      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter 
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
          procedure, public :: get_Phot4k_CtrCF_CovCF_BoundReflec   =>  &
                               get_Phot4k_CtrCF_CovCF_BoundReflec_Sub 
      end type Photon_Emitter
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: get_Phot4k_CtrCF_CovCF_BoundReflec_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      contains  
!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_BoundReflec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, theta, phi, w_ini, v, fn, nu, index_nu
      real(mcp) :: p0(1: 3), p2(1: 3), cos_BigTheta, r1, sin_theta, cos_theta, &
                   cos_Theta1, sin_Theta1, chi_phi, mu0, sin_phi, cos_phi, Chi
      integer :: i
  
!~~~ The incident photons are monochromatic, the energy is E_ini, which is a constant.
      this%Phot4k_CtrCF_ini(1) = this%E_ini
!~~~ Set the inital direction of the photon, which is along the z-axis, i.e., p=(0, 0, 1).
      cos_theta = one 
      sin_theta = zero 
      phi = twopi * ranmar()
      sin_phi = dsin( phi )
      cos_phi = dcos( phi )

      this%Vector_of_Momentum_ini(1) = sin_theta * cos_phi
      this%Vector_of_Momentum_ini(2) = sin_theta * sin_phi
      this%Vector_of_Momentum_ini(3) = cos_theta 
      this%Phot4k_CtrCF_ini(2: 4) = this%Phot4k_CtrCF_ini(1) * this%Vector_of_Momentum_ini(1:3)

      !this%Phot4k_CovCF = this%Phot4k_CtrCF 
 

      this%w_ini_em = one
      !this%w_ini = this%w_ini_em
      !this%w_ini0 = this%w_ini_em

      this%z_tau = zero 
  
      end subroutine get_Phot4k_CtrCF_CovCF_BoundReflec_Sub

      end module PhotonEmitter





