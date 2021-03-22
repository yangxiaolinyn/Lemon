      module PhotonEmitter
      use Basic_Variables_And_Methods
      implicit none
      integer, parameter :: N_wt1 = 2000

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

      contains  
          procedure, public :: get_Phot4k_CtrCF_CovCF_Reflection  =>  get_Phot4k_CtrCF_CovCF_Reflection_Sub
          procedure, public :: max_plankexp    =>  max_plankexp_fn
          procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn 

      end type Photon_Emitter
 
      private :: get_Phot4k_CtrCF_CovCF_Reflection_Sub
      private :: max_plankexp_fn
      private :: scatter_photon_blackbody_fn 

      contains
 
!*******************************************************************************************************
      function scatter_photon_blackbody_fn( this, maxn, T_s, w_ini_em )
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: scatter_photon_blackbody_fn
      real(mcp) :: fn, maxfn, v, xn
      real(mcp),intent(in) :: maxn, T_s
      real(mcp),intent(out) :: w_ini_em
      !real(mcp) ::  Tb = 2.D-3
      integer :: n = 2
 
      !maxfn = (maxn)**2/( dexp(maxn) - one )        

      do
          !v = ranmar()*(26.D-3 - 0.1D-3) + 0.1D-3
          v = ranmar() * ( this%E_up1 - this%E_low1 ) + this%E_low1
          fn = one/( dexp(v/T_s) - one )*(v/T_s)**2
          if( ranmar() <= fn/this%f_max )exit
      enddo 
      scatter_photon_blackbody_fn = v
      w_ini_em = one 
      end function scatter_photon_blackbody_fn

!*******************************************************************************************************
      function max_plankexp_fn( this, n )
!******************************************************************************************************* 
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: max_plankexp_fn
      real(mcp) :: xStart = 1.D-3, xn, xn1, fn, fn1 
      integer, intent(in) :: n
      
      xn = xStart 
      do
          fn = n*( one - dexp(-xn) )
          xn1 = xn
          xn = fn 
          if(dabs(xn-xn1) < 1.D-10)exit
      enddo 
      max_plankexp_fn = xn
      end function max_plankexp_fn
  

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Reflection_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: phi0, v, nu, index_nu
      real(mcp) :: r1, sin_theta, cos_theta, mu0, sin_phi, cos_phi, sin_mu0, s1 
      integer :: i
  
      mu0 = 0.2D0
      sin_mu0 = dsqrt( one - mu0**2 )
      cos_theta = -mu0 !one - two * ranmar() !-0.8D0
      sin_theta = dsqrt( one - cos_theta**2 )
      phi0 = zero !twopi * ranmar()
      sin_phi = dsin( phi0 )
      cos_phi = dcos( phi0 )

      this%phi_ini = phi0
      this%cosphi_ini = cos_phi
      this%sinphi_ini = sin_phi
      this%cos2phi_ini = dcos( two*phi0 )
      this%sin2phi_ini = dsin( two*phi0 )
  
      this%Vector_of_Momentum(1) = sin_theta * cos_phi
      this%Vector_of_Momentum(2) = sin_theta * sin_phi
      this%Vector_of_Momentum(3) = cos_theta
      !this%Phot4k_CtrCF(2: 4) = this%Phot4k_CtrCF(1) * this%Vector_of_Momentum(1: 3)
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%z_tau = zero 
      this%x = - this%z_tau * dsqrt( one - mu0**2 )
      this%y = zero
      this%z = - this%z_tau !this%tau_max ! - s1 * mu0
      this%r = zero !dsqrt( this%z * this%z + this%y * this%y + this%x * this%x )
      this%Vector_of_position(1) = zero !this%x
      this%Vector_of_position(2) = zero !this%y
      this%Vector_of_position(3) = zero !this%z 
  
      this%w_ini_em = one 
      !write(*, *)'ff=', this%w_ini_em, ( one + cos_Theta1 ** 2 ) 
      
      end subroutine get_Phot4k_CtrCF_CovCF_Reflection_Sub


      end module PhotonEmitter





