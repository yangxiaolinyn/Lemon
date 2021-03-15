      module PhotonEmitterBB
      use Basic_Variables_And_Methods
      implicit none
      integer, parameter :: N_wt1 = 2000

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
          real(mcp) :: weights(1: N_wt1)
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: dnu
          real(mcp) :: Z_max

      contains  
          procedure, public :: get_Phot4k_CtrCF_CovCF_IQ_only  =>  get_Phot4k_CtrCF_CovCF_IQ_only_Sub 
          procedure, public :: max_plankexp    =>  max_plankexp_fn
          procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn
          !procedure, public :: Set_Emin_Emax  =>  Set_Emin_Emax_Sub 

      end type Photon_Emitter_BB
 
      private :: get_Phot4k_CtrCF_CovCF_IQ_only_Sub 
      private :: max_plankexp_fn
      private :: scatter_photon_blackbody_fn
      !private :: Set_Emin_Emax_Sub 

      contains
 
!*******************************************************************************************************
      function scatter_photon_blackbody_fn( this, maxn, T_s, w_ini_em )
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
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
      class(Photon_Emitter_BB) :: this
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
      subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: v,  index_nu, phi
      real(mcp) :: sin_theta, cos_theta, C_tau
      !integer :: i
 
      this%ln_nu2 = 22.D0
      this%ln_nu1 = 11.D0
      index_nu = this%ln_nu1 + ( this%ln_nu2 - this%ln_nu1 ) * 0.5D0
      !write(*,*)'emitter=', index_nu,10.D0**( index_nu )
      v = 10.D0**( index_nu ) * h_ev * 1.D-6
      this%Phot4k_CtrCF(1) = v
      !theta = pi / two * ranmar()
      cos_theta = ranmar() !one - two * ranmar() !dcos( pi /two * ranmar() ) 
      sin_theta = dsqrt( one - cos_theta**2 )
      phi = twopi * ranmar()
      !this%Phot4k_CtrCF(2) = this%Phot4k_CtrCF(1) * sin_theta * DCOS(phi)
      !this%Phot4k_CtrCF(3) = this%Phot4k_CtrCF(1) * sin_theta * DSIN(phi)
      !this%Phot4k_CtrCF(4) = this%Phot4k_CtrCF(1) * cos_theta
    
      this%Vector_of_Momentum(1) = sin_theta * DCOS(phi)
      this%Vector_of_Momentum(2) = sin_theta * DSIN(phi)
      this%Vector_of_Momentum(3) = cos_theta

      this%Phot4k_CtrCF(2: 4) = this%Phot4k_CtrCF(1) * this%Vector_of_Momentum(1: 3)

      !this%Phot4k_CovCF = this%Phot4k_CtrCF
      !this%Phot4k_CovCF(1) = - this%Phot4k_CovCF(1) 

      this%z_tau = this%tau_max ! * ranmar()
      !C_tau = one - dexp( - this%tau_max )
      !this%z_tau = - dlog( one - C_tau * ranmar() )
      this%w_ini_em = one !cos_theta 
          !cos_theta !cos_theta**2 * ( this%z_tau + 100.705927D0 )
      
      end subroutine get_Phot4k_CtrCF_CovCF_IQ_only_Sub

      end module PhotonEmitterBB





