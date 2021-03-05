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
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: dnu 

      contains  
          procedure, public :: get_Phot4k_CtrCF_CovCF_Reflection  =>  get_Phot4k_CtrCF_CovCF_Reflection_Sub
          procedure, public :: max_plankexp    =>  max_plankexp_fn
          procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn 
          procedure, public :: Set_Emin_Emax  =>  Set_Emin_Emax_Sub 
          procedure, public :: get_Phot4k_CtrCF_CovCF_BoundReflec   =>  &
                               get_Phot4k_CtrCF_CovCF_BoundReflec_Sub 
      end type Photon_Emitter_BB
 
      private :: get_Phot4k_CtrCF_CovCF_Reflection_Sub
      private :: max_plankexp_fn
      private :: scatter_photon_blackbody_fn 
      private :: Set_Emin_Emax_Sub 
      private :: get_Phot4k_CtrCF_CovCF_BoundReflec_Sub

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
      subroutine get_Phot4k_CtrCF_CovCF_Reflection_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: phi0, v, nu, index_nu
      real(mcp) :: r1, sin_theta, cos_theta, mu0, sin_phi, cos_phi, sin_mu0, s1 
      integer :: i
  
      mu0 = 0.5D0
      sin_mu0 = dsqrt( one - mu0**2 )
      cos_theta = -0.5D0 !one - two * ranmar() !-0.8D0
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

!*******************************************************************************************************
      subroutine Set_Emin_Emax_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: maxn, f1, a2, b2, dE, dnu
      integer :: i, istat
 
      this%E_max = this%max_plankexp(2) * this%T_s
      this%f_max = (this%E_max / this%T_s)**2/( dexp(this%E_max / this%T_s) - one ) 
      this%E_low1 = this%E_max * 0.9D0
      DO
          f1 = (this%E_low1 / this%T_s)**2/( dexp(this%E_low1 / this%T_s) - one )
          if(f1 < this%f_max * 1.D-8)exit
          this%E_low1 = this%E_low1 * 0.9D0
      ENDDO
      this%E_up1 = this%E_max * 1.1D0
      this%E_up1 = 10**(two) * this%T_s
      this%E_low1 = 10**(-two) * this%T_s
      DO
          f1 = (this%E_up1 / this%T_s)**2/( dexp(this%E_up1 / this%T_s) - one )
          if(f1 < this%f_max * 1.D-8)exit
          this%E_up1 = this%E_up1 * 1.1D0
      ENDDO
      if(.false.)then
          this%ln_nu2 = 22.D0
          this%ln_nu1 = 11.D0
          this%nu_low = 10.D0**(11.D0)
          this%nu_up = 10.D0**22.D0
          this%E_low1 = this%nu_low * h_ev * 1.D-6
          this%E_up1 = this%nu_up * h_ev * 1.D-6
      else
          !this%nu_low = this%E_low1*1.D6/h_ev
          !this%nu_up = this%E_up1*1.D6/h_ev
          !this%ln_nu2 = dlog10( this%nu_up )
          !this%ln_nu1 = dlog10( this%nu_low )
          !this%E_low1 = 10.D0**(-5.D0) * mec2
          !this%E_up1 = 10.D0**(-2.5D0) * mec2
          this%nu_low = this%E_low1*1.D6/h_ev
          this%nu_up = this%E_up1*1.D6/h_ev
          this%ln_nu2 = dlog10( this%nu_up )
          this%ln_nu1 = dlog10( this%nu_low )
      endif
      this%dnu = ( this%ln_nu2 - this%ln_nu1 ) / N_wt1
      !write(*, *)'mm=', this%E_low1*1.D6/h_ev, this%E_max*1.D6/h_ev, this%E_up1*1.D6/h_ev
      !write(*, *)'mm=', this%E_low1, this%E_max, this%E_up1
      !stop   
      
      end subroutine Set_Emin_Emax_Sub
  

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_BoundReflec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: maxn, theta, phi, w_ini, v, fn, nu, index_nu
      real(mcp) :: p0(1: 3), p2(1: 3), cos_BigTheta, r1, sin_theta, cos_theta, &
                   cos_Theta1, sin_Theta1, chi_phi, mu0, sin_phi, cos_phi, Chi
      integer :: i

      !do
          !index_nu = this%ln_nu1 + ranmar()*( this%ln_nu2 - this%ln_nu1 )
          !v = 10.D0**( index_nu ) * h_ev * 1.D-6
          !v = !ranmar() * ( this%E_up1 - this%E_low1 ) + this%E_low1
          !fn = one/( dexp(v/this%T_s) - one )*(v/this%T_s)**2
          !if( ranmar() <= fn/this%f_max )exit
      !enddo 
      index_nu = this%ln_nu1 + ( this%ln_nu2 - this%ln_nu1 ) * ranmar()
      !write(*,*)'emitter=', index_nu, this%E_low1, 10.D0**( index_nu ) * h_ev * 1.D-6, this%E_up1
      !write(*,*)'emitter=', one / ( this%ln_nu2 - this%ln_nu1 ), this%ln_nu2, this%ln_nu1
      v = 10.D0 **( index_nu ) * h_ev * 1.D-6
      this%Phot4k_CtrCF(1) = v 
      cos_theta = ranmar()
      sin_theta = dsqrt( one - cos_theta**2 )
      phi = twopi * ranmar()
      sin_phi = dsin( phi )
      cos_phi = dcos( phi )

      this%Vector_of_Momentum(1) = sin_theta * cos_phi
      this%Vector_of_Momentum(2) = sin_theta * sin_phi
      this%Vector_of_Momentum(3) = cos_theta 
      this%Phot4k_CtrCF(2: 4) = this%Phot4k_CtrCF(1) * this%Vector_of_Momentum(1:3)

      this%Phot4k_CovCF = this%Phot4k_CtrCF 

      this%w_ini_em = ( v / this%T_s )**3 / ( dexp( v / this%T_s ) - one )! * &
                      !(cos_theta + two * cos_theta**2)

      this%z_tau = this%z_max ! this%tau_max * ranmar()
      !write(*,*)'emitter=', 
      !this%medium_case = 1 
      !this%Psi_I = one
      !this%Psi_Q = one
      !this%Psi_U = one
      !this%Psi_V = one

      this%f4_CF(1) = zero
      this%f4_CF(2) = sin_phi
      this%f4_CF(3) = - cos_phi
      this%f4_CF(4) = zero
      
      end subroutine get_Phot4k_CtrCF_CovCF_BoundReflec_Sub

      end module PhotonEmitterBB





