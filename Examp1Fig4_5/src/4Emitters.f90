      module PhotonEmitter 
      use Basic_Variables_And_Methods
      use Integrations
      use CrossSection

      implicit none
          real(mcp), parameter :: c1 = 20.D0
          real(mcp), parameter :: c2 = 2.D0**(11.D0/6.D0)
          real(mcp), parameter :: c3 = 2.D0**(11.D0/12.D0) * 8.D0
          real(mcp), parameter :: c4 = 2.D0**(11.D0/12.D0)
          real(mcp), parameter :: prob1 = c1 / ( c1 + c2 + c3 )
          real(mcp), parameter :: prob2 = ( c1 + c2 ) / ( c1 + c2 + c3 )
          !integer, parameter :: N_wt = 3000

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter
          real(mcp) :: R_in
          real(mcp) :: R_out
          real(mcp) :: w_ini_em
          real(mcp) :: mag_B = 1.D0
          real(mcp) :: n_e 
          real(mcp) :: Big_Theta_e = 100.D0 
          real(mcp) :: j_Enu_theta
          !real(mcp) :: weights_y(1: N_wt) = zero
          real(mcp) :: Optical_Depth_absorption
          real(mcp), private :: Theta_e
          real(mcp) :: n_obs(1: 3)
          !real(mcp) :: sin_theta_obs
          !real(mcp) :: cos_theta_obs
          real(mcp) :: sin_phi_obs
          real(mcp) :: cos_phi_obs
          real(mcp) :: T_e
          real(mcp) :: phi_obs
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2

      contains 
          procedure, public :: get_Phot4k_CtrCF_CovCF    =>  get_Phot4k_CtrCF_CovCF_Sub
          procedure, public :: j_nu_theta_emissity_sampling  =>  j_nu_theta_emissity_sampling_Sub
          procedure, public :: j_theta_nu_emissity   =>   j_theta_nu_emissity_fn
          procedure, public :: j_theta_E_nu_emissity   =>   j_theta_E_nu_emissity_fn
          procedure, public :: get_Phot4k    =>  get_Phot4k_Sub 

      end type Photon_Emitter

      private :: get_Phot4k_CtrCF_CovCF_Sub
      private :: j_nu_theta_emissity_sampling_Sub
      private :: j_theta_nu_emissity_fn
      private :: get_Phot4k_Sub
      private :: j_theta_E_nu_emissity_fn 

      contains

!*******************************************************************************************************
      real(mcp) function chisquare(n)
!*******************************************************************************************************
      implicit none
      integer, intent(in) :: n
      integer :: i
     !***************************************************** 

      chisquare = zero
      i = 0
      do 
          chisquare = chisquare + GAUSSIAN1()**2
          i = i + 1
          if ( i == n) exit
      end do
      chisquare = dsqrt( chisquare / two )
      return
      end function chisquare

 
!*******************************************************************************************************
      function scatter_photon_blackbody(maxn, w_ini_em)
!*******************************************************************************************************
      implicit none
      real(mcp) :: scatter_photon_blackbody
      real(mcp) :: fn, maxfn, v, xn
      real(mcp),intent(in) :: maxn
      real(mcp),intent(out) :: w_ini_em
      real(mcp) ::  Tb = 2.D-3
      integer :: n = 2
 
      maxfn = (maxn)**2/( dexp(maxn) - one )        

      do
          v = ranmar()*(26.D-3 - 0.1D-3) + 0.1D-3 !(26D0/0.1D0)**ranmar()*0.1D0
          fn = one/( dexp(v/Tb) - one )*(v/Tb)**2
          if ( ranmar() <= fn/maxfn ) exit
      enddo 
      scatter_photon_blackbody = v
      !w_ini_em = one/( dexp(v/Tb) - one )*(v/Tb)**2 / maxfn
      w_ini_em = one !one/( dexp(v/Tb) - one ) * two * planck_h / Cv**2
      !log10_dv_em = ( dlog10(26.D-3 * 1.D6 / h_ev) - dlog10(0.1D-3 * 1.D6 / h_ev) ) / 500.D0
      !dv_em =  (26.D-3 - 0.1D-3)/ 1000.D0
      end function scatter_photon_blackbody

!*******************************************************************************************************
      function max_plankexp(n)
!******************************************************************************************************* 
      implicit none
      real(mcp) :: max_plankexp
      real(mcp) :: xStart = 1.D-3, xn, xn1, fn, fn1
      real(mcp) ::  Tb = 2.D-3
      integer, intent(in) :: n
      
      xn = xStart 
      do
          fn = n*( one - dexp(-xn) )
          xn1 = xn
          xn = fn 
          if(dabs(xn-xn1) < 1.D-10)exit
      enddo 
      max_plankexp = xn
      end function max_plankexp

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, theta, phi, w_ini, log10_dv, dv
      integer :: i

      i = 2
      maxn = max_plankexp(i)
      this%Phot4k_CtrCF(1) = scatter_photon_blackbody(maxn, w_ini)
      this%w_ini_em = w_ini
      !this%dv_em = dv
      theta = pi * ranmar()
      phi = twopi * ranmar()
      this%Phot4k_CtrCF(2) = this%Phot4k_CtrCF(1) * DSIN(theta) * DCOS(phi)
      this%Phot4k_CtrCF(3) = this%Phot4k_CtrCF(1) * DSIN(theta) * DSIN(phi)
      this%Phot4k_CtrCF(4) = this%Phot4k_CtrCF(1) * DCOS(theta)

      this%Phot4k_CovCF = this%Phot4k_CtrCF
      this%Phot4k_CovCF(1) = - this%Phot4k_CovCF(1)
      
      end subroutine get_Phot4k_CtrCF_CovCF_Sub

!*******************************************************************************************************
      subroutine j_nu_theta_emissity_sampling_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      !real(mcp), intent(in) :: T_e
      real(mcp) :: nu, phi
      real(mcp) :: nu0, nus, Big_Theta
      real(mcp) :: r1, r2, y1, a1, b1 
      real(mcp) :: del_ln_nu, nu1, nu2, del_nu
      real(mcp) :: a2, b2, del_ab, ln_nu
      integer :: i, n

      Big_Theta = this%T_e / mec2
      nu0 = two / nine * electron_Charge * this%mag_B / &
                 ( twopi * electron_mass * Cv) * Big_Theta**2 
 
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      ln_nu = this%ln_nu1 + ranmar() * ( this%ln_nu2 - this%ln_nu1 )  
      nu = 10.D0**( ln_nu ) 
      phi = ranmar() * twopi

      this%Phot4k_CtrCF_ini(1) = nu * h_ev * 1.D-6
      this%Vector_of_Momentum_ini(1) = this%sin_theta_obs * DCOS(phi)
      this%Vector_of_Momentum_ini(2) = this%sin_theta_obs * DSIN(phi)
      this%Vector_of_Momentum_ini(3) = this%cos_theta_obs
      this%Phot4k_CtrCF_ini(2:4) = this%Phot4k_CtrCF_ini(1) * &
                                   this%Vector_of_Momentum_ini(1: 3) 

      this%j_Enu_theta = this%j_theta_nu_emissity( nu, this%cos_theta_obs, this%sin_theta_obs )
      this%E_ini = this%Phot4k_CtrCF_ini(1) 
      this%w_ini_em = this%j_Enu_theta * nu
      this%w_ini = this%w_ini_em

      this%Optical_Depth_absorption = this%j_Enu_theta / &
          ( two * planck_h * ( this%E_ini * 1.D6 * erg_of_one_ev / planck_h )**3 / Cv**2 / &
                       ( dexp( this%E_ini / this%T_e ) - one ) ) 

      this%Phot4k_CovCF_ini = this%Phot4k_CtrCF_ini
      this%Phot4k_CovCF_ini(1) = - this%Phot4k_CovCF_ini(1)   
      end subroutine j_nu_theta_emissity_sampling_Sub
 
!*******************************************************************************************************
      real(mcp) function j_theta_nu_emissity_fn( this, nu, costheta, sintheta )
!******************************************************************************************************* 
      implicit none 
      class(Photon_Emitter) :: this
      real(mcp), intent(in) :: nu, costheta, sintheta
      real(mcp) :: nu0, nus, Big_Theta, K2_Theta
      real(mcp) :: X, coef_A
      integer :: i, n

      Big_Theta = this%T_e / mec2
      nu0 = two / nine * electron_Charge * this%mag_B / &
                 ( twopi * electron_mass * Cv) * this%Big_Theta_e**2
      K2_Theta = bessk(2, one / Big_Theta)
      nus = nu0 * sintheta
      coef_A = sqrt2 * pi * electron_Charge**2 * this%n_e * nus / three / Cv / K2_Theta 
      X = nu / nus
       
      j_theta_nu_emissity_fn = coef_A * ( dsqrt(X) + c4 * &
                       X**( one / 6.D0 ) )**2 * dexp( - X**(one/three) )
      !write(*, *)'sss====', i, j, j_nu, nu, costheta, sintheta
      end function j_theta_nu_emissity_fn

!*******************************************************************************************************
      real(mcp) function j_theta_E_nu_emissity_fn( this, E_nu, costheta, sintheta )
!******************************************************************************************************* 
      implicit none 
      class(Photon_Emitter) :: this
      real(mcp), intent(in) :: E_nu, costheta, sintheta
      real(mcp) :: nu0, nus, Big_Theta, K2_Theta
      real(mcp) :: X, coef_A, nu
      integer :: i, n

      Big_Theta = this%T_e / mec2
      nu0 = two / nine * electron_Charge * this%mag_B / &
                 ( twopi * electron_mass * Cv) * this%Big_Theta_e**2
      K2_Theta = bessk(2, one / Big_Theta)
      nus = nu0 * sintheta
      coef_A = sqrt2 * pi * electron_Charge**2 * this%n_e * nus / three / Cv / K2_Theta
      nu = E_nu / h_ev * 1.D6 
      X = nu / nus
       
      j_theta_E_nu_emissity_fn = coef_A * ( dsqrt(X) + c4 * &
                       X**( one / 6.D0 ) )**2 * dexp( - X**(one/three) )
      !write(*, *)'sss====', nu, E_nu
      end function j_theta_E_nu_emissity_fn

!*******************************************************************************************************
      subroutine get_Phot4k_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, theta, phi, w_ini, log10_dv, dv
      integer :: i

      i = 2
      maxn = max_plankexp(i)
      this%Phot4k_CtrCF(1) = scatter_photon_blackbody(maxn, w_ini)
      this%w_ini_em = w_ini
      !this%dv_em = dv
      theta = pi * ranmar()
      phi = twopi * ranmar()
      this%Phot4k_CtrCF(2) = this%Phot4k_CtrCF(1) * DSIN(theta) * DCOS(phi)
      this%Phot4k_CtrCF(3) = this%Phot4k_CtrCF(1) * DSIN(theta) * DSIN(phi)
      this%Phot4k_CtrCF(4) = this%Phot4k_CtrCF(1) * DCOS(theta)

      this%Phot4k_CovCF = this%Phot4k_CtrCF
      this%Phot4k_CovCF(1) = - this%Phot4k_CovCF(1)
      
      end subroutine get_Phot4k_Sub

      end module PhotonEmitter





