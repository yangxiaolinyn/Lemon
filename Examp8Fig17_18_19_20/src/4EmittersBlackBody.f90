      module PhotonEmitterBB
      use Basic_Variables_And_Methods
      implicit none 
      real(mcp), parameter :: nu2MeV = h_ev * 1.D-6

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon_Emitter_BB 
          real(mcp) :: w_ini_em
          real(mcp) :: w_ini0
          real(mcp) :: T_s
          real(mcp) :: max_nu
          real(mcp) :: E_low1
          real(mcp) :: E_up1
          real(mcp) :: E_max
          real(mcp) :: f_max 
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: y_em1
          real(mcp) :: y_em2
          real(mcp) :: dy_em

      contains   
          procedure, public :: max_plankexp    =>  max_plankexp_fn 
          procedure, public :: Set_Emin_Emax  =>  Set_Emin_Emax_Sub 
          procedure, public :: get_Phot4k_CtrCF_CovCF   =>  &
                               get_Phot4k_CtrCF_CovCF_Sub 
      end type Photon_Emitter_BB
  
      private :: max_plankexp_fn 
      private :: Set_Emin_Emax_Sub 
      private :: get_Phot4k_CtrCF_CovCF_Sub

      contains 

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
      subroutine Set_Emin_Emax_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: maxn, f1, a2, b2, dE 
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
  
      this%nu_low = this%E_low1 * 1.D6 / h_ev
      this%nu_up = this%E_up1 * 1.D6 / h_ev
      this%y_em2 = dlog10( this%nu_up )
      this%y_em1 = dlog10( this%nu_low )
      this%dy_em = this%y_em2 - this%y_em1
  
      end subroutine Set_Emin_Emax_Sub
  

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter_BB) :: this
      real(mcp) :: phi, v, index_nu
      real(mcp) :: sin_theta, cos_theta, sin_phi, cos_phi  
 
      index_nu = this%y_em1 + this%dy_em * ranmar()
 
      v = 10.D0 **( index_nu ) * nu2MeV 
      this%E_ini = v

      this%Phot4k_CtrCF_ini(1) = v 
      cos_theta = ranmar()
      sin_theta = dsqrt( one - cos_theta**2 )
      phi = twopi * ranmar()
      sin_phi = dsin( phi )
      cos_phi = dcos( phi )

      this%Vector_of_Momentum_ini(1) = sin_theta * cos_phi
      this%Vector_of_Momentum_ini(2) = sin_theta * sin_phi
      this%Vector_of_Momentum_ini(3) = cos_theta 
      this%Phot4k_CtrCF_ini(2: 4) = this%Phot4k_CtrCF_ini(1) * &
                                    this%Vector_of_Momentum_ini(1:3)

      !this%Phot4k_CovCF_ini = this%Phot4k_CtrCF_ini 

      this%w_ini_em = ( v / this%T_s )**3 / ( dexp( v / this%T_s ) - one ) * &
                      (cos_theta + two * cos_theta**2)

      this%w_ini = this%w_ini_em 
      this%w_ini0 = this%w_ini_em

      this%z_tau = this%z_max  
      
      end subroutine get_Phot4k_CtrCF_CovCF_Sub

      end module PhotonEmitterBB





