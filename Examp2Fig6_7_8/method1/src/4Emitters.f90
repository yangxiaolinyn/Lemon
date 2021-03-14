      module PhotonEmitter 
      use Basic_Variables_And_Methods
      USE MPI 
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
          integer :: num_process, my_ID
          character*80 :: CrossSectFileName
          real(mcp) :: Scattered_Phot4k_In_Elec(1: 4), Temp_Matrix_1X3(1: 3)
        

      contains 
          procedure, public :: get_Phot4k_CtrCF_CovCF    =>  get_Phot4k_CtrCF_CovCF_Sub
          procedure, public :: max_plankexp    =>  max_plankexp_fn
          procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn
          procedure, public :: Set_Emin_Emax  =>  Set_Emin_Emax_Sub
          procedure, public :: Emitter_A_Photon   =>    Emitter_A_Photon_Sub 
      end type Photon_Emitter

      private :: get_Phot4k_CtrCF_CovCF_Sub
      private :: max_plankexp_fn
      private :: scatter_photon_blackbody_fn
      private :: Set_Emin_Emax_Sub 
      private :: Emitter_A_Photon_Sub
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
      subroutine Set_Emin_Emax_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, f1, a2, b2, dE, dnu
      integer :: i, istat
 
      this%E_max = this%max_plankexp(2) * this%T_s
      this%f_max = (this%E_max / this%T_s)**2/( dexp(this%E_max / this%T_s) - one ) 
      this%E_low1 = this%E_max * 0.9D0
      DO
          f1 = (this%E_low1 / this%T_s)**2/( dexp(this%E_low1 / this%T_s) - one )
          if(f1 < this%f_max * 1.D-7)exit
          this%E_low1 = this%E_low1 * 0.9D0
      ENDDO
      this%E_up1 = this%E_max * 1.1D0
      DO
          f1 = (this%E_up1 / this%T_s)**2/( dexp(this%E_up1 / this%T_s) - one )
          if(f1 < this%f_max * 1.D-7)exit
          this%E_up1 = this%E_up1 * 1.1D0
      ENDDO
      if( .false. )then
          this%ln_nu2 = 22.D0
          this%ln_nu1 = 11.D0
          this%nu_low = 10.D0**(11.D0)
          this%nu_up = 10.D0**22.D0
          this%E_low1 = this%nu_low * h_ev * 1.D-6
          this%E_up1 = this%nu_up * h_ev * 1.D-6
      else
          this%nu_low = this%E_low1*1.D6/h_ev
          this%nu_up = this%E_up1*1.D6/h_ev
          this%ln_nu2 = dlog10( this%nu_up )
          this%ln_nu1 = dlog10( this%nu_low )
      endif
      this%dnu = ( this%ln_nu2 - this%ln_nu1 ) / N_wt1
      !write(*, *)'mm=', this%E_low1*1.D6/h_ev, this%E_max*1.D6/h_ev, this%E_up1*1.D6/h_ev
      !write(*, *)'mm=', this%E_low1, this%E_max, this%E_up1
 
      end subroutine Set_Emin_Emax_Sub
 

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, theta, phi, w_ini, v, fn, nu, index_nu
      integer :: i
  
      index_nu = this%ln_nu1 + ranmar()*( this%ln_nu2 - this%ln_nu1 )
 
      v = 10.D0**( index_nu ) * h_ev * 1.D-6
      this%Phot4k_CtrCF_ini(1) = v

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%E_ini = v
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      theta = pi * ranmar()
      phi = twopi * ranmar()
 
      this%Vector_of_Momentum_ini(1) = DSIN(theta) * DCOS(phi)
      this%Vector_of_Momentum_ini(2) = DSIN(theta) * DSIN(phi)
      this%Vector_of_Momentum_ini(3) = DCOS(theta)

      this%Phot4k_CtrCF_ini(2:4) = this%Vector_of_Momentum_ini * this%E_ini

      this%Phot4k_CovCF_ini = this%Phot4k_CtrCF_ini
      this%Phot4k_CovCF_ini(1) = - this%Phot4k_CovCF_ini(1)

      this%w_ini_em = ( v / this%T_s )**3 / ( dexp( v / this%T_s ) - one )! * dsin(theta)
      this%w_ini = this%w_ini_em
 
      !write(*,*)'ss3===', this%w_ini, this%w_ini_em, dexp( v / this%T_s ), v
      !if(this%w_ini_em == zero)stop
      
      end subroutine get_Phot4k_CtrCF_CovCF_Sub


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon_Emitter) :: this
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r = zero  
      this%mucos = dcos( this%theta )
      this%musin = dsin( this%theta )
      this%phi = zero!twopi * ranmar()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%x = this%r * this%musin * dcos( this%phi )
      this%y = this%r * this%musin * dsin( this%phi )
      this%z = this%r * this%mucos
      this%Vector_of_position(1) = this%x
      this%Vector_of_position(2) = this%y
      this%Vector_of_position(3) = this%z
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      this%r_ini = zero 
      !this%mucos_ini  = ranmar()
      !this%musin_ini  = dsqrt( one - this%mucos_ini**2 )
      !this%phi_ini    = twopi * ranmar()
      !this%x_ini  = this%r_ini * this%musin_ini * dcos( this%phi_ini )
      !this%y_ini  = this%r_ini * this%musin_ini * dcos( this%phi_ini )
      !this%z_ini  = this%r_ini * this%mucos_ini
      this%Vector_of_position_ini = zero  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL this%get_Phot4k_CtrCF_CovCF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE Emitter_A_Photon_Sub

      end module PhotonEmitter





