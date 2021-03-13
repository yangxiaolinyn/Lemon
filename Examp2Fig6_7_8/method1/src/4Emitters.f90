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
          real(mcp) :: weights(1: N_wt1)
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          real(mcp) :: ln_nu1
          real(mcp) :: ln_nu2
          real(mcp) :: dnu
        

      contains 
          procedure, public :: get_Phot4k_CtrCF_CovCF    =>  get_Phot4k_CtrCF_CovCF_Sub
          procedure, public :: max_plankexp    =>  max_plankexp_fn
          procedure, public :: scatter_photon_blackbody  =>  scatter_photon_blackbody_fn
          procedure, public :: Set_Emin_Emax  =>  Set_Emin_Emax_Sub
          procedure, public :: Integration3   =>  Integration3_fn
          !procedure, public :: get_Phot4k_BL    =>  get_Phot4k_BL_Sub
          !procedure, public :: get_Phot_Covariant4k_BL    =>  get_Phot_Covariant4k_BL_Sub
          !procedure, public :: Set_Phot4k_BL    =>  Set_Phot4k_BL_Sub
          !procedure, public :: get_Phot4k_CovBL =>  get_Phot4k_CovBL_Sub 
          !procedure, public :: Get_Phot4k_inCF_From_4k_inBL => Get_Phot4k_inCF_From_4k_inBL_Sub
          !procedure, public :: temp1

      end type Photon_Emitter

      private :: get_Phot4k_CtrCF_CovCF_Sub
      private :: max_plankexp_fn
      private :: scatter_photon_blackbody_fn
      private :: Set_Emin_Emax_Sub
      private :: Integration3_fn
      !private :: Get_Phot4k_inCF_From_4k_inBL_Sub
      !private :: get_Phot4k_CovBL_Sub

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
      if(.false.)then
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
      write(*, *)'mm=', this%E_low1*1.D6/h_ev, this%E_max*1.D6/h_ev, this%E_up1*1.D6/h_ev
      write(*, *)'mm=', this%E_low1, this%E_max, this%E_up1
      !stop
      open(unit=11, file = './data/weights.txt', status = "old", &
                                    action = "read", iostat = istat)  
      if(istat == 0)then
          write(unit = *, fmt = *)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'  
          write(unit = *, fmt = *)' Now reading the weights from the data file..... ' 
          write(unit = *, fmt = *)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
          do i = 1, N_wt1
              read(unit = 11, fmt = *)this%weights(i)
          enddo
      else 
          open(unit=12, file = './data/weights.txt', status = "replace", &
                                   action = "write", iostat = istat) 
          if(istat == 0)then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the value of weigts into a data File......'  
              write(unit = *, fmt = *)'************************************************************'
              do i = 1, N_wt1 
                  b2 = 10.D0**( this%ln_nu1 + this%dnu * i ) * h_ev * 1.D-6 / this%T_s
                  a2 = 10.D0**( this%ln_nu1 + this%dnu * ( i - 1 ) ) * h_ev * 1.D-6 / this%T_s  
                  this%weights(i) = this%Integration3( a2, b2, 1700000 ) 
                  write(unit = 12, fmt = *)this%weights(i)
                  write(unit = *, fmt = *)i, this%weights(i)
              enddo
          else
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)' wrongss !!!......'  
              write(unit = *, fmt = *)'************************************************************'
              stop
          endif
      endif 
      close(unit=11)
      close(unit=12)
      !write(*, *)'mm=', this%E_low, this%E_max, this%E_up
      
      end subroutine Set_Emin_Emax_Sub

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration3_fn( this, a, b, n )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none 
      class(Photon_Emitter) :: this
      real(mcp), intent(in) :: a, b
      real(mcp) :: Del_y, ints, y, Del_ints, Del_ints_0, Del_ints_n
      integer, intent(in) :: n
      integer :: i
      real(mcp), parameter :: c1 = 20.D0
      real(mcp), parameter :: c2 = 2.D0**(11.D0/6.D0)
      real(mcp), parameter :: c3 = 2.D0**(11.D0/12.D0) * 8.D0
      real(mcp), parameter :: c4 = 2.D0**(11.D0/12.D0)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          ints = zero  
          i = 0 
          Del_y = ( b - a ) / (n + 1) 
          y = a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Del_ints = y**2/( dexp(y) - one )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              ints = ints + Del_ints    
              If ( i == 0 ) then 
                  Del_ints_0 = Del_ints
              Endif  
              If ( i > n ) then
                  Del_ints_n = Del_ints
                  exit
              Endif  
              i = i + 1  
              y = Del_y + y
              !write(*,*)'sd', y, Del_y, ints
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          Integration3_fn = ( ints - Del_ints_0 / two - Del_ints_n / two )*Del_y
      end function Integration3_fn

!*******************************************************************************************************
      subroutine get_Phot4k_CtrCF_CovCF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_Emitter) :: this
      real(mcp) :: maxn, theta, phi, w_ini, v, fn, nu, index_nu
      integer :: i

      !do
          !index_nu = this%ln_nu1 + ranmar()*( this%ln_nu2 - this%ln_nu1 )
          !v = 10.D0**( index_nu ) * h_ev * 1.D-6
          !v = !ranmar() * ( this%E_up1 - this%E_low1 ) + this%E_low1
          !fn = one/( dexp(v/this%T_s) - one )*(v/this%T_s)**2
          !if( ranmar() <= fn/this%f_max )exit
      !enddo 
      index_nu = this%ln_nu1 + ranmar()*( this%ln_nu2 - this%ln_nu1 )
      !write(*,*)'emitter=', index_nu,10.D0**( index_nu )
      !write(*, *)'mm1=', this%ln_nu1, this%ln_nu2 * 4.D8
      v = 10.D0**( index_nu ) * h_ev * 1.D-6
      this%Phot4k_CtrCF(1) = v
      theta = pi * ranmar()
      phi = twopi * ranmar()
      this%Phot4k_CtrCF(2) = this%Phot4k_CtrCF(1) * DSIN(theta) * DCOS(phi)
      this%Phot4k_CtrCF(3) = this%Phot4k_CtrCF(1) * DSIN(theta) * DSIN(phi)
      this%Phot4k_CtrCF(4) = this%Phot4k_CtrCF(1) * DCOS(theta)

      this%Phot4k_CovCF = this%Phot4k_CtrCF
      this%Phot4k_CovCF(1) = - this%Phot4k_CovCF(1)
      Do i = 1, N_wt1
          If( this%ln_nu1 + this%dnu * i > index_nu )then 
              this%w_ini_em = this%weights(i) 
              exit
          endif
      ENDDO
      
      end subroutine get_Phot4k_CtrCF_CovCF_Sub



      end module PhotonEmitter





