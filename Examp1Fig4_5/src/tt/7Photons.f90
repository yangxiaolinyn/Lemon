      module Photons 
      use ScatterPhoton
      use CrossSection 
      implicit none
      integer, parameter :: N_sigma = 1000
      !integer, parameter :: N_H3 = 4000
      !real(mcp), dimension(0:N_sigma) :: sigmaaTeE
      !real(mcp), dimension(0:N_H3) :: H3Array

      type, public, extends(ScatPhoton) :: Photon
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE
          real(mcp), dimension(0:400, 0:400) :: disk_image = zero, disk_image_Recv
          integer(kind=8) :: effect_number
          integer(kind=8) :: scatter_times
          real(mcp) :: n_e
          real(mcp) :: time_arrive_observer
          real(mcp) :: r_obs
          real(mcp) :: R_in
          real(mcp) :: R_out
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(0:100) :: delta_pds(0:100)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          logical :: Go2infinity
          logical :: fall2BH
          logical :: At_outer_Shell
          logical :: At_inner_Shell
          logical :: direct_escaped
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: A_normal
          real(mcp) :: Sigma_Max
          real(mcp) :: p_out
          real(mcp) :: p_boundary1
          real(mcp) :: p_boundary2
          real(mcp) :: p_boundary3
          real(mcp) :: p_maxs
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          logical :: fall2BHs, escapteds
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          real(mcp) :: Red_Shift_g    ! Red_Shift_g = Phto_E_CF_end / Phto_E_CF_ini
          real(mcp) :: Phot_E_CF_ini  ! Phto_E_CF_ini = P_{\mu} * U^{\mu}(ini)
          real(mcp) :: Phot_E_CF_end  ! Phto_E_CF_end = P_{\mu} * U^{\mu}(end)
          real(mcp) :: logE_low
          real(mcp) :: logE_up
          real(mcp) :: dindexE
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          !real(mcp) :: K1_pw      !pw means Penrose & Walker constant, its real part is K1
          !real(mcp) :: K2_pw      !pw means Penrose & Walker constant, its image part is K2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: frequency_v
          real(mcp) :: v_L_v_i(1:500)
          real(mcp) :: v_L_v_i_ET(1:500, 1:200)
          real(mcp) :: Optical_Depth_scatter
          !real(mcp) :: log10_dv_ini
          !real(mcp) :: log10_dv_p
          !real(mcp) :: dv_ini
          !real(mcp) :: dv_p
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: R_gc ! Radius of the electron gas cloud
          real(mcp) :: T_e ! = 1.D-3  ! Temparature of the electron gas cloud
          !real(mcp) :: R_out, R_in

      contains 
!*******************************************************************************************************
          procedure, public :: get_Tau                           =>  get_Tau_fn
          procedure, public :: p_scatter_distance
          procedure, public :: sigma_fn
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Set_Cross_Section                 => Set_Cross_Section_sub
          !procedure, public :: Get_Phot4k_CF_From_ContraVariant4k_BL     => &
          !                     Get_Phot4k_CF_From_ContraVariant4k_BL_Sub
          !procedure, public :: Get_Phot4k_CF_From_CoVariant4k_BL         => &
          !                     Get_Phot4k_CF_From_CoVariant4k_BL_Sub
          !procedure, public :: Contra_Variant_Vector_inner_Product       => &
          !                     Contra_Variant_Vector_inner_Product_Sub
          !procedure, public :: get_f4_CtrBL_From_f4_CF     =>  get_f4_CtrBL_From_f4_CF_Sub
          !procedure, public :: get_f4_CovBL_From_f4_CovCF  =>  get_f4_CovBL_From_f4_CovCF_Sub
          !procedure, public :: Get_Ctrf4_CF_From_Ctrf4_BL  => Get_Ctrf4_CF_From_Ctrf4_BL_Sub
          !procedure, public :: Get_Ctr4Vec_CF_From_Ctr4Vec_BL  => Get_Ctr4Vec_CF_From_Ctr4Vec_BL_Sub
          !procedure, public :: Get_Cov4Vec_CF_From_Cov4Vec_BL  => Get_Cov4Vec_CF_From_Cov4Vec_BL_Sub
          procedure, public :: Calc_Phot_Polarization_At_Observor  => &
                                Calc_Phot_Polarization_At_Observor_Sub
!*******************************************************************************************************
          !procedure, public :: Set_Elec4U_CtrBL_At_ini         => Set_Elec4U_CtrBL_At_ini_Sub
          procedure, public :: Get_bias_parameter              => Get_bias_parameter_Sub
!*******************************************************************************************************
          procedure, public :: Get_Optical_Depth_At_p  =>      Get_Optical_Depth_At_p_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Get_Max_Value_of_Sigma  =>  Get_Max_Value_of_Sigma_Sub
          procedure, public :: r_p2
          procedure, public :: get_Tau2                           =>  get_Tau2_fn
          procedure, public :: Get_scatter_distance
          procedure, public :: Get_scatter_distance2
          procedure, public :: Calc_Phot_Informations_At_Observor =>  &
                                Calc_Phot_Informations_At_Observor_Sub
          procedure, public :: Get_Max_Value_of_Sigma_2zones => &
                                Get_Max_Value_of_Sigma_2zones_Sub
          procedure, public :: Calc_Phot_Informations_At_Observor_2zones => &
                                Calc_Phot_Informations_At_Observor_2zones_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon
 
      private :: get_Tau_fn, Set_Cross_Section_sub
      !private :: Set_Phot_BLCoordinates_At_p_Sub
      !private :: Get_Phot4k_CF_From_ContraVariant4k_BL_Sub
      !private :: Get_Phot4k_CF_From_CoVariant4k_BL_Sub
      !private :: Contra_Variant_Vector_inner_Product_Sub
      !private :: get_f4_CtrBL_From_f4_CF_Sub
      !private :: Get_Ctrf4_CF_From_Ctrf4_BL_Sub
      !private :: Get_Ctr4Vec_CF_From_Ctr4Vec_BL_Sub
      !private :: Get_Cov4Vec_CF_From_Cov4Vec_BL_Sub
      private :: Calc_Phot_Polarization_At_Observor_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !private :: p_scatter_distance
      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !private :: Set_Elec4U_CtrBL_At_ini_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Get_Optical_Depth_At_p_Sub
      private :: Get_bias_parameter_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 

!*******************************************************************************************************
      !subroutine get_Phot4k_CF_Sub(this)
!*******************************************************************************************************
      !implicit none
      !class(Photon) :: this 

      !call Matrix_Multiplication14X44_Sub( this%Phot4k_CF, &
      !                  this%Lorentz_Matrix, this%Phot4k_LNRF )
       !write(*,*)'Hello My friend'
      !end subroutine get_Phot4k_CF_Sub
!*******************************************************************************************************
      real(mcp) function p_scatter_distance2( this, T_e )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: T_e
      !TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr
      logical :: p_lt_pmax
      integer(kind=8) :: i
 
      call this%Set_Photon_Tetrad_In_CF()
      !func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( R_max**2+ this%aspin**2*this%mucos**2 ) * this%E
      func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( this%R_out**2 + &
                                  this%aspin**2 ) * this%E_ini
      !write(*,*)'p_scatter_distance, func1_tau_max_value === ', func1_tau_max_value
      p_max = p_total(this%Phot4k_LNRF_ini(2), this%lambda_ini, this%q_ini, &
                       this%musin_ini, this%mucos_ini, this%aspin, this%r_ini, one)
      this%p_maxs = p_max
      tempA = one - DEXP(-p_max*func1_tau_max_value)
      !write(*,*)'p_total === ', p_max
      i = 0 
      DO
        Do
          p = zero
          p_lt_pmax = .TRUE.
          do
              r1 = ranmar()
              p1 = - DLOG( one - r1*tempA ) / func1_tau_max_value
              If ( p1 < zero )then
                  write(*,*)' p1 < 0', p1, one - r1*tempA, tempA, func1_tau_max_value, p_max
                  stop
              endif
              p = p + p1
              rp = radius(p,this%Phot4k_LNRF_ini(2), this%lambda_ini, &
                               this%q_ini, this%aspin, this%r_ini, one, sign_pr)  
              !write(*,*)'p, rp, R_max, p_max=',p ,func1_tau_max_value, tempA !rp, R_max, p_max
              !if (p >= p_max) then
              if ( rp >= this%R_out .OR. rp <= this%R_in ) then
                  p_lt_pmax = .FALSE.
                  exit
              endif 
              Sigma_atP = this%get_Tau(p, T_e)
              !write(*,*)'p_scatter_distance, p, Sigma_P === ', p, Sigma_atP, func1_tau_max_value, &
              !                    Sigma_atP/func1_tau_max_value
              rprobability = Sigma_atP / func1_tau_max_value
              r3 = ranmar()
              !write(*,*)'p_scatter_distance, p, Sigma_P === ',r3,rprobability
              !write(*,*)' r1 r2 === ',r3,rprobability
              if ( r3 <=  rprobability ) exit 
          enddo
          if (p_lt_pmax) exit
          i = i + 1
          !write(*,*)'i ======= ', i, p_max
        Enddo 
        !write(*,*)'rr probobility===',dexp( -( one/tempA - one ) * func1_tau_max_value * p )
        If ( ranmar() <= dexp( -( one/tempA - one ) * func1_tau_max_value * p ) ) exit
      ENDDO
      p_scatter_distance2 = p
      end function p_scatter_distance2

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance( this, T_e )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: T_e
      !TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i
 
      !func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( R_max**2+ this%aspin**2*this%mucos**2 ) * this%E
      !func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( (this%R_out*1.D0)**2 + &
      !                            this%aspin**2 ) * this%E_ini
      !write(*,*)'p_scatter_distance, func1_tau_max_value === ', func1_tau_max_value
      !p_max = p_total(this%Phot4k_LNRF_ini(2), this%lambda_ini, this%q_ini, &
      !                 this%musin_ini, this%mucos_ini, this%aspin, this%r_ini, one) 
      p_max = this%R_in * 1.2D0
      this%p_maxs = p_max 
      call this%Get_Max_Value_of_Sigma( T_e ) 
      !call this%Get_bias_parameter( T_e ) ! Through this call, Sigma_Max is obtained!
      func1_tau_max_value = this%Sigma_Max * 1.D0
      p_out = this%p_out/1.D0
      !p_out = dlog( 100.D0 / 99.D0 ) / this%Sigma_Max
      !p_out = p_max 
      tempA = one - DEXP( - p_out * func1_tau_max_value )
      this%A_normal = one / tempA
      !write(*,*)'p_total === ', p_out, func1_tau_max_value, this%A_normal
      i = 0 
      loop_stop = .False. 
        !Do
          p = zero
          p_lt_pmax = .TRUE.
          do
              r1 = ranmar()
              !p1 = - DLOG( one - r1*tempA ) / func1_tau_max_value
              p1 = - DLOG( one - r1 ) / func1_tau_max_value
              If ( p1 < zero )then
                  write(*,*)' p1 < 0', p1, one - r1*tempA, tempA, func1_tau_max_value, p_max
                  stop
              endif
              p = p + p1
              !write(*,*)' p1 < 0',this%p_out, p, p1
              !rp = radius(p,this%Phot4k_LNRF_ini(2), this%lambda_ini, &
              !                 this%q_ini, this%aspin, this%r_ini, one, sign_pr)  
              rp = this%r_p2( this%Vector_of_position_ini, this%Vector_of_Momentum_ini, p )
              !if ( isnan( dsqrt(rp - 1.D8) + 12.D0 ) ) write(*,*)'sdfsdf===',rp
              if ( isnan(rp) ) then
                  write(*,*)'rp is NaN', rp
                  stop
              endif
              !write(*,*)'p, rp, R_max, p_max=',p ,func1_tau_max_value, tempA !rp, R_max, p_max
              !if (p >= p_max) then
              !write(*,*)' p2 < 0',rp, this%R_in, func1_tau_max_value
              if ( rp > this%R_in )then
                  this%At_outer_Shell = .True.
                  p_lt_pmax = .FALSE.
                  loop_stop = .True.
                  exit
              endif

              !If ( rp < this%R_in ) then
              !    this%At_inner_Shell = .TRUE.
              !    p_lt_pmax = .FALSE.
              !    loop_stop = .True.
              !    exit
              !endif
              !write(*,*)' p3 < 0',loop_stop
              Sigma_atP = this%get_Tau2(p, T_e)
              !write(*,*)'p_scatter_distance, p, Sigma_P === ', p, Sigma_atP, func1_tau_max_value, &
              !                    Sigma_atP/func1_tau_max_value
              rprobability = Sigma_atP / func1_tau_max_value
              r3 = ranmar()
              !write(*,*)'p_scatter_distance, p, Sigma_P === ',r3,rprobability
              !write(*,*)' p2 < 0',rp, this%R_in
              !write(*,*)' r1 r2 === ',r3,rprobability
              if ( r3 <=  rprobability ) exit 
              i = i + 1
              if (i > 1.D3) write(*,*)i,r3,rprobability,rp,this%R_in 
          enddo
          !if (p_lt_pmax) exit
          !if (loop_stop) exit
          !write(*,*)'i ======= ', i, p_max,p
        !write(*,*)'rr probobility===',dexp( -( one/tempA - one ) * func1_tau_max_value * p ),&
        !     ( one/tempA - one ) * func1_tau_max_value * p, &
        !     ( one/tempA - one ),func1_tau_max_value*p
        !If ( ranmar() <= dexp( -( one/tempA - one ) * func1_tau_max_value * p ) ) exit
        !If ( loop_stop ) exit
      !ENDDO
      Get_scatter_distance = p
      end function Get_scatter_distance


!*******************************************************************************************************
      real(mcp) function Get_scatter_distance2( this, T_e )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: T_e
      !TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out, r_times_p
      real(mcp) :: p_c, p_out1, p_1, p_2, delta
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i, path_cases
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  First ditermine which zone the photon is included in.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - r_times_p + dsqrt( r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )
      delta = r_times_p**2 - ( this%r_ini**2 - this%R_in**2 )
      If ( ( this%R_in < this%r_ini .AND. this%r_ini < this%R_out ) .OR. &
           ( this%r_ini == this%R_out .AND. r_times_p < zero )      .OR. &
           ( this%At_inner_Shell .AND. r_times_p > zero )            &
      ) then
          !if( this%r_ini == this%R_out .AND. r_times_p < zero )write(*,*)'The first time scattering!'
          If ( (delta < zero) .or. this%At_inner_Shell ) then 
              this%p_maxs = p_out1
              path_cases = 1
              this%p_boundary1 = p_out1 
              call this%Get_Max_Value_of_Sigma_2zones( T_e, p_out1, path_cases ) 
              func1_tau_max_value = this%Sigma_Max * 1.D0
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              i = 0    
              p = zero
              Do
                  r1 = ranmar() 
                  p1 = - DLOG( one - r1 ) / func1_tau_max_value
                  If ( p1 < zero )then
                      write(*,*)' p1 < 0', p1,r1, func1_tau_max_value, p_max
                      stop
                  endif 
                  p = p + p1 
                  if ( p > p_out1 )then
                      this%At_outer_Shell = .True. 
                      exit
                  endif 
                  Sigma_atP = this%get_Tau2(p, T_e) 
                  rprobability = Sigma_atP / func1_tau_max_value
                  r3 = ranmar() 
                  if ( r3 <=  rprobability ) then
                      this%At_inner_Shell = .False. 
                      exit 
                  endif
                  i = i + 1
                  if (i > 1.D3) write(*,*)i,r3,rprobability,rp,this%R_in 
              Enddo
              Get_scatter_distance2 = p
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
              p_1 = - r_times_p - dsqrt( delta )
              p_2 = - r_times_p + dsqrt( delta ) 
              this%p_maxs = p_1
              path_cases = 2
              this%p_boundary1 = p_1 
              this%p_boundary2 = p_2 
              this%p_boundary3 = p_out1 
              call this%Get_Max_Value_of_Sigma_2zones( T_e, p_1, path_cases ) 
              func1_tau_max_value = this%Sigma_Max * 1.D0
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              i = 0    
              p = zero
              Do
                  r1 = ranmar() 
                  p1 = - DLOG( one - r1 ) / func1_tau_max_value
                  If ( p1 < zero )then
                      write(*,*)' p1 < 0', p1,r1, func1_tau_max_value, p_max
                      stop
                  endif 
                  p = p + p1 
                  if ( p > p_1 )then
                      this%At_inner_Shell = .True. 
                      Get_scatter_distance2 = p_1
                      exit
                  endif 
                  Sigma_atP = this%get_Tau2(p, T_e) 
                  rprobability = Sigma_atP / func1_tau_max_value
                  r3 = ranmar() 
                  if ( r3 <=  rprobability ) then
                      this%At_inner_Shell = .False. 
                      Get_scatter_distance2 = p
                      exit
                  endif
                  i = i + 1
                  if (i > 1.D3) write(*,*)i,r3,rprobability,rp,this%R_in 
              Enddo
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
      else If ( ( this%r_ini < this%R_in ) .OR. &
                !( this%r_ini == this%R_in .AND. r_times_p < zero )      &
                 ( this%At_inner_Shell .AND. r_times_p < zero )    &
      ) then 
          p_1 = - r_times_p + dsqrt( delta ) 
          this%p_maxs = p_1
          path_cases = 3
          this%p_boundary1 = p_1
          this%p_boundary2 = p_out1
          call this%Get_Max_Value_of_Sigma_2zones( T_e, p_1, path_cases ) 
          func1_tau_max_value = this%Sigma_Max * 1.D0
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          i = 0    
          p = zero
          Do
              r1 = ranmar() 
              p1 = - DLOG( one - r1 ) / func1_tau_max_value
              If ( p1 < zero )then
                  write(*,*)' p1 < 0', p1,r1, func1_tau_max_value, p_max
                  stop
              endif 
              p = p + p1 
              if ( p > p_1 )then
                  this%At_inner_Shell = .True. 
                  Get_scatter_distance2 = p_1
                  exit
              endif 
              Sigma_atP = this%get_Tau2(p, T_e) 
              rprobability = Sigma_atP / func1_tau_max_value
              r3 = ranmar() 
              if ( r3 <=  rprobability ) then
                  this%At_inner_Shell = .False. 
                  Get_scatter_distance2 = p
                  exit
              endif
              i = i + 1
              if (i > 1.D3) write(*,*)i,r3,rprobability,rp,this%R_in 
          Enddo
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !Get_scatter_distance2 = p
       !write(*,*)this%Sigma_Max,Get_scatter_distance2
      end function Get_scatter_distance2

!*******************************************************************************************************
      real(mcp) function p_scatter_distance( this, T_e )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: T_e
      !TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i
 
      !func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( R_max**2+ this%aspin**2*this%mucos**2 ) * this%E
      !func1_tau_max_value = this%n_e * Sigma_T * rg_SUN * ( (this%R_out*1.D0)**2 + &
      !                            this%aspin**2 ) * this%E_ini
      !write(*,*)'p_scatter_distance, func1_tau_max_value === ', func1_tau_max_value
      p_max = p_total(this%Phot4k_LNRF_ini(2), this%lambda_ini, this%q_ini, &
                       this%musin_ini, this%mucos_ini, this%aspin, this%r_ini, one)
      this%p_maxs = p_max
      call this%Get_bias_parameter( T_e )
      func1_tau_max_value = this%Sigma_Max * 1.D0
      p_out = this%p_out/1.D0
      !p_out = dlog( 100.D0 / 99.D0 ) / this%Sigma_Max
      !p_out = p_max 
      tempA = one - DEXP( - p_out * func1_tau_max_value )
      this%A_normal = one / tempA
      !write(*,*)'p_total === ', p_out, func1_tau_max_value, this%A_normal
      i = 0 
      loop_stop = .False.
      DO
        Do
          p = zero
          p_lt_pmax = .TRUE.
          do
              r1 = ranmar()
              p1 = - DLOG( one - r1*tempA ) / func1_tau_max_value
              If ( p1 < zero )then
                  write(*,*)' p1 < 0', p1, one - r1*tempA, tempA, func1_tau_max_value, p_max
                  stop
              endif
              p = p + p1
              !write(*,*)' p1 < 0',this%p_out, p, p1
              rp = radius(p,this%Phot4k_LNRF_ini(2), this%lambda_ini, &
                               this%q_ini, this%aspin, this%r_ini, one, sign_pr)  
              !write(*,*)'p, rp, R_max, p_max=',p ,func1_tau_max_value, tempA !rp, R_max, p_max
              !if (p >= p_max) then
              !write(*,*)' p2 < 0',rp, this%R_out
              if ( rp > this%R_out )then
                  this%At_outer_Shell = .True.
                  p_lt_pmax = .FALSE.
                  loop_stop = .True.
                  exit
              endif
              If ( rp < this%R_in ) then
                  this%At_inner_Shell = .TRUE.
                  p_lt_pmax = .FALSE.
                  loop_stop = .True.
                  exit
              endif
              !write(*,*)' p3 < 0',loop_stop
              Sigma_atP = this%get_Tau(p, T_e)
              !write(*,*)'p_scatter_distance, p, Sigma_P === ', p, Sigma_atP, func1_tau_max_value, &
              !                    Sigma_atP/func1_tau_max_value
              rprobability = Sigma_atP / func1_tau_max_value
              r3 = ranmar()
              !write(*,*)'p_scatter_distance, p, Sigma_P === ',r3,rprobability
              !write(*,*)' r1 r2 === ',r3,rprobability
              if ( r3 <=  rprobability ) exit 
          enddo
          !if (p_lt_pmax) exit
          !if (loop_stop) exit
          i = i + 1
          !write(*,*)'i ======= ', i, p_max
          exit
        Enddo 
        !write(*,*)'rr probobility===',dexp( -( one/tempA - one ) * func1_tau_max_value * p ),&
        !     ( one/tempA - one ) * func1_tau_max_value * p, &
        !     ( one/tempA - one ),func1_tau_max_value*p
        If ( ranmar() <= dexp( -( one/tempA - one ) * func1_tau_max_value * p ) ) exit
        If ( loop_stop ) exit
      ENDDO
      p_scatter_distance = p
      end function p_scatter_distance

!*******************************************************************************************************
      real(mcp) function get_Tau_fn(this, samp, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: samp, T_e
      real(mcp) :: Sigmap, mup, musinp, rp, phip, timep, &
                 sigp, P_mu_Multipl_U_mu, PhotE_In_Elec_CF
   
      !CALL YNOGK(p,this%Phot4k_LNRF, this%lambda, this%q, this%musin, this%mucos, &
      !                 this%aspin, this%r, one, rp, mup, phip, timep, sigp) 
      !write(*,*)'Ra mu ,22222 === ', rp, mup, phip, timep, sigp 

      CALL this%Get_Phot_BLCoor_At_samp( samp )
      Sigmap = this%r_samp**2 + ( this%aspin * this%mucos_samp )**2
      !write(*,*)'Sigma === ', Sigmap, this%r_p, this%mucos_p, this%musin_p, this%phi_p 
        
      CALL this%Get_Phot4k_CovBL_At_samp( samp )
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^Set Parameters of The Electrons^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CALL this%Set_Kerr_Metric_At_samp() 
      CALL this%Set_Elec4U_CtrBL_At_samp()
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
      PhotE_In_Elec_CF = Vector4D_Inner_Product(this%Phot4k_CovBL_At_samp, &
                                                 this%Elec4U_CtrBL_At_samp)

      !write(unit = *, fmt = *)'***!!!!!!!!!!!!We are hererer You Dare!!!****************' 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'Phot4k_BLp1 === ', this%Phot4k_CovBL_At_samp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'Elec4U_BLp2 === ', this%Elec4U_CtrBL_At_samp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'PhotE_In_Elec_CF === ', PhotE_In_Elec_CF
      !write(unit = *, fmt = *)'************************************************************'
      get_Tau_fn = this%n_e * this%sigma_fn(T_e, DABS(PhotE_In_Elec_CF)) &
                       * DABS(PhotE_In_Elec_CF) / this%E_ini * Sigmap * rg_SUN
      end function get_Tau_fn
 

!*******************************************************************************************************
      real(mcp) function get_Tau2_fn(this, samp, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: samp, T_e
      real(mcp) :: Sigmap, mup, musinp, rp, phip, timep, &
                 sigp, P_mu_Multipl_U_mu, PhotE_In_Elec_CF
   
      !CALL YNOGK(p,this%Phot4k_LNRF, this%lambda, this%q, this%musin, this%mucos, &
      !                 this%aspin, this%r, one, rp, mup, phip, timep, sigp) 
      !write(*,*)'Ra mu ,22222 === ', rp, mup, phip, timep, sigp 

      !CALL this%Get_Phot_BLCoor_At_samp( samp )
      !Sigmap = this%r_samp**2 + ( this%aspin * this%mucos_samp )**2
      !write(*,*)'Sigma === ', Sigmap, this%r_p, this%mucos_p, this%musin_p, this%phi_p 
        
      !CALL this%Get_Phot4k_CovBL_At_samp( samp )
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^Set Parameters of The Electrons^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !CALL this%Set_Kerr_Metric_At_samp() 
      !CALL this%Set_Elec4U_CtrBL_At_samp()
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
      !PhotE_In_Elec_CF = this%Vector_Dot_Product(this%Phot4k_CovBL_At_samp, &
      !                                           this%Elec4U_CtrBL_At_samp)
      PhotE_In_Elec_CF = this%E_ini
 
      !write(unit = *, fmt = *)'***!!!!!!!!!!!!We are hererer You Dare!!!****************' 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'Phot4k_BLp1 === ', this%Phot4k_CovBL_At_samp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'Elec4U_BLp2 === ', this%Elec4U_CtrBL_At_samp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'PhotE_In_Elec_CF === ', PhotE_In_Elec_CF
      !write(unit = *, fmt = *)'************************************************************'
      get_Tau2_fn = this%n_e * this%sigma_fn(T_e, DABS(PhotE_In_Elec_CF))! &
      !                 * DABS(PhotE_In_Elec_CF) / this%E_ini * Sigmap * rg_SUN
      end function get_Tau2_fn
 
!*******************************************************************************************************
      subroutine Set_Cross_Section_sub(this, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: T_e
      real(mcp) :: temp, Ephoton
      real(mcp) :: dindexE
      integer(kind = 8) :: i, j, k, N, istat

      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma)
      open(unit=18, file = './data/SigmaArray.dat', status = "old", &
                           action = "read", iostat = istat)
      if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... '
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)'************************************************************' 
          do i = 0, N_sigma
              read(unit = 18, fmt = *)this%sigmaaTeE(i)
          enddo
      else
          open(unit=19, file = './data/SigmaArray.dat', status = "replace", &
                                   action = "write", iostat = istat)
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....'
              write(unit = *, fmt = *)'************************************************************' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE(i) = sigma_a(T_e, Ephoton)
                  write(unit = 19, fmt = *)this%sigmaaTeE(i)
                  write(unit = *, fmt = *)i
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File are Open Failed. The code have to Stop!!'
              Stop
          endif
          close(unit=19)
      endif
      close(unit=18)
      !i = 0
      !do 
      !    Ephoton = 10**( this%logE_low + this%dindexE * i )
           !this%sigmaaTeE(i) = sigma_a(T_e, Ephoton)
      !    sigmaaTeE(i) = sigma_a(T_e, Ephoton)
      !    i = i + 1
      ! !   write(*,*)i, N_sigma
      !    if ( i > N_sigma ) exit
      !enddo
      end subroutine Set_Cross_Section_sub


!*******************************************************************************************************
      function sigma_fn(this, T_e, Ephoton)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      real(mcp) :: sigma_fn, Ei, Ei1
      real(mcp), intent(in) :: T_e, Ephoton
      integer(kind = 8) :: i
      real(mcp) :: Delta_p, Theta_p, R_p
      real(mcp) :: Ep, sign_pth
      !integer, intent(in) :: n
!*******************************************************************************************************
    
      Ep = Ephoton
      if ( dlog10( Ep ) < this%logE_low ) then
          Ep = 10.D0**this%logE_low
      endif
      i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      !write(unit = *, fmt = *)'ttsfsfds!!',Ephoton, dlog10( Ephoton ), this%logE_low,this%dindexE, i
      if (i<0) then
          write(*,*)'sigma_fn = =', dlog10( Ephoton ), Ephoton, this%logE_low, this%dindexE,&
                ( dlog10( Ephoton ) - this%logE_low ) / this%dindexE!,&
          write(unit = *, fmt = *)'************************************************************' 
          write(*,*)'Phot4k_BLp1 === ', this%Phot4k_CovBL_At_samp
          write(unit = *, fmt = *)'************************************************************' 
          write(*,*)'Elec4U_BLp2 === ', this%Elec4U_CtrBL_At_samp
          write(unit = *, fmt = *)'************************************************************'  
          associate( rp => this%r_samp, &
                 a => this%aspin, &
                 mucosp => this%mucos_samp, &
                 musinp => this%musin_samp, &
                 l => this%lambda_ini, &
                 q => this%q_ini )
          Delta_p = rp*rp - two*rp + a**2
          R_p = rp**4 - ( q + l**2 - a**2 )*rp**2 + two*( q + (l - a)**2 )*rp - a**2*q
          if (l .NE. zero) then
              Theta_p = q + ( a*mucosp )**2 - ( l * mucosp / musinp )**2
          else
              Theta_p = q + ( a*mucosp )**2
          endif
          write(unit = *, fmt = *)'************************************************************'  
          write(*,*)'R_p1  = ', R_p, Theta_p, this%p_samp
          write(unit = *, fmt = *)'************************************************************'  
          write(*,*)'R_p2  = ',  this%r_samp,  this%mucos_samp, this%musin_samp
          write(unit = *, fmt = *)'************************************************************'  
          write(*,*)'R_p3  = ',  this%lambda_ini, this%q_ini, l ,q
          write(unit = *, fmt = *)'************************************************************' 
          write(*,*)'R_p4  = ',  this%E_ini, this%p_maxs
          write(unit = *, fmt = *)'************************************************************'
          end associate
          Delta_p = mucos(this%p_samp, this%Phot4k_LNRF_ini(4), &
                                 this%Phot4k_LNRF_ini(3), &
                                 this%lambda_ini, &
                                 this%q_ini,  &
                                 this%musin_ini, &
                                 this%mucos_ini, &
                                 this%aspin,-one, sign_pth)
          write(*,*)'mucos_samp aaaa = ',  Delta_p, this%p_samp
          write(unit = *, fmt = *)'************************************************************' 
      endif 

      Ei = 10**( this%logE_low + this%dindexE * i )
      Ei1 = 10**( this%logE_low + this%dindexE * (i + 1) )
      sigma_fn = ( this%sigmaaTeE(i + 1) - this%sigmaaTeE(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE(i) 
      end function sigma_fn

!*******************************************************************************************************
      subroutine Calc_Phot_Polarization_At_Observor_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      !real(mcp), intent(in) :: a, b 
      real(mcp), parameter :: al0 = -50.D0, be0 = -50.D0
      real(mcp) :: gama_S, beta_T, f_the, f_phi, Sin_The_Obs, Cos_The_Obs
      real(mcp) :: mup, mup0, p_max, p_temp, sign_p_th
      real(mcp) :: Q_obs, U_obs, Psi_Obs, rdp
      real(mcp) :: alpha, beta, delta_th = two / 100.D0, mu_obs
      !real(mcp), dimension(0:300, 0:300) :: disk_image = zero
      real(mcp) :: delta_al = -two*al0 / 400.D0, delta_be = -two*be0 / 400.D0
      real(mcp) :: sign_pth
      integer :: i, j, h, k, mu_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !ph1%Phot4k_CF! = sph%Scattered_Phot4k_CF
      !ph1%Phot4k_CovCF! = sph%Scattered_Phot4k_CovCF
      !ph1%f4_CF = sph%f4_scat_CF
      !ph1%f4_CovCF = sph%f4_scat_CovCF
      !ph1%delta_pd = sph%delta_pd_scat
      !call ph1%get_Phot4k_LNRF()
      !write(unit = *, fmt = *)'************************************************************'  
      !write(unit = *, fmt = *)'************************************************************' 
      !CALL ph1%get_lambdaq_At_p()
      !write(*,*)'elambda q  ==', pe1%lambda, pe1%q
      !write(*,*)'elambda q  ==', ph1%lambda, ph1%q
      !**************************************************************
      !CALL ph1%get_Phot4k_BL()  
      !CALL ph1%get_Phot4k_CovBL()  
      !ph1%E = DABS( ph1%Phot4k_CovBL(1) )    !E is a constant, just like the lambda and q
      !CALL ph1%get_f4_CtrBL_From_f4_CF()     !Now f4_CtrBL has been obtained.
      !CALL ph1%get_f4_CovBL_From_f4_CovCF()  !Now f4_CovBL has been obtained.
      !******************************************************************************************
      p_max = p_total(this%Phot4k_LNRF_ini(2), this%lambda_ini, &
                      this%q_ini, this%musin_ini, this%mucos_ini, &
                      this%aspin, this%r_ini, one)
      Cos_The_Obs = mucos(p_max, this%Phot4k_LNRF_ini(4), &
                                 this%Phot4k_LNRF_ini(3), &
                                 this%lambda_ini, this%q_ini, &
                                 this%musin_ini, this%mucos_ini, &
                                 this%aspin, one, sign_pth)
      Sin_The_Obs = DSQRT( one - Cos_The_Obs**2 )
      !p_temp = p_max - p_max * 1.D-5 
      !mup = mucos(p_temp, this%Phot4k_LNRF_ini(4), &
      !                    this%Phot4k_LNRF_ini(3), &
      !                    this%lambda_ini, this%q_ini, &
      !                    this%musin_ini, this%mucos_ini, &
      !                    this%aspin, one, sign_pth)
      !rdp = radius(p_max, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                    this%q_ini, this%aspin, this%r_ini, one)
      !sign_p_th = - DSIGN(one, Cos_The_Obs - mup )
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**mup === ', p_max, mup, mup0, sign_p_th, mup - mup0, rdp
      !write(unit = *, fmt = *)'************************************************************' 
      !******************************************************************************************
      gama_S = this%lambda_ini / Sin_The_Obs - this%aspin * Sin_The_obs
      beta_T = sign_pth * DSQRT( this%q_ini + ( this%aspin * Cos_The_Obs )**2 - &
                                ( this%lambda_ini * Cos_The_Obs / Sin_The_Obs )**2 )
      f_the =   ( this%K2_pw_ini * gama_S - this%K1_pw_ini * &
                          beta_T ) / ( gama_S**2 + beta_T**2 )
      f_phi = - ( this%K1_pw_ini * gama_S + this%K2_pw_ini * &
                          beta_T ) / ( gama_S**2 + beta_T**2 )
      Psi_Obs = DATAN( DABS( f_the / f_phi ) )
      If ( f_the * f_phi < zero ) then
          Psi_Obs = Pi - Psi_Obs
      Endif
      Q_obs = this%delta_pd * DCOS( two*Psi_Obs )
      U_obs = this%delta_pd * DSIN( two*Psi_Obs )
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'** Q U_obs === ', Q_obs, U_obs, this%delta_pd
      !write(unit = *, fmt = *)'************************************************************'
      !******************************************************************************************
      mu_obs = dcos( 30.D0 * pi / 180.D0  )
      If ( mu_obs - delta_th / two < Cos_The_Obs .AND. &
           mu_obs + delta_th / two > Cos_The_Obs ) Then
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          alpha = - this%lambda_ini / dsqrt( one - Cos_The_Obs**2 )
          beta = sign_pth * dsqrt( this%q_ini + Cos_The_Obs**2 * &
                                 ( this%aspin**2 - alpha**2 ) )
          i = (alpha - al0) / delta_al
          j = (beta - be0) / delta_be
          If (i < 0 .or. j <0) then
              write(*,*)'jjjii111 = = ', i, j, alpha, beta
              write(*,*)'jjjii222 = = ', i, j, this%lambda_ini, Cos_The_Obs
          endif
          this%disk_image(i,j) = this%disk_image(i,j) + this%delta_pd
          this%effect_number = this%effect_number + 1
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************'
      !write(unit = *, fmt = *)'** alpha beta === ', i, j, beta, alpha, delta_al, delta_be
      !write(unit = *, fmt = *)'************************************************************'
      !stop
      !******************************************************************************************
      Endif
      mu_i = floor( (one - Cos_The_Obs ) / delta_th )
      this%delta_pds( mu_i ) = this%delta_pd + this%delta_pds( mu_i )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          this%frequency_v = DABS( this%Phot4k_CovBL_ini(1) ) * 1.D6 / h_ev
          k = floor( dlog10( this%frequency_v / 1.D15 ) / 0.016D0 )
          !write(*, *)'ssss = = ', k,dlog10( this%frequency_v ), this%frequency_v, this%Phot4k_CovBL_ini(1)
          !stop
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * this%frequency_v**3
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      return
      end subroutine Calc_Phot_Polarization_At_Observor_Sub

!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      !real(mcp), intent(in) :: a, b 
      real(mcp), parameter :: al0 = -50.D0, be0 = -50.D0
      real(mcp) :: gama_S, beta_T, f_the, f_phi, Sin_The_Obs, Cos_The_Obs
      real(mcp) :: mup, mup0, p_max, p_temp, sign_p_th
      real(mcp) :: Q_obs, U_obs, Psi_Obs, rdp
      real(mcp) :: alpha, beta, delta_th = two / 100.D0, mu_obs
      !real(mcp), dimension(0:300, 0:300) :: disk_image = zero
      real(mcp) :: delta_al = -two*al0 / 400.D0, delta_be = -two*be0 / 400.D0
      real(mcp) :: sign_pth, r_times_p, r_ini2, p_ini_obs
      integer :: i, j, h, k, mu_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !ph1%Phot4k_CF! = sph%Scattered_Phot4k_CF
      !ph1%Phot4k_CovCF! = sph%Scattered_Phot4k_CovCF
      !ph1%f4_CF = sph%f4_scat_CF
      !ph1%f4_CovCF = sph%f4_scat_CovCF
      !ph1%delta_pd = sph%delta_pd_scat
      !call ph1%get_Phot4k_LNRF()
      !write(unit = *, fmt = *)'************************************************************'  
      !write(unit = *, fmt = *)'************************************************************' 
      !CALL ph1%get_lambdaq_At_p()
      !write(*,*)'elambda q  ==', pe1%lambda, pe1%q
      !write(*,*)'elambda q  ==', ph1%lambda, ph1%q
      !**************************************************************
      !CALL ph1%get_Phot4k_BL()  
      !CALL ph1%get_Phot4k_CovBL()  
      !ph1%E = DABS( ph1%Phot4k_CovBL(1) )    !E is a constant, just like the lambda and q
      !CALL ph1%get_f4_CtrBL_From_f4_CF()     !Now f4_CtrBL has been obtained.
      !CALL ph1%get_f4_CovBL_From_f4_CovCF()  !Now f4_CovBL has been obtained.
      !******************************************************************************************
      !p_max = p_total(this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                this%q_ini, this%musin_ini, this%mucos_ini, &
      !                this%aspin, this%r_ini, one)
      !Cos_The_Obs = mucos(p_max, this%Phot4k_LNRF_ini(4), &
      !                           this%Phot4k_LNRF_ini(3), &
      !                           this%lambda_ini, this%q_ini, &
      !                           this%musin_ini, this%mucos_ini, &
      !                           this%aspin, one, sign_pth)
      !Sin_The_Obs = DSQRT( one - Cos_The_Obs**2 )
      !p_temp = p_max - p_max * 1.D-5 
      !mup = mucos(p_temp, this%Phot4k_LNRF_ini(4), &
      !                    this%Phot4k_LNRF_ini(3), &
      !                    this%lambda_ini, this%q_ini, &
      !                    this%musin_ini, this%mucos_ini, &
      !                    this%aspin, one, sign_pth)
      !rdp = radius(p_max, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                    this%q_ini, this%aspin, this%r_ini, one)
      !sign_p_th = - DSIGN(one, Cos_The_Obs - mup )
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**mup === ', p_max, mup, mup0, sign_p_th, mup - mup0, rdp
      !write(unit = *, fmt = *)'************************************************************' 
      !******************************************************************************************
      !gama_S = this%lambda_ini / Sin_The_Obs - this%aspin * Sin_The_obs
      !beta_T = sign_pth * DSQRT( this%q_ini + ( this%aspin * Cos_The_Obs )**2 - &
      !                          ( this%lambda_ini * Cos_The_Obs / Sin_The_Obs )**2 )
      !f_the =   ( this%K2_pw_ini * gama_S - this%K1_pw_ini * &
      !                    beta_T ) / ( gama_S**2 + beta_T**2 )
      !f_phi = - ( this%K1_pw_ini * gama_S + this%K2_pw_ini * &
      !                    beta_T ) / ( gama_S**2 + beta_T**2 )
      !Psi_Obs = DATAN( DABS( f_the / f_phi ) )
      !If ( f_the * f_phi < zero ) then
      !    Psi_Obs = Pi - Psi_Obs
      !Endif
      !Q_obs = this%delta_pd * DCOS( two*Psi_Obs )
      !U_obs = this%delta_pd * DSIN( two*Psi_Obs )
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'** Q U_obs === ', Q_obs, U_obs, this%delta_pd
      !write(unit = *, fmt = *)'************************************************************'
      !******************************************************************************************
      !mu_obs = dcos( 30.D0 * pi / 180.D0  )
      !If ( mu_obs - delta_th / two < Cos_The_Obs .AND. &
      !     mu_obs + delta_th / two > Cos_The_Obs ) Then
      If ( .false. ) then
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          alpha = - this%lambda_ini / dsqrt( one - Cos_The_Obs**2 )
          beta = sign_pth * dsqrt( this%q_ini + Cos_The_Obs**2 * &
                                 ( this%aspin**2 - alpha**2 ) )
          i = (alpha - al0) / delta_al
          j = (beta - be0) / delta_be
          If (i < 0 .or. j <0) then
              write(*,*)'jjjii111 = = ', i, j, alpha, beta
              write(*,*)'jjjii222 = = ', i, j, this%lambda_ini, Cos_The_Obs
          endif
          this%disk_image(i,j) = this%disk_image(i,j) + this%delta_pd
          this%effect_number = this%effect_number + 1
      !******************************************************************************************
      !write(unit = *, fmt = *)'************************************************************'
      !write(unit = *, fmt = *)'** alpha beta === ', i, j, beta, alpha, delta_al, delta_be
      !write(unit = *, fmt = *)'************************************************************'
      !stop
      !******************************************************************************************
      Endif
      !mu_i = floor( (one - Cos_The_Obs ) / delta_th )
      !this%delta_pds( mu_i ) = this%delta_pd + this%delta_pds( mu_i )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      r_times_p = Vector3D_Inner_Product( this%vector_of_position_ini,&
                                 this%Vector_of_Momentum_ini )
      r_ini2 = Vector3D_Inner_Product( this%vector_of_position_ini,&
                                       this%vector_of_position_ini )
      p_ini_obs = - r_times_p + dsqrt( r_times_p**2 + this%r_obs**2 - r_ini2 )
      this%time_arrive_observer = this%time_arrive_observer + p_ini_obs / Cv
      !write(*,*)'times == ',this%time_arrive_observer
      If ( .Not. this%At_outer_Shell ) then
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          k = floor( dlog10( this%frequency_v / 1.D15 ) / 0.016D0 )
          !this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * this%frequency_v * h_ev * &
          !                  dexp( - this%Optical_Depth_scatter )
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * &
                            dexp( - this%Optical_Depth_scatter )
          !write(*,*)'w_ini = ', k,dexp( - this%Optical_Depth_scatter )
      Else
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          k = floor( dlog10( this%frequency_v / 1.D15 ) / 0.016D0 )
          !write(*, *)'ssss = = ', k,dlog10( this%frequency_v ), this%frequency_v, this%Phot4k_CovBL_ini(1)
          !stop
          !this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * this%frequency_v * h_ev
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini
          !write(*,*)'w_ini = ', this%w_ini
      Endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      return
      end subroutine Calc_Phot_Informations_At_Observor_Sub


!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this 
      !real(mcp), intent(in) :: a, b 
      real(mcp), parameter :: al0 = -50.D0, be0 = -50.D0
      real(mcp) :: gama_S, beta_T, f_the, f_phi, Sin_The_Obs, Cos_The_Obs
      real(mcp) :: mup, mup0, p_max, p_temp, sign_p_th
      real(mcp) :: Q_obs, U_obs, Psi_Obs, rdp
      real(mcp) :: alpha, beta, delta_th = two / 100.D0, mu_obs
      !real(mcp), dimension(0:300, 0:300) :: disk_image = zero
      real(mcp) :: delta_al = -two*al0 / 400.D0, delta_be = -two*be0 / 400.D0
      real(mcp) :: sign_pth, r_times_p, r_ini2, p_ini_obs
      integer :: i, j, h, k, mu_i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !r_times_p = Vector3D_Inner_Product( this%vector_of_position_ini,&
      !                           this%Vector_of_Momentum_of_photon_ini )
      !r_ini2 = Vector3D_Inner_Product( this%vector_of_position_ini,&
      !                                 this%vector_of_position_ini )
      !p_ini_obs = - r_times_p + dsqrt( r_times_p**2 + this%r_obs**2 - r_ini2 )
      !this%time_arrive_observer = this%time_arrive_observer + p_ini_obs / Cv
      !write(*,*)'times == ',this%time_arrive_observer
      If ( this%direct_escaped ) then
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          k = floor( dlog10( this%frequency_v / 1.D15 ) / 0.016D0 )
          !this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * this%frequency_v * h_ev * &
          !                  dexp( - this%Optical_Depth_scatter )
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * &
                            dexp( - this%Optical_Depth_scatter )
          !write(*,*)'w_ini = ', k,dexp( - this%Optical_Depth_scatter ) 
      Endif
      If ( this%At_outer_Shell ) then 
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          k = floor( dlog10( this%frequency_v / 1.D15 ) / 0.016D0 )
          !write(*, *)'ssss = = ', k,dlog10( this%frequency_v ), this%frequency_v, this%Phot4k_CovBL_ini(1)
          !stop
          !this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * this%frequency_v * h_ev
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini
          !write(*,*)'w_ini = ', this%w_ini
      Endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      return
      end subroutine Calc_Phot_Informations_At_Observor_2zones_Sub
!*******************************************************************************************************
!*******************************************************************************************************
!**  SubRoutines In terms of the Electron that sacttering or interacting with the photon  **************
!*******************************************************************************************************
!*******************************************************************************************************

!*******************************************************************************************************
      !subroutine Set_Elec4U_CtrBL_At_ini_Sub(this)
!*******************************************************************************************************
      !class(Photon) :: this
      !real(mcp) :: Omega_K, Omega
   
      !associate( r => this%r_ini, &
      !!           theta => this%theta_ini )
      !Omega_K = one / ( this%aspin + r**1.5D0 ) !Kepler velocity
      !Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
      !              ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega_ini
      !this%Elec4U_CtrBL_At_ini(1) = one / DSQRT( - ( this%g_tt_ini + two * this%g_tp_ini * Omega + &
      !                                               this%g_pp_ini * Omega**2 ) )
      !end associate
      
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'****** Set_Elec4U_BL === ',this%g_tt,this%g_tp ,this%g_pp ,&
      !                                     this%g_up_tt,this%g_up_tp ,this%g_up_pp ,'***************'
      !write(unit = *, fmt = *)'************************************************************' 
  
      !this%Elec4U_CtrBL_At_ini(2) = zero
      !this%Elec4U_CtrBL_At_ini(3) = zero
      !this%Elec4U_CtrBL_At_ini(4) = this%Elec4U_CtrBL_At_ini(1) * Omega
      !end subroutine Set_Elec4U_CtrBL_At_ini_Sub

!*******************************************************************************************************
!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!~* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!~* 
!~* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!*******************************************************************************************************
      subroutine Get_Optical_Depth_At_p_Sub(this, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      integer, parameter :: n = 200
      real(mcp), intent(in) :: T_e
      real(mcp) :: alpha_nu, E0
      real(mcp) :: Elec4U_CtrBL(1:4), Omega_K, Omega, p, rp, theta, mucosp, musinp, &
                phip, tp, sign_pr, sign_pth, a, q, l, R_p, Theta_p, Phot4k_CovBL(1:4), &
                PhotE_In_ComvFrame, Delta_p, Sigmap, Tau_nu, Del_p, Del_tau0, Del_tau_n
      integer :: i
      real(mcp) :: Tau_nu_i(0:n-1), Del_Tau_nu, Integral_jvlv3, Del_Integral_jvlv3, &
               Integral_jvlv3_0, Integral_jvlv3_n, jvlv3, kappa = 8.D-46, v_frequency, &
               Tau_scat, Del_Tau_scat, Del_Tau_scat_n, Del_Tau_scat_0
  
      Tau_nu = zero
      Tau_scat = zero
      Tau_nu_i = zero 
      !Integral_jvlv3 = zero
      Del_p = this%p_scattering / (n + 1)
      i = n
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          p = Del_p*i
          rp = radius(p,this%Phot4k_LNRF_ini(2), &
                             this%lambda_ini, &
                             this%q_ini, &
                             this%aspin, &
                             this%r_ini, one, sign_pr)
          mucosp = mucos(p, this%Phot4k_LNRF_ini(4), &
                                 this%Phot4k_LNRF_ini(3), &
                                 this%lambda_ini, &
                                 this%q_ini,  &
                                 this%musin_ini, &
                                 this%mucos_ini, &
                                 this%aspin, one, sign_pth)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          theta = DACOS( mucosp )
          musinp = DSQRT( one - mucosp**2 )
          phip = zero
          tp = zero 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          this%r = rp
          this%mucos = mucosp
          this%musin = musinp
          call this%Set_Kerr_Metric()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
          Omega_K = one / ( this%aspin + rp**1.5D0 ) !Kepler velocity
          Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
                    ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega
          Elec4U_CtrBL(1) = one / DSQRT( - ( this%g_tt + two * this%g_tp * Omega + &
                                                     this%g_pp * Omega**2 ) ) 
  
          Elec4U_CtrBL(2) = zero
          Elec4U_CtrBL(3) = zero
          Elec4U_CtrBL(4) = Elec4U_CtrBL(1) * Omega  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Get_Phot4k_CovBL_At_samp( samp )  
          a = this%aspin 
          l = this%lambda_ini
          q = this%q_ini

          Delta_p = rp**2 - two*rp + a**2
          R_p = rp**4 - ( q + l**2 - a**2 )*rp**2 + two*( q + (l - a)**2 )*rp - a**2*q
          if (l .NE. zero) then
              Theta_p = q + ( a*mucosp )**2 - ( l * mucosp / musinp )**2
          else
              Theta_p = q + ( a*mucosp )**2
          endif 

          If (R_p < zero  ) then
              !write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'R_p ff === ', rdp - rp , R_p , DSQRT(dabs(R_p)), rp
              !write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'r_ini  P_r === ', this%r_ini, this%Phot4k_LNRF_ini(2), r_tp1, r_tp2
              write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'r_tp1, r_tp2 === ',&
              !r_tp1**4 - ( q + l**2 - a**2 )*r_tp1**2 + two*( q + (l - a)**2 )*r_tp1 - a**2*q, & 
              !r_tp2**4 - ( q + l**2 - a**2 )*r_tp2**2 + two*( q + (l - a)**2 )*r_tp2 - a**2*q
              write(*,*)'rdp sssss === ',R_p
              R_p = - R_p
          endif
          If (Theta_p < zero )then! .and. dabs(Theta_p)<1.D-9) then 
              write(*,*)'rdp sssss === ',Theta_p!this%Phot4k_LNRF_ini, &
              !                           this%Phot4k_LNRF_ini(3)/this%Phot4k_LNRF_ini(1)
               Theta_p = -Theta_p
          endif 
          Phot4k_CovBL(1) = - one
          Phot4k_CovBL(2) = sign_pr * DSQRT(R_p) / Delta_p
          !this%Phot4k_CovBL_At_samp(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
          Phot4k_CovBL(3) = sign_pth * DSQRT(Theta_p)
          Phot4k_CovBL(4) = l
          Phot4k_CovBL = Phot4k_CovBL !* this%E_ini 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Set_Kerr_Metric_At_samp() 
          !CALL this%Set_Elec4U_CtrBL_At_samp()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          E0 = - Vector4D_Inner_Product(Phot4k_CovBL, Elec4U_CtrBL)
          v_frequency = E0 * this%E_ini * 1.D6 / h_ev ! the unit of E0 is MeV, so 1.D6 is used to eV.
          Sigmap = ( rp**2 + ( a * mucosp )**2 )*rg_SUN
          !alpha_nu = AbsorptionCoefficient(E0, this%n_e)
          alpha_nu = (this%n_e)*sigma_T
          Del_Tau_nu = alpha_nu * E0 * Sigmap 
          !write(*,*)'11111111', Phot4k_CovBL, Elec4U_CtrBL, E0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
            !get_Tau_fn = this%n_e * this%sigma_fn(T_e, DABS(PhotE_In_Elec_CF)) &
            !           * DABS(PhotE_In_Elec_CF) * Sigmap * rg_SUN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Tau_nu = Tau_nu + Del_Tau_nu 
          Tau_scat = Tau_scat + Del_Tau_scat
          !jvlv3 = this%n_e**2 / ( E0 * this%E_ini * 1.D3 ) / dsqrt( this%T_e*1.D3 ) * &
          !        dexp( - E0 * this%E_ini / this%T_e ) / v_frequency**3 !(E0 * this%E_ini)**3 
          If (i /= 0 .and. i /= n) then
              !Tau_nu_i(i) = ( Tau_nu - Del_tau_n / two - Del_Tau_nu / two ) * Del_p
              !jvlv3 = this%n_e**2 / E0 / dsqrt( T_e )*dexp( - E0 / T_e ) / E0**3
              !Del_Integral_jvlv3 = jvlv3 * dexp( - Tau_nu_i(i) ) * E0 * Sigmap 
              !Integral_jvlv3 = Integral_jvlv3 + Del_Integral_jvlv3
          Else if ( i == n ) then
              Del_tau_n = Del_Tau_nu
              Del_Tau_scat_n = Del_Tau_scat
              !Integral_jvlv3_n = jvlv3 * E0 * Sigmap
              !Integral_jvlv3 = Integral_jvlv3 + Integral_jvlv3_n
          Else if ( i == 0 ) then
              Del_tau0 = Del_Tau_nu 
              Del_Tau_scat_0 = Del_Tau_scat
          Endif
          i = i - 1
          If ( i < 0 ) exit
          !write(*,*)'aasdsff = ',jvlv3, Del_Integral_jvlv3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Enddo
      Tau_nu = ( Tau_nu - Del_tau0 / two - Del_tau_n / two )*Del_p
      Tau_scat = ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p
      !Tau_nu_i(0) = Tau_nu
      !Integral_jvlv3_0 = jvlv3 * dexp( - Tau_nu_i(0) ) * E0 * Sigmap
      !Integral_jvlv3 = ( Integral_jvlv3 - Integral_jvlv3_n / two + &
      !                   Integral_jvlv3_0 / two ) * Del_p * Kappa
 
      this%w_p = this%w_ini * dexp( - Tau_nu ) / this%A_normal * &
                 dexp( (this%A_normal - one)*Tau_scat ) ! + Integral_jvlv3 
      !write(*,*)'sstau_nu1 = ', this%w_ini, this%w_p, this%A_normal,Tau_scat,Tau_nu
      !write(*,*)'sstau_nu2 = ', (this%A_normal - one)*Tau_scat, dexp( (this%A_normal - one)*Tau_scat )
      
      end subroutine Get_Optical_Depth_At_p_Sub
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


!*******************************************************************************************************
      subroutine Get_bias_parameter_Sub(this, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      integer, parameter :: n = 200
      real(mcp), intent(in) :: T_e
      real(mcp) :: alpha_nu, E0
      real(mcp) :: Elec4U_CtrBL(1:4), Omega_K, Omega, p, rp, theta, mucosp, musinp, &
                phip, tp, sign_pr, sign_pth, a, q, l, R_p, Theta_p, Phot4k_CovBL(1:4), &
                PhotE_In_ComvFrame, Delta_p, Sigmap, Tau_nu, Del_p, Del_tau0, Del_tau_n
      integer :: i
      real(mcp) :: Tau_nu_i(0:n-1), Del_Tau_nu, Integral_jvlv3, Del_Integral_jvlv3, &
               Integral_jvlv3_0, Integral_jvlv3_n, jvlv3, kappa = 8.D-46, v_frequency, &
               Tau_scat, Del_Tau_scat, Del_Tau_scat_n, Del_Tau_scat_0
  
      Tau_nu = zero
      Tau_scat = zero
      Tau_nu_i = zero 
      !Integral_jvlv3 = zero
      Del_p = this%p_maxs / (n + 1)
      this%Sigma_Max = zero
      i = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          p = Del_p*i
          rp = radius(p,this%Phot4k_LNRF_ini(2), &
                             this%lambda_ini, &
                             this%q_ini, &
                             this%aspin, &
                             this%r_ini, one, sign_pr)
          mucosp = mucos(p, this%Phot4k_LNRF_ini(4), &
                                 this%Phot4k_LNRF_ini(3), &
                                 this%lambda_ini, &
                                 this%q_ini,  &
                                 this%musin_ini, &
                                 this%mucos_ini, &
                                 this%aspin, one, sign_pth)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          theta = DACOS( mucosp )
          musinp = DSQRT( one - mucosp**2 )
          phip = zero
          tp = zero 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          this%r = rp
          this%mucos = mucosp
          this%musin = musinp
          call this%Set_Kerr_Metric()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   
          Omega_K = one / ( this%aspin + rp**1.5D0 ) !Kepler velocity
          Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
                    ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega
          Elec4U_CtrBL(1) = one / DSQRT( - ( this%g_tt + two * this%g_tp * Omega + &
                                                     this%g_pp * Omega**2 ) ) 
  
          Elec4U_CtrBL(2) = zero
          Elec4U_CtrBL(3) = zero
          Elec4U_CtrBL(4) = Elec4U_CtrBL(1) * Omega  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Get_Phot4k_CovBL_At_samp( samp )  
          a = this%aspin 
          l = this%lambda_ini
          q = this%q_ini

          Delta_p = rp**2 - two*rp + a**2
          R_p = rp**4 - ( q + l**2 - a**2 )*rp**2 + two*( q + (l - a)**2 )*rp - a**2*q
          if (l .NE. zero) then
              Theta_p = q + ( a*mucosp )**2 - ( l * mucosp / musinp )**2
          else
              Theta_p = q + ( a*mucosp )**2
          endif 

          If (R_p < zero  ) then
              !write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'R_p ff === ', rdp - rp , R_p , DSQRT(dabs(R_p)), rp
              !write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'r_ini  P_r === ', this%r_ini, this%Phot4k_LNRF_ini(2), r_tp1, r_tp2
              write(unit = *, fmt = *)'************************************************************'
              !write(*,*)'r_tp1, r_tp2 === ',&
              !r_tp1**4 - ( q + l**2 - a**2 )*r_tp1**2 + two*( q + (l - a)**2 )*r_tp1 - a**2*q, & 
              !r_tp2**4 - ( q + l**2 - a**2 )*r_tp2**2 + two*( q + (l - a)**2 )*r_tp2 - a**2*q
              write(*,*)'rdp sssss === ',R_p
              R_p = - R_p
          endif
          If (Theta_p < zero )then! .and. dabs(Theta_p)<1.D-9) then 
              write(*,*)'rdp sssss === ',Theta_p!this%Phot4k_LNRF_ini, &
              !                           this%Phot4k_LNRF_ini(3)/this%Phot4k_LNRF_ini(1)
               Theta_p = -Theta_p
          endif 
          Phot4k_CovBL(1) = - one
          Phot4k_CovBL(2) = sign_pr * DSQRT(R_p) / Delta_p
          !this%Phot4k_CovBL_At_samp(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
          Phot4k_CovBL(3) = sign_pth * DSQRT(Theta_p)
          Phot4k_CovBL(4) = l
          !Phot4k_CovBL = Phot4k_CovBL * this%E_ini 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Set_Kerr_Metric_At_samp() 
          !CALL this%Set_Elec4U_CtrBL_At_samp()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          E0 = - Vector4D_Inner_Product(Phot4k_CovBL, Elec4U_CtrBL) 
          Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
          !alpha_nu = (this%n_e)*sigma_T
          !Del_Tau_nu = alpha_nu * E0 * Sigmap 
          !E0 = E0 * this%E_ini
          !write(*,*)'sspp==',sign_pr,DSQRT(R_p), Delta_p, rp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !write(*,*)'ssssssss====',T_e, E0
          !write(*,*)'222222', Phot4k_CovBL, Elec4U_CtrBL, E0
          Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
          If( this%Sigma_Max < Del_Tau_scat )then
              this%Sigma_Max = Del_Tau_scat
          endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Tau_scat = Tau_scat + Del_Tau_scat    
          If ( i == 0 ) then 
              Del_Tau_scat_0 = Del_Tau_scat
          Endif
          i = i + 1
          If( rp > this%R_out*1.1D0 .or. rp < this%R_in*0.99D0 )then
              Del_Tau_scat_n = Del_Tau_scat
              this%p_out = p
              exit
          endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Enddo 
      Tau_scat = ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p 
  
      !write(*,*)'Sigma_Max = ', this%Sigma_Max, this%p_maxs, this%p_out, Del_p
      !this%Sigma_Max = Tau_scat
      
      end subroutine Get_bias_parameter_Sub 


!*******************************************************************************************************
      subroutine Get_Max_Value_of_Sigma_Sub(this, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      integer, parameter :: n = 500
      real(mcp), intent(in) :: T_e
      real(mcp) :: alpha_nu, E0
      real(mcp) :: Elec4U_CtrBL(1:4), Omega_K, Omega, p, rp, theta, mucosp, musinp, &
                phip, tp, sign_pr, sign_pth, a, q, l, R_p, Theta_p, Phot4k_CovBL(1:4), &
                PhotE_In_ComvFrame, Delta_p, Sigmap, Tau_nu, Del_p, Del_tau0, Del_tau_n
      integer :: i
      real(mcp) :: Tau_nu_i(0:n-1), Del_Tau_nu, Integral_jvlv3, Del_Integral_jvlv3, &
               Integral_jvlv3_0, Integral_jvlv3_n, jvlv3, kappa = 8.D-46, v_frequency, &
               Tau_scat, Del_Tau_scat, Del_Tau_scat_n, Del_Tau_scat_0
      real(mcp) :: r_ini
  
      Tau_nu = zero
      Tau_scat = zero
      Tau_nu_i = zero 
      !Integral_jvlv3 = zero
      Del_p = this%p_maxs / (n + 1)
      this%Sigma_Max = zero
      r_ini = this%r_ini!dsqrt( this%x_ini**2 + this%y_ini**2 + this%z_ini**2 )
      i = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          p = Del_p*i
          !this%x_p = this%x_ini + p * this%Vector_of_Momentum_of_photon_ini(1)
          !this%y_p = this%y_ini + p * this%Vector_of_Momentum_of_photon_ini(2)
          !this%z_p = this%z_ini + p * this%Vector_of_Momentum_of_photon_ini(3)
          !rp = dsqrt( this%x_p**2 + this%y_p**2 + this%z_p**2 )
          rp = this%r_p2( this%Vector_of_position_ini, this%Vector_of_Momentum_ini, p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !Omega_K = one / ( this%aspin + rp**1.5D0 ) !Kepler velocity
          !Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
          !          ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega
          !Elec4U_CtrBL(1) = one / DSQRT( - ( this%g_tt + two * this%g_tp * Omega + &
          !                                           this%g_pp * Omega**2 ) ) 
  
          !Elec4U_CtrBL(2) = zero
          !Elec4U_CtrBL(3) = zero
          !Elec4U_CtrBL(4) = Elec4U_CtrBL(1) * Omega  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Get_Phot4k_CovBL_At_samp( samp )   
          !Phot4k_CovBL(1) = - one
          !Phot4k_CovBL(2) = sign_pr * DSQRT(R_p) / Delta_p
          !this%Phot4k_CovBL_At_samp(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
          !Phot4k_CovBL(3) = sign_pth * DSQRT(Theta_p)
          !Phot4k_CovBL(4) = l 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !CALL this%Set_Kerr_Metric_At_samp() 
          !CALL this%Set_Elec4U_CtrBL_At_samp()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !E0 = - this%Vector_Dot_Product(Phot4k_CovBL, Elec4U_CtrBL)
          E0 = one
          !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
          Sigmap = one 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !write(*,*)'ssssssss====',E0,this%E_ini,this%Sigma_Max
          !write(*,*)'222222', Phot4k_CovBL, Elec4U_CtrBL, E0
          !write(*,*)'222222', this%n_e , this%sigma_fn(T_e, E0*this%E_ini) , E0 , Sigmap
          Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
          If( this%Sigma_Max < Del_Tau_scat )then
              this%Sigma_Max = Del_Tau_scat
          endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Tau_scat = Tau_scat + Del_Tau_scat    
          If ( i == 0 ) then 
              Del_Tau_scat_0 = Del_Tau_scat
          Endif
          i = i + 1
!write(*,*)'pp=',r_ini,this%R_in,rp
          If( r_ini <= this%R_in )then
!write(*,*)'sss=',r_ini,this%R_in,rp
              If ( rp > this%R_in * 1.001D0 ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  this%p_out = p
                  exit
              Endif  
          endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !If( r_ini < this%R_in )then
          !    If ( r_p > this%R_in * 1.001D0 ) then
          !        Del_Tau_scat_n = Del_Tau_scat
          !        this%p_out = p
          !        exit
          !    Endif 
          !Else If( r_ini < this%R_out )then
          !    If ( r_p > this%R_out * 1.001D0 .or. r_p < this%R_in * 0.999D0 ) then
          !        Del_Tau_scat_n = Del_Tau_scat
          !        this%p_out = p
          !        exit
          !    Endif 
          !endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Enddo 
      this%Optical_Depth_scatter = ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p 
  
      !write(*,*)'Sigma_Max = ', this%Optical_Depth_scatter, &
      !        this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * this%p_out
      !this%Sigma_Max = Tau_scat
      
      end subroutine Get_Max_Value_of_Sigma_Sub
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      function r_p2( this, r_ini, vector_of_momentum, p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: r_ini(1:3), vector_of_momentum(1:3), p
      real(mcp) :: vector_p(1:3), r_p2
 
      vector_p = r_ini + p * vector_of_momentum
      this%vector_of_position = vector_p
      r_p2 = dsqrt( vector_p(1)**2 + vector_p(2)**2 + vector_p(3)**2 )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function r_p2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!*******************************************************************************************************
      subroutine Get_Max_Value_of_Sigma_2zones_Sub(this, T_e, p_max, path_cases)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      integer, parameter :: n = 500
      real(mcp), intent(in) :: T_e, p_max
      integer(kind=8) :: path_cases
      real(mcp) :: alpha_nu, E0
      real(mcp) :: Elec4U_CtrBL(1:4), Omega_K, Omega, p, rp, theta, mucosp, musinp, &
                phip, tp, sign_pr, sign_pth, a, q, l, R_p, Theta_p, Phot4k_CovBL(1:4), &
                PhotE_In_ComvFrame, Delta_p, Sigmap, Tau_nu, Del_p, Del_tau0, Del_tau_n
      integer :: i
      real(mcp) :: Tau_nu_i(0:n-1), Del_Tau_nu, Integral_jvlv3, Del_Integral_jvlv3, &
               Integral_jvlv3_0, Integral_jvlv3_n, jvlv3, kappa = 8.D-46, v_frequency, &
               Tau_scat, Del_Tau_scat, Del_Tau_scat_n, Del_Tau_scat_0
      real(mcp) :: r_ini, p_length
   
      this%Sigma_Max = zero 
      select case( path_cases ) 
      case(1)
          Tau_scat = zero  
          i = 0
          Del_p = this%p_boundary1 / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one
              !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
              Sigmap = one 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
              If( this%Sigma_Max < Del_Tau_scat )then
                  this%Sigma_Max = Del_Tau_scat
              endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > this%p_boundary1 ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = ( Tau_scat - Del_Tau_scat_0 / &
                                       two - Del_Tau_scat_n / two )*Del_p 
      case(2)
          Tau_scat = zero  
          i = 0
          Del_p = this%p_boundary1 / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one
              !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
              Sigmap = one 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
              If( this%Sigma_Max < Del_Tau_scat )then
                  this%Sigma_Max = Del_Tau_scat
              endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > this%p_boundary1 ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = ( Tau_scat - Del_Tau_scat_0 / &
                                       two - Del_Tau_scat_n / two )*Del_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Tau_scat = zero  
          i = 0
          p_length = ( this%p_boundary2 - this%p_boundary1 )
          Del_p = p_length / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one
              !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
              Sigmap = one  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > p_length ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
              ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Tau_scat = zero  
          i = 0
          p_length = ( this%p_boundary2 - this%p_boundary3 )
          Del_p = p_length / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one  
              Sigmap = one  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > p_length ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
              ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p
      case(3)
          Tau_scat = zero  
          i = 0
          Del_p = this%p_boundary1 / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one
              !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
              Sigmap = one 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
              If( this%Sigma_Max < Del_Tau_scat )then
                  this%Sigma_Max = Del_Tau_scat
              endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > this%p_boundary1 ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = ( Tau_scat - Del_Tau_scat_0 / &
                                       two - Del_Tau_scat_n / two )*Del_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Tau_scat = zero  
          i = 0
          p_length = ( this%p_boundary2 - this%p_boundary1 )
          Del_p = p_length / (n + 1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
              E0 = one
              !Sigmap = ( rp**2 + ( a * mucosp )**2 ) * rg_SUN 
              Sigmap = one  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e * this%sigma_fn(T_e, E0*this%E_ini) * E0 * Sigmap
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif 
              If ( p > p_length ) then
                  Del_Tau_scat_n = Del_Tau_scat
                  !this%p_out = p
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
              ( Tau_scat - Del_Tau_scat_0 / two - Del_Tau_scat_n / two )*Del_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end select 
      end subroutine Get_Max_Value_of_Sigma_2zones_Sub
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Photons





