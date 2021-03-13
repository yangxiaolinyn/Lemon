      module Photons 
      use ScatDistance
      !use CrossSection 
      implicit none 

      type, public, extends(Photon_With_ScatDistance) :: Photon
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(1:500)
          real(mcp) :: v_L_v_i_ET(1:500, 1:1000)
          real(mcp) :: nu_low
          real(mcp) :: nu_up
          !real(mcp) :: Optical_Depth_scatter 

      contains 
!*******************************************************************************************************
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !procedure, public :: Set_Cross_Section                 => Set_Cross_Section_sub  
          !procedure, public :: Get_bias_parameter              => Get_bias_parameter_Sub 
          !procedure, public :: Get_Optical_Depth_At_p  =>      Get_Optical_Depth_At_p_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          procedure, public :: r_p2
          !procedure, public :: Get_Max_Value_of_Sigma_2zones => &
          !                      Get_Max_Value_of_Sigma_2zones_Sub
          !procedure, public :: Get_Max_Value_of_Sigma_2zones2 => &
          !                      Get_Max_Value_of_Sigma_2zones2_Sub
          procedure, public :: Calc_Phot_Informations_At_Observor_2zones => &
                               Calc_Phot_Informations_At_Observor_2zones_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Calc_Phot_Informations_At_Observor_2zones_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains   
 
 
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub(this, cases)
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      integer, intent(in) :: cases  
      real(mcp), parameter :: al0 = -50.D0, be0 = -50.D0
      real(mcp) :: gama_S, beta_T, f_the, f_phi, Sin_The_Obs, Cos_The_Obs
      real(mcp) :: mup, mup0, p_max, p_temp, sign_p_th
      real(mcp) :: Q_obs, U_obs, Psi_Obs, rdp
      real(mcp) :: alpha, beta, delta_th = two / 100.D0, mu_obs 
      real(mcp) :: delta_al = -two*al0 / 400.D0, delta_be = -two*be0 / 400.D0
      real(mcp) :: sign_pth, r_times_p, r_ini2, p_ini_obs, index_i, a, b, del_x
      real(mcp) :: index_i1, a1, b1, del_x1, angle
      integer :: i, j, h, k, mu_i, N_low, N_low1

      a = -10.D0
      b = 1.D0
      del_x = (b-a) / 1000.D0
      N_low = -floor(a/del_x)
      p_ini_obs = - this%r_times_p + dsqrt( this%r_times_p**2 + this%r_obs**2 - this%r_ini**2 )
      this%time_arrive_observer = this%time_travel + p_ini_obs / Cv 
      If( cases == 1 )then
          a1 = dlog10( this%nu_low )
          b1 = dlog10( this%nu_up * 4.D8 ) !dlog10( this%nu_up )
          del_x1 = (b1 - a1) / 500.D0
          N_low1 =  floor( a1 / del_x1 )
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          index_i1 = dlog10( this%frequency_v / this%nu_low )
           !write(*, *)'ss=', this%nu_low, this%frequency_v, this%nu_up * 4.D8
           !write(*, *)'ss=', a1, b1
          if( this%frequency_v >= this%nu_up * 4.D8 .or. this%frequency_v <= this%nu_low )then 
              return
          endif
          !if( index_i1 >= b1 .or. index_i1 <= a1 )then 
          !    return
          !endif 
          k = floor( index_i1 / del_x1 ) + 1
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * &
                            dexp( - this%Optical_Depth_scatter )
          !If( DABS( this%Phot4k_CovCF_ini(1) ) < 10.D-3 .OR. &
          !    DABS( this%Phot4k_CovCF_ini(1) ) > 20.D-3 )return
          !index_i = dlog10( this%time_arrive_observer - this%t_standard ) 
          !if( index_i > b .or. index_i < a )then
              !write(*,*)'dddd===', index_i
              !return
          !endif
          !i = floor( index_i / del_x ) + 2 + N_low
          !If (i > 1000) then
              !i=1000
              !write(*,*)'sss2=', i, 
              !return
          !endif
          !this%v_L_v_i_ET(k, i) = this%v_L_v_i_ET(k, i) + this%w_ini * 0.01D0* &
          !                  dexp( - this%Optical_Depth_scatter )! * DABS( this%Phot4k_CovCF_ini(1) )
          !write(*,*)'i,k1',i, k, this%time_travel,  p_ini_obs / Cv  
          !this%time_arrive_observer = zero
      else If ( cases == 2 ) then 
          a1 = dlog10( this%nu_low )
          b1 = dlog10( this%nu_up ) !dlog10( this%nu_up )
          del_x1 = (b1 - a1) / 500.D0
          N_low1 =  floor( a1 / del_x1 )
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
          index_i1 = dlog10( this%frequency_v / this%nu_low )
          !write(*, *)'ss1=', this%nu_low, this%frequency_v, this%nu_up * 4.D8
          !write(*, *)'ss2=', this%w_ini
          !write(*, *)'ss=', a1, b1
          if( this%frequency_v >= this%nu_up .or. this%frequency_v <= this%nu_low )then 
              return
          endif
          k = floor( index_i1 / del_x1 ) + 1
          this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini
          !this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini 
          If( DABS( this%Phot4k_CovCF_ini(1) ) < 10.D-3 .OR. &
              DABS( this%Phot4k_CovCF_ini(1) ) > 20.D-3 )return
          index_i = dlog10( this%time_arrive_observer - this%t_standard ) 
          !index_i = dlog10( this%time_arrive_observer ) 
          if( index_i > b .or. index_i < a )then
              !write(*,*)'dddd===', index_i
              return
          endif
          !i = floor( index_i / del_x ) + 2 + N_low
          If (i > 1000) then
              i=1000
              !write(*,*)'sss2=', i, 
              return
          endif
          !this%v_L_v_i_ET(k, i) = this%v_L_v_i_ET(k, i) + this%w_ini! * this%Phot4k_CovCF_ini(1) 
      Endif 
      return
      end subroutine Calc_Phot_Informations_At_Observor_2zones_Sub
!******************************************************************************************************* 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      function r_p2( this, r_ini, vector_of_momentum, p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      class(Photon) :: this
      real(mcp), intent(in) :: r_ini(1:3), vector_of_momentum(1:3), p
      real(mcp) :: vector_p(1:3), r_p2
 
      vector_p = r_ini + p * vector_of_momentum
      this%vector_of_position = vector_p
      !write(*,*)'ddd==', vector_p, p, vector_of_momentum, r_ini
      r_p2 = dsqrt( vector_p(1)**2 + vector_p(2)**2 + vector_p(3)**2 )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function r_p2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!*******************************************************************************************************
 
      end module Photons





