      module ScatDistance
      use ScatterPhoton
      use CrossSection 
      implicit none 
      integer, parameter :: N_sigma = 2000 

      type, public, extends(ScatPhoton) :: Photon_With_ScatDistance
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_85
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_230
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_50
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_400
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_FST
          integer(kind=8) :: effect_number
          integer(kind=8) :: scatter_times
          logical :: mymethod
          real(mcp) :: NormalA
          real(mcp) :: n_e_in
          real(mcp) :: n_e_out
          real(mcp) :: n_e
          real(mcp) :: n_e0
          real(mcp) :: r_times_p
          real(mcp) :: T_e
          real(mcp) :: time_arrive_observer
          real(mcp) :: time_travel
          real(mcp) :: t_standard
          real(mcp) :: r_obs
          !real(mcp) :: R_in
          !real(mcp) :: R_out
          real(mcp) :: Important_Sampling_Const
          logical :: test_it = .FALSE.
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
          real(mcp) :: frequency_v 
          real(mcp) :: Optical_Depth_scatter 

      contains 
!******************************************************************************************************* 
          procedure, public :: sigma_fn
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !procedure, public :: Set_Cross_Section => Set_Cross_Section_sub  
          procedure, public :: Set_Cross_Section_Array_Whth_Te => Set_Cross_Section_Array_Whth_Te_sub
          !procedure, public :: Set_Cross_Section_Te => Set_Cross_Section_Te_sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          procedure, public :: get_Tau2                           =>  get_Tau2_fn
          procedure, public :: Get_scatter_distance
          procedure, public :: Get_scatter_distance2 
          procedure, public :: n_e_p => n_e_p_fn
          procedure, public :: get_Tau2_At_out_zone => get_Tau2_At_out_zone_fn
          procedure, public :: Integration_Of_OpticalDepth
          procedure, public :: Get_Max_Value_of_Sigma_2zones3 => &
                                Get_Max_Value_of_Sigma_2zones3_sub
          procedure, public :: RandomNum2p_case2   =>   RandomNum2p_case2_fn
          procedure, public :: Integrals   =>    Integrals_fn
          procedure, public :: RandomNum2p_case3   =>   RandomNum2p_case3_fn
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: RandNum2p_case1_Finit_Space   =>  &
                               RandNum2p_case1_Finit_Space_fn
          procedure, public :: RandNum2p_case2_Finit_Space   =>  &
                               RandNum2p_case2_Finit_Space_fn
          procedure, public :: RandNum2p_case3_Finit_Space   =>  &
                               RandNum2p_case3_Finit_Space_fn
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !private :: Set_Cross_Section_sub  
      !private :: Set_Cross_Section_Te_sub
      private :: Set_Cross_Section_Array_Whth_Te_sub 
      private :: n_e_p_fn
      private :: get_Tau2_At_out_zone_fn
      private :: RandomNum2p_case3_fn
      private :: RandomNum2p_case2_fn
      private :: Integrals_fn
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: RandNum2p_case1_Finit_Space_fn
      private :: RandNum2p_case2_Finit_Space_fn
      private :: RandNum2p_case3_Finit_Space_fn
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 
          

!*******************************************************************************************************
      real(mcp) function Integrals_fn(this, p1, p2, T_e, coef)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: p1, p2, T_e, coef 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, p)
      Integrals_fn = coef * dlog( ( p2 + this%r_times_p + dsqrt( p2**2 + this%r_ini**2 + &
                      two * p2 * this%r_times_p ) ) / &
                    ( p1 + this%r_times_p + dsqrt( p1**2 + this%r_ini**2 + &
                      two * p1 * this%r_times_p ) ) ) 
      end function Integrals_fn

!*******************************************************************************************************
      real(mcp) function RandomNum2p_case2_fn(this, r, p1, p2, p3, T_e, NormalA)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: r, p1, p2, p3, T_e
      real(mcp), intent(out) :: NormalA
      real(mcp) :: PI0 
      real(mcp) :: w1, w2, t, x, s1, s2, y1, y2, coef, k2, I_infty, C2 
      real(mcp) :: Sigma_I, Sigma_III, tau, tau1, tau2, tau3, eta, ksi
!*******************************************************************************************************

      Sigma_I = this%n_e_in * this%sigma_fn(T_e, this%E_ini, zero)
      !Sigma_III = this%n_e0 * this%sigma_fn(T_e, this%E_ini, zero) * &
      !                               ( this%R_in / this%R_out )**1.5D0
      If( this%r_ini > dabs(this%r_times_p) )then  
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = this%Integrals(zero, p1, T_e, coef)
          tau2 = (p2 - p1) * Sigma_I
          tau3 = this%Integrals(p2, p3, T_e, coef)
          NormalA = one! - dexp( - tau1 - tau2 - tau3 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              ksi = dexp( eta / coef )
              RandomNum2p_case2_fn = ( this%r_times_p - this%r_ini + &
                            ksi**2*( this%r_times_p + this%r_ini ) - &
                            two*ksi*this%r_times_p ) / two / ksi 
          else if( eta < tau1 + tau2 )then
              RandomNum2p_case2_fn = p1 + ( eta - tau1 ) / Sigma_I
          else 
              ksi = dexp( (eta - tau1 - tau2 ) / coef )
              C2 = p2 + this%r_times_p + dsqrt( p2**2 + this%r_ini**2 + two*this%r_times_p*p2 )
              RandomNum2p_case2_fn = - this%r_times_p + &
                     ( C2 * ksi - ( this%r_ini**2 - this%r_times_p**2 ) / C2 / ksi ) / two 
          !else
          !    RandomNum2p_case2_fn = p3 + ( eta - tau1 - tau2 - tau3 ) / Sigma_III
          end if 
      else 
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero) 
          tau1 = coef * dlog( this%r_ini / this%R_in )
          tau2 = (p2 - p1) * Sigma_I
          tau3 = coef * dlog( this%R_out / this%R_in )
          NormalA = one! - dexp( - tau1 - tau2 - tau3 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              !PI0 = one / dsqrt( this%r_ini ) + eta / coef
              !RandomNum2p_case2_fn = this%r_ini - one / PI0**2
              RandomNum2p_case2_fn = this%r_ini * ( one - dexp( - eta / coef ) )
          else if( eta < tau1 + tau2 )then
              RandomNum2p_case2_fn = p1 + ( eta - tau1 ) / Sigma_I
          else if( eta < tau1 + tau2 + tau3 )then 
              !PI0 = one / dsqrt( this%r_ini ) - ( eta - tau1 - tau2 ) / coef
              !RandomNum2p_case2_fn = one / PI0**2 - this%r_ini + p2
              RandomNum2p_case2_fn = this%R_in * dexp( ( eta - tau1 - tau2 ) / coef ) - this%R_in + p2
          else
              RandomNum2p_case2_fn = p3 + ( eta - tau1 - tau2 - tau3 ) / Sigma_III
          endif
      endif
      end function RandomNum2p_case2_fn
   
!*******************************************************************************************************
      real(mcp) function RandomNum2p_case3_fn(this, r, p1, p2, T_e, NormalA)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: r, p1, p2, T_e
      real(mcp), intent(out) :: NormalA
      real(mcp) :: PI0, C2
      real(mcp) :: a, w1, w2, t, x, s1, s2, y1, y2, coef, k2, I_infty 
      real(mcp) :: Sigma_I, Sigma_III, tau, tau1, tau2, tau3, eta, ksi
!*******************************************************************************************************

      Sigma_I = this%n_e_in * this%sigma_fn(T_e, this%E_ini, zero)
      !Sigma_III = this%n_e0 * this%sigma_fn(T_e, this%E_ini, zero) * &
      !                               ( this%R_in / this%R_out )**1.5D0
      If( this%r_ini > dabs(this%r_times_p) )then  
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = p1 * Sigma_I
          tau2 = this%Integrals(p1, p2, T_e, coef)
          NormalA = one! - dexp( - tau1 - tau2 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              RandomNum2p_case3_fn = eta / Sigma_I
          else 
              ksi = dexp( (eta - tau1) / coef )    
              C2 = p1 + this%r_times_p + dsqrt( p1**2 + this%r_ini**2 + two*this%r_times_p*p1 )
              RandomNum2p_case3_fn = - this%r_times_p + &
                     ( C2 * ksi - ( this%r_ini**2 - this%r_times_p**2 ) / C2 / ksi ) / two
          end if 
      else 
          !coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero) * &
          !       dsqrt( this%R_in ) * two 
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = p1 * Sigma_I
          tau2 = coef * dlog( this%R_out / this%R_in )!( one / dsqrt(this%R_in) - one / dsqrt(this%R_out) )
          NormalA = one! - dexp( - tau1 - tau2 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              RandomNum2p_case3_fn = eta / Sigma_I
          else if( eta < tau1 + tau2 )then
              !PI0 = one / dsqrt( this%R_in ) - ( eta - tau1 ) / coef
              !RandomNum2p_case3_fn = !one / PI0**2 - this%r_ini + p1
              RandomNum2p_case3_fn = this%R_in * dexp( (eta - tau1)/coef ) - this%R_in + p1
          else
              RandomNum2p_case3_fn = p2 + ( eta - tau1 - tau2 ) / Sigma_III
          endif
          !write(*,*)'ssdf=', one - r * NormalA, NormalA, tau1, tau2, coef, Sigma_I
      endif 
      end function RandomNum2p_case3_fn

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance( this, T_e )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: T_e
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out, eta, Sigma_I
      real(mcp) :: p_c, p_out1, p_1, p_2, delta
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i, path_cases
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  First ditermine which zone the Photon_With_ScatDistance is included in.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !path_cases = 10
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) ) 
      !call this%Get_Max_Value_of_Sigma_2zones3( T_e, p_out1, path_cases )  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( this%r_ini /= zero )then
          Sigma_I = this%n_e * this%sigma_fn(T_e, this%E_ini, zero) + &
                this%r_times_p / this%r_ini * this%Important_Sampling_Const
      else
          Sigma_I = this%n_e * this%sigma_fn(T_e, this%E_ini, zero) + &
                    this%Important_Sampling_Const
      endif
      Sigma_I = dabs(Sigma_I)
      this%Optical_Depth_scatter = p_out1 * Sigma_I
      r1 = ranmar()
      p = - dlog( one - r1 ) / Sigma_I
      if ( p > p_out1 )then
          this%At_outer_Shell = .True. 
      endif 
      Get_scatter_distance = p 
      if( isnan(p) )then
          write(*,*)'sdf==1', Sigma_I, this%n_e * this%sigma_fn(T_e, this%E_ini, zero),&
                this%r_times_p / this%r_ini , this%Important_Sampling_Const
          stop
      endif
      If (Get_scatter_distance < zero) then
          write(*,*)'sdf==2',Get_scatter_distance,Sigma_I, this%n_e * this%sigma_fn(T_e, this%E_ini, zero),&
                this%r_times_p / this%r_ini , this%Important_Sampling_Const
          stop
      endif
      If (Get_scatter_distance == zero) then
          write(*,*)'ends=',Get_scatter_distance, path_cases, eta, p, delta,this%At_inner_Shell
          write(*,*)'ends2=', this%r_ini, this%R_in, this%R_out
          stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      end function Get_scatter_distance
 
!*******************************************************************************************************
      real(mcp) function RandNum2p_case1_Finit_Space_fn(this, r, p1, T_e, NormalA)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: r, p1, T_e
      real(mcp), intent(out) :: NormalA
      real(mcp) :: w1, w2, t, x, s1, s2, y1, y2, coef, k2, I_infty, PI0  
      real(mcp) :: Sigma_I, Sigma_III, tau, tau1, tau2, tau3, eta, ksi
!*******************************************************************************************************

      !Sigma_I = this%n_e_in * this%sigma_fn(T_e, this%E_ini, zero) 
      If( this%r_ini > dabs(this%r_times_p) )then  
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = this%Integrals(zero, p1, T_e, coef) 
          NormalA = one - dexp( - tau1 )
          eta = - dlog( one - r * NormalA ) 
          ksi = dexp( eta / coef )
          RandNum2p_case1_Finit_Space_fn = ( this%r_times_p - this%r_ini + &
                            ksi**2*( this%r_times_p + this%r_ini ) - &
                            two*ksi*this%r_times_p ) / two / ksi  
      else 
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero) * &
                 dsqrt( this%R_in ) * two
          tau1 = coef * ( one / dsqrt(this%r_ini) - one / dsqrt(this%R_out) ) 
          NormalA = one - dexp( - tau1 )
          eta = - dlog( one - r * NormalA ) 
          PI0 = one / dsqrt( this%r_ini ) - eta / coef
          RandNum2p_case1_Finit_Space_fn = one / PI0**2 - this%r_ini
      endif
      end function RandNum2p_case1_Finit_Space_fn
    
!*******************************************************************************************************
      real(mcp) function RandNum2p_case2_Finit_Space_fn(this, r, p1, p2, p3, T_e, NormalA)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: r, p1, p2, p3, T_e
      real(mcp), intent(out) :: NormalA
      real(mcp) :: PI0, C2
      real(mcp) :: w1, w2, t, x, s1, s2, y1, y2, coef, k2, I_infty 
      real(mcp) :: Sigma_I, Sigma_III, tau, tau1, tau2, tau3, eta, ksi
!*******************************************************************************************************

      Sigma_I = this%n_e_in * this%sigma_fn(T_e, this%E_ini, zero) 
      If( this%r_ini > dabs(this%r_times_p) )then  
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = this%Integrals(zero, p1, T_e, coef)
          tau2 = (p2 - p1) * Sigma_I
          tau3 = this%Integrals(p2, p3, T_e, coef)
          NormalA = one - dexp( - tau1 - tau2 - tau3 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              ksi = dexp( eta / coef )
              RandNum2p_case2_Finit_Space_fn = ( this%r_times_p - this%r_ini + &
                            ksi**2*( this%r_times_p + this%r_ini ) - &
                            two*ksi*this%r_times_p ) / two / ksi 
          else if( eta < tau1 + tau2 )then
              RandNum2p_case2_Finit_Space_fn = p1 + ( eta - tau1 ) / Sigma_I
          else 
              ksi = dexp( (eta - tau1 - tau2 ) / coef ) 
              C2 = p2 + this%r_times_p + dsqrt( p2**2 + this%r_ini**2 + two*this%r_times_p*p2 )
              RandNum2p_case2_Finit_Space_fn = - this%r_times_p + &
                     ( C2 * ksi - ( this%r_ini**2 - this%r_times_p**2 ) / C2 / ksi ) / two
          end if 
      else 
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero) * &
                 dsqrt( this%R_in ) * two
          tau1 = coef * ( one / dsqrt(this%R_in) - one / dsqrt(this%r_ini) )
          tau2 = (p2 - p1) * Sigma_I
          tau3 = coef * ( one / dsqrt(this%R_in) - one / dsqrt(this%R_out) )
          NormalA = one - dexp( - tau1 - tau2 - tau3 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              PI0 = one / dsqrt( this%r_ini ) + eta / coef
              RandNum2p_case2_Finit_Space_fn =this%r_ini - one / PI0**2
          else if( eta < tau1 + tau2 )then
              RandNum2p_case2_Finit_Space_fn = p1 + ( eta - tau1 ) / Sigma_I
          else if( eta < tau1 + tau2 + tau3 )then 
              PI0 = one / dsqrt( this%r_ini ) - ( eta - tau1 - tau2 ) / coef
              RandNum2p_case2_Finit_Space_fn = one / PI0**2 - this%r_ini + p2 
          endif
      endif
      end function RandNum2p_case2_Finit_Space_fn
   
!*******************************************************************************************************
      real(mcp) function RandNum2p_case3_Finit_Space_fn(this, r, p1, p2, T_e, NormalA)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp), intent(in) :: r, p1, p2, T_e
      real(mcp), intent(out) :: NormalA 
      real(mcp) :: a, w1, w2, t, x, s1, s2, y1, y2, coef, k2, I_infty, PI0
      real(mcp) :: Sigma_I, tau, tau1, tau2, tau3, eta, ksi, C2
!*******************************************************************************************************

      Sigma_I = this%n_e_in * this%sigma_fn(T_e, this%E_ini, zero) 
      If( this%r_ini > dabs(this%r_times_p) )then  
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero)
          tau1 = p1 * Sigma_I
          tau2 = this%Integrals(p1, p2, T_e, coef)
          NormalA = one - dexp( - tau1 - tau2 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              RandNum2p_case3_Finit_Space_fn = eta / Sigma_I
          else 
              ksi = dexp( (eta - tau1) / coef )  
              C2 = p1 + this%r_times_p + dsqrt( p1**2 + this%r_ini**2 + two*this%r_times_p*p1 )
              RandNum2p_case3_Finit_Space_fn = - this%r_times_p + &
                     ( C2 * ksi - ( this%r_ini**2 - this%r_times_p**2 ) / C2 / ksi ) / two
          end if 
      else 
          coef = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, zero) 
          tau1 = p1 * Sigma_I
          tau2 = coef * dlog( this%R_out / this%R_in )
          NormalA = one - dexp( - tau1 - tau2 )
          eta = - dlog( one - r * NormalA )
          If( eta < tau1 )then
              RandNum2p_case3_Finit_Space_fn = eta / Sigma_I
          else if( eta < tau1 + tau2 )then
              !PI0 = one / dsqrt( this%R_in ) - ( eta - tau1 ) / coef
              !RandNum2p_case3_Finit_Space_fn = one / PI0**2 - this%r_ini + p1 
              RandNum2p_case3_Finit_Space_fn = this%R_in * ( dexp( (eta - tau1)/coef ) - one ) + p1
          endif 
      endif 
      end function RandNum2p_case3_Finit_Space_fn


!*******************************************************************************************************
      real(mcp) function Get_scatter_distance2( this, T_e1 )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: T_e1 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_out1, Sigma_I
      integer(kind=8) :: i, path_cases
   
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )
      !call this%Get_Max_Value_of_Sigma_2zones3( T_e1, p_out1, path_cases )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Sigma_I = this%n_e * this%sigma_fn(T_e1, this%E_ini, zero)
      this%Optical_Depth_scatter = p_out1 * Sigma_I
      this%NormalA = one - dexp( - p_out1 * Sigma_I )
      r1 = ranmar()
      Get_scatter_distance2 = - dlog( one - r1 * this%NormalA ) / Sigma_I
      !write(*,*)'sdf1==', Get_scatter_distance2, this%NormalA, Sigma_I
      !stop 
      !Get_scatter_distance2 = this%RandNum2p_case1_Finit_Space(r1, p_out1, T_e1, this%NormalA)
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance2, p_out1, Sigma_I, this%R_out
          write(*,*)'sdf3==',  this%r_times_p, this%r_ini, this%R_out
      endif
      if(p_out1 < zero .or. this%r_ini > this%R_out)then
          write(*,*)'sdf1==',Get_scatter_distance2, r1, this%NormalA,  Sigma_I, this%r_times_p, &
          this%r_ini, this%R_out 
      endif
      if( Get_scatter_distance2 > p_out1 )then
          write(*,*)'I111',Get_scatter_distance2, p_out1, Get_scatter_distance2 - p_out1
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance2 < zero) then
          write(*,*)'sdf==',Get_scatter_distance2, r1, this%NormalA,  Sigma_I, this%r_times_p, &
          this%r_ini, this%R_out
          stop
      endif
      If (Get_scatter_distance2 == zero) then
          write(*,*)'ends=',Get_scatter_distance2, r1, this%NormalA,  Sigma_I
          stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function Get_scatter_distance2
  


!*******************************************************************************************************
      real(mcp) function Get_scatter_distance3( this, T_e1 )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: T_e1 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_out1, Sigma_I
      integer(kind=8) :: i, path_cases
   
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )
      !call this%Get_Max_Value_of_Sigma_2zones3( T_e1, p_out1, path_cases )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Sigma_I = this%n_e * this%sigma_fn(T_e1, this%E_ini, zero)
      this%Optical_Depth_scatter = p_out1 * Sigma_I
      this%NormalA = one - dexp( - p_out1 * Sigma_I )
      r1 = ranmar()
      Get_scatter_distance3 = - dlog( one - r1 * this%NormalA ) / Sigma_I
      !write(*,*)'sdf1==', Get_scatter_distance2, this%NormalA, Sigma_I
      !stop 
      !Get_scatter_distance2 = this%RandNum2p_case1_Finit_Space(r1, p_out1, T_e1, this%NormalA)
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance3, p_out1, Sigma_I, this%R_out
          write(*,*)'sdf3==',  this%r_times_p, this%r_ini, this%R_out
      endif
      if(p_out1 < zero .or. this%r_ini > this%R_out)then
          write(*,*)'sdf1==',Get_scatter_distance3, r1, this%NormalA,  Sigma_I, this%r_times_p, &
          this%r_ini, this%R_out 
      endif
      if( Get_scatter_distance3 > p_out1 )then
          write(*,*)'I111',Get_scatter_distance3, p_out1, Get_scatter_distance3 - p_out1
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance3 < zero) then
          write(*,*)'sdf==',Get_scatter_distance3, r1, this%NormalA,  Sigma_I, this%r_times_p, &
          this%r_ini, this%R_out
          stop
      endif
      If (Get_scatter_distance3 == zero) then
          write(*,*)'ends=',Get_scatter_distance3, r1, this%NormalA,  Sigma_I
          stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function Get_scatter_distance3
  
!*******************************************************************************************************
      real(mcp) function get_Tau2_fn(this, samp, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: samp, T_e
      real(mcp) :: Sigmap, mup, musinp, rp, phip, timep, &
                 sigp, P_mu_Multipl_U_mu, PhotE_In_Elec_CF, p
    
      PhotE_In_Elec_CF = this%E_ini
   
      get_Tau2_fn = this%n_e_in * this%sigma_fn(T_e, DABS(PhotE_In_Elec_CF), p)! &
      !                 * DABS(PhotE_In_Elec_CF) / this%E_ini * Sigmap * rg_SUN
      end function get_Tau2_fn

!*******************************************************************************************************
      real(mcp) function get_Tau2_At_out_zone_fn(this, samp, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: samp, T_e 
      
      get_Tau2_At_out_zone_fn = this%n_e_p( samp ) * this%sigma_fn(T_e, DABS(this%E_ini), samp)! &
      !                 * DABS(PhotE_In_Elec_CF) / this%E_ini * Sigmap * rg_SUN
      end function get_Tau2_At_out_zone_fn
 

!*******************************************************************************************************
      real(mcp) function n_e_p_fn(this, p)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: p
   
      n_e_p_fn = this%n_e0 * this%R_in / dsqrt( this%r_ini**2 + p**2 + two*p*this%r_times_p ) 
      end function n_e_p_fn
  
 
!*******************************************************************************************************
      subroutine Set_Cross_Section_Array_Whth_Te_sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      !real(mcp) :: T_e
      real(mcp) :: temp, Ephoton
      real(mcp) :: dindexE, NorN
      integer(kind = 8) :: i, j, k, N, istat
  
      open(unit=18, file = this%CrossSectFileName, status = "old", &
                           action = "read", iostat = istat)
 
      if (istat == 0) then
              if( this%my_ID == this%num_process - 1 )then
                  write(unit = *, fmt = *)'************************************************************'    
                  write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... ' 
                  write(unit = *, fmt = *)'************************************************************'
              endif  
          do i = 0, N_sigma
              read(unit = 18, fmt = *)this%sigmaaTeE_FST(i)
          enddo
      else
          open(unit=19, file = this%CrossSectFileName, status = "replace", &
                                   action = "write", iostat = istat)
          !call Set_xi_wi_all()
          call gauleg_x_w( -one, one, x1000, w1000, 1000 )
          call gaulag( x0la362, w0la362, 362, zero ) 
 
          if (istat == 0) then
              if( this%my_ID == this%num_process - 1 )then
                  write(unit = *, fmt = *)'************************************************************'  
                  write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
                  write(unit = *, fmt = *)'************************************************************'
              endif 
              do i = 0, N_sigma 
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_FST(i) = gama_Integration( this%T_e, Ephoton, &
                                   x1000, w1000, 1000, x0la362, w0la362, 362 )
 
                  if( mod(i, 100)==0 )this%sigmaaTeE_400(i) = sigma_a( this%T_e, Ephoton )
                  write(unit = 19, fmt = *)this%sigmaaTeE_FST(i)
                  if( mod(i, 100)==0 )write(unit = *, fmt = *)i, 'Sigma_Te = ', &
                                     this%sigmaaTeE_FST(i), this%sigmaaTeE_400(i)
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File Open Failed. The code have to Stop!!', istat
              Stop
          endif
          close(unit=19)
          write(unit = *, fmt = *)trim('The  ')//trim(this%CrossSectFileName)//&
                                  trim(' File has been successfully created!!!')     
      endif
      close(unit=18) 
      end subroutine Set_Cross_Section_Array_Whth_Te_sub


!*******************************************************************************************************
      function sigma_fn(this, T_e, Ephoton, p)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp) :: sigma_fn, Ei, Ei1
      real(mcp), dimension(0:N_sigma) :: sigmaaTeE
      real(mcp), intent(in) :: T_e, Ephoton, p
      integer(kind = 8) :: i
      real(mcp) :: Delta_p, Theta_p, R_p
      real(mcp) :: Ep, sign_pth
      !integer, intent(in) :: n
!*******************************************************************************************************
    
      Ep = Ephoton
      i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      if ( dlog10( Ep ) < this%logE_low ) then
          Ep = 10.D0**this%logE_low
          i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      endif
      if ( dlog10( Ep ) > this%logE_up ) then
          Ep = 10.D0**this%logE_up
          i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE ) - 1
      endif
      !write(unit = *, fmt = *)'ttsfsfds!!',Ephoton, dlog10( Ephoton ), this%logE_low,this%dindexE, i
      if (i<0) then
          write(*,*)'sigma_fn = =', dlog10( Ephoton ), Ephoton, this%logE_low, this%dindexE,&
                ( dlog10( Ephoton ) - this%logE_low ) / this%dindexE!,&
          write(unit = *, fmt = *)'************************************************************'  
          write(unit = *, fmt = *)'************************************************************' 
      endif 

      Ei = 10**( this%logE_low + this%dindexE * i )
      Ei1 = 10**( this%logE_low + this%dindexE * (i + 1) )
      If ( T_e == 4.D0 * mec2 ) then  
          sigma_fn = ( this%sigmaaTeE_400(i + 1) - this%sigmaaTeE_400(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_400(i)  
      Else If (T_e == 85.D-3) then 
          sigma_fn = ( this%sigmaaTeE_85(i + 1) - this%sigmaaTeE_85(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_85(i)  
      Else If (T_e == 230.D-3) then 
          sigma_fn = ( this%sigmaaTeE_230(i + 1) - this%sigmaaTeE_230(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_230(i)  
      Endif
      !sigma_fn = ( sigmaaTeE(i + 1) - sigmaaTeE(i) ) / &
      !           ( Ei1 - Ei ) * (Ep - Ei) + sigmaaTeE(i)  
      end function sigma_fn
 
!*******************************************************************************************************
!*******************************************************************************************************
!**  SubRoutines In terms of the Electron that sacttering or interacting with the Photon_With_ScatDistance  **************
!*******************************************************************************************************
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      real(mcp) function Integration_Of_OpticalDepth( this, a, b, n, T_e, p_bias )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: a, b, T_e, p_bias
      real(mcp) :: Tau_scat, p_length, Del_p, p, Del_Tau_scat, &
                   Del_Tau_scat_0, Del_Tau_scat_n
      integer, intent(in) :: n
      integer :: i

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          Tau_scat = zero  
          i = 0 
          Del_p = ( b - a ) / (n + 1) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              p = Del_p*i 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Del_Tau_scat = this%n_e_p( p + p_bias ) * this%sigma_fn(T_e, this%E_ini, p)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              Tau_scat = Tau_scat + Del_Tau_scat    
              If ( i == 0 ) then 
                  Del_Tau_scat_0 = Del_Tau_scat
              Endif  
              If ( i > n ) then
                  Del_Tau_scat_n = Del_Tau_scat 
                  exit
              Endif  
              i = i + 1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Enddo 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          Integration_Of_OpticalDepth = ( Tau_scat - Del_Tau_scat_0 / &
                                     two - Del_Tau_scat_n / two )*Del_p
      end function Integration_Of_OpticalDepth


!*******************************************************************************************************
      subroutine Get_Max_Value_of_Sigma_2zones3_Sub(this, T_e1, p_max, path_cases)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      integer :: n = 50
      real(mcp), intent(in) :: T_e1, p_max
      integer(kind=8) :: path_cases
      real(mcp) :: alpha_nu, Tau_scat, Del_Tau_scat, Del_Tau_scat_0, Del_Tau_scat_n
      real(mcp) :: p, Del_p, E0, Sigmap, temp_I
      integer :: i 
      real(mcp) :: r_ini, p_length, T_e, p_bias
   
      this%Sigma_Max = zero 
      T_e = 50.D-3 
      select case( path_cases ) 
      case(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          this%Sigma_Max = this%n_e_p( zero ) * this%sigma_fn(T_e, this%E_ini, zero)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          n = 2000   
          !this%Optical_Depth_scatter = this%Integration_Of_OpticalDepth(zero, &
          !                              this%p_boundary1, n, T_e, zero)
          !write(*,*)'sss1111=',this%Optical_Depth_scatter  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
          this%Optical_Depth_scatter = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, p) * &
            dlog( ( this%p_boundary1 + this%r_times_p + dsqrt( this%p_boundary1**2 + this%r_ini**2 + &
                      two * this%p_boundary1 * this%r_times_p ) ) / (this%r_ini + this%r_times_p) )
          !write(*,*)'sss3333=',this%Optical_Depth_scatter  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      case(2) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          this%Sigma_Max = this%n_e_p( this%p_boundary1 ) * &
                           this%sigma_fn(T_e, this%E_ini, this%p_boundary1)  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          n = 2000   
          !this%Optical_Depth_scatter = this%Integration_Of_OpticalDepth(zero, &
          !                                       this%p_boundary1, n, T_e, zero) 
          this%Optical_Depth_scatter = this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, p) * &
            dlog( ( this%p_boundary1 + this%r_times_p + dsqrt( this%p_boundary1**2 + this%r_ini**2 + &
                      two * this%p_boundary1 * this%r_times_p ) ) / (this%r_ini + this%r_times_p) ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          p_length = ( this%p_boundary2 - this%p_boundary1 )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
                       this%n_e_in * this%sigma_fn(T_e, this%E_ini, p_length) * p_length
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          p_length = ( this%p_boundary3 - this%p_boundary2 )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !temp_I = this%Integration_Of_OpticalDepth(this%p_boundary2, &
          !               this%p_boundary3, n, T_e, this%p_boundary2)
          !this%Optical_Depth_scatter =  this%Optical_Depth_scatter + temp_I
          !write(*,*)'sdf1=',temp_I
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
                          this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, p) * &
            dlog( ( this%p_boundary3 + this%r_times_p + dsqrt( this%p_boundary3**2 + this%r_ini**2 + &
                      two * this%p_boundary3 * this%r_times_p ) ) / &
                    ( this%p_boundary2 + this%r_times_p + dsqrt( this%p_boundary2**2 + this%r_ini**2 + &
                      two * this%p_boundary2 * this%r_times_p ) )  &
                )  
          !write(*,*)'sdf2=',temp_I
      case(3) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
          this%Sigma_Max = this%n_e_in * this%sigma_fn(T_e, this%E_ini, p)
          this%Optical_Depth_scatter = this%n_e_in * this%sigma_fn(T_e, this%E_ini, p) * this%p_boundary1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
          p_length = ( this%p_boundary2 - this%p_boundary1 ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !this%Optical_Depth_scatter =  this%Optical_Depth_scatter + &
          !    this%Integration_Of_OpticalDepth(this%p_boundary1, this%p_boundary2, n, &
          !                                      T_e, this%p_boundary1) 
          this%Optical_Depth_scatter = this%Optical_Depth_scatter + &
                          this%n_e0 * this%R_in * this%sigma_fn(T_e, this%E_ini, p) * &
            dlog( ( this%p_boundary2 + this%r_times_p + dsqrt( this%p_boundary2**2 + this%r_ini**2 + &
                      two * this%p_boundary2 * this%r_times_p ) ) / &
                    ( this%p_boundary1 + this%r_times_p + dsqrt( this%p_boundary1**2 + this%r_ini**2 + &
                      two * this%p_boundary1 * this%r_times_p ) )  &
                ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end select 
      end subroutine Get_Max_Value_of_Sigma_2zones3_Sub

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module ScatDistance





