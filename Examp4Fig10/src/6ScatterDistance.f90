      module ScatDistance_FlatSP
      use ScatterPhoton 
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type, public, extends(ScatPhoton) :: Photon_With_ScatDistance_FlatSP 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_85
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_230
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_50
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_400
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_FST
          integer(kind=8) :: effect_number
          integer(kind=8) :: scatter_times
          logical :: mymethod
          real(mcp) :: NormalA
          real(mcp) :: w_ini0
          real(mcp) :: n_e_in
          real(mcp) :: n_e_out
          real(mcp) :: n_e
          real(mcp) :: n_e0
          real(mcp) :: r_times_p
          real(mcp) :: T_e
          real(mcp) :: time_arrive_observer
          real(mcp) :: time_travel
          real(mcp) :: t_standard
          !real(mcp) :: r_obs
          real(mcp) :: R_in
          real(mcp) :: R_out
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
          integer :: cases, InterSection_Cases
          real(mcp) :: Sigma_I0
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          real(mcp) :: Z_max
          real(mcp) :: ratio
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      contains 
!*******************************************************************************************************   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Get_scatter_distance_IQ    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance_FlatSP
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains  
!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_IQ( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA, Temp_log
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_out1, Sigma_I, tau1, eta
      integer(kind=8) :: i, path_cases
    
      if( this%mu_zp_ini > zero )then

          p_out1 =  this%z_tau / this%mu_zp_ini
          this%InterSection_Cases = 1
          this%Optical_Depth_scatter = p_out1 
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA ) 
          Get_scatter_distance_IQ = this%z_tau + Temp_log * this%mu_zp_ini
          !write(*, *)'ff11=', Get_scatter_distance4, p_out1, this%Z_max, this%z_ini, &
          !             this%mu_zp_ini, dexp( - p_out1 * Sigma_I )

      else if( this%mu_zp_ini < zero )then

          p_out1 = - ( this%tau_max - this%z_tau ) / this%mu_zp_ini
          this%InterSection_Cases = - 1
          this%Optical_Depth_scatter = p_out1 
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA ) 
          !Get_scatter_distance_IQ = - dlog( one - r1 * this%NormalA ) 
          Get_scatter_distance_IQ = this%z_tau - Temp_log * dabs( this%mu_zp_ini )

      else if( this%mu_zp_ini == zero )then

          !p_out1 = - ( this%tau_max - this%z_tau ) / this%mu_zp_ini
          this%InterSection_Cases = - 2
          !this%Optical_Depth_scatter = p_out1 
          this%NormalA = one! - dexp( - p_out1 )
          r1 = ranmar()
          Get_scatter_distance_IQ = this%z_tau !- dlog( one - r1 * this%NormalA ) 

      endif
     
      if( isnan(Get_scatter_distance_IQ) ) then
           write(*,*)'ends=',Get_scatter_distance_IQ, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%mu_zp_ini
           stop
      endif 

      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_IQ, p_out1, Sigma_I, this%R_out
          write(*,*)'sdf3==',  this%r_times_p, this%r_ini, this%R_out
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf12==',Get_scatter_distance_IQ, this%InterSection_Cases, p_out1, this%z_ini,&
           this%mu_zp_ini
       endif 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_IQ < zero) then
          write(*,*)'sdf==',Get_scatter_distance_IQ, p_out1, Sigma_I, this%Z_max, this%z_ini, &
              this%mu_zp_ini, this%InterSection_Cases
          stop
      endif
      If (Get_scatter_distance_IQ == zero) then
          write(*,*)'ends=',Get_scatter_distance_IQ, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau
          stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function Get_scatter_distance_IQ

 
 
!*******************************************************************************************************
 
      end module ScatDistance_FlatSP





