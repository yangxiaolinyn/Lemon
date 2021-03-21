      module ScatDistance_FlatSP
      use ScatterPhoton
      !use CrossSection 
      implicit none 

      type, public, extends(ScatPhoton) :: Photon_With_ScatDistance_FlatSP
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
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
          real(mcp) :: Important_Sampling_Const
          logical :: test_it = .FALSE.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(0:100) :: delta_pds(0:100)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
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
          real(mcp) :: Optical_Depth_scatter 
          integer :: cases, InterSection_Cases
          !real(mcp) :: Sigma_I0
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          !real(mcp) :: Z_max
          !real(mcp) :: ratio
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      contains 
!*******************************************************************************************************   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Get_scatter_distance_IQ    
          procedure, public :: Get_scatter_distance_tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance_FlatSP
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 
           
 

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_tau( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, dp
      real(mcp) :: r1, r2, r3
      real(mcp) :: sign_pr, Temp_log
      real(mcp) :: p_out1, Sigma_I, tau1, eta
      integer(kind=8) :: i, path_cases
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      if( this%Vector_of_Momentum_ini(3) > zero )then

          p_out1 = this%z_tau / this%Vector_of_Momentum_ini(3)
          this%InterSection_Cases = 1
          this%Optical_Depth_scatter = p_out1 
          this%NormalA = one - dexp( - p_out1  )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA )
          Get_scatter_distance_tau = this%z_tau + Temp_log * this%Vector_of_Momentum_ini(3)   

      else if( this%Vector_of_Momentum_ini(3) < zero )then

          p_out1 = - ( this%tau_max - this%z_tau ) / this%Vector_of_Momentum_ini(3)
          this%InterSection_Cases = - 1
          this%Optical_Depth_scatter = p_out1 
          this%NormalA = one  - dexp( - p_out1 )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA )
          Get_scatter_distance_tau = this%z_tau + Temp_log * this%Vector_of_Momentum_ini(3)  
          !write(*, *)'ff22=',  Get_scatter_distance_tau, r1, this%NormalA

      else if( this%Vector_of_Momentum_ini(3) == zero )then
 
          this%InterSection_Cases = - 3
          this%Optical_Depth_scatter = Infinity 
          this%NormalA = one 
          r1 = ranmar()
          Get_scatter_distance_tau = this%z_tau

      endif 
  
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_tau, p_out1, Sigma_I 
          write(*,*)'sdf3==',  this%r_times_p, this%r_ini 
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf12==',Get_scatter_distance_tau, this%InterSection_Cases, p_out1, this%z_ini,&
           this%Vector_of_Momentum_ini(3), Sigma_I
       endif
      !if( Get_scatter_distance_tau > p_out1 )then
      !    write(*,*)'I111',Get_scatter_distance_tau, p_out1, Get_scatter_distance_tau - p_out1
      !endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_tau < zero) then
          write(*,*)'sdfss==',Get_scatter_distance_tau, p_out1, Sigma_I , this%z_ini, &
              this%Vector_of_Momentum_ini(3), this%InterSection_Cases
          stop
      endif 
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end function Get_scatter_distance_tau
    

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_IQ( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA, Temp_log
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_out1, Sigma_I, tau1, eta
      integer(kind=8) :: i, path_cases
    
      if( this%Vector_of_Momentum_ini(3) > zero )then

          p_out1 =  this%z_tau / this%Vector_of_Momentum_ini(3)
          this%InterSection_Cases = 1
          this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(1)!p_out1 
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA ) 
          Get_scatter_distance_IQ = this%z_tau + Temp_log * this%Vector_of_Momentum_ini(3)  
          !write(*, *)'ff11=', Get_scatter_distance4, p_out1 , this%z_ini, &
          !             this%Vector_of_Momentum_ini(3), dexp( - p_out1 * Sigma_I )

      else if( this%Vector_of_Momentum_ini(3) < zero )then

          p_out1 = - ( this%tau_max - this%z_tau ) / this%Vector_of_Momentum_ini(3)
          this%InterSection_Cases = - 1
          !this%Optical_Depth_scatter = p_out1 
          this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(1)!p_out1 
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          Temp_log = dlog( one - r1 * this%NormalA ) 
          !Get_scatter_distance_IQ = - dlog( one - r1 * this%NormalA ) 
          Get_scatter_distance_IQ = this%z_tau - Temp_log * dabs( this%Vector_of_Momentum_ini(3) )

      else if( this%Vector_of_Momentum_ini(3) == zero )then

          !p_out1 = - ( this%tau_max - this%z_tau ) / this%Vector_of_Momentum_ini(3)
          this%InterSection_Cases = - 2
          !this%Optical_Depth_scatter = p_out1 
          this%NormalA = one! - dexp( - p_out1 )
          r1 = ranmar()
          Get_scatter_distance_IQ = this%z_tau !- dlog( one - r1 * this%NormalA ) 

      endif
     
      if( isnan(Get_scatter_distance_IQ) ) then
           write(*,*)'ends=',Get_scatter_distance_IQ, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini(3)
           stop
      endif 

      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_IQ, p_out1, Sigma_I 
          write(*,*)'sdf3==',  this%r_times_p, this%r_ini 
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf12==',Get_scatter_distance_IQ, this%InterSection_Cases, p_out1, this%z_ini,&
           this%Vector_of_Momentum_ini(3) 
       endif 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_IQ < zero) then
          write(*,*)'sdf==',Get_scatter_distance_IQ, p_out1, Sigma_I , this%z_ini, &
              this%Vector_of_Momentum_ini(3), this%InterSection_Cases
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





