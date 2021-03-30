      module ScatDistance_FlatSP
      use ScatterPhoton_KN 
      use CrossSection   
      implicit none 

      type, public, extends(ScatPhoton_KN) :: Photon_With_ScatDistance_FlatSP
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
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
          real(mcp) :: Important_Sampling_Const
          logical :: test_it = .FALSE.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !real(mcp), dimension(0:100) :: delta_pds(0:100)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
          real(mcp) :: A_normal
          real(mcp) :: Sigma_Max
          !real(mcp) :: p_out
          !real(mcp) :: p_boundary1
          !real(mcp) :: p_boundary2
          !real(mcp) :: p_boundary3
          !real(mcp) :: p_maxs
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          logical :: fall2BHs, escapteds
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: Optical_Depth_scatter 
          integer :: cases, InterSection_Cases 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          character*80 :: CrossSectFileName 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      contains 
!*******************************************************************************************************   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Set_Cross_Section_3Te   =>   Set_Cross_Section_3Te_sub
          procedure, public :: sigma_fn
          procedure, public :: sigma_KNs
          procedure, public :: Get_scatter_distance_BoundReflec
          procedure, public :: Get_scatter_distance_BoundReflec1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance_FlatSP
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Set_Cross_Section_3Te_sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains   
!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_BoundReflec( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp), intent(in) :: T_e1  
      real(mcp) :: r1, r2 
      real(mcp) :: p_out1, eta, delta_z
      integer(kind=8) :: i, path_cases
     
      if( this%Vector_of_Momentum_ini(3) > zero )then
          this%InterSection_Cases = 1  
          p_out1 = this%z_tau / this%Vector_of_Momentum_ini(3) * this%ne_times_Sigma_a
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          delta_z = dlog( one - r1 * this%NormalA ) * this%Vector_of_Momentum_ini(3) / &
                       this%ne_times_Sigma_a
          Get_scatter_distance_BoundReflec = this%z_tau + delta_z
          !if( this%E_ini <= 1.D-2 )then
          !    this%CROS_absorption = dabs( delta_z ) * this%n_e1 * PhotoElect_CSect( this%E_ini )
          !else
          !    this%CROS_absorption = zero
          !endif
      else 
          this%InterSection_Cases = 2  
          p_out1 = - ( this%z_max - this%z_tau ) / this%Vector_of_Momentum_ini(3) * this%ne_times_Sigma_a
          this%NormalA = one - dexp( - p_out1 ) 
          r1 = ranmar()
          !eta = - dlog( one - r1 * this%NormalA ) 
          delta_z = dlog( one - r1 * this%NormalA ) * this%Vector_of_Momentum_ini(3) / &
                       this%ne_times_Sigma_a
          Get_scatter_distance_BoundReflec = this%z_tau + delta_z
          !if( this%E_ini <= 1.D-2 )then
          !    this%CROS_absorption = dabs( delta_z ) * this%n_e1 * PhotoElect_CSect( this%E_ini )
          !else
          !    this%CROS_absorption = zero
          !endif
      endif
      !write(*,*)'sdf==', p_out1, this%NormalA, this%z_max, this%z_tau
     
      if( isnan(Get_scatter_distance_BoundReflec) ) then
           write(*,*)'ends=', this%medium_case, Get_scatter_distance_BoundReflec, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
           stop
      endif 
  
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_BoundReflec, p_out1
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf123==',Get_scatter_distance_BoundReflec, this%InterSection_Cases, p_out1, this%z_ini,&
           this%Vector_of_Momentum_ini(3), this%NormalA, this%medium_case
          stop
       endif 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_BoundReflec < zero) then
          !write(*,*)'sdf==',Get_scatter_distance_BoundReflec, p_out1 , this%Z_max, this%z_ini, &
          !    this%Vector_of_Momentum_ini(3), this%InterSection_Cases
          !stop
      endif
      If (Get_scatter_distance_BoundReflec == zero) then
          write(*,*)'endsf=', this%medium_case, Get_scatter_distance_BoundReflec, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
          stop
      endif
      return 
      end function Get_scatter_distance_BoundReflec
!******************************************************************************************************* 



!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_BoundReflec1( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp), intent(in) :: T_e1  
      real(mcp) :: r1, r2, tau1, tau2
      real(mcp) :: p_out1, Sigma_E1, Sigma_E2, eta
      integer(kind=8) :: i, path_cases
    
      Sigma_E1 = this%sigma_fn( this%E_ini ) * this%n_e1
      Sigma_E2 = this%sigma_KNs( this%E_ini ) * this%n_e2
      if( this%Vector_of_Momentum_ini(3) > zero )then
          this%InterSection_Cases = 1 
          if(this%medium_case == 1)then
              p_out1 = this%z_tau / this%Vector_of_Momentum_ini(3) * Sigma_E1
              this%Optical_Depth_scatter = p_out1
              this%NormalA = one - dexp( - p_out1 )
              r1 = ranmar()
              Get_scatter_distance_BoundReflec1 = this%z_tau + dlog( one - r1 * this%NormalA ) * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E1
          else if(this%medium_case == 2)then
              tau1 = this%z_max / this%Vector_of_Momentum_ini(3) * Sigma_E1
              tau2 = ( this%z_tau - this%z_max ) / this%Vector_of_Momentum_ini(3) * Sigma_E2
              p_out1 = tau1 + tau2
              this%Optical_Depth_scatter = p_out1
              this%NormalA = one - dexp( - p_out1 )
              r1 = ranmar()
              eta = - dlog( one - r1 * this%NormalA )
              if(eta <= tau2)then
                  Get_scatter_distance_BoundReflec1 = this%z_tau - eta * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E2
              else
                  Get_scatter_distance_BoundReflec1 = this%z_max - (eta - tau2) * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E1
                  this%medium_case = 1
              endif
          endif
          !write(*,*)'ends1=', this%medium_case, Sigma_E2, this%z_tau, this%z_max, eta, &
          !                  tau1, tau2, Get_scatter_distance_BoundReflec1
      else 
          this%InterSection_Cases = 2 
          !write(*,*)'ends2=', this%medium_case, this%z_tau
          if(this%medium_case == 1)then
              tau1 = - ( this%z_max - this%z_tau ) / this%Vector_of_Momentum_ini(3) * Sigma_E1 
              p_out1 = tau1
              this%Optical_Depth_scatter = p_out1
              this%NormalA = one ! - dexp( - p_out1 )
              r1 = ranmar()
              eta = - dlog( one - r1 * this%NormalA )
              if(eta <= tau1)then
                  Get_scatter_distance_BoundReflec1 = this%z_tau - eta * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E1
              else
                  Get_scatter_distance_BoundReflec1 = this%z_max - (eta - tau1) * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E2
                  this%medium_case = 2
              endif 
          else if(this%medium_case == 2)then 
              !tau1 = this%z_max / this%Vector_of_Momentum_ini(3) * Sigma_E1
              !tau2 = ( this%z_tau - this%z_max ) / this%Vector_of_Momentum_ini(3) * Sigma_E2
              !p_out1 = tau1 + tau2
              !this%Optical_Depth_scatter = p_out1
              this%NormalA = one! - dexp( - p_out1 )
              r1 = ranmar()
              eta = - dlog( one - r1 * this%NormalA )
              Get_scatter_distance_BoundReflec1 = this%z_tau - eta * &
                                     this%Vector_of_Momentum_ini(3) / Sigma_E2 
          endif 
      endif
     
      if( isnan(Get_scatter_distance_BoundReflec1) ) then
           write(*,*)'ends=', this%medium_case, Get_scatter_distance_BoundReflec1, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
           stop
      endif 
  
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_BoundReflec1, p_out1
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf123==',Get_scatter_distance_BoundReflec1, this%InterSection_Cases, p_out1, this%z_ini,&
           this%Vector_of_Momentum_ini(3), this%NormalA, tau1, tau2, this%medium_case
          stop
       endif 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_BoundReflec1 < zero) then
          !write(*,*)'sdf==',Get_scatter_distance_BoundReflec1, p_out1 , this%Z_max, this%z_ini, &
          !    this%Vector_of_Momentum_ini(3), this%InterSection_Cases
          !stop
      endif
      If (Get_scatter_distance_BoundReflec1 == zero) then
          write(*,*)'endsf=', this%medium_case, Get_scatter_distance_BoundReflec1, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
          stop
      endif
      return 
      end function Get_scatter_distance_BoundReflec1
!******************************************************************************************************* 


      
!*******************************************************************************************************
      subroutine Set_Cross_Section_3Te_sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp) :: T_e
      real(mcp) :: temp, Ephoton
      real(mcp) :: dindexE
      integer(kind = 8) :: i, j, k, N, istat

      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma)
      !write(*,*)'crossection=', this%logE_up, this%logE_low
      open(unit=18, file = this%CrossSectFileName, status = "old", &
                           action = "read", iostat = istat)
      !write(*, *)'fsd', istat
      if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... ' 
              write(unit = *, fmt = *)'************************************************************' 
          do i = 0, N_sigma
              read(unit = 18, fmt = *)this%sigmaaTeE_FST(i)
          enddo
      else
          open(unit=19, file = this%CrossSectFileName, status = "replace", &
                                   action = "write", iostat = istat)
  
          call gauleg_x_w( -one, one, x1000, w1000, 1000 )
          call gaulag( x0la362, w0la362, 362, zero ) 
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_FST(i) = gama_Integration( this%T_e, Ephoton, &
                                   x1000, w1000, 1000, x0la362, w0la362, 362 ) 
                  if(i == N_sigma)this%sigmaaTeE_400(i) = sigma_a( this%T_e, Ephoton )
                  write(unit = 19, fmt = *)this%sigmaaTeE_FST(i)
                  write(unit = *, fmt = *)i, this%sigmaaTeE_FST(i), this%sigmaaTeE_400(i)
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File are Open Failed. The code have to Stop!!'
              Stop
          endif
          close(unit=19)
      endif
      close(unit=18)      
      end subroutine Set_Cross_Section_3Te_sub

!*******************************************************************************************************
      function sigma_fn(this, Ephoton)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      real(mcp) :: sigma_fn, Ei, Ei1
      real(mcp), dimension(0:N_sigma) :: sigmaaTeE
      real(mcp), intent(in) :: Ephoton
      integer(kind = 8) :: i
      real(mcp) :: Delta_p, Theta_p, R_p
      real(mcp) :: Ep, sign_pth
      !integer, intent(in) :: n
!*******************************************************************************************************
    
      Ep = Ephoton
      i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      !write(*, *)'f2 = ', i, this%logE_low, this%logE_up, this%dindexE 
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
          write(*,fmt="(' ', A13, 5ES20.10)")'sigma_fn = =', dlog10( Ephoton ), Ephoton, &
                this%logE_low, this%dindexE, ( dlog10( Ephoton ) - this%logE_low ) / this%dindexE!,&
          write(unit = *, fmt = *)'************************************************************'  
          write(unit = *, fmt = *)'************************************************************' 
      endif 

      Ei = 10**( this%logE_low + this%dindexE * i )
      Ei1 = 10**( this%logE_low + this%dindexE * (i + 1) )
      sigma_fn = ( this%sigmaaTeE_FST(i + 1) - this%sigmaaTeE_FST(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_FST(i)   
      end function sigma_fn
  
!*******************************************************************************************************
 

!*******************************************************************************************************
      function sigma_KNs( this, Ephoton )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp) :: sigma_fn, Ei, Ei1
      !real(mcp), dimension(0:N_sigma) :: sigmaaTeE
      real(mcp), intent(in) :: Ephoton  
      real(mcp) :: sigma_KNs, epsi, sigma_KN1
      !integer, intent(in) :: n
!*******************************************************************************************************
     
      epsi = two * Ephoton / mec2 
      if ( epsi > 8.19D-4 ) then
          sigma_KN1 = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
      else 
          sigma_KN1 = one - epsi
      end if 
      !write(*,*)'sigma_fn = =', sigma_KN1 * Sigma_T, sigma_KN1, Sigma_T
      sigma_KNs = sigma_KN1 * Sigma_T
      return
      end function sigma_KNs
!*******************************************************************************************************


      end module ScatDistance_FlatSP





