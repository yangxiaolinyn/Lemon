      module ScatDistance_FlatSP
      use ScatterPhoton_KN 
      use CrossSection 
      !use CrossSection 
      implicit none 

      type, public, extends(ScatPhoton_KN) :: Photon_With_ScatDistance_FlatSP
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_400 
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
          procedure, public :: Get_scatter_distance_tau
          procedure, public :: Set_Cross_Section   =>   Set_Cross_Section_sub
          procedure, public :: sigma_fn
          procedure, public :: Get_scatter_distance_BoundReflec
          procedure, public :: sigma_KN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance_FlatSP
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Set_Cross_Section_sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 
           

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_BoundReflec( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp  
      real(mcp) :: r1,r2, rprobability, r3, tempA 
      real(mcp) :: p_out1, Sigma_I, tau1, eta
      integer(kind=8) :: i, path_cases
    
      this%sigma_KNs = this%sigma_KN( this%E_ini ) * this%n_e
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if( this%Vector_of_Momentum_ini(3) > zero )then

          p_out1 = this%z_tau / this%Vector_of_Momentum_ini(3) * this%sigma_KNs
          this%InterSection_Cases = 1
          this%Optical_Depth_scatter = p_out1
          this%NormalA = one - dexp( - p_out1 )
          r1 = ranmar()
          Get_scatter_distance_BoundReflec = this%z_tau + dlog( one - r1 * this%NormalA ) * &
                                     this%Vector_of_Momentum_ini(3) / this%sigma_KNs

      else if( this%Vector_of_Momentum_ini(3) < zero )then

          !p_out1 = - ( this%tau_max - this%z_tau ) / this%Vector_of_Momentum_ini(3) * this%sigma_KNs
          this%InterSection_Cases = - 1
          !this%Optical_Depth_scatter = p_out1 
          !this%NormalA = one - dexp( - p_out1 )
          this%NormalA = one !- dexp( - p_out1 )
          r1 = ranmar() 
          Get_scatter_distance_BoundReflec = this%z_tau + dlog( one - r1 * this%NormalA ) * &
                                     this%Vector_of_Momentum_ini(3) / this%sigma_KNs
          !write(*, *)'ffss11111=', this%z_tau, Get_scatter_distance_BoundReflec
      else if( this%Vector_of_Momentum_ini(3) == zero )then

          p_out1 = 1.D100
          this%InterSection_Cases = - 2
          this%Optical_Depth_scatter = p_out1 
          this%NormalA = one 
          r1 = ranmar()
          !Get_scatter_distance_BoundReflec = - dlog( one - r1 * this%NormalA ) 
          Get_scatter_distance_BoundReflec = this%z_tau 
      endif
      !write(*,*)'endss=', p_out1
     
      if( isnan(Get_scatter_distance_BoundReflec) ) then
           write(*,*)'endss=',Get_scatter_distance_BoundReflec, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
           stop
      endif 
  
      if( this%test_it )then
          write(*,*)'sdf2==', r1, this%NormalA, Get_scatter_distance_BoundReflec, p_out1 
      endif
      if(p_out1 < zero )then
          write(*,*)'sdf12==',Get_scatter_distance_BoundReflec, this%InterSection_Cases, &
             p_out1, this%z_tau, this%Vector_of_Momentum_ini(3), this%sigma_KNs
       endif 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If (Get_scatter_distance_BoundReflec < zero) then
          !write(*,*)'sdf==',Get_scatter_distance_BoundReflec, p_out1, Sigma_I, this%Z_max, this%z_ini, &
          !    this%Vector_of_Momentum_ini(3), this%InterSection_Cases
          !stop
      endif
      If (Get_scatter_distance_BoundReflec == zero) then
          write(*,*)'endsf=',Get_scatter_distance_BoundReflec, r1, this%NormalA,  &
             this%InterSection_Cases, p_out1, this%z_tau, this%Vector_of_Momentum_ini
          write(*,*)'endsfs=', this%tau_max, this%z_tau, this%Vector_of_Momentum_ini(3), this%sigma_KNs
          write(*,*)'**********************************************************************'
      endif
      return 
      end function Get_scatter_distance_BoundReflec
!*******************************************************************************************************

 

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance_tau( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, Temp_log
      real(mcp) :: p_out1, Sigma_I, tau1, eta
      integer(kind=8) :: i, path_cases
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      if( this%Vector_of_Momentum_ini(3) > zero )then

          p_out1 = this%z_tau / this%Vector_of_Momentum_ini(3) * this%sigma_KNs
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
      
!*******************************************************************************************************
      subroutine Set_Cross_Section_sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp) :: T_e
      real(mcp) :: temp, Ephoton
      real(mcp) :: dindexE
      integer(kind = 8) :: i, j, k, N, istat

      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma)
      write(*,*)'crossection=', this%logE_up, this%logE_low
      open(unit=18, file = this%CrossSectFileName, status = "old", &
                           action = "read", iostat = istat)
      write(*, *)'fsd', istat
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
          call Set_xi_wi_all()
          !T_e = 4.D0 * mec2
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_FST(i) = gama_Integration( this%T_e, Ephoton, &
                                   x1000, w1000, 1000, x0la1000, w0la1000, 362 ) 
                  !this%sigmaaTeE_400(i) = sigma_a( this%T_e, Ephoton )
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
      end subroutine Set_Cross_Section_sub

!*******************************************************************************************************
      function sigma_fn(this, T_e, Ephoton, p)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
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
 

      sigma_fn = ( this%sigmaaTeE_FST(i + 1) - this%sigmaaTeE_FST(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_FST(i)  
 
      !sigma_fn = ( sigmaaTeE(i + 1) - sigmaaTeE(i) ) / &
      !           ( Ei1 - Ei ) * (Ep - Ei) + sigmaaTeE(i)  
      end function sigma_fn
  

!*******************************************************************************************************
      function sigma_KN( this, Ephoton )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance_FlatSP) :: this
      !real(mcp) :: sigma_fn, Ei, Ei1
      !real(mcp), dimension(0:N_sigma) :: sigmaaTeE
      real(mcp), intent(in) :: Ephoton  
      real(mcp) :: sigma_KN, epsi, sig, alpha
      real(mcp) :: Ep, sign_pth
      !integer, intent(in) :: n
!*******************************************************************************************************
     
      epsi = two * Ephoton / mec2
      alpha = Ephoton / mec2 
      if ( epsi > 8.19D-4 ) then
          sigma_KN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
             + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four) * sigma_T
          sig = three * sigma_T / eight / alpha * ( ( one - two *(one + alpha)/alpha**2 )*&
             dlog(two*alpha + one) + 0.5D0 + four / alpha - one /two/(two*alpha + one)**2 ) 
          !write(*, *)'fffs= ', sigma_KN, sig
      else 
          sigma_KN = ( one - epsi ) * sigma_T
      end if 
      !sigma_KN = one 
      return
      end function sigma_KN
!*******************************************************************************************************
 
      end module ScatDistance_FlatSP





