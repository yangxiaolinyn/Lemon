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
          real(mcp) :: frequency_v 
          real(mcp) :: Optical_Depth_scatter 

      contains 
!******************************************************************************************************* 
          procedure, public :: sigma_fn    
          procedure, public :: Get_scatter_distance2  
          procedure, public :: Set_Cross_Section_Array_Whth_Te => &
                               Set_Cross_Section_Array_Whth_Te_sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Set_Cross_Section_Array_Whth_Te_sub  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains 

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance2( this )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_out1, Sigma_I
      integer(kind=8) :: i, path_cases
   
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) ) 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Sigma_I = this%n_e * this%sigma_fn( this%E_ini )
      this%Optical_Depth_scatter = p_out1 * Sigma_I
      this%NormalA = one - dexp( - p_out1 * Sigma_I )
      r1 = ranmar()
      Get_scatter_distance2 = - dlog( one - r1 * this%NormalA ) / Sigma_I
      !write(*,*)'sdf1==', Get_scatter_distance2, this%NormalA, Sigma_I 
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
          write(*,*)'sdfss==',Get_scatter_distance2, r1, this%NormalA,  Sigma_I, this%r_times_p, &
          this%r_ini, this%R_out, p_out1
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
      function sigma_fn(this, Ephoton)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this 
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
         ! write(*, *)'f1 = ', i, this%logE_low, this%logE_up, this%dindexE 
          Ep = 10.D0**this%logE_low
          i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      endif
      if ( dlog10( Ep ) > this%logE_up ) then
          !write(*, *)'f2 = ', i, this%logE_low, this%logE_up, dlog10( Ep ), this%dindexE 
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
      end function sigma_fn


!*******************************************************************************************************
 
      end module ScatDistance




