      module ScatDistance
      use ScatterPhoton
      use CrossSection 
      implicit none 
      integer, parameter :: N_sigma = 1000 

      type, public, extends(ScatPhoton) :: Photon_With_ScatDistance
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_85
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_230
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_50
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_100
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE
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
          real(mcp) :: frequency_v 
          real(mcp) :: Optical_Depth_scatter
          real(mcp) :: Optical_Depth_absorption
          real(mcp) :: j_Enu_theta

      contains 
!******************************************************************************************************* 
          procedure, public :: sigma_fn
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Set_Cross_Section   => Set_Cross_Section_sub  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          procedure, public :: get_Tau2                           =>  get_Tau2_fn
          procedure, public :: Get_scatter_distance
          procedure, public :: Get_scatter_distance2 
          procedure, public :: Set_Cross_Section_3Te => Set_Cross_Section_3Te_sub
          procedure, public :: n_e_p => n_e_p_fn
          procedure, public :: get_Tau2_At_out_zone => get_Tau2_At_out_zone_fn
          procedure, public :: Integration_Of_OpticalDepth
          procedure, public :: Get_Max_Value_of_Sigma_2zones3 => &
                                Get_Max_Value_of_Sigma_2zones3_sub
          !procedure, public :: RandomNum2p_case2   =>   RandomNum2p_case2_fn
          !procedure, public :: Integrals   =>    Integrals_fn
          !procedure, public :: RandomNum2p_case3   =>   RandomNum2p_case3_fn
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !procedure, public :: RandNum2p_case1_Finit_Space   =>  &
          !                      RandNum2p_case1_Finit_Space_fn
          !procedure, public :: RandNum2p_case2_Finit_Space   =>  &
          !                      RandNum2p_case2_Finit_Space_fn
          !procedure, public :: RandNum2p_case3_Finit_Space   =>  &
          !                      RandNum2p_case3_Finit_Space_fn
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_With_ScatDistance
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private ::  Set_Cross_Section_sub  
      private :: Set_Cross_Section_3Te_sub 
      private :: n_e_p_fn
      private :: get_Tau2_At_out_zone_fn
      !private :: RandomNum2p_case3_fn
      !private :: RandomNum2p_case2_fn
      !private :: Integrals_fn
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !private :: RandNum2p_case1_Finit_Space_fn
      !private :: RandNum2p_case2_Finit_Space_fn
      !private :: RandNum2p_case3_Finit_Space_fn
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains  

!*******************************************************************************************************
      real(mcp) function Get_scatter_distance( this, T_e1 )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: T_e1 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out, eta
      real(mcp) :: p_c, p_out1, p_1, p_2, delta
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i, path_cases
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  First ditermine which zone the Photon_With_ScatDistance is included in.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) ) 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !Emitter%j_theta_E_nu_emissity( T_e1, this%E_ini, zero, one )
      this%Optical_Depth_absorption = this%j_Enu_theta / &
          ( two * planck_h * (this%E_ini*1.D6*erg_of_one_ev / planck_h )**3 / Cv**2 / &
                       ( dexp(this%E_ini/51.1D0) - one ) ) * p_out1
      !write(*,*)'dist2===', this%j_Enu_theta, &
      !     ( two*(this%E_ini*1.D6*erg_of_one_ev)**3/Cv**2 / ( dexp(this%E_ini/51.1D0) - one ) ), p_out1
      this%Optical_Depth_scatter = this%n_e * this%sigma_fn(T_e1, this%E_ini, one) * p_out1
1123  continue
      r1 = ranmar() 
      Get_scatter_distance = - dlog( one - r1 ) / &
                               ( this%n_e * this%sigma_fn(T_e1, this%E_ini, one) )
      !write(*,*)'dist2===', Get_scatter_distance, p_out1,  - dlog( one - r1 ), r1
      if( Get_scatter_distance > p_out1 )then
          this%At_outer_Shell = .true.
      else
          goto 1123
      endif
      If (Get_scatter_distance < zero) then
          write(*,*)'sdf==',Get_scatter_distance,this%n_e , this%sigma_fn(T_e1, this%E_ini, one), &
              - dlog( one - r1 )
          stop
      endif
      If (Get_scatter_distance == zero) then
          write(*,*)'ends=',Get_scatter_distance, this%n_e * this%sigma_fn(T_e1, this%E_ini, one)
          write(*,*)'ends2=', this%r_ini, this%R_out
          stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      end function Get_scatter_distance
   
 
!*******************************************************************************************************
      real(mcp) function Get_scatter_distance2( this, T_e1 )
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp), intent(in) :: T_e1
      real(mcp) :: T_e 
      real(mcp) :: p, p1, p_max, func1_tau_max_value, rp, rtp 
      real(mcp) :: temp,temp2, Sigma_atP, dp
      real(mcp) :: r1,r2, rprobability, r3, tempA
      real(mcp) :: sign_pr, p_out 
      real(mcp) :: p_c, p_out1, p_1, p_2, delta
      logical :: p_lt_pmax, loop_stop
      integer(kind=8) :: i, path_cases
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  First ditermine which zone the Photon is included in.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - this%r_times_p + dsqrt( this%r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )  
          T_e = 100.D0 * mec2
              this%p_maxs = p_out1
              path_cases = 1
              this%p_boundary1 = p_out1 
              !call this%Get_Max_Value_of_Sigma_2zones3( T_e, p_out1, path_cases )
              this%Optical_Depth_scatter = this%n_e * this%sigma_fn(T_e1, this%E_ini, p) * p_out1
              this%NormalA = one - dexp( - this%Optical_Depth_scatter )
              !write(*,*)'cass1==',this%n_e, this%sigma_fn(T_e1, this%E_ini, p), p_out1, &
              !           this%Optical_Depth_scatter
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              r1 = ranmar() 
              Get_scatter_distance2 = - dlog( one - r1 * this%NormalA ) / &
                              this%n_e * this%sigma_fn(T_e1, this%E_ini, p)
              if( Get_scatter_distance2 > p_out1 )then
                  write(*,*)'I111',Get_scatter_distance2, p_out1, Get_scatter_distance2 - p_out1
              endif
              !this%w_ini = this%w_ini * tempA
              !write(*,*)'I am in zone 12',   
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      If (Get_scatter_distance2 < zero) then
          write(*,*)'sdf==',Get_scatter_distance2,p,p_1 
          stop
      endif
      If (Get_scatter_distance2 == zero) then
          write(*,*)'ends=',Get_scatter_distance2, r1, this%NormalA,  - dlog( one - r1 * this%NormalA ) , &
                              this%n_e * this%sigma_fn(T_e1, this%E_ini, p)
          write(*,*)'tts=', p_out1, this%r_ini
          !stop
      endif
      return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !Get_scatter_distance2 = p
       !write(*,*)this%Sigma_Max,Get_scatter_distance2
      end function Get_scatter_distance2
  

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
      subroutine Set_Cross_Section_sub(this, T_e)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
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
      subroutine Set_Cross_Section_3Te_sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_With_ScatDistance) :: this
      real(mcp) :: T_e
      real(mcp) :: temp, Ephoton
      real(mcp) :: dindexE
      integer(kind = 8) :: i, j, k, N, istat

      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma)
      write(*,*)'sssssss111=', this%logE_up, this%logE_low, this%dindexE
      open(unit=18, file = './data/SigmaArray1000.dat', status = "old", &
                           action = "read", iostat = istat)
      if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... ' 
              write(unit = *, fmt = *)'************************************************************' 
          do i = 0, N_sigma
              read(unit = 18, fmt = *)this%sigmaaTeE_100(i)
          enddo
      else
          open(unit=19, file = './data/SigmaArray1000.dat', status = "replace", &
                                   action = "write", iostat = istat)
          T_e = 100.D0 * mec2
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_100(i) = sigma_a(T_e, Ephoton)
                  write(unit = 19, fmt = *)this%sigmaaTeE_100(i)
                  write(unit = *, fmt = *)i,this%sigmaaTeE_100(i)
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File are Open Failed. The code have to Stop!!'
              Stop
          endif
          close(unit=19)
      endif
      close(unit=18)     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      open(unit=20, file = './data/SigmaArray50.dat', status = "old", &
                           action = "read", iostat = istat)
      if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... ' 
              write(unit = *, fmt = *)'************************************************************' 
          do i = 0, N_sigma
              read(unit = 20, fmt = *)this%sigmaaTeE_50(i)
          enddo
      else
          open(unit=21, file = './data/SigmaArray50.dat', status = "replace", &
                                   action = "write", iostat = istat)
          T_e = 50.D-3
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_50(i) = sigma_a(T_e, Ephoton)
                  write(unit = 21, fmt = *)this%sigmaaTeE_50(i)
                  write(unit = *, fmt = *)i
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File are Open Failed. The code have to Stop!!'
              Stop
          endif
          close(unit=21)
      endif
      close(unit=20)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      open(unit=22, file = './data/SigmaArray230.dat', status = "old", &
                           action = "read", iostat = istat)
      if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)' Now reading the SigmaArray data from file..... ' 
              write(unit = *, fmt = *)'************************************************************' 
          do i = 0, N_sigma
              read(unit = 22, fmt = *)this%sigmaaTeE_230(i)
          enddo
      else
          open(unit=23, file = './data/SigmaArray230.dat', status = "replace", &
                                   action = "write", iostat = istat)
          T_e = 230.D-3
          if (istat == 0) then
              write(unit = *, fmt = *)'************************************************************'  
              write(unit = *, fmt = *)'Now Writting the SigmaArray data to the file.....' 
              write(unit = *, fmt = *)'************************************************************' 
              do i = 0, N_sigma
                  write(unit = *, fmt = *)'herer'
                  Ephoton = 10**( this%logE_low + this%dindexE * i )
                  this%sigmaaTeE_230(i) = sigma_a(T_e, Ephoton)
                  write(unit = 23, fmt = *)this%sigmaaTeE_230(i)
                  write(unit = *, fmt = *)i
              enddo
          else
              write(unit = *, fmt = *)'The SigmaArray File are Open Failed. The code have to Stop!!'
              Stop
          endif
          close(unit=23)
      endif
      close(unit=22) 
      end subroutine Set_Cross_Section_3Te_sub

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
      !write(*,*)'2222 === ',Ephoton, 10.D0**this%logE_up
      Ep = Ephoton
      if ( dlog10( Ep ) < this%logE_low ) then
          Ep = 10.D0**this%logE_low
      endif
      i = floor( ( dlog10( Ep ) - this%logE_low ) / this%dindexE )
      if( i >= 1000 )i=999
      !write(unit = *, fmt = *)'sigma=',Ephoton, dlog10( Ephoton ), this%logE_low,this%dindexE, i
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
          !write(*,*)'sigma2==', i, T_e
      If ( T_e == 100.D0 .or. T_e == 100.D0 * mec2 ) then   
          sigma_fn = ( this%sigmaaTeE_100(i + 1) - this%sigmaaTeE_100(i) ) / &
                 ( Ei1 - Ei ) * (Ep - Ei) + this%sigmaaTeE_100(i) 
         ! write(*,*)'sigma2==', i, T_e, sigma_fn 
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
      real(mcp) :: r_ini, p_length, p_bias
   
          !this%Sigma_Max = zero  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !this%Sigma_Max = this%n_e_p( zero ) * this%sigma_fn(T_e1, this%E_ini, zero)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !n = 2000    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
          this%Optical_Depth_scatter = this%n_e0 * this%sigma_fn(T_e1, this%E_ini, p) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      end subroutine Get_Max_Value_of_Sigma_2zones3_Sub

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module ScatDistance




 
