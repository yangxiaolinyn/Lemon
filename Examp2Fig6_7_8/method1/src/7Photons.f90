      module Photons 
      use ScatDistance
      !use CrossSection 
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photon_With_ScatDistance) :: Photon 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(1:500)
          real(mcp) :: v_L_v_i_ET(1:500, 1:1000)
          !real(mcp) :: nu_low
          !real(mcp) :: nu_up 

      contains 
!*******************************************************************************************************  
          procedure, public :: r_p2 
          procedure, public :: Calc_Phot_Informations_At_Observor_2zones => &
                               Calc_Phot_Informations_At_Observor_2zones_Sub
          procedure, public :: Determine_P_Of_Scatt_Site_And_Quantities_At_p => &
                               Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub 
          procedure, public :: Set_InI_Conditions_For_Next_Scattering  =>  &
                               Set_InI_Conditions_For_Next_Scattering_Sub
          procedure, public :: Determine_Next_Scattering_Site  =>  &
                               Determine_Next_Scattering_Site_Sub
          procedure, public :: Photon_Electron_Scattering   =>   &
                               Photon_Electron_Scattering_Sub
          procedure, public :: Set_Initial_Parameters_And_Conditions   =>   &
                               Set_Initial_Parameters_And_Conditions_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Calc_Phot_Informations_At_Observor_2zones_Sub
      private :: Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub 
      private :: Set_InI_Conditions_For_Next_Scattering_Sub
      private :: Determine_Next_Scattering_Site_Sub
      private :: Photon_Electron_Scattering_Sub
      private :: Set_Initial_Parameters_And_Conditions_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains
!************************************************************************************ 
      SUBROUTINE Set_Initial_Parameters_And_Conditions_Sub( this, &
                            tau, T_e, T_s, n_e, CrossSec_filename ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: tau, T_e, T_s, n_e
      character*80, intent(in) :: CrossSec_filename
      integer cases
      real(mcp) :: E_up, E_low
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL this%Set_Emin_Emax() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !E_low = Emitter%E_low1 !1.D-5
      E_up = 2.D1  !this%E_up1 !1.D1 
      this%T_e = T_e 

      this%logE_low = DLOG10( this%E_low1 )
      this%logE_up = DLOG10( E_up )
      this%dindexE = ( this%logE_up - this%logE_low )/dfloat(N_sigma) !Here indexE means a
          ! new variable y and E = 10^y, and logE_up =y2, logE_low = y1, dindexE = dy, 
          ! E is the energy of the photon.

      this%n_e = n_e
      this%R_out = tau / Sigma_T / this%n_e
      this%CrossSectFileName = CrossSec_filename
      !CALL this%Set_Cross_Section_Te( CrossSec_filename )
      call this%Set_Cross_Section_Array_Whth_Te( )
      this%effect_number = 0

      RETURN
      END SUBROUTINE Set_Initial_Parameters_And_Conditions_Sub


!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub( this ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this  
      integer cases 
   
      this%p_scattering = this%Get_scatter_distance2( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      cases = 1
      Call this%Calc_Phot_Informations_At_Observor_2zones( cases ) 
      this%w_ini = this%w_ini * this%NormalA  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
      this%r_p = this%r_p2( this%Vector_of_position_ini, &
                      this%Vector_of_Momentum_ini, this%p_scattering )
      this%Vector_of_position_p = this%Vector_of_position
  
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini
      this%Phot4k_CovCF_At_p = this%Phot4k_CovCF_ini 
      !write(unit = *, fmt = *)'************************************************************'

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub
 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this, T_e, phot, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Phot4k_CtrCF_ini = sphot%Scattered_Phot4k_CF
      this%Phot4k_CovCF_ini = sphot%Scattered_Phot4k_CovCF
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / this%Phot4k_CtrCF_ini(1)
   
      this%E_ini = DABS( this%Phot4k_CovCF_ini(1) ) 
      !write(*,*)'5555==', this%E_ini, this%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r_ini = this%r_p
      this%vector_of_position_ini = this%vector_of_position_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site_Sub( this, T_e, phot, sphot, Scatter_Times )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: T_e
      integer(kind = 8), intent(IN) :: Scatter_Times
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !write(*,*)'6666==', this%E_ini, this%Phot4k_CovCF_ini(1)  
      this%p_scattering = this%Get_scatter_distance2( )  
      !write(*,*)'7777==', this%E_ini, this%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !this%direct_escaped = .True.
      cases = 1
      Call this%Calc_Phot_Informations_At_Observor_2zones( cases )
      !write(*,*)'ss3=',this%w_ini
      this%w_ini = this%w_ini * this%NormalA
      !write(*,*)'ss4=',this%w_ini
      !this%direct_escaped = .false. 
      !write(*,*)'8888==', this%E_ini, this%Phot4k_CovCF_ini(1)
          !write(*,*)'666=', this%time_travel 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !If ( this%At_outer_Shell ) RETURN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%time_travel = this%time_travel + this%p_scattering / Cv
          !write(*,*)'101=', this%time_travel,  this%p_scattering/ CV
      this%r_p = this%r_p2( this%Vector_of_position_ini, &
                 this%Vector_of_Momentum_ini, this%p_scattering ) 
      if( this%r_p > this%R_out )then
          this%test_it = .true.
          this%p_scattering = this%Get_scatter_distance2( )  
          write(*,*)'101=', this%time_travel,  this%p_scattering/ CV
          write(*,*)'102=', this%r_p, this%p_scattering, this%NormalA
          write(*,*)'103=', dsqrt( this%Vector_of_position_ini(1)**2+&
                   this%Vector_of_position_ini(2)**2+this%Vector_of_position_ini(3)**2 ), this%r_ini, &
               dsqrt( this%Vector_of_Momentum_ini(1)**2+this%Vector_of_Momentum_ini(2)**2+&
               this%Vector_of_Momentum_ini(3)**2 ), this%r_p2( this%Vector_of_position_ini, &
                 this%Vector_of_Momentum_ini, this%p_scattering ) 
          stop
      endif
      ! the vector of position also obtained.
      this%Vector_of_position_p = this%Vector_of_position 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !write(unit = *, fmt = *)'************************************************************'
      this%Phot4k_CtrCF_At_p = this%Phot4k_CtrCF_ini
      this%Phot4k_CovCF_At_p = this%Phot4k_CovCF_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering_Sub( this, T_e, sphot )
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      sphot%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Set_Photon_Tetrad_In_CF() 
      CALL sphot%Get_gama_mu_phi_Of_Scat_Elec(T_e) 
      CALL sphot%Set_Elec_Tetrad_In_CF()
      CALL sphot%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Compton_Scattering_WithOut_Polarizations()  
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub( this, cases )
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





