      module Photons 
      use ScatDistance
      !use CrossSection 
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photon_With_ScatDistance) :: Photon 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(0: Num_y)  
          real(mcp) :: y1, y2, dy

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
                            tau, T_e, T_s, n_e, y1, y2, CrossSec_filename ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: tau, T_e, T_s, n_e, y1, y2
      character*80, intent(in) :: CrossSec_filename
      integer cases
      real(mcp) :: E_up, E_low
      integer :: ierr
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      this%T_s = T_s
      CALL this%Set_Emin_Emax()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !E_low = Emitter%E_low1 !1.D-5
      E_up = 3.D1  !this%E_up1 !1.D1 
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

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
      If( this%my_id == this%num_process-1 )then  
          call this%Set_Cross_Section_Array_Whth_Te( )
      endif 
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
      CALL MPI_BCAST( this%sigmaaTeE_FST, N_sigma + 1, MPI_DOUBLE_PRECISION, &
                        this%num_process-1, MPI_COMM_WORLD, ierr ) 

      this%effect_number = 0
 
      this%y1 = y1
      this%y2 = y2
      this%dy = ( this%y2 - this%y1 ) / Num_y

      RETURN
      END SUBROUTINE Set_Initial_Parameters_And_Conditions_Sub


!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p_Sub( this ) 
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this  
  
      this%p_scattering = this%Get_scatter_distance2( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      Call this%Calc_Phot_Informations_At_Observor_2zones( ) 
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
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub( this )
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Phot4k_CtrCF_ini = this%Scattered_Phot4k_CF
      this%Phot4k_CovCF_ini = this%Scattered_Phot4k_CovCF
      this%Vector_of_Momentum_ini(1:3) = this%Phot4k_CtrCF_ini(2:4) / this%Phot4k_CtrCF_ini(1)
 
      this%E_ini = DABS( this%Phot4k_CtrCF_ini(1) ) 
      !write(*,*)'5555==', this%E_ini, this%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%r_ini = this%r_p
      this%vector_of_position_ini = this%vector_of_position_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site_Sub( this )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      class(Photon) :: this 

      this%p_scattering = this%Get_scatter_distance2( ) 
      Call this%Calc_Phot_Informations_At_Observor_2zones( ) 
      this%w_ini = this%w_ini * this%NormalA 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
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
      SUBROUTINE Photon_Electron_Scattering_Sub( this, T_e )
!************************************************************************************
      IMPLICIT NONE
      class(Photon) :: this
      REAL(mcp), INTENT(IN) :: T_e 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      this%Phot4k_CtrCF = this%Phot4k_CtrCF_At_p 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL this%Set_Photon_Tetrad_In_CF( this%Phot4k_CtrCF_At_p(2:4) / &
                                          this%Phot4k_CtrCF_At_p(1) ) 
      CALL this%Get_gama_mu_phi_Of_Scat_Elec(T_e) 
      CALL this%Set_Elec_Tetrad_In_CF()
      CALL this%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !CALL this%Compton_Scattering_WithOut_Polarizations() 
      CALL this%Compton_Scattering_With_Zero_QU() 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      !integer, intent(in) :: cases   
      real(mcp) :: index_i1, a1, b1, del_x1, angle
      integer :: i, j, h, k, mu_i, N_low1, i_E, i_1
  
          this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev 
          i_1 = floor( ( dlog10(this%frequency_v) - this%y1 ) / (this%dy/two) )
          if(mod(i_1, 2)==0)then
              i_E = i_1 / 2
          else
              i_E = (i_1 + 1) / 2
          endif 
          if( i_E > Num_y .or. i_E < 0 )return  
          this%v_L_v_i( i_E ) = this%v_L_v_i( i_E ) + this%w_ini * &
                            dexp( - this%Optical_Depth_scatter ) 
 
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





