!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !USE constants
      USE RandUtils
      !USE PhotonEmitterBB
      USE Photons_FlatSP
      USE MPI
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Set_initial_parameter_values( Phot, Emitter, tau )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter
      TYPE(Photon_FlatSP), INTENT(INOUT) :: Phot 
      REAL(mcp), INTENT(INOUT) :: tau
      integer :: i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!| Set Initial conditions for the Emitter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      Emitter%tau_max = tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set Initial conditions for the Photon                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%tau_max = tau   
 
      CALL Set_psi_phi_chi_zera_array() 
      phot%effect_number = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      phot%delta_pds = zero
      phot%v_L_v_i = zero
      phot%d_theta = one / Num_PolDeg
      phot%d_phi =  twopi / Num_Phi 
      phot%mu_estimat(1) = one / 100.D0 * 5.D0
      phot%mu_estimat(2) = one / 100.D0 * 30.D0
      phot%mu_estimat(3) = one / 100.D0 * 60.D0
      phot%mu_estimat(4) = one / 100.D0 * 85.D0
      do i = 0, Num_Phi
          phot%phi_estimates( i ) = phot%d_phi * i
      enddo
   
      RETURN
      END SUBROUTINE Set_initial_parameter_values


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Calculate_The_Diffuse_Reflection_of_Chandra( Phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      TYPE(Photon_FlatSP), INTENT(INOUT) :: Phot 
      INTEGER :: i, j
      REAL(mcp) :: mu, mu0, phi, phi0, I_r, I_l, U, V, IQUV_in(1: 3), IQU(1: 3), &
                   Irluv_in(1: 3), Irluv(1: 3)
      real(mcp), dimension(1:3, 1:3) :: QS_Matrix, QS_Matrix_rluv

      !If ( myid == 0 ) then
          IQUV_in(1) = one
          IQUV_in(2) = one / four
          IQUV_in(3) = one / four
          !Irluv_in(1) = (IQUV_in(1) + IQUV_in(2)) / two
          !Irluv_in(2) = (IQUV_in(1) - IQUV_in(2)) / two
          !Irluv_in(3) = IQUV_in(3)
          mu0 = 0.8D0
          phi0 = zero
          100 FORMAT(' ', 4F20.13) 
          200 FORMAT(' ', 4F20.13) 
          open(unit=9, file='./spectrum/phi_dist/Iquv10.txt', status="replace") 
          open(unit=10, file='./spectrum/phi_dist/Iquv30.txt', status="replace") 
          open(unit=11, file='./spectrum/phi_dist/Iquv60.txt', status="replace") 
          open(unit=12, file='./spectrum/phi_dist/Iquv80.txt', status="replace") 

          mu = one / 100.D0 * 5.D0
          do i = 0, 200 
              phi = twopi / 200.D0 * i   
              CALL Phot%Diffuse_Reflection_Ir_Il_U_Setting( mu, mu0, phi, phi0, I_r, I_l, U, V )
              !write(unit = 9, fmt = 100)I_r, I_l, U, V

              CALL Phot%Diffuse_Reflection_QS_Matrix_of_IQUV_Setting( mu, mu0, phi, phi0, QS_Matrix ) 
              CALL Matrix_Multiplication33X31_Sub( QS_Matrix, IQUV_in, IQU )
              write(unit = 9, fmt = 200)IQU(1), IQU(2), IQU(3), V 
          enddo

          mu = one / 100.D0 * 30.D0
          do i = 0, 200
              phi = twopi / 200.D0 * i  
              CALL Phot%Diffuse_Reflection_Ir_Il_U_Setting( mu, mu0, phi, phi0, I_r, I_l, U, V )

              CALL Phot%Diffuse_Reflection_QS_Matrix_of_IQUV_Setting( mu, mu0, phi, phi0, QS_Matrix ) 
              CALL Matrix_Multiplication33X31_Sub( QS_Matrix, IQUV_in, IQU )
              write(unit = 10, fmt = 200)IQU(1), IQU(2), IQU(3), V 
          enddo 

          mu = one / 100.D0 * 60.D0
          do i = 0, 200
              phi = twopi / 200.D0 * i  
              CALL Phot%Diffuse_Reflection_Ir_Il_U_Setting( mu, mu0, phi, phi0, I_r, I_l, U, V )

              CALL Phot%Diffuse_Reflection_QS_Matrix_of_IQUV_Setting( mu, mu0, phi, phi0, QS_Matrix ) 
              CALL Matrix_Multiplication33X31_Sub( QS_Matrix, IQUV_in, IQU )
              write(unit = 11, fmt = 200)IQU(1), IQU(2), IQU(3), V 
          enddo 

          mu = one / 100.D0 * 85.D0
          do i = 0, 200
              phi = twopi / 200.D0 * i  
              CALL Phot%Diffuse_Reflection_Ir_Il_U_Setting( mu, mu0, phi, phi0, I_r, I_l, U, V )

              CALL Phot%Diffuse_Reflection_QS_Matrix_of_IQUV_Setting( mu, mu0, phi, phi0, QS_Matrix ) 
              CALL Matrix_Multiplication33X31_Sub( QS_Matrix, IQUV_in, IQU )
              write(unit = 12, fmt = 200)IQU(1), IQU(2), IQU(3), V 
          enddo 
 
          close(unit=9) 
          close(unit=10)
          close(unit=11) 
          close(unit=12)
          !stop
      !endif 
      RETURN
      END SUBROUTINE Calculate_The_Diffuse_Reflection_of_Chandra

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon( Emitter )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL Emitter%get_Phot4k_CtrCF_CovCF_Reflection()  

      RETURN
      END SUBROUTINE Emitter_A_Photon

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(IN) :: Emitter
      TYPE(Photon_FlatSP), INTENT(INOUT) :: Phot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      phot%r_ini      = Emitter%r
      !phot%theta_ini  = Emitter%theta
      !phot%mucos_ini  = Emitter%mucos
      !phot%musin_ini  = Emitter%musin
      !phot%phi_ini    = Emitter%phi 
      phot%x_ini  = Emitter%x
      phot%y_ini  = Emitter%y
      phot%z_ini  = Emitter%z
      phot%z_tau  = Emitter%z_tau
      phot%Vector_of_Momentum_ini = Emitter%Vector_of_Momentum
      phot%Vector_of_position_ini = Emitter%Vector_of_position
      !phot%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF 
      !phot%Phot4k_CovCF_ini = Emitter%Phot4k_CovCF 
      phot%cosphi_ini = Emitter%cosphi_ini
      phot%sinphi_ini = Emitter%sinphi_ini
      phot%cos2phi_ini = Emitter%cos2phi_ini
      phot%sin2phi_ini = Emitter%sin2phi_ini
      phot%phi_ini = Emitter%phi_ini 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%w_ini = Emitter%w_ini_em 
      phot%w_ini0 = Emitter%w_ini_em 
      phot%Psi_I = one
      phot%Psi_Q = one / 4.D0
      phot%Psi_U = one / 4.D0
      phot%Psi_V = one
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
      RETURN
      END SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon

!************************************************************************************ 
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p( phot ) 
!************************************************************************************
      IMPLICIT NONE 
      TYPE(Photon_FlatSP), INTENT(INOUT) :: phot 
      integer cases
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      phot%z_tau = phot%Get_scatter_distance_tau( )   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL phot%Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi() 
      phot%Psi_I = phot%Psi_I * phot%NormalA  
      phot%Psi_Q = phot%Psi_Q * phot%NormalA  
      phot%Psi_U = phot%Psi_U * phot%NormalA  
      phot%Psi_V = phot%Psi_V * phot%NormalA  
      !write(*,*)'ss2=', phot%z_tau, phot%p_scattering, phot%Vector_of_Momentum_ini(3) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      phot%Vector_of_Momentum_p = phot%Vector_of_Momentum_ini
      !phot%Vector_of_position_p = phot%Vector_of_position_ini
      !write(unit = *, fmt = *)'************************************************************'

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p

!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE( phot, sphot )
!************************************************************************************
      IMPLICIT NONE 
      TYPE(Photon_FlatSP), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer :: mu_i
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Vector_of_Momentum_ini = phot%Vector_of_Momentum_p
      !sphot%Vector_of_position_ini = phot%Vector_of_position_p
      sphot%Psi_I = phot%Psi_I
      sphot%Psi_Q = phot%Psi_Q
      sphot%Psi_U = phot%Psi_U
      sphot%Psi_V = phot%Psi_V
      sphot%cosphi_ini = phot%cosphi_ini
      sphot%sinphi_ini = phot%sinphi_ini
      sphot%cos2phi_ini = phot%cos2phi_ini
      sphot%sin2phi_ini = phot%sin2phi_ini
      sphot%phi_ini = phot%phi_ini
      !sphot%Vector_of_Momentum_ini(3) = phot%Vector_of_Momentum_ini(3) 
      !sphot%z_tau = phot%z_tau
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !write(*,*)'5555==', phot%E_ini, phot%Phot4k_CovCF_ini(1)   
   If( .true. )then
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(1)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(2)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(3)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(4)
   else
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(1) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(2) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(3) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(4) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(4)
   endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection() 
      !write(*,*)'6666==', phot%E_ini, phot%Phot4k_CovCF_ini(1)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering( phot, sphot )
!************************************************************************************
      IMPLICIT NONE 
      TYPE(Photon_FlatSP), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%Vector_of_Momentum_ini(3) = sphot%Vector_of_Momentum_ini(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      if( isnan( phot%Vector_of_Momentum_ini(3) ) )write(*, *)'mms=', phot%Phot4k_CtrCF_ini 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      !phot%phi_estimat = sphot%phi_estimat

      phot%Psi_I = sphot%Psi_I
      phot%Psi_Q = sphot%Psi_Q
      phot%Psi_U = sphot%Psi_U
      phot%Psi_V = sphot%Psi_V
      phot%cosphi_ini = sphot%cosphi_ini
      phot%sinphi_ini = sphot%sinphi_ini
      phot%cos2phi_ini = sphot%cos2phi_ini
      phot%sin2phi_ini = sphot%sin2phi_ini
      phot%phi_ini = sphot%phi_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !phot%E_ini = DABS( phot%Phot4k_CovCF_ini(1) ) 
      !write(*,*)'5555==', phot%E_ini, phot%Phot4k_CovCF_ini(1) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site( phot, sphot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE 
      TYPE(Photon_FlatSP), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      phot%z_tau = phot%Get_scatter_distance_tau( ) 
      !write(*,*)'7777==', phot%z_tau
      !write(*,*)'7777==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat()
      !CALL phot%Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi() 
      !phot%w_ini = phot%w_ini * phot%NormalA   
      phot%Psi_I = phot%Psi_I * phot%NormalA
      phot%Psi_Q = phot%Psi_Q * phot%NormalA  
      phot%Psi_U = phot%Psi_U * phot%NormalA  
      phot%Psi_V = phot%Psi_V * phot%NormalA   
      !write(*,*)'ss3=',phot%w_ini, phot%NormalA, phot%r_one_hvlmec2_one_cosE    
      phot%Vector_of_Momentum_p = phot%Vector_of_Momentum_ini  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering( phot, sphot )
!************************************************************************************
      IMPLICIT NONE 
      TYPE(Photon_FlatSP), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Vector_of_Momentum_ini = phot%Vector_of_Momentum_p
      !sphot%Vector_of_position_ini = phot%Vector_of_position_p
      sphot%Psi_I = phot%Psi_I
      sphot%Psi_Q = phot%Psi_Q
      sphot%Psi_U = phot%Psi_U
      sphot%Psi_V = phot%Psi_V
      sphot%cosphi_ini = phot%cosphi_ini
      sphot%sinphi_ini = phot%sinphi_ini
      sphot%cos2phi_ini = phot%cos2phi_ini
      sphot%sin2phi_ini = phot%sin2phi_ini
      sphot%phi_ini = phot%phi_ini
      !sphot%Vector_of_Momentum_ini(3) = phot%Vector_of_Momentum_ini(3) 
      !sphot%z_tau = phot%z_tau   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2( sphot%mu_estimat(1) ) 
   If( .true. )then
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(1)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(2)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(3)
      CALL phot%Calc_Phot_Inform_At_Observer_with_mu_phi_Given(4)
   else
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(1) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(2) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(3) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection2(4) 
      phot%f_IQUV_estimat = sphot%f_IQUV_estimat 
      phot%phi_estimat = sphot%phi_estimat
      CALL phot%Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat(4)
   endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Tompson_Scat_With_Polarized_Diffuse_Reflection() 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!**************************************************************************************
    SUBROUTINE mimick_of_ph_Emerge_DiffuseReflec( Total_Phot_Num, tau )
!************************************************************************************** 
    implicit none
    real(mcp), intent(inout) :: tau
    real(mcp) :: E, E_low, E_up  
    integer(kind = 8) :: Num_Photons
    integer(kind = 8), intent(in) :: Total_Phot_Num 
    type(Photon_Emitter) :: Emitter
    type(Photon_FlatSP) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE)   
    !real(mcp) :: I_Recv(0 : Num_PolDeg), Q_Recv(0 : Num_PolDeg) 
    real(mcp) :: IQUV10_Recv(1:4, 0 : Num_phi), IQUV60_Recv(1:4, 0 : Num_phi), &
                 IQUV30_Recv(1:4, 0 : Num_phi), IQUV80_Recv(1:4, 0 : Num_phi) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer :: i,j 
!--------------MPI---------------------------------------------------------------
    integer np, myid, ierr
    integer :: namelen
    character*(MPI_MAX_PROCESSOR_NAME) processor_name
    integer(kind = 8) :: mydutyPhot
!===========================================================================
    call MPI_INIT (ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, np, ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
    call MPI_GET_PROCESSOR_NAME(processor_name, namelen, ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write(*,*)'MyId Is:  ', np, myid, namelen
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    if( myid == np-1 ) then
        mydutyphot = Total_Phot_Num / np + &
                   mod( Total_Phot_Num, np )
    else
        mydutyphot = Total_Phot_Num / np
    endif 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL Set_initial_parameter_values( Phot, Emitter, tau ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    If( myid == np - 1 )then
        call Calculate_The_Diffuse_Reflection_of_Chandra( Phot )
    endif
 
    Num_Photons = 0   
    sphot%mu_estimat = phot%mu_estimat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0   
        CALL Emitter_A_Photon( Emitter )
        CALL Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )
        CALL Determine_P_Of_Scatt_Site_And_Quantities_At_p( phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        CALL FIRST_SCATTERING_OF_PHOT_ELCE( phot, sphot )
  !write(*,*)'****************************************************************************'
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1
            !write(*,*)'ss2===', phot%scatter_times, phot%z_tau!, phot%E_ini, phot%nu_up*h_ev/1.D6 
            CALL Set_InI_Conditions_For_Next_Scattering( phot, sphot )   
            CALL Determine_Next_Scattering_Site( phot, sphot )
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
            !if( phot%Optical_Depth_scatter >= 1.D2 )exit
             !if( phot%I_IQ / phot%w_ini0 <= 1.D-6 )exit
              if( phot%Psi_I <= 1.D-4 .or. phot%z_tau > 100.D0 )exit
            !if( phot%scatter_times > 100 )exit
            !if( phot%z_tau > 50.D0 .or. dabs(phot%I_IQ) <= 1.D-5) exit
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CALL Photon_Electron_Scattering( phot, sphot ) 
            !if(phot%scatter_times > 10)stop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
        If ( mod(Num_Photons, 1000)==0 .and. myid == np-1 ) then 
        write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'***** The', Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ', myid, mydutyphot 
        write(unit = *, fmt = *)'*************************************************************************'
        endif
    If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enddo  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          send_num = ( Num_phi + 1 ) * 4
          send_tag = 1  
          call MPI_SEND( phot%PolarArrayIQUV10, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayI to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          send_tag = 2  
          call MPI_SEND( phot%PolarArrayIQUV30, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayQ to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          send_tag = 3
          call MPI_SEND( phot%PolarArrayIQUV60, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayU to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          send_tag = 4
          call MPI_SEND( phot%PolarArrayIQUV80, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayV to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      else 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  recv_num = ( Num_phi + 1 ) * 4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( IQUV10_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                1, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives I data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( IQUV30_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                2, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives Q data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( IQUV60_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                3, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives U data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( IQUV80_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                4, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives U data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  phot%PolarArrayIQUV10 = phot%PolarArrayIQUV10 + IQUV10_Recv
                  phot%PolarArrayIQUV30 = phot%PolarArrayIQUV30 + IQUV30_Recv
                  phot%PolarArrayIQUV60 = phot%PolarArrayIQUV60 + IQUV60_Recv
                  phot%PolarArrayIQUV80 = phot%PolarArrayIQUV80 + IQUV80_Recv  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              enddo
          endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!' ,phot%times_counter 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !do i = 0, Num_PolDeg
              !phot%PolarArrayd( i ) = DSQRT( phot%PolarArrayQ(i)**2 + phot%PolarArrayU(i)**2 &
              !                               ) / phot%PolarArrayI(i)
          !enddo 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          open(unit=10, file='./spectrum/phi_dist/IQUV10_mu0=08_nozeroQUV.txt', status="replace")  
          open(unit=11, file='./spectrum/phi_dist/IQUV30_mu0=08_nozeroQUV.txt', status="replace")  
          open(unit=12, file='./spectrum/phi_dist/IQUV60_mu0=08_nozeroQUV.txt', status="replace")  
          open(unit=13, file='./spectrum/phi_dist/IQUV80_mu0=08_nozeroQUV.txt', status="replace") 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          do j = 0, Num_phi
              write(unit = 10, fmt = 100)phot%PolarArrayIQUV10(1:4, j) 
              write(unit = 11, fmt = 100)phot%PolarArrayIQUV30(1:4, j) 
              write(unit = 12, fmt = 100)phot%PolarArrayIQUV60(1:4, j) 
              write(unit = 13, fmt = 100)phot%PolarArrayIQUV80(1:4, j) 
              100 FORMAT(' ', 4F22.12)
          enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          close(unit=10)  
          close(unit=11)  
          close(unit=12) 
          close(unit=13)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    call MPI_FINALIZE ( ierr ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE mimick_of_ph_Emerge_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    
