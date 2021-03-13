  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE constants
      USE RandUtils
      USE PhotonEmitter
      USE Photons
      USE MPI
      !USE ScatterPhoton
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon( Emitter )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set the Emitter's BL Coordinate and four Velocity ComFrameU4_BL, from which one      ||
!| Can compute the physical velocity of the Emitter with respect to the LNRF reference. ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      Emitter%r = zero!Emitter%R_out
      Emitter%theta = zero! * ranmar()
      Emitter%mucos = dcos( Emitter%theta )
      Emitter%musin = dsin( Emitter%theta )
      Emitter%phi = zero!twopi * ranmar()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Emitter%x = Emitter%r * Emitter%musin * dcos( Emitter%phi )
      Emitter%y = Emitter%r * Emitter%musin * dsin( Emitter%phi )
      Emitter%z = Emitter%r * Emitter%mucos
      Emitter%Vector_of_position(1) = Emitter%x
      Emitter%Vector_of_position(2) = Emitter%y
      Emitter%Vector_of_position(3) = Emitter%z
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL Emitter%get_Phot4k_CtrCF_CovCF() 
  
      Emitter%Vector_of_Momentum(1:3) = Emitter%Phot4k_CtrCF(2:4) / Emitter%Phot4k_CtrCF(1)  

      RETURN
      END SUBROUTINE Emitter_A_Photon

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(IN) :: Emitter
      TYPE(Photon), INTENT(INOUT) :: Phot 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      phot%r_ini      = Emitter%r
      phot%theta_ini  = Emitter%theta
      phot%mucos_ini  = Emitter%mucos
      phot%musin_ini  = Emitter%musin
      phot%phi_ini    = Emitter%phi 
      phot%x_ini  = Emitter%x
      phot%y_ini  = Emitter%y
      phot%z_ini  = Emitter%z
      phot%Vector_of_Momentum_ini = Emitter%Vector_of_Momentum
      phot%Vector_of_position_ini = Emitter%Vector_of_position
      phot%Phot4k_CtrCF_ini = Emitter%Phot4k_CtrCF 
      phot%Phot4k_CovCF_ini = Emitter%Phot4k_CovCF 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%E_ini = DABS( Emitter%Phot4k_CovCF(1) )  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%w_ini = Emitter%w_ini_em 
  
      RETURN
      END SUBROUTINE Transmit_Data_And_Parameters_From_Emitter2Photon

!************************************************************************************
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p( T_e, phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!************************************************************************************
      IMPLICIT NONE
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(Photon), INTENT(INOUT) :: phot 
      integer cases
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
      !write(*,*)'1111==', phot%E_ini, phot%r_ini 
      phot%p_scattering = phot%Get_scatter_distance2( T_e )  
      !write(*,*)'2222==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !phot%direct_escaped = .True.
      cases = 1
      Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
      !write(*,*)'ss1=',phot%w_ini
      phot%w_ini = phot%w_ini * phot%NormalA
      !write(*,*)'ss2=',phot%w_ini 
      !phot%direct_escaped = .false.
      !write(*,*)'3333=', phot%E_ini, phot%Phot4k_CovCF_ini(1)
      !write(*,*)'111=', phot%time_travel 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !If ( phot%At_outer_Shell ) RETURN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%time_travel = phot%time_travel + phot%p_scattering / Cv 
          !write(*,*)'555=', phot%time_travel,  phot%p_scattering / CV
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                      phot%Vector_of_Momentum_ini, phot%p_scattering )
      phot%Vector_of_position_p = phot%Vector_of_position
  
      phot%Phot4k_CtrCF_At_p = phot%Phot4k_CtrCF_ini
      phot%Phot4k_CovCF_At_p = phot%Phot4k_CovCF_ini 
      !write(unit = *, fmt = *)'************************************************************'

      RETURN
      END SUBROUTINE Determine_P_Of_Scatt_Site_And_Quantities_At_p

!************************************************************************************
      SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE( T_e, phot, sphot )
!************************************************************************************
      IMPLICIT NONE
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      sphot%Phot4k_CtrCF = phot%Phot4k_CtrCF_At_p 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Photon_Tetrad_In_CF()
      CALL sphot%Get_gama_mu_phi_Of_Scat_Elec(T_e)
      CALL sphot%Set_Elec_Tetrad_In_CF()
      CALL sphot%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      CALL sphot%Compton_Scattering_WithOut_Polarizations()      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE FIRST_SCATTERING_OF_PHOT_ELCE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!************************************************************************************
      SUBROUTINE Set_InI_Conditions_For_Next_Scattering( T_e, phot, sphot )
!************************************************************************************
      IMPLICIT NONE
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Now we transform All quantities that in the Electron Comoving Frame into the LNRF 
!   Frame and to the BL coordinates, wiht which then we can calculate the corresponding 
!   quantities, and take them as the initial conditions for the next scattering. Note 
!   that affer the first scattering, the Polarisation of the photon is not zero.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      phot%Phot4k_CtrCF_ini = sphot%Scattered_Phot4k_CF
      phot%Phot4k_CovCF_ini = sphot%Scattered_Phot4k_CovCF
      phot%Vector_of_Momentum_ini(1:3) = phot%Phot4k_CtrCF_ini(2:4) / phot%Phot4k_CtrCF_ini(1)
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      phot%E_ini = DABS( phot%Phot4k_CovCF_ini(1) ) 
      !write(*,*)'5555==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%r_ini      = phot%r_p
      phot%vector_of_position_ini = phot%vector_of_position_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END SUBROUTINE Set_InI_Conditions_For_Next_Scattering
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Determine_Next_Scattering_Site( T_e, phot, sphot, Scatter_Times )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      REAL(mcp), INTENT(IN) :: T_e
      integer(kind = 8), intent(IN) :: Scatter_Times
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
      integer cases
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !write(*,*)'6666==', phot%E_ini, phot%Phot4k_CovCF_ini(1)  
      phot%p_scattering = phot%Get_scatter_distance2( T_e )  
      !write(*,*)'7777==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !phot%direct_escaped = .True.
      cases = 1
      Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
      !write(*,*)'ss3=',phot%w_ini
      phot%w_ini = phot%w_ini * phot%NormalA
      !write(*,*)'ss4=',phot%w_ini
      !phot%direct_escaped = .false. 
      !write(*,*)'8888==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
          !write(*,*)'666=', phot%time_travel 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !If ( phot%At_outer_Shell ) RETURN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%time_travel = phot%time_travel + phot%p_scattering / Cv
          !write(*,*)'101=', phot%time_travel,  phot%p_scattering/ CV
      phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                 phot%Vector_of_Momentum_ini, phot%p_scattering ) 
      if( phot%r_p > phot%R_out )then
          phot%test_it = .true.
          phot%p_scattering = phot%Get_scatter_distance2( T_e )  
          write(*,*)'101=', phot%time_travel,  phot%p_scattering/ CV
          write(*,*)'102=', phot%r_p, phot%p_scattering, phot%NormalA
          write(*,*)'103=', dsqrt( phot%Vector_of_position_ini(1)**2+&
                   phot%Vector_of_position_ini(2)**2+phot%Vector_of_position_ini(3)**2 ), phot%r_ini, &
               dsqrt( phot%Vector_of_Momentum_ini(1)**2+phot%Vector_of_Momentum_ini(2)**2+&
               phot%Vector_of_Momentum_ini(3)**2 ), phot%r_p2( phot%Vector_of_position_ini, &
                 phot%Vector_of_Momentum_ini, phot%p_scattering ) 
          stop
      endif
      ! the vector of position also obtained.
      phot%Vector_of_position_p = phot%Vector_of_position 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !write(unit = *, fmt = *)'************************************************************'
      phot%Phot4k_CtrCF_At_p = phot%Phot4k_CtrCF_ini
      phot%Phot4k_CovCF_At_p = phot%Phot4k_CovCF_ini
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Determine_Next_Scattering_Site
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!************************************************************************************
      SUBROUTINE Photon_Electron_Scattering( T_e, phot, sphot )
!************************************************************************************
      IMPLICIT NONE
      REAL(mcp), INTENT(IN) :: T_e
      TYPE(Photon), INTENT(INOUT) :: phot
      TYPE(ScatPhoton), INTENT(INOUT) :: sphot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      sphot%Phot4k_CtrCF = phot%Phot4k_CtrCF_At_p 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Set_Photon_Tetrad_In_CF()
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Get_gama_mu_phi_Of_Scat_Elec(T_e)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL sphot%Set_Elec_Tetrad_In_CF()
      CALL sphot%Set_Phot4k_In_Elec_CF() 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CALL sphot%Compton_Scattering_WithOut_Polarizations()  
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END SUBROUTINE Photon_Electron_Scattering
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!**************************************************************************************
    SUBROUTINE mimick_of_photon_with_finity_zone( Total_Phot_Num, tau )
!**************************************************************************************
    !use Basic_Variables_And_Methods
    !use PhotonEmitter
    !use Photons
    !use ScatterPhoton
    !use SubroutineFunction
    implicit none
    integer(kind = 8), intent(in) :: Total_Phot_Num
    real(mcp), intent(in) :: tau
    real(mcp) :: T_e, E, dt, t = 0.D0, p_length, E_low, E_up
    real(mcp), parameter :: n_e1 = 1.6D17
    integer(kind = 8) :: scatter_Times = 0
    integer(kind = 8) :: Num_Photons
    integer(kind = 4), parameter :: n_sample = 1000000
    real(mcp) :: R_out, R_in, rea, img
    logical :: At_outer_Shell, At_inner_Shell
    type(Photon_Emitter) :: Emitter
    type(Photon) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE) 
    integer :: effect_Photon_Number_Sent 
    integer :: effect_Photon_Number_Recv, cases
    real(mcp), dimension(0:400, 0:400) :: disk_image_Sent, disk_image_Recv
    real(mcp), dimension(1:500, 1:1000) :: ET_array_Recv, ET_array_Sent
    real(mcp) :: v_L_v_Sent(1:500), v_L_v_Recv(1:500), temp_ET
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Complex*16 :: r1,r2,r3,r4,s1,s2,s3,temp(1:4),temp1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    integer i,j,del,reals
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
    call InitRandom( myid )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    if( myid == np-1 ) then
        mydutyphot = Total_Phot_Num / np + &
                   mod( Total_Phot_Num, np )
    else
        mydutyphot = Total_Phot_Num / np
    endif
    !stop 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!| Set Initial conditions for the Emitter
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Emitter%R_out = R_out 
    !Emitter%aspin = a_spin
    Emitter%T_s = 1.D-8 * mec2
    CALL Emitter%Set_Emin_Emax()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Set Initial conditions for the Photon                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    E_low = Emitter%E_low1 !1.D-5
    E_up = 2.D1!this%E_up1 !1.D1
    phot%nu_low = Emitter%nu_low
    phot%nu_up = Emitter%nu_up
    phot%T_e = 4.D0 * mec2
    T_e = phot%T_e
    phot%logE_low = DLOG10(E_low)
    phot%logE_up = DLOG10(E_up)
    phot%n_e = 1.D17
    phot%R_out = tau / Sigma_T / phot%n_e
    !phot%aspin = a_spin
    write(*,*)'ddds=',phot%R_out 
    !CALL phot%Set_Cross_Section( T_e )
    CALL phot%Set_Cross_Section_3Te()
    phot%effect_number = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Num_Photons = 0
    phot%delta_pds = zero
    phot%v_L_v_i = zero
    phot%r_obs = R_out!2.63D0 * Cv *1.D0
    phot%t_standard =  phot%R_out / Cv  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !****************************************************
    Do
    !****************************************************
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0
        phot%At_outer_Shell = .false.
        phot%At_inner_Shell = .false.
        phot%time_arrive_observer = zero
        phot%time_travel = zero
        CALL Emitter_A_Photon( Emitter )
        CALL Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )
        CALL Determine_P_Of_Scatt_Site_And_Quantities_At_p( T_e, phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            If ( phot%At_outer_Shell ) then
                cases = 2
                !Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
                phot%At_outer_Shell = .False.
                phot%Go2infinity = .TRUE.
                goto 112 
            endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CALL FIRST_SCATTERING_OF_PHOT_ELCE( T_e, phot, sphot )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !scatter_Times = Scatter_Times + 1
            phot%scatter_times = phot%scatter_times + 1
            !write(*,*)'ss2===', phot%scatter_times, phot%p_scattering, phot%E_ini, phot%nu_up*h_ev/1.D6
            !if( phot%scatter_times > 1.D5)write(*,*)'ss===',phot%scatter_times,phot%p_scattering, phot%R_in
            !If ( myid == np-1 ) write(*,*)'Scatter_Times = ', phot%W_ini, phot%w_p, Scatter_Times
            !If ( myid == np-1 ) write(*,*)'Scatter_Times = ', phot%dv_ini, phot%dv_p, &
            !sphot%Scal_Factor_after, sphot%Scal_Factor_before, Scatter_Times
            CALL Set_InI_Conditions_For_Next_Scattering( T_e, phot, sphot )   
            CALL Determine_Next_Scattering_Site( T_e, phot, sphot, Scatter_Times )
           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
            If ( phot%At_outer_Shell ) then
                cases = 2
                !Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
                phot%At_outer_Shell = .False.
                phot%Go2infinity = .TRUE. 
                Exit 
            endif
            !write(*,*)'ss=',phot%w_ini
            if( phot%w_ini <= 1.D-15 )exit
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CALL Photon_Electron_Scattering( T_e, phot, sphot )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

112     continue
        If ( mod(Num_photons,1000000)==0 ) then
        !write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'*************************************************************************'
        !write(unit = *, fmt = *)'*****                                                               *****' 
        !write(unit = *, fmt = *)'*****                                                               *****' 
        write(unit = *, fmt = *)'***** The',Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ',myid, mydutyphot, R_in
        !write(unit = *, fmt = *)'***** The Photon goes to infinity is', phot%Go2infinity
        !write(unit = *, fmt = *)'***** The Photon has fall into BH is', phot%fall2BH
        !write(unit = *, fmt = *)'*****                                                               *****' 
        !write(unit = *, fmt = *)'*****                                                               *****' 
        !write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'*************************************************************************'
        endif
    If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enddo 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      If ( myid /= np-1 ) then
          send_num = 1
          send_tag = 0
          !effect_Photon_Number_Sent = phot%effect_number
          !call MPI_SEND(effect_Photon_Number_Sent,send_num,MPI_INTEGER8,np-1,send_tag,MPI_COMM_WORLD,ierr)
          !write(*,*)'Processor ',myid,' has send Effective Photon Number:', effect_Photon_Number_Sent

          !send_num = 500*1000
          !send_tag = 1 
          !ET_array_Sent = phot%v_L_v_i_ET
          !call MPI_SEND(ET_array_Sent ,send_num,MPI_DOUBLE_PRECISION,np-1,send_tag,MPI_COMM_WORLD,ierr)
          !write(*,*)'Processor ',myid,' has send disk_image:'
          !write(*,*)'ppp===', ET_array_Sent
 
          send_num = 500
          send_tag = 2
          v_L_v_Sent = phot%v_L_v_i
          call MPI_SEND(v_L_v_Sent,send_num,MPI_DOUBLE_PRECISION,np-1,send_tag,MPI_COMM_WORLD,ierr)
          write(*,*)'Processor ',myid,' has send v_L_v_Sent'
 
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      else 
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !call MPI_RECV(effect_Photon_Number_Recv,1,MPI_INTEGER8,RECV_SOURCE,&
                  !              0,MPI_COMM_WORLD,status,ierr)
                  !write(*,*)'master Processor ',myid,' has receved :',effect_Photon_Number_Recv, &
                  !          'Photons from Processor:', RECV_SOURCE
                  !phot%effect_number = phot%effect_number + effect_Photon_Number_Recv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !recv_num = 500*1000
                  !call MPI_RECV(ET_array_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                  !              1, MPI_COMM_WORLD, status, ierr) 
                  !write(*,*)'master Processor ',myid,' ET_array from Processor:', RECV_SOURCE
                  !write(*,*)'sdf==', ET_array_Recv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !do i = 0, 300
                  !    do j = 0, 300
                  !        phot%v_L_v_i_ET = phot%v_L_v_i_ET + ET_array_Recv
                  !    enddo 
                  !enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  recv_num = 500
                  call MPI_RECV(v_L_v_Recv,recv_num,MPI_DOUBLE_PRECISION,RECV_SOURCE,&
                                2,MPI_COMM_WORLD,status,ierr) 
                  !write(*,*)'sdf500==', v_L_v_Recv
                  write(*,*)'master Processor ',myid,' has receved v_L_v Processor:', RECV_SOURCE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  !do j = 1, 500
                      !phot%v_L_v_i(j) = phot%v_L_v_i(j) + v_L_v_Recv(j)
                  !enddo 
                  phot%v_L_v_i = phot%v_L_v_i + v_L_v_Recv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              enddo
          endif
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!'
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          !open(unit=15, file='tdiskg.txt', status="replace") 
          !open(unit=16, file='./image/ET2.txt', status="replace") 
          open(unit=17, file='./spectrum/vLv_finity_tau=004.txt', status="replace")
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          !do i = 1, 1000
          !    temp_ET = zero
          !    do j = 1,500
          !        temp_ET = temp_ET + phot%v_L_v_i_ET(j,i)
          !    enddo
          !    write(unit = 16, fmt = *)temp_ET
          !enddo
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          !do j = 0,100
          !    write(unit = 16, fmt = *)phot%delta_pds(j)
          !enddo
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          do j = 1,500
              write(unit = 17, fmt = *)phot%v_L_v_i(j)
          enddo
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          !close(unit=15) 
          !close(unit=16) 
          close(unit=17) 
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      endif  
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
!===========================================================
    call MPI_FINALIZE ( ierr )
!===========================================================
    END SUBROUTINE mimick_of_photon_with_finity_zone


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|        ||
!|        ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|        ||
!|        ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~














 
