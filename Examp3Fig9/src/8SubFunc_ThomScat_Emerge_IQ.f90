!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Method_Of_FLST_ThomScat_Emerge_IQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !USE constants
      USE RandUtils
      USE Photons_FlatSP
      USE MPI
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS

!**************************************************************************************
    SUBROUTINE mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau )
!************************************************************************************** 
    implicit none
    real(mcp), intent(inout) :: tau
    real(mcp) :: E, E_low, E_up  
    integer(kind = 8) :: Num_Photons
    integer(kind = 8), intent(in) :: Total_Phot_Num 
    type(Photon_Emitter_BB) :: Emitter
    type(Photon_FlatSP) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE)   
    real(mcp) :: I_Recv(0 : Num_PolDeg), Q_Recv(0 : Num_PolDeg)  
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
    CALL phot%Set_initial_parameter_values( Emitter, tau )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !CALL Set_xi_wi_all()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !CALL phot%Set_Cross_Section_3Te()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !CALL Set_psi_phi_chi_zera_array()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Num_Photons = 0 
    !phot%v_L_v_i = zero   
    phot%tau_max = tau
    !phot%Z_max = tau / Sigma_T / phot%n_e
    !phot%d_theta = pi / two / 100.D0 
    phot%PolarArrayI = zero
    phot%PolarArrayQ = zero   
    phot%PolarArrayU = zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0   
        CALL phot%Emitter_A_Photon( Emitter )
        CALL phot%Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter )
        CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        CALL phot%FIRST_SCATTERING_OF_PHOT_ELCE( sphot )
  !write(*,*)'****************************************************************************'
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1
            !write(*,*)'ss2===', phot%scatter_times, phot%p_scattering!, phot%E_ini, phot%nu_up*h_ev/1.D6 
            CALL phot%Set_InI_Conditions_For_Next_Scattering( sphot )   
            CALL phot%Determine_Next_Scattering_Site( sphot )
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
            !if(phot%I_IQ < 0.D0)write(*,*)'w_ini =', phot%w_ini, phot%I_IQ!, phot%Optical_Depth_scatter,&
                   !phot%z_tau
             !write(*,*)'w_ini =',phot%w_ini, phot%I_IQ, phot%Q_IQ , phot%z_tau, phot%scatter_times !,&
              !      phot%z_tau
             !if( phot%scatter_times >= 2.D3 )exit
             !if( phot%I_IQ / phot%w_ini0 <= 1.D-6 )exit
              if( phot%I_IQ <= 1.D-25 )exit
            !if( phot%scatter_times > 100 )exit
            !if( phot%z_tau > 50.D0 .or. dabs(phot%I_IQ) <= 1.D-5) exit
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CALL phot%Photon_Electron_Scattering( sphot )
            !if( dabs( sphot%Scattered_Phot4k_CF(1) ) >= 1.D-1 * mec2 )exit
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
        If ( mod(Num_Photons, 5000)==0 ) then 
        write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'***** The', Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ',myid, mydutyphot,  Num_PolDeg   
        write(unit = *, fmt = *)'*************************************************************************'
        endif
    If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enddo  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          send_num = Num_PolDeg + 1
          send_tag = 1  
          call MPI_SEND( phot%PolarArrayI, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayI to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
          send_tag = 2  
          call MPI_SEND( phot%PolarArrayQ, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayQ to Processor:', np-1  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      else 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  recv_num = Num_PolDeg + 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( I_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                1, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives I data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( Q_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                2, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives Q data from Processor:', RECV_SOURCE   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                  phot%PolarArrayI = phot%PolarArrayI + I_Recv
                  phot%PolarArrayQ = phot%PolarArrayQ + Q_Recv   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              enddo
          endif 
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!'
          do i = 0, Num_PolDeg
              phot%PolarArrayd( i ) = DSQRT( phot%PolarArrayQ(i)**2 + phot%PolarArrayU(i)**2 &
                                             ) / phot%PolarArrayI(i) 
          enddo 
          open(unit=10, file='./spectrum/EIQ_I_tau=50.txt', status="replace") 
          open(unit=11, file='./spectrum/EIQ_Q_tau=50.txt', status="replace")     
          do j = 0, Num_PolDeg
              write(unit = 10, fmt = *)phot%PolarArrayI(j)
              write(unit = 11, fmt = *)phot%PolarArrayQ(j) 
          enddo 
          close(unit=10) 
          close(unit=11)   
      endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    call MPI_FINALIZE ( ierr ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE mimick_of_ph_finity_zone_Emerge_IQ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END MODULE Method_Of_FLST_ThomScat_Emerge_IQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    
