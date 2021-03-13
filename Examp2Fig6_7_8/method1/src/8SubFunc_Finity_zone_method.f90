!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    USE RandUtils 
    USE Photons
    USE MPI 
    IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CONTAINS 
!**************************************************************************************
    SUBROUTINE mimick_of_photon_with_finity_zone( Total_Phot_Num, tau, T_e, T_s, n_e )
!************************************************************************************** 
    implicit none
    integer(kind = 8), intent(in) :: Total_Phot_Num
    real(mcp), intent(in) :: tau, T_e, T_s, n_e
    !real(mcp) :: T_e, E, dt, t = 0.D0, p_length, E_low, E_up 
    integer(kind = 8) :: scatter_Times = 0
    integer(kind = 8) :: Num_Photons 
    !real(mcp) :: R_out, R_in, rea, img
    logical :: At_outer_Shell, At_inner_Shell
    !type(Photon_Emitter) :: Emitter
    type(Photon) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE) 
    integer :: effect_Photon_Number_Sent 
    integer :: effect_Photon_Number_Recv, cases
    real(mcp), dimension(0:400, 0:400) :: disk_image_Sent, disk_image_Recv
    real(mcp), dimension(1:500, 1:1000) :: ET_array_Recv, ET_array_Sent
    real(mcp) :: v_L_v_Sent(1:500), v_L_v_Recv(1:500), temp_ET 
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Num_Photons = 0 
    phot%v_L_v_i = zero
    call phot%Set_Initial_Parameters_And_Conditions( tau, T_e, T_s, n_e ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0
        phot%At_outer_Shell = .false.
        phot%At_inner_Shell = .false.
        phot%time_arrive_observer = zero
        phot%time_travel = zero
        CALL phot%Emitter_A_Photon( ) 
        CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( T_e, phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            If ( phot%At_outer_Shell ) then
                cases = 2
                !Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
                phot%At_outer_Shell = .False.
                phot%Go2infinity = .TRUE.
                goto 112 
            endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CALL phot%FIRST_SCATTERING_OF_PHOT_ELCE( T_e, phot, sphot )
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
            CALL phot%Set_InI_Conditions_For_Next_Scattering( T_e, phot, sphot )   
            CALL phot%Determine_Next_Scattering_Site( T_e, phot, sphot, Scatter_Times )
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
            CALL phot%Photon_Electron_Scattering( T_e, phot, sphot )
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
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ',myid, mydutyphot 
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
  
      else  
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
      END MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 














 
