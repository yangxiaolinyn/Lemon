!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    USE RandUtils 
    USE Photons
    IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CONTAINS 
!**************************************************************************************
    SUBROUTINE mimick_of_photon_with_finity_zone( Total_Phot_Num, tau, &
                           T_e, T_s, n_e, y1, y2, CrossSec_filename, MCResults )
!************************************************************************************** 
    implicit none
    integer(kind = 8), intent(in) :: Total_Phot_Num
    real(mcp), intent(in) :: tau, T_e, T_s, n_e, y1, y2
    character*80, intent(in) :: CrossSec_filename, MCResults
    integer(kind = 8) :: scatter_Times = 0
    integer(kind = 8) :: Num_Photons 
    logical :: At_outer_Shell, At_inner_Shell 
    type(Photon) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE) 
    integer :: effect_Photon_Number_Sent 
    integer :: effect_Photon_Number_Recv, cases  
    real(mcp) :: v_L_v_Sent(0: Num_y), v_L_v_Recv(0: Num_y), temp_ET 
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
    phot%my_ID = myid
    phot%num_process = np
    call phot%Set_Initial_Parameters_And_Conditions( tau, T_e, &
                                   T_s, n_e, y1, y2, CrossSec_filename ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0   
        CALL phot%Emitter_A_Photon( ) 
        CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( )    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        CALL phot%Photon_Electron_Scattering( phot%T_e )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1
            !write(*,*)'ss2===', phot%scatter_times, phot%w_ini, phot%w_ini_em
            !if( phot%scatter_times > 1.D5)write(*,*)'ss===',phot%scatter_times,phot%p_scattering, phot%R_in
            !If ( myid == np-1 ) write(*,*)'Scatter_Times = ', phot%W_ini, phot%w_p, Scatter_Times
            !If ( myid == np-1 ) write(*,*)'Scatter_Times = ', phot%dv_ini, phot%dv_p, &
            !sphot%Scal_Factor_after, sphot%Scal_Factor_before, Scatter_Times
            CALL phot%Set_InI_Conditions_For_Next_Scattering( ) 
            CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( )  
            if( phot%w_ini / phot%w_ini_em <= 1.D-15 )exit 
            !if( phot%scatter_times > 100 )stop
            CALL phot%Photon_Electron_Scattering( phot%T_e )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
        If ( mod(Num_photons,20000)==0 .and. myid == np-1 ) then 
        write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'***** The',Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ',myid, mydutyphot   
        write(unit = *, fmt = *)'*************************************************************************'
        endif
        If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enddo 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      If ( myid /= np-1 ) then 

          send_num = Num_y + 1
          send_tag = 2
          v_L_v_Sent = phot%v_L_v_i
          call MPI_SEND(v_L_v_Sent,send_num,MPI_DOUBLE_PRECISION,np-1,send_tag,MPI_COMM_WORLD,ierr)
          write(*,*)'Processor ',myid,' has send v_L_v_Sent'
  
      else  
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2   
                  recv_num = Num_y + 1
                  call MPI_RECV(v_L_v_Recv,recv_num,MPI_DOUBLE_PRECISION,RECV_SOURCE,&
                                2,MPI_COMM_WORLD,status,ierr)  
                  write(*,*)'master Processor ',myid,' has receved v_L_v Processor:', RECV_SOURCE 
                  phot%v_L_v_i = phot%v_L_v_i + v_L_v_Recv 
              enddo
          endif 

          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!'
 
          open(unit=17, file = MCResults, status="replace")  

          do j = 0, Num_y
              write(unit = 17, fmt = *)phot%v_L_v_i(j)
          enddo
          close(unit=17)  
      endif  
!===========================================================
      call MPI_FINALIZE ( ierr )
!===========================================================
      END SUBROUTINE mimick_of_photon_with_finity_zone
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      END MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 














 
