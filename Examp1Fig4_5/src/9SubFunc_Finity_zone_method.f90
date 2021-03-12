!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !USE constants
      USE RandUtils
      USE SemiAnalyMethod
      !USE PhotonEmitter
      !USE Photons
      USE MPI
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS 
!**************************************************************************************
    SUBROUTINE mimick_of_photon_with_finity_zone( Total_Phot_Num )
!************************************************************************************** 
    implicit none 
    real(mcp) :: T_e, E, dt, t = 0.D0, p_length, E_low, E_up, Tb = 2.D-3
    real(mcp), parameter :: n_e1 = 1.6D17 
    integer(kind = 8) :: Num_Photons
    integer(kind = 8), intent(in) :: Total_Phot_Num  
    !real(mcp) :: R_out, R_in, rea, img, E_max = 1.D15*h_ev*1.D-6 
    !type(Photon_Emitter) :: Emitter
    type(Photon) :: phot 
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE) 
    !integer :: effect_Photon_Number_Sent 
    !integer :: effect_Photon_Number_Recv, cases
    !real(mcp), dimension(0:400, 0:400) :: disk_image_Sent, disk_image_Recv
    !real(mcp), dimension(1:500, 1:1000) :: ET_array_Recv, ET_array_Sent
    real(mcp) :: v_L_v_Sent(0:500), v_L_v_Recv(0:500) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
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
!| Set Initial conditions for the Photon                     !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Num_Photons = 0 
    phot%v_L_v_i = zero 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL phot%Set_initial_parameter_values( ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( myid == np-1 ) then
        !CALL Emitter_A_Photon2( Emitter, phot )  
    endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !****************************************************
    Do
    !****************************************************
        Num_Photons = Num_Photons + 1  
        CALL phot%Emitter_A_Photon( )   
        CALL phot%Calc_Phot_Informations_At_Observor_2zones( ) 
 
        If ( mod(Num_photons,5000000)==0 .and. myid==np-1 ) then 
        !write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = *)'***** The',Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ', mydutyphot 
        write(unit = *, fmt = *)'*************************************************************************'
        endif
        If( Num_Photons > mydutyphot )EXIT 
    Enddo  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then
 
          send_num = 500
          send_tag = 2
          v_L_v_Sent = phot%v_L_v_i
          call MPI_SEND(v_L_v_Sent, send_num, MPI_DOUBLE_PRECISION, np-1, &
                                    send_tag, MPI_COMM_WORLD, ierr)
          write(*,*)'MPI_Processor ', myid, ' has send v_L_v_Sent to Master Processor!'
  
      else  
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2 
                  recv_num = 500
                  call MPI_RECV(v_L_v_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE,&
                                2, MPI_COMM_WORLD, status, ierr) 
                  !write(*,*)'sdf500==', v_L_v_Recv
                  write(*,*)'master Processor ',myid,' has receved v_L_v Processor:', RECV_SOURCE  
                  phot%v_L_v_i = phot%v_L_v_i + v_L_v_Recv 
              enddo
          endif 
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!' 
          open(unit=17, file='./spectrum/vLv_MC_2.txt', status="replace")  
          do j = 0, 499
              write(unit = 17, fmt = *)phot%v_L_v_i(j)
          enddo 
          close(unit=17)  
      endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!===========================================================
    call MPI_FINALIZE ( ierr )
!===========================================================
    END SUBROUTINE mimick_of_photon_with_finity_zone  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END MODULE Statistial_Method_Of_Finity_zone
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~














 
