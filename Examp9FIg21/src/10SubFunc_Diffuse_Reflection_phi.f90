!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE constants
      USE RandUtils 
      USE Photons_FlatSP 
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS  
!*********************************************************************** 
    SUBROUTINE mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, &
                        tau, T_bb, T_elec, E1_scat, E2_scat, mu_estis, &
                        CrossSec_filename, Hxm_analyticResult )
!*********************************************************************** 
    implicit none 
    real(mcp), intent(inout) :: tau, T_bb, T_elec, E1_scat, E2_scat, mu_estis(1: 8)
    integer(kind = 8), intent(in) :: Total_Phot_Num
    character*80, intent(inout) :: CrossSec_filename, Hxm_analyticResult
    real(mcp) :: E, E_low, E_up  
    integer(kind = 8) :: Num_Photons 
    type(Photon_Emitter_BB) :: Emitter
    type(Photon_FlatSP) :: phot
    type(ScatPhoton_KN) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, &
              status(MPI_STATUS_SIZE), send_num2, recv_num2  
    real(mcp) :: Imu_Recv(1: 8, 0: 6, 0 : vL_sc_up) 
    real(mcp) :: E0, Es, mu2, g1Emu2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    integer :: i,j, i1, j1, k1
!--------------MPI---------------------------------------- 
    integer np, myid, ierr
    integer :: namelen
    character*(MPI_MAX_PROCESSOR_NAME) processor_name
    integer(kind = 8) :: mydutyPhot
!=========================================================== 
    call MPI_INIT (ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, np, ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
    call MPI_GET_PROCESSOR_NAME(processor_name, namelen, ierr)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    write(*,*)'MyId Is:  ', np, myid, namelen
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    if( myid == np-1 ) then
        mydutyphot = Total_Phot_Num / np + &
                   mod( Total_Phot_Num, np )
    else
        mydutyphot = Total_Phot_Num / np
    endif 

    phot%myid = myid
    phot%num_np = np 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL phot%Set_Initial_Values_For_Photon_Parameters( T_elec, &
                          T_bb, tau, E1_scat, E2_scat, mu_estis, CrossSec_filename )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~Calculate the analytical results of Hua Xin-Min of 1992~~~~~~~~~~
    if( myid == np-1 )then
        E0 = mec2
        do i1 = 1, 8
            mu2 = phot%mu_esti(i1)
            do k1 = 0, vL_sc_up
                Es = k1 * phot%dy
                call phot%Hxm_analytical_results( E0, Es, mu2, g1Emu2 )
                Imu_Recv(i1, 1, k1) = g1Emu2
            enddo
        enddo 
        open(unit=10, file = Hxm_analyticResult, status="replace") 
        do j = 0, vL_sc_up 
            write(unit = 10, fmt = 200)Imu_Recv(1: 8, 1,  j)   
        enddo
        close(unit=10)
    endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Num_Photons = 0   
    sphot%mu_estimat = phot%mu_estimat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0   
        CALL phot%Generate_A_Photon( Emitter )    
        CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( )    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        CALL phot%Photon_Electron_Scattering( )  
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1 
            if( phot%scatter_times >= 2 )exit   
            CALL phot%Set_InI_Conditions_For_Next_Scattering( )    
            !write(*, *)'times2 = ', phot%scatter_times
            CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( )  
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            CALL phot%Photon_Electron_Scattering( )  
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
        If ( mod(Num_Photons, 50000)==0 .and. myid == np-1 ) then  
        write(unit = *, fmt = *)'***** The', Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = *)'***** My Duty Photon Number is: ', myid, mydutyphot 
        write(unit = *, fmt = *)'************************************************************'
        endif
    If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Enddo  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then  
          send_num = ( vL_sc_up + 1 ) * 7 * 8
          send_tag = 1  
          call MPI_SEND( phot%PolArrImu, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayI to Processor:', np-1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      else 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  recv_num = ( vL_sc_up + 1 ) * 7 * 8
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  call MPI_RECV( Imu_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                1, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives I data from Processor:', RECV_SOURCE  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  phot%PolArrImu = phot%PolArrImu + Imu_Recv   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              enddo
          endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!'  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          open(unit= 9, file='./spectrum/Imu0.dat', status="replace")   
          open(unit=10, file='./spectrum/Imu1.dat', status="replace")  
          open(unit=11, file='./spectrum/Imu2.dat', status="replace")  
          open(unit=12, file='./spectrum/Imu3.dat', status="replace") 
          open(unit=13, file='./spectrum/Imu6.dat', status="replace")    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          do j = 0, vL_sc_up
              write(unit = 9, fmt = 200)phot%PolArrImu(1: 8, 0,  j)  
              write(unit = 10, fmt = 200)phot%PolArrImu(1: 8, 1,  j)  
              write(unit = 11, fmt = 200)phot%PolArrImu(1: 8, 2,  j) 
              write(unit = 12, fmt = 200)phot%PolArrImu(1: 8, 3,  j) 
              write(unit = 13, fmt = 200)phot%PolArrImu(1: 8, 6,  j)  
          enddo 
          200 FORMAT(' ', ES15.6, '  ', ES15.6, '  ', ES15.6, '  ', ES15.6 ,'  ', &
                          ES15.6 ,'  ', ES15.6 ,'  ', ES15.6, '  ', ES15.6 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          close(unit=9)
          close(unit=10)
          close(unit=11)
          close(unit=12)
          close(unit=13)
      endif    
    call MPI_FINALIZE ( ierr ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    END SUBROUTINE mimick_of_ph_Slab_BoundReflc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    END MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    
