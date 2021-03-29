!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      USE RandUtils 
      USE Photons_FlatSP
      USE MPI
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS  
!**************************************************************************************
    SUBROUTINE mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, T_elec, &
                      CrossSec_filename, Spentrum_filename, S_in, alp, gama1, gama2, vy1, vy2, E_ini )
!************************************************************************************** 
    implicit none 
    real(mcp), intent(in) :: T_elec, S_in(1: 4), alp, gama1, gama2, vy1, vy2, E_ini
    integer(kind = 8), intent(in) :: Total_Phot_Num
    character*80, intent(inout) :: CrossSec_filename, Spentrum_filename
    real(mcp) :: E, E_low, E_up  
    integer(kind = 8) :: Num_Photons 
    type(Photon_Emitter) :: Emitter
    type(Photon_FlatSP) :: phot
    type(ScatPhoton_KN) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE), send_num2, recv_num2
    real(mcp) :: IQUV_Recv(1: 4, 0: 6, 1: Num_mu, 0: vL_sc_up)
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
        
        write(*, *)'S_im = ', S_in
    else
        mydutyphot = Total_Phot_Num / np
    endif 

        write(*, fmt=*)'ffsdfa='
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL phot%Set_Initial_Values_For_Photon_Parameters( T_elec, Emitter, CrossSec_filename, &
                      S_in, alp, gama1, gama2, vy1, vy2, E_ini )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Num_Photons = 0   
    sphot%mu_estimat = phot%mu_estimat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0  

     !write(*,*)'**********************************************************' 
        CALL phot%Generate_A_Photon( Emitter ) 
        !CALL Emitter_A_Photon( Emitter )
        !CALL Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )   
        !CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
        CALL phot%FIRST_SCATTERING_OF_PHOT_ELCE( sphot )
        !CALL phot%Photon_Electron_Scattering( sphot )
        !write(*,fmt=*)'ssff222===', Num_Photons

    if( myid == np-1 ) then
        If ( mod(Num_Photons, 5000000)==0 ) then 
        write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = "(' ', A14, I10, A40, I10, A40)")'***** The', Num_Photons,&
              '  th Photons have been scattered', phot%scatter_times, &
                         '  times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = 120)'***** My Duty Photon Number is: ', myid, mydutyphot 
        write(unit = *, fmt = *)'*************************************************************************'
        endif
    endif

        If( Num_Photons > mydutyphot )EXIT
        cycle
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1 
            !CALL Set_InI_Conditions_For_Next_Scattering( phot, sphot )   
            CALL phot%Set_InI_Conditions_For_Next_Scattering( sphot )   
            !CALL Determine_Next_Scattering_Site( phot, sphot ) 
            !write(*,fmt="(' ', A8, 3ES15.5, I5)")'ssff222===', phot%Phot4k_CtrCF_ini(1), &
            !        mec2 * 10.D0, phot%w_ini, phot%scatter_times
            !write(*,fmt=*)'ssff222===', phot%Phot4k_CtrCF_ini(1), mec2 * 10.D0, phot%w_ini, phot%scatter_times
            CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p() 
            !if( phot%Phot4k_CtrCF_ini(1) > mec2 * 1.D1)exit 
            if( phot%w_ini / phot%w_ini0 <= 1.D-40)exit
            !write(*,fmt=*)'ssff222===', T_bb * 1.1D4,  mec2 * 8.D-1
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
            !if( phot%scatter_times >= 20 )exit  
            !if(   phot%Psi_I <= 1.D-1 )exit 
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CALL phot%Photon_Electron_Scattering( sphot ) 
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !write(*,*)'************************************* *****',phot%scatter_times, phot%medium_case
  
        If ( mod(Num_Photons, 20000)==0 ) then 
        write(unit = *, fmt = *)'*************************************************************************' 
        write(unit = *, fmt = "(' ', A14, I5, A50, I5, A50)")'***** The', Num_Photons,&
              'th Photons have been scattered', phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
        write(unit = *, fmt = 120)'***** My Duty Photon Number is: ', myid, mydutyphot 
        write(unit = *, fmt = *)'*************************************************************************'
        endif
    If( Num_Photons > mydutyphot )EXIT
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Enddo  
    120 format(' ', A35, 2I10)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          send_num  = ( vL_sc_up + 1 ) * 4 * 7 * Num_mu
          send_tag = 1  
          call MPI_SEND( phot%PolArrIQUV, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayI to Processor:', np-1
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      else 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  recv_num = ( vL_sc_up + 1 ) * 4 * 7 * Num_mu
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call MPI_RECV( IQUV_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                1, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives I data from Processor:', RECV_SOURCE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  phot%PolArrIQUV = phot%PolArrIQUV + IQUV_Recv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              enddo
          endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          open(unit=13, file = Spentrum_filename, status="replace") 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          do i = 1, Num_mu
              do j = 0, vL_sc_up 
                  write(unit = 13, fmt = 300)phot%PolArrIQUV(1: 4, 6, i, j)! * 1.D60
              enddo
          enddo
          100 FORMAT(' ', '    ', ES15.6, '    ', ES15.6, '    ', ES15.6, '    ', ES15.6, '    ', ES15.6)
          200 FORMAT(' ', ES15.6, '   ', ES15.6, '   ', ES15.6, '   ', ES15.6 ,'   ', &
                          ES15.6 ,'   ', ES15.6 ,'   ', ES15.6 )
          300 FORMAT(' ', '   ', ES15.6, '   ', ES15.6, '   ', ES15.6, '   ', ES15.6)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          close(unit=13)   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    call MPI_FINALIZE ( ierr ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE mimick_of_ph_Slab_BoundReflc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END MODULE Method_Of_FLST_DiffuseReflec
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    
