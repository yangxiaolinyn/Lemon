!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE Method_Of_FLST_ThomScat_Emerge_IQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      USE RandUtils 
      USE PhotonModule
      USE MPI
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Calculate_The_Diffuse_Reflection_of_Chandra( Phot, I0, Q0, U0, V0, ChandraFNphi )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE  
      TYPE(Photons), INTENT(INOUT) :: Phot 
      REAL(mcp), INTENT(IN) :: I0, Q0, U0, V0 
      character*80, intent(in) :: ChandraFNphi
      INTEGER :: i, j
      REAL(mcp) :: mu, mu0, phi, phi0, I_r, I_l, U, V, IQUV_in(1: 3), IQU(1: 3), &
                   Irluv_in(1: 3), Irluv(1: 3)
      real(mcp), dimension(1:3, 1:3) :: QS_Matrix, QS_Matrix_rluv
  
          IQUV_in(1) = I0
          IQUV_in(2) = Q0
          IQUV_in(3) = U0 

          mu0 = dabs( Phot%cos_Theta0 )
          phi0 = zero 
          200 FORMAT(' ', 4F20.13) 
          open(unit=9, file = ChandraFNphi, status="replace")  

          do i = 1, Phot%N_mu_esti
              mu = Phot%mu_estimat(i)
              do j = 0, 200 
                  phi = twopi / 200.D0 * j
                  CALL Phot%Diffuse_Reflection_Ir_Il_U_Setting( mu, mu0, phi, phi0, I_r, I_l, U, V ) 

                  CALL Phot%Diffuse_Reflection_QS_Matrix_of_IQUV_Setting( mu, mu0, phi, phi0, QS_Matrix ) 
                  CALL Matrix_Multiplication33X31_Sub( QS_Matrix, IQUV_in, IQU )
                  write(unit = 9, fmt = 200)IQU(1), IQU(2), IQU(3), V
              enddo
          enddo

          close(unit=9)  
          !stop
 
      RETURN
      END SUBROUTINE Calculate_The_Diffuse_Reflection_of_Chandra


!**************************************************************************************
    SUBROUTINE Mimick_Photon_Diffuse_Transfer( Total_Phot_Num, tau, epsilon_crit, mu0, phi0, &
           I0, Q0, U0, V0, mu_esti, Num_mu_esti, MCResultsFNphi, ChandraFNphi )
!************************************************************************************** 
    implicit none
    real(mcp), intent(inout) :: tau, epsilon_crit, mu0, phi0, I0, Q0, U0, V0, mu_esti(1: Num_mu_esti)
    integer, intent(inout) :: Num_mu_esti
    character*80, intent(in) :: MCResultsFNphi, ChandraFNphi
    real(mcp) :: E, E_low, E_up  
    integer(kind = 8) :: Num_Photons
    integer(kind = 8), intent(in) :: Total_Phot_Num
    type(Photons) :: phot
    type(ScatPhoton) :: sphot
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE)    
    real(mcp) :: IQUV_Recv(1: 4, 1: Num_mu_esti, 0: Num_Phi)
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
    call InitRandom( myid )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL phot%Set_initial_parameter_values( tau, mu0, phi0, I0, &
                               Q0, U0, V0, mu_esti, Num_mu_esti ) 
!~~~~~~Calculate the semi-analytical results of Chandrasekha~~~~~~~~~~~~~~~~~~~~~
    If( myid == np - 1 )then
        call Calculate_The_Diffuse_Reflection_of_Chandra( Phot, &
                                   I0, Q0, U0, V0, ChandraFNphi )
    endif
    Num_Photons = 0   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        Num_Photons = Num_Photons + 1 
        phot%scatter_times = 0   
        CALL phot%Emitter_A_Photon( ) 
        CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( ) 
        CALL phot%Photon_Electron_Scattering( )  
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Scattering_loop: Do
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            phot%scatter_times = phot%scatter_times + 1  
            !write(unit = *, fmt = *)'******* ', phot%scatter_times,  &
            !         phot%Psi_I, phot%z_tau
            CALL phot%Determine_P_Of_Scatt_Site_And_Quantities_At_p( ) 
            if( phot%Psi_I <= epsilon_crit .or. phot%z_tau > tau )exit
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
            !if( phot%Optical_Depth_scatter >= 1.D2 )exit
            !if( phot%I_IQ / phot%w_ini0 <= 1.D-6 )exit
            !if( phot%scatter_times > 100 )exit
            !if( phot%z_tau > 50.D0 .or. dabs(phot%I_IQ) <= 1.D-5) exit 
            !if(phot%scatter_times > 20)stop
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            CALL phot%Photon_Electron_Scattering( )
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END DO Scattering_loop
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
        If ( mod(Num_Photons, 500) == 0 .and. myid == np-1 ) then  
            write(unit = *, fmt = *)'***** The', Num_Photons,'th Photons have been scattered', &
                                  phot%scatter_times, &
                         'times and Escapted from the region !!!!!!!'      
            write(unit = *, fmt = *)'***** My Duty Photon Number is: ',myid, mydutyphot 
            write(unit = *, fmt = *)'*****  tau === ', phot%z_tau, phot%Psi_I
            write(unit = *, fmt = *)'*************************************************************************'
        endif
        If( Num_Photons > mydutyphot )EXIT 
    Enddo  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( myid /= np-1 ) then  
          send_num = (Num_Phi + 1) * Num_mu_esti * 4
          send_tag = 1  
          call MPI_SEND( phot%PolArrayIQUV_Phi, send_num, MPI_DOUBLE_PRECISION, np-1, &
                        send_tag, MPI_COMM_WORLD, ierr)
          write(*, *)'Processor ', myid, ' has send PolarArrayI to Processor:', np-1   
      else 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if (np-2 >= 0) then
              do RECV_SOURCE = 0, np-2  

                  recv_num = (Num_Phi + 1) * Num_mu_esti * 4 
                  call MPI_RECV( IQUV_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE, &
                                1, MPI_COMM_WORLD, status, ierr) 
                  write(*,*)'master Processor ', myid,' Receives I data from Processor:', RECV_SOURCE    
                  phot%PolArrayIQUV_Phi = phot%PolArrayIQUV_Phi + IQUV_Recv   

              enddo
          endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          open(unit=9, file = MCResultsFNphi, status="replace")  
          do i = 1, phot%N_mu_esti
              do j = 0, Num_phi
                  write(unit = 9, fmt = 100)phot%PolArrayIQUV_Phi(1: 4, i, j)
              enddo
          enddo
          100 FORMAT(' ', '   ', 1F22.12, '   ', 1F22.12, '   ', 1F22.12, '   ', 1F22.12)
   
          close(unit=9)      
      endif   

      deallocate( phot%PolArrayIQUV_Phi, stat = ierr )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    call MPI_FINALIZE ( ierr ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END SUBROUTINE Mimick_Photon_Diffuse_Transfer
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END MODULE Method_Of_FLST_ThomScat_Emerge_IQ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

    
