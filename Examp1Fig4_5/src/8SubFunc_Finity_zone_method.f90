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
      SUBROUTINE Emitter_A_Photon1( Emitter )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter 
      REAL(mcp) :: j_nu, D, L, R_sphere, r_xy, r_max, dr, dtheta, theta_xy
      REAL(mcp) :: dnu, nu_low, nu_up, v_Lv(0: 500)
      REAL(mcp) :: alpha, nu, costheta, sintheta
      integer(kind=8) :: i, j, k, Nr, Nt, Nn

      v_Lv = zero
      R_sphere = one
      Nr = 300
      Nt = 300
      Nn = 500
      D = one
      L = 1.D4
      r_max = D * R_sphere / dsqrt( L**2 - R_sphere**2 )
      dr = r_max / Nr
      dtheta = twopi / Nt
      nu_low = 8.D0
      nu_up = 15.D0
      dnu = ( nu_up - nu_low ) / Nn
      DO i = 0, Nr
          r_xy = dr * i
          write(*, *)'sss====', i, r_xy
          DO j = 0, Nt
              Do k = 0, Nn
                  nu = 10**( k * dnu + nu_low )
                  theta_xy = dtheta * j
                  costheta = r_xy * dsin( theta_xy ) / dsqrt( r_xy**2 + D**2 )
                  sintheta = dsqrt( one - costheta**2 ) 
                  v_Lv( k ) = v_Lv( k ) + nu * Emitter%j_theta_nu_emissity&
                              ( nu, costheta, sintheta ) * two * &
                              dsqrt( R_sphere**2 - ( L*r_xy )**2 / ( D**2 + r_xy**2 ) )
                  !write(*, *)'sss====', i, j, j_nu, nu, costheta, sintheta
              enddo
          ENDDO
      ENDDO
      open(unit=16, file='./image/v_Lv1.txt', status="replace") 
      do k = 0,Nn
          write(unit = 16, fmt = *)v_Lv( k )
      enddo 
      RETURN
      END SUBROUTINE Emitter_A_Photon1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Emitter_A_Photon2( Emitter, phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter
      TYPE(Photon), INTENT(INOUT) :: Phot 
      REAL(mcp) :: j_nu, D, L, R_sphere, r_xy, r_max, dr, dtheta, theta_xy
      REAL(mcp) :: dnu, nu_low, nu_up, v_Lv(0: 500), Length_0
      REAL(mcp) :: alpha, nu, costheta, sintheta, alpha_nu, T_e
      integer(kind=8) :: i, j, k, Nr, Nt, Nn

      v_Lv = zero
      R_sphere = one
      Nr = 300
      Nt = 300
      Nn = 500
      D = one
      L = 1.D10
      T_e = 100.D0 * mec2
      r_max = D * R_sphere / dsqrt( L**2 - R_sphere**2 )
      dr = r_max / Nr
      dtheta = twopi / Nt
      nu_low = 8.D0
      nu_up = 15.D0
      dnu = ( nu_up - nu_low ) / Nn
      DO i = 0, Nr
          r_xy = dr * i
          write(*, *)'sss====', i, r_xy
          DO j = 0, Nt
              Do k = 0, Nn
                  nu = 10**( k * dnu + nu_low )
                  theta_xy = dtheta * j
                  costheta = r_xy * dsin( theta_xy ) / dsqrt( r_xy**2 + D**2 )
                  sintheta = dsqrt( one - costheta**2 ) 
                  Length_0 = two * dsqrt( R_sphere**2 - ( L*r_xy )**2 / ( D**2 + r_xy**2 ) ) 
                  alpha_nu = Emitter%j_theta_nu_emissity( nu, costheta, sintheta ) / &
                             ( two*planck_h*nu**3/Cv**2 / ( dexp(h_ev*nu*1.D-6/51.1D0) - one ) ) 
                  if( Length_0 * alpha_nu <= 1.D-10 )then
                      v_Lv( k ) = v_Lv( k ) + nu * Emitter%j_theta_nu_emissity &
                              ( nu, costheta, sintheta ) * Length_0  
                  else
                      v_Lv( k ) = v_Lv( k ) + nu * Emitter%j_theta_nu_emissity &
                              ( nu, costheta, sintheta ) / alpha_nu * &
                              ( one - dexp( - Length_0 * alpha_nu ) )  
                  endif 
                  !write(*, *)'sss====', Emitter%n_e * phot%sigma_fn( phot%T_e, nu*h_ev*1.D-6 ),alpha_nu, &
                  !       alpha_nu - Emitter%n_e * phot%sigma_fn( phot%T_e, nu*h_ev*1.D-6 )
                  !       Emitter%j_theta_nu_emissity( one, nu, costheta, sintheta )
              enddo
          ENDDO
      ENDDO
      open(unit=16, file='./spectrum/vLv_Integral2.txt', status="replace") 
      do k = 0,Nn
          write(unit = 16, fmt = *)v_Lv( k )
      enddo 
      RETURN
      END SUBROUTINE Emitter_A_Photon2

  
  

!**************************************************************************************
    SUBROUTINE mimick_of_photon_with_finity_zone( Total_Phot_Num )
!************************************************************************************** 
    implicit none 
    real(mcp) :: T_e, E, dt, t = 0.D0, p_length, E_low, E_up, Tb = 2.D-3
    real(mcp), parameter :: n_e1 = 1.6D17
    integer(kind = 8) :: scatter_Times = 0
    integer(kind = 8) :: Num_Photons
    integer(kind = 8), intent(in) :: Total_Phot_Num 
    integer(kind = 4), parameter :: n_sample = 1000000
    real(mcp) :: R_out, R_in, rea, img, E_max = 1.D15*h_ev*1.D-6
    logical :: At_outer_Shell, At_inner_Shell
    type(Photon_Emitter) :: Emitter
    type(Photon) :: phot 
    integer :: send_num, recv_num, send_tag, RECV_SOURCE, status(MPI_STATUS_SIZE) 
    integer :: effect_Photon_Number_Sent 
    integer :: effect_Photon_Number_Recv, cases
    real(mcp), dimension(0:400, 0:400) :: disk_image_Sent, disk_image_Recv
    real(mcp), dimension(1:500, 1:1000) :: ET_array_Recv, ET_array_Sent
    real(mcp) :: v_L_v_Sent(0:500), v_L_v_Recv(0:500), temp_ET
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
        !write(unit = *, fmt = *)'ggg==', Num_Photons 
        CALL phot%Emitter_A_Photon( ) 
        !CALL Transmit_Data_And_Parameters_From_Emitter2Photon( Emitter, Phot )  
        CALL phot%Calc_Phot_Informations_At_Observor_2zones( ) 

112     continue
        If ( mod(Num_photons,5000000)==0 .and. myid==np-1 ) then
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  recv_num = 500
                  call MPI_RECV(v_L_v_Recv, recv_num, MPI_DOUBLE_PRECISION, RECV_SOURCE,&
                                2, MPI_COMM_WORLD, status, ierr) 
                  !write(*,*)'sdf500==', v_L_v_Recv
                  write(*,*)'master Processor ',myid,' has receved v_L_v Processor:', RECV_SOURCE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                  phot%v_L_v_i = phot%v_L_v_i + v_L_v_Recv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              enddo
          endif
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          write(*,*)'There are', phot%effect_number, 'of total', Total_Phot_Num, 'photons',&
                        'arrive at the plate of observer!!'
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
          open(unit=17, file='./spectrum/vLv_MC_2.txt', status="replace") 
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          do j = 0, 499
              write(unit = 17, fmt = *)phot%v_L_v_i(j)
          enddo
!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
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














 
