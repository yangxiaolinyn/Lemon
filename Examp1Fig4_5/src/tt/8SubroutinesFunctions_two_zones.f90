!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE SubroutineFunction
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE constants
      USE RandUtils
      USE PhotonEmitter
      USE Photons
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
      !Emitter%ComFrameU4_BL(1) = one
      !Emitter%ComFrameU4_BL(2) = zero
      !Emitter%ComFrameU4_BL(3) = zero
      !Emitter%ComFrameU4_BL(4) = one / ( Emitter%aspin + Emitter%r**1.5D0 )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Here we begin to set the several transformation Matrices, by which we can complete     ||
!| Transformation of the Photon's four velocity, Phot4k_CtrCF(1:4), from Comoving Frame   ||
!| to the BL coordinate, i.e. Phot4k_CtrBL and Phot4k_CovBL.                              ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Now we generate a Photon in the CF(Comoving Frame), i.e. we get the photon's four      ||
!| Velocity: Phot4k_CtrCF, then transform it into the LNRF reference: Phot4k_LNRF, into   ||
!| BL coordinate: Phot4k_CtrBL, Phot4k_CovBL, they are the Contravariant four Velocity    ||
!| and the Covariant four Velocity, respectively. With Phot4k_LNRF, we can calculate the  ||
!| constants of the geodesic motion: lambda and q.                                        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !Do
          CALL Emitter%get_Phot4k_CtrCF_CovCF()
          !write(*,*)'ssss=',Vector3D_Inner_Product( Emitter%Phot4k_CtrCF(2:4), Emitter%Vector_of_position )
          !If ( Vector3D_Inner_Product( Emitter%Phot4k_CtrCF(2:4), Emitter%Vector_of_position ) < zero ) exit
      !Enddo 
  
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
      If( phot%mymethod )then
          phot%p_scattering = phot%Get_scatter_distance2( T_e ) 
      else
          phot%p_scattering = phot%Get_scatter_distance( T_e ) 
      endif 
      !write(*,*)'2222==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !phot%direct_escaped = .True.
      cases = 1
      Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
      phot%w_ini = phot%w_ini * phot%NormalA
      !phot%direct_escaped = .false.
      !write(*,*)'3333=', phot%E_ini, phot%Phot4k_CovCF_ini(1)
      !write(*,*)'111=', phot%time_travel 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      If ( phot%At_outer_Shell ) RETURN
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
      If( phot%mymethod )then
          phot%p_scattering = phot%Get_scatter_distance2( T_e ) 
      else
          phot%p_scattering = phot%Get_scatter_distance( T_e ) 
      endif 
      !write(*,*)'7777==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !phot%direct_escaped = .True.
      cases = 1
      Call phot%Calc_Phot_Informations_At_Observor_2zones( cases )
      phot%w_ini = phot%w_ini * phot%NormalA
      !phot%direct_escaped = .false. 
      !write(*,*)'8888==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
          !write(*,*)'666=', phot%time_travel 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      If ( phot%At_outer_Shell ) RETURN 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      If ( phot%At_inner_Shell .and. phot%r_ini > phot%R_in ) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          phot%time_travel = phot%time_travel + phot%p_boundary1 / Cv
          !write(*,*)'777=', phot%time_travel,  phot%p_boundary1/ CV
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          phot%r_ini      = phot%R_in 
          phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                      phot%Vector_of_Momentum_ini, phot%p_scattering )
          phot%Vector_of_position_ini = phot%Vector_of_position
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          If( phot%mymethod )then
              phot%p_scattering = phot%Get_scatter_distance2( T_e ) 
          else
              phot%p_scattering = phot%Get_scatter_distance( T_e ) 
          endif 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          If ( phot%At_inner_Shell ) then
              phot%time_travel = phot%time_travel + phot%p_boundary1 / Cv
              !write(*,*)'888=', phot%time_travel,  phot%p_boundary1/ CV
              phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                      phot%Vector_of_Momentum_ini, phot%p_scattering )
              phot%Vector_of_position_ini = phot%Vector_of_position
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              If( phot%mymethod )then
                  phot%p_scattering = phot%Get_scatter_distance2( T_e ) 
              else
                  phot%p_scattering = phot%Get_scatter_distance( T_e ) 
              endif 
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              If ( phot%At_outer_Shell ) RETURN
          endif
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      If ( phot%At_inner_Shell .and. phot%r_ini < phot%R_in ) then
          phot%time_travel = phot%time_travel + phot%p_boundary1 / Cv
          !write(*,*)'999=', phot%time_travel,  phot%p_boundary1/ CV
          phot%r_ini      = phot%R_in  
          phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                      phot%Vector_of_Momentum_ini, phot%p_scattering )
          phot%Vector_of_position_ini = phot%Vector_of_position
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          If( phot%mymethod )then
              phot%p_scattering = phot%Get_scatter_distance2( T_e ) 
          else
              phot%p_scattering = phot%Get_scatter_distance( T_e ) 
          endif 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          If ( phot%At_outer_Shell ) RETURN
      endif
      !write(*,*)'9999==', phot%E_ini, phot%Phot4k_CovCF_ini(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      phot%time_travel = phot%time_travel + phot%p_scattering / Cv
          !write(*,*)'101=', phot%time_travel,  phot%p_scattering/ CV
      phot%r_p = phot%r_p2( phot%Vector_of_position_ini, &
                 phot%Vector_of_Momentum_ini, phot%p_scattering ) ! After this call
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|        ||
!|        ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END MODULE SubroutineFunction
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|        ||
!|        ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~














 
