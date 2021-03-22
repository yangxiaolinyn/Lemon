!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_ThomScat_Emerge_IQ  
    implicit none
    real(mcp) :: theta_obs, tau, epsilon_crit, mu0, phi0, I0, Q0, U0, V0
    real(mcp) :: mu_esti(1: 4), phi_esti(1: 1)
    integer :: Num_mu_esti, Num_phi_esti
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, mu0s, MCResultsFNphi, ChandraFNphi
  
    tau = 0.2D0 ! The total optical depth of the plane-parallel atmosphere from the up 
                ! boundary to the bottom boundary.
    epsilon_crit = 1.D-4   ! The threshold value to terminate the scattering sequence.
    mu0 = 0.8D0  ! Notice that mu0 must be negative !!!
    phi0 = zero    ! the phi angle of incident beam.
    Total_Phot_Num = 2.D5
    !~~~~~~~~~~~~The Stokes Parameters of incident beam!~~~~~~~~~~~~~~~~~~~~~~ 
    I0 = one
    Q0 = zero
    U0 = zero
    V0 = one
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Num_mu_esti = 4 ! If you change the value of Num_mu_esti, you should modify the 
                    ! array mu_esti and set its initial values as well. Here we set to 
                    ! it to be 4.
    !mu_esti(1) = 5.D0 / 100.D0
    !mu_esti(2) = 30.D0 / 100.D0
    !mu_esti(3) = 60.D0 / 100.D0
    !mu_esti(4) = 85.D0 / 100.D0 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Num_phi_esti = 1 ! If you change the value of Num_mu_esti, you should modify the 
                    ! array mu_esti and set its initial values as well. Here we set to 
                    ! it to be 4.
    phi_esti(1) = zero
    !phi_esti(2) = pi / two
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(mu0s, "(f8.4)")mu0
    MCResultsFNphi = trim('./spectrum/IQUVphi_mu0=')//trim(adjustl(mu0s))//trim('.dat')
    ChandraFNphi = trim('./spectrum/ChandraIQUV_phi_mu0=')//trim(adjustl(mu0s))//trim('.dat') 

    CALL Mimick_Photon_Diffuse_Transfer( Total_Phot_Num, tau, epsilon_crit, mu0, &
                  phi0, I0, Q0, U0, V0, phi_esti, Num_phi_esti, &
                  MCResultsFNphi, ChandraFNphi )      

    END PROGRAM main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








 
