!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_ThomScat_Emerge_IQ  
    implicit none
    real(mcp) :: theta_obs, tau, mu0, phi0 
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, mu0s, MCResultsFNphi0180, MCResultsFNphi90
  
    tau = 200.D0   ! tau is a critical value to terminate the scattering sequence, i.e., 
                   ! As z_tau > tau, the sequence is terminated!
    mu0 = - 0.2D0  ! Notice that mu0 must be negative !!!
    phi0 = zero
    Total_Phot_Num = 2.D5

    write(mu0s, "(f8.4)")mu0
    MCResultsFNphi0180 = trim('./spectrum/IQUVphi=0180_mu0=')//trim(adjustl(mu0s))//trim('.dat')
    MCResultsFNphi90 = trim('./spectrum/IQUVphi=90_mu0=')//trim(adjustl(mu0s))//trim('.dat') 

    CALL Mimick_Photon_Diffuse_Transfer( Total_Phot_Num, tau, mu0, &
                  phi0, MCResultsFNphi0180, MCResultsFNphi90 )      

    END PROGRAM main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








 
