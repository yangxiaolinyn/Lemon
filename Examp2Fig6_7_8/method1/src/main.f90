!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Statistial_Method_Of_Finity_zone
    implicit none
    real(mcp) :: tau
    integer(kind = 8) :: Total_Phot_Num 
    real(mcp) :: the_obs, phi_obs, y1, y2, T_e, ne, Rout, T_s
    character*80 :: SemiAnalyResults, MCResults, tobs, Te, nes, &
                 CrossSec_filename, taus
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call InitRandom( )
 
    tau = 1.D-1
    T_e = 4.D0 * mec2
    T_s = 1.D-8 * mec2
    ne = 1.D17

    y1 = 11.D0
    y2 = 22.D0  ! \nu = 10^y
  
    write(Te, "(f8.3)")T_e
    write(taus, "(f8.4)")tau
    CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(taus))//trim('.txt') ! the total scattering 
    MCResults = trim('./spectrum/MCRTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(taus))//trim('.txt')


    Total_Phot_Num = 1D7
    call mimick_of_photon_with_finity_zone( Total_Phot_Num, tau, T_e, &
                                T_s, ne, y1, y2, CrossSec_filename, MCResults )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END PROGRAM main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








 
