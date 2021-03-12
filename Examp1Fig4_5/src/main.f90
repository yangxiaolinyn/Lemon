    program main  
    use Statistial_Method_Of_Finity_zone
    implicit none  
    integer(kind = 8) :: Total_Phot_Num 
    real(mcp) :: the_obs, phi_obs, y1, y2, Te, Rout
    character*80 :: SemiAnalyResults, tobs, pobs
 
    the_obs = 90.D0
    phi_obs = zero
    y1 = 8.D0
    y2 = 15.D0
    Te = 100.D0 * mec2
    Rout = one

        write(tobs, "(f8.3)")the_obs
        write(pobs, "(f8.5)")phi_obs
        SemiAnalyResults = trim('./spectrum/tobs=')//trim(adjustl(tobs))//&
                      trim('pobs=')//trim(adjustl(pobs))//trim('.txt') 
  
    Total_Phot_Num = 1D10
    call mimick_of_photon_with_finity_zone( Total_Phot_Num, &
                     the_obs, phi_obs, y1, y2, Te, Rout, SemiAnalyResults )  

    END PROGRAM main








 
