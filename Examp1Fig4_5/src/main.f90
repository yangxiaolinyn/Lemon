    program main  
    use Statistial_Method_Of_Finity_zone
    implicit none  
    integer(kind = 8) :: Total_Phot_Num     
    type(Photon) :: phot 
    real(mcp) :: the_obs, phi_obs, y1, y2, Te, Rout
 
    the_obs = 90.D0
    phi_obs = zero
    y1 = 8.D0
    y2 = 15.D0
    Te = 100.D0 * mec2
    Rout = one

    CALL phot%Set_initial_parameter_values( the_obs, phi_obs, y1, y2, Te, Rout ) 
    CALL Semi_Analytical_Calculations2( phot )

    Total_Phot_Num = 1D10
    call mimick_of_photon_with_finity_zone( Total_Phot_Num, phot )  

    END PROGRAM main








 
