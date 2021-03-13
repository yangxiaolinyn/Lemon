    program main
    use Statistial_Method_Of_Finity_zone
    implicit none
    real(mcp) :: tau
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
 
    tau = 1.D-4 

    Total_Phot_Num = 4D8
    call mimick_of_photon_with_finity_zone( Total_Phot_Num, tau)  

    END PROGRAM main








 
