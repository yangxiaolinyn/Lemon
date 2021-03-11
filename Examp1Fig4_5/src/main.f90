    program main 
    !use Statistial_Method_Of_Infinity_zone 
    use Statistial_Method_Of_Finity_zone
    implicit none
    real(mcp) :: a_spin
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    TYPE(Photon_Emitter) :: Emitter 
    TYPE(Photon) :: Phot 
 
    call InitRandom()
    methods_cases = 1
    if( methods_cases == 1 )then
        Total_Phot_Num = 1D10
        call mimick_of_photon_with_finity_zone( Total_Phot_Num )
    else
        !Total_Phot_Num = 2D8
        !call mimick_of_photon_with_infinity_zone(a_spin, Total_Phot_Num)
    endif

    END PROGRAM main








 
