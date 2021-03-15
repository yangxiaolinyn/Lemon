!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    use Method_Of_FLST_ThomScat_Emerge_IQ 
    implicit none
    real(mcp) :: tau
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array
 
    call InitRandom() 
!~~~~~~~~~~Optical depth~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = 10.D0 
!~~~~~~~~~~Number of photons need be traced!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Total_Phot_Num = 3.D6

    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau )   

    END PROGRAM main








 
