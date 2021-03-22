!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_ThomScat_Emerge_IQ 
    !use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: theta_obs, tau
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    theta_obs = 1.D0  
    tau = 400.D0 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Total_Phot_Num = 2.D5
    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau )  
    !call mimick_of_ph_Emerge_DiffuseReflec( Total_Phot_Num, tau )   

    END PROGRAM main








 
