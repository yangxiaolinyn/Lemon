    program main  
    use Method_Of_FLST_ThomScat_Emerge_IQ 
    implicit none
    real(mcp) :: theta_obs, tau
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    !theta_obs = 1.D0
    call InitRandom()
    !methods_cases = 7
    !tau = 8.D9
    tau = 10.D0
    !write(*,*)'m=', 1.D11, 10.D11
    !CrossSec_filename = './data/SigmaArrayTe=4mec2.dat'
    !CrossSec_filename = './data/SigmaArrayTe=4mec2_Seq.dat'
    !filenameH3Array = './data/H3Array.dat'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    Total_Phot_Num = 3.D6
    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau )   

    END PROGRAM main








 
