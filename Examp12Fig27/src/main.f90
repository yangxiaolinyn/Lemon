!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_ThomScat_Emerge_IQ  
    implicit none
    real(mcp) :: tau 
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: AnalyticResult_Fig27_left, AnalyticResult_Fig27_right, &
                    MCResultsFile
    logical :: Fig27_left_panel

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    call InitRandom()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    AnalyticResult_Fig27_left = './figures/I_Q.dat'
    AnalyticResult_Fig27_right = './figures/QUV.dat'
    Fig27_left_panel = .true.
    Fig27_left_panel = .false.
    if( Fig27_left_panel )then 
        tau = 10.D0
        MCResultsFile = './figures/MC_IQ.dat'
        Total_Phot_Num = 2.0D7
    else
        tau = 3.2D0
        MCResultsFile = './figures/MC_QUV.dat'
        Total_Phot_Num = 4.0D7
    endif
    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau, &
              AnalyticResult_Fig27_left, AnalyticResult_Fig27_right, &
              MCResultsFile, Fig27_left_panel )   

    END PROGRAM main








 
