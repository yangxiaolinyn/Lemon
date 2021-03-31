!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_ThomScat_Emerge_IQ  
    implicit none
    real(mcp) :: tau 
    real(mcp) :: jIQUV(1: 4), alpIQUV(1: 4), rhoQUV(1: 3), IQUV0(1: 4), alpha
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

        jIQUV(1) = 2.1D0
        jIQUV(2) = 1.2D0
        jIQUV(3) = 0.D0
        jIQUV(4) = 0.D0

        alpIQUV(1) = 2.D0
        alpIQUV(2) = 1.5D0
        alpIQUV(3) = 0.D0
        alpIQUV(4) = 0.D0
 
        rhoQUV(1) = 0.D0
        rhoQUV(2) = 0.D0
        rhoQUV(3) = 0.D0

        IQUV0(1) = 0.D0
        IQUV0(2) = 0.D0
        IQUV0(3) = 0.D0
        IQUV0(4) = 0.D0 

        alpha = 3.D0

        MCResultsFile = './figures/MC_IQ.dat'
        Total_Phot_Num = 2.0D7
    else
        tau = 3.2D0

        jIQUV(1) = 0.D0
        jIQUV(2) = -10.71710 
        jIQUV(3) = 9.033D0 
        jIQUV(4) = 9.0587D0 

        alpIQUV(1) = 0.D0
        alpIQUV(2) = 0.D0
        alpIQUV(3) = 0.D0
        alpIQUV(4) = 0.D0
 
        rhoQUV(1) = 7.5D0
        rhoQUV(2) = 3.4D0
        rhoQUV(3) = 7.2D0

        IQUV0(1) = 0.D0
        IQUV0(2) = 0.D0
        IQUV0(3) = 0.D0
        IQUV0(4) = 0.D0 

        alpha = 60.D0

        MCResultsFile = './figures/MC_QUV.dat'
        Total_Phot_Num = 4.0D7
    endif
    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, tau, &
              jIQUV, alpIQUV, rhoQUV, IQUV0, alpha, &
              AnalyticResult_Fig27_left, AnalyticResult_Fig27_right, &
              MCResultsFile, Fig27_left_panel )   

    END PROGRAM main








 
