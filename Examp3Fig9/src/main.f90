!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    use Method_Of_FLST_ThomScat_Emerge_IQ 
    implicit none
    real(mcp) :: tau
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: MCResultsFileNameI, MCResultsFileNameQ, taus
 
    call InitRandom() 
!~~~~~~~~~~Optical depth~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = 10.D0 
!~~~~~~~~~~Number of photons need be traced!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Total_Phot_Num = 3.D6
 
    write(taus, "(f8.4)")tau
    MCResultsFileNameI = trim('./spectrum/EIQ_I_tau=')//trim(adjustl(taus))//trim('.dat')
    MCResultsFileNameQ = trim('./spectrum/EIQ_Q_tau=')//trim(adjustl(taus))//trim('.dat')
    call mimick_of_ph_finity_zone_Emerge_IQ( Total_Phot_Num, &
                 tau, MCResultsFileNameI, MCResultsFileNameQ )   

    END PROGRAM main








 
