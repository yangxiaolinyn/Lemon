    program main   
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec
    real(mcp) :: E1_scat, E2_scat, y_obs1, y_obs2
    real(mcp) :: mu_esti(1: 4), sin_esti(1: 4)
    integer :: methods_cases, Num_mu_esti
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: Savefilename, CrossSec_filename, filenameH3Array, Te, Ta, E1, E2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    filenameH3Array = './data/H3Array.dat' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methods_cases = 5
    E1_scat = 1.D-7 ! In unit of MeV. E1_scat, E2_scat define the energy range 
    E2_scat = 1.D1  ! of scattering. The total scattering coefficient: Sigma_a(T_e, E_p) also
                    ! define on this interval, i.e., E_p \in (E1_scat, E2_scat). E_p is the
                    ! photon energy.
    y_obs1 = -1.D0  ! where \nu_obs = 10^{y_obs}, or y_obs = log10( \nu_obs )
    y_obs2 = 5.D0   ! where y_obs1 and y_obs2 determine the range of observational frequency. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods_cases == 1)then
!~~~~~~~Initial parameters for Fig.17 of our paper ~~~~~~~~~~~~~~~~~~~~~~~~~ 
        tau = 0.15D0            ! Optical depth of the atmosphere.
        T_bb = 1.D-5 * mec2     ! The source Blackbody's temperature, In unit of MeV
        T_elec = 1.7D0 * mec2   ! The electron's temperature. In unit of MeV 
    elseif(methods_cases == 2)then
!~~~~~~~Initial parameters for Fig.18 of our paper ~~~~~~~~~~~~~~~~~~~~~~~~~ 
        tau = 0.01D0 
        T_bb = 1.D-5 * mec2  
        T_elec = 1.59D0 * mec2 
    elseif(methods_cases == 3)then
!~~~~~~~Initial parameters for Fig.19 of our paper ~~~~~~~~~~~~~~~~~~~~~~~~~ 
        tau = 0.06D0 
        T_bb = 1.D-5 * mec2 
        T_elec = 0.188D0 * mec2 
    elseif(methods_cases == 4)then
!~~~~~~~Initial parameters for left panel of Fig.20 of our paper ~~~~~~~~~~~~~~~~~~~~~~~~~ 
        tau = 0.1D0 
        T_bb = 1.D-5 * mec2  
        T_elec = 0.5D0 * mec2 
    elseif(methods_cases == 5)then
!~~~~~~~Initial parameters for right panel of Fig.20 of our paper ~~~~~~~~~~~~~~~~~~~~~~~~~ 
        tau = 0.1D0 
        T_bb = 1.D-5 * mec2     
        T_elec = 1.D0 * mec2   
    endif
  

!~~~~~~~~Set the cosine and sine of observational angle~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( methods_cases == 1 .or. methods_cases == 2 .or. methods_cases == 3 )then
        Num_mu_esti = 3
        mu_esti(1) = 0.1D0 
        mu_esti(2) = 0.5D0 
        mu_esti(3) = 0.9D0
        mu_esti(4) = 0.D0
    else if( methods_cases == 4 .or. methods_cases == 5 )then
        Num_mu_esti = 4
        mu_esti(1) = 1.0D0 
        mu_esti(2) = 0.5D0 
        mu_esti(3) = 0.1D0
        mu_esti(4) = -0.5D0 
    endif
    sin_esti = dsqrt( one - mu_esti**2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    Total_Phot_Num = 1.D9
    write(Te, "(ES10.3)")T_elec
    write(Ta, "(ES10.3)")tau
    write(E1, "(ES10.3)")E1_scat
    write(E2, "(ES10.3)")E2_scat
    CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
               trim('tau=')//trim(adjustl(Ta))//trim('E1=')//&
               trim(adjustl(E1))//trim('E2')//trim(adjustl(E2))//trim('.dat') !
    Savefilename = trim('./spectrum/dataTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(Ta))//trim('.dat') !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ _DiffuseReflec
    !Total_Phot_Num = 1.2D8
    !call mimick_of_ph_FST_Slab_BoundReflc(Total_Phot_Num, tau, T_bb, T_elec, CrossSec_filename) 
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, &
                      T_elec, E1_scat, E2_scat, y_obs1, y_obs2, &
                      mu_esti, sin_esti, Num_mu_esti, &
                      CrossSec_filename, Savefilename )

    END PROGRAM main








 
