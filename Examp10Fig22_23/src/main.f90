    program main  
    !use Method_Of_FLST_ThomScat_Emerge_IQ 
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec
    real(mcp) :: E1_scat, E2_scat, y_obs1, y_obs2
    real(mcp) :: mu_esti(1: 4), sin_esti(1: 4), Terminate_Tolerence
    integer :: methods_cases, TerminateTime, Num_mu_esti
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: Savefilename, CrossSec_filename, &
                   Spentrum_filename, filenameH3Array, Te, Ta, E1, E2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom() 
    filenameH3Array = './data/H3Array.dat'
    !CrossSec_filename = './data/SigmaArrayT_e=56kev.txt' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methods_cases = 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    E1_scat = 1.D-7 ! In unit of MeV. E1_scat, E2_scat define the energy range 
    E2_scat = 1.D1  ! of scattering. The total scattering coefficient: Sigma_a(T_e, E_p) also
                    ! define on this interval, i.e., E_p \in (E1_scat, E2_scat). E_p is the
                    ! photon energy.
    y_obs1 = -1.D0  ! where \nu_obs = 10^{y_obs}, or y_obs = log10( \nu_obs )
    y_obs2 = 5.D0   ! where y_obs1 and y_obs2 determine the range of observational frequency. 
    TerminateTime = 5 ! The Number of scattering times to terminate the scattering sequence.
    Terminate_Tolerence = 1.D-40 ! The Tolerence value to terminate the scattering sequence.
                 ! And Terminate_Tolerence = w / w_ini, w_ini and w are the initial weight and 
                 ! weight, and w will decrease as scattering sequence increase.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods_cases == 1)then  
        tau = 0.05D0 
        T_bb = 10.D0 * 1.D-6    ! In unit of MeV
        T_elec = 352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=352kev.txt' 
        Spentrum_filename = './spectrum/IQUV352_8.txt'
        Total_Phot_Num = 1.0D9
        mu_esti(1) = 0.11D0  
    elseif(methods_cases == 2)then
        tau = 0.5D0 
        T_bb = 10.D0 * 1.D-6     ! In unit of MeV
        T_elec = 56.D0 * 1.D-3 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=56kev1.txt' 
        Spentrum_filename = './spectrum/IQUV56_8.txt'
        Total_Phot_Num = 5.0D8
        mu_esti(1) = 0.5D0  
    endif
    Num_mu_esti = 1
    sin_esti = dsqrt( one - mu_esti**2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

    Total_Phot_Num = 5.D8
    write(Te, "(ES10.3)")T_elec
    write(Ta, "(ES10.3)")tau
    write(E1, "(ES10.3)")E1_scat
    write(E2, "(ES10.3)")E2_scat
    CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
               trim('tau=')//trim(adjustl(Ta))//trim('E1=')//&
               trim(adjustl(E1))//trim('E2=')//trim(adjustl(E2))//trim('.dat') !
    Savefilename = trim('./spectrum/dataTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(Ta))//trim('.dat') !
 
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, T_elec, &
                    E1_scat, E2_scat, y_obs1, y_obs2, mu_esti, sin_esti, &
                    Num_mu_esti, Terminate_Tolerence, &
                    CrossSec_filename, Spentrum_filename )

    END PROGRAM main








 
