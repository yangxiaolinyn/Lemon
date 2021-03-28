!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec, E1_scat, E2_scat, mu_estis(1: 8)
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: Savefilename, Hxm_analyticResult, CrossSec_filename, &
                       filenameH3Array, Te, Ta, E1, E2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom() 
    filenameH3Array = './data/H3Array.dat'  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = 1.05D100 
    !T_bb = 10.D0 * 1.D-6     ! In unit of MeV    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    E1_scat = 1.D-7 ! In unit of MeV. E1_scat, E2_scat define the energy range 
    E2_scat = 1.D1  ! of scattering. The total scattering coefficient: Sigma_a(T_e, E_p) also
                    ! define on this interval, i.e., E_p \in (E1_scat, E2_scat). E_p is the
                    ! photon energy.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mu_estis(1) = 0.05D0
    mu_estis(2) = 0.1D0
    mu_estis(3) = 0.15D0
    mu_estis(4) = 0.25D0
    mu_estis(5) = 0.50D0
    mu_estis(6) = 0.75D0 
    mu_estis(7) = 0.95D0
    mu_estis(8) = 1.D0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    write(Te, "(ES10.3)")T_elec 
    write(E1, "(ES10.3)")E1_scat
    write(E2, "(ES10.3)")E2_scat
    CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
               trim('E1=')//trim(adjustl(E1))//trim('E2=')//trim(adjustl(E2))//trim('.dat')
    !Savefilename = trim('./spectrum/dataTe=')//trim(adjustl(Te))//trim('.dat')
    Hxm_analyticResult = './spectrum/HxmImu1.dat'

    Total_Phot_Num = 2.D8
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, &
                    T_elec, E1_scat, E2_scat, mu_estis, CrossSec_filename, &
                    Hxm_analyticResult )

    END PROGRAM main








 
