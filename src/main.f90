    program main  
    !use Method_Of_FLST_ThomScat_Emerge_IQ 
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom() 
    filenameH3Array = './data/H3Array.dat'
    !CrossSec_filename = './data/SigmaArrayT_e=56kev.txt' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methods_cases = 5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods_cases == 1)then
        tau = 0.15D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 1.7D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=352kev1.txt' 
    elseif(methods_cases == 2)then
        tau = 0.01D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 1.59D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=352kev2.txt' 
    elseif(methods_cases == 3)then
        tau = 0.06D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 0.188D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=352kev3.txt' 
    elseif(methods_cases == 4)then
        tau = 0.1D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 0.5D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=250.txt' 
        Spentrum_filename = './spectrum/IQUV250.txt'
    elseif(methods_cases == 5)then
        tau = 0.05D0 
        T_bb = 10.D0 * 1.D-6    ! In unit of MeV
        T_elec = 352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=352kev.txt' 
        Spentrum_filename = './spectrum/IQUV352_5.txt'
        Total_Phot_Num = 1.0D6
    elseif(methods_cases == 6)then
        tau = 0.5D0 
        T_bb = 10.D0 * 1.D-6     ! In unit of MeV
        T_elec = 56.D0 * 1.D-3 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=56kev1.txt' 
        Spentrum_filename = './spectrum/IQUV56_5.txt'
        Total_Phot_Num = 4.0D8
    endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ _DiffuseReflec
    !Total_Phot_Num = 1.2D8
    !call mimick_of_ph_FST_Slab_BoundReflc(Total_Phot_Num, tau, T_bb, T_elec, CrossSec_filename) 
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, T_elec, &
                    CrossSec_filename, Spentrum_filename )

    END PROGRAM main








 
