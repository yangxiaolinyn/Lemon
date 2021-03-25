    program main   
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: Savefilename, CrossSec_filename, filenameH3Array, Te, Ta

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    filenameH3Array = './data/H3Array.dat' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methods_cases = 5
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods_cases == 1)then
        tau = 0.15D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 1.7D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV 
    elseif(methods_cases == 2)then
        tau = 0.01D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 1.59D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV 
    elseif(methods_cases == 3)then
        tau = 0.06D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 0.188D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV 
    elseif(methods_cases == 4)then
        tau = 0.1D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 0.5D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV 
    elseif(methods_cases == 5)then
        tau = 0.1D0 
        T_bb = 1.D-5 * mec2     ! In unit of MeV
        T_elec = 1.D0 * mec2 !352.D0 * 1.D-3 !0.11D0 * mec2   ! In unit of MeV 
    endif

    Total_Phot_Num = 1.D9
    write(Te, "(f8.3)")T_elec
    write(Ta, "(f8.5)")tau
    CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(Ta))//trim('.txt') !
    Savefilename = trim('./spectrum/dataTe=')//trim(adjustl(Te))//&
                      trim('tau=')//trim(adjustl(Ta))//trim('.txt') !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ _DiffuseReflec
    !Total_Phot_Num = 1.2D8
    !call mimick_of_ph_FST_Slab_BoundReflc(Total_Phot_Num, tau, T_bb, T_elec, CrossSec_filename) 
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, &
                            T_elec, CrossSec_filename, Savefilename, methods_cases )

    END PROGRAM main








 
