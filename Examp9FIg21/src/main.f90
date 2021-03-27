!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
    use Method_Of_FLST_DiffuseReflec
    implicit none
    real(mcp) :: tau, T_bb, T_elec
    integer :: methods_cases
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom()
    methods_cases = 1  
    filenameH3Array = './data/H3Array.dat'
    CrossSec_filename = './data/SigmaArrayT_e=56kev.txt'  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tau = 1.05D100 
    T_bb = 10.D0 * 1.D-6     ! In unit of MeV
    T_elec = 352.D0 * 1.D-3    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Total_Phot_Num = 1.D9 
    call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, tau, T_bb, &
                                        T_elec, CrossSec_filename )

    END PROGRAM main








 
