!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use Method_Of_FLST_DiffuseReflec
    use BCS_simulations
    implicit none
    real(mcp) :: tau, T_bb, T_elec, gam ,Theta_e, temp_v1, freq_s, J_I, J_I1, J_I2
    real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS), S_in(1: 4), alp, &
                 gama1, gama2, vy1, vy2, E_ini, theta, phis
    integer :: methods_cases, i
    integer(kind = 8) :: Total_Phot_Num 
    character*80 :: image_filename, Spentrum_filename, CrossSec_filename, &
              CrossSec_filename_Seq, filenameH3Array, filename

    TYPE(BCS_photons) :: BCS_phot

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom() 
    filenameH3Array = './data/H3Array.dat'
    !CrossSec_filename = './data/SigmaArrayT_e=56kev.txt' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    methods_cases = 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(methods_cases == 1)then 
        !call BCS_phot%BCS_analytical_formula( Theta_e )
    elseif(methods_cases == 2)then
        Theta_e = 10.D0
        S_in(1) = one
        S_in(2) = one * zero!one / dsqrt(3.D0)
        S_in(3) = one * zero!one / dsqrt(3.D0)
        S_in(4) = one !* zero!one / dsqrt(3.D0)

        theta = 30.D0
        phis = 60.D0
        S_in(2) = dsin(theta*dtor) * dcos(phis*dtor)!one / dsqrt(3.D0)
        S_in(3) = dsin(theta*dtor) * dsin(phis*dtor)!one / dsqrt(3.D0)
        S_in(4) = dcos(theta*dtor) !* zero!one / dsqrt(3.D0)
        filename = './spectrum/bcs_V_alp=3.txt'

        alp = 3.D0
        gama1 = 60.D0
        gama2 = 10.D9
        vy1 = 1.D0
        vy2 = 7.D0
        E_ini = 2.5D-11 * mec2

        BCS_phot%alp = alp
        BCS_phot%gama1 = gama1
        BCS_phot%gama2 = gama2
        BCS_phot%vy1 = vy1
        BCS_phot%vy2 = vy2
        BCS_phot%E_ini = E_ini

        call BCS_phot%BCS_analytical_formula( Theta_e, S_in, filename )

        !tau = 0.05D0 
        !T_bb = 10.D0 * 1.D-6     ! In unit of MeV
        T_elec = 50.D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=100kev.txt' 
        Spentrum_filename = './spectrum/I_Te=100_V_alp=3.txt'

        !write(*, fmt=*)'ffsdfa='

        Total_Phot_Num = 1.1D10
        call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, T_elec, &
                    CrossSec_filename, Spentrum_filename, S_in, alp, gama1, gama2, vy1, vy2, E_ini )
        !write(*, fmt=*)'ffsdfa=', S_in
    endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

    END PROGRAM main








 
