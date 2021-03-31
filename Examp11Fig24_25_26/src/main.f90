!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    use Method_Of_FLST_DiffuseReflec
    use BCS_simulations
    implicit none
    real(mcp) :: T_elec, Theta_e 
    real(mcp) :: S_in(1: 4), alp, &
                 gama1, gama2, vy1, vy2, E_ini, theta, phis
    real(mcp) :: mu_obs
    integer :: input_cases, i
    integer(kind = 8) :: Total_Phot_Num 
    character*80 ::  Spentrum_filename, CrossSec_filename, &
              filename, alps, Te, S_in_case
    logical :: case_powerlaw
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    TYPE(BCS_photons) :: BCS_phot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
    call InitRandom()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
   
    vy1 = 1.D0   ! vy is the index of the observed energy: E_obs, i.e.
    vy2 = 7.D0   ! (E_obs / E_ini) = 10^vy, E_ini is the energy of incident beam.
    E_ini = 2.5D-11 * mec2  ! the energy of incident beam.
!~~~~~~~The cosine of the polar angle of the observer~~~~~~~~~~~~~~~~
    mu_obs = dcos( 85.D0 * dtor )

    case_powerlaw = .true.
    if( case_powerlaw )then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

!~~~~~~~~Set the normalized incident Stokes parameters~~~~~~~~~~~~~~~~~
        input_cases = 1
        select case ( input_cases )
        case(1)
            theta = 30.D0
            phis = 60.D0
            S_in(1) = one
            S_in(2) = dsin(theta*dtor) * dcos(phis*dtor)!one / dsqrt(3.D0)
            S_in(3) = dsin(theta*dtor) * dsin(phis*dtor)!one / dsqrt(3.D0)
            S_in(4) = dcos(theta*dtor) !* zero!one / dsqrt(3.D0)
            S_in_case = "_fig26"
        case(2)
            S_in(1) = one
            S_in(2) = one
            S_in(3) = zero
            S_in(4) = zero
            S_in_case = "_fig25_1"
        case(3)
            S_in(1) = one
            S_in(2) = -one
            S_in(3) = zero
            S_in(4) = zero
            S_in_case = "_fig25_2"
        case(4)
            S_in(1) = one
            S_in(2) = zero
            S_in(3) = one
            S_in(4) = zero
            S_in_case = "_fig25_3"
        case(5)
            S_in(1) = one
            S_in(2) = zero
            S_in(3) = zero
            S_in(4) = one
            S_in_case = "_fig25_4"
        end select 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        alp = 3.D0   ! the power law index of electron gas.
        gama1 = 60.D0  ! the Lorentz factor gamma of the power law electron gas is 
        gama2 = 10.D9  ! distributed between gama1 and gama2, i.e. gama1 < gamma < gama2.

!~~~~~~~Set initial paramters for the semi-analytic formulae of BCS 1970~~~~~~~~~~~~~~
        BCS_phot%alp = alp
        BCS_phot%gama1 = gama1
        BCS_phot%gama2 = gama2
        BCS_phot%vy1 = vy1
        BCS_phot%vy2 = vy2
        BCS_phot%E_ini = E_ini
!~~ N_coef is the normalization factor for the power law distributed electron gas.
        BCS_phot%N_coef = ( BCS_phot%alp - one ) / &
                   ( BCS_phot%gama1**(one - BCS_phot%alp) - &
                     BCS_phot%gama2**(one - BCS_phot%alp) )
 
 
        write(alps, "(f8.4)")alp
        if( input_cases == 1 )then
            Spentrum_filename = trim('./spectrum/Fig26/dataPowerlaw_alp=')//&
                       trim(adjustl(alps))//trim(S_in_case)//trim('.dat') 
            filename = trim('./spectrum/Fig26/BCS_Powerlaw_alp=')//&
                       trim(adjustl(alps))//trim(S_in_case)//trim('.dat')
        else
            Spentrum_filename = trim('./spectrum/Fig25/dataPowerlaw_alp=')//&
                       trim(adjustl(alps))//trim(S_in_case)//trim('.dat') 
            filename = trim('./spectrum/Fig25/BCS_Powerlaw_alp=')//&
                       trim(adjustl(alps))//trim(S_in_case)//trim('.dat')
        endif

!~~~~~~~ Implement the calculations and the results it saved in file: filename.~~~~~~~~~
        call BCS_phot%BCS_analytical_formula_Powerlaw( S_in, mu_obs, filename )

        !T_elec = 50.D0 * mec2   ! In unit of MeV 
        !write(Te, "(ES10.3)")T_elec
        !CrossSec_filename = trim('./data/SigArrTe=')//trim(adjustl(Te))//&
        !       trim('_alp=')//trim(adjustl(alps))//trim('.dat') !
  
        Total_Phot_Num = 1.1D9
        call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, T_elec, &
                    CrossSec_filename, Spentrum_filename, S_in, alp, & 
                    gama1, gama2, vy1, vy2, E_ini, mu_obs, case_powerlaw ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    else
!~~~~~~~~Set the normalized incident Stokes parameters~~~~~~~~~~~~~~~~~
        input_cases = 4
        select case ( input_cases ) 
        case(1)
            S_in(1) = one
            S_in(2) = one
            S_in(3) = zero
            S_in(4) = zero
            S_in_case = "_fig24_1"
        case(2)
            S_in(1) = one
            S_in(2) = -one
            S_in(3) = zero
            S_in(4) = zero
            S_in_case = "_fig24_2"
        case(3)
            S_in(1) = one
            S_in(2) = zero
            S_in(3) = one
            S_in(4) = zero
            S_in_case = "_fig24_3"
        case(4)
            S_in(1) = one
            S_in(2) = zero
            S_in(3) = zero
            S_in(4) = one
            S_in_case = "_fig24_4"
        end select 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        T_elec = 100.D0 * mec2   ! In unit of MeV
        Theta_e = 100.D0  ! Theta_e = T_e / mec2, which is the 
                          ! dimensionless temperature of the electron gas.

        write(Te, "(f8.4)")Theta_e
        Spentrum_filename = trim('./spectrum/Fig24/dataHotElectron_Te=')//&
                   trim(adjustl(Te))//trim(S_in_case)//trim('.dat') 
        filename = trim('./spectrum/Fig24/BCS_HotElectron_Te=')//&
                   trim(adjustl(Te))//trim(S_in_case)//trim('.dat')

        call BCS_phot%BCS_analytical_formula_HotElectron( Theta_e, S_in, filename )
 
        Total_Phot_Num = 1.2D9

        call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, T_elec, &
                    CrossSec_filename, Spentrum_filename, S_in, alp, & 
                    gama1, gama2, vy1, vy2, E_ini, mu_obs, case_powerlaw ) 
    endif

    END PROGRAM main








 
