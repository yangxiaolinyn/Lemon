!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    program main
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    use Method_Of_FLST_DiffuseReflec
    use BCS_simulations
    implicit none
    real(mcp) :: tau, T_bb, T_elec, gam ,Theta_e, temp_v1, freq_s, J_I, J_I1, J_I2
    real(mcp) :: Intensity(0: N_BCS), log_fs1f(0, N_BCS), S_in(1: 4), alp, &
                 gama1, gama2, vy1, vy2, E_ini, theta, phis
    real(mcp) :: mu_obs
    integer :: methods_cases, i
    integer(kind = 8) :: Total_Phot_Num 
    character*80 ::  Spentrum_filename, CrossSec_filename, &
              filename 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    TYPE(BCS_photons) :: BCS_phot
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
        call InitRandom()  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
        Theta_e = 10.D0  ! Theta_e = T_e / mec2, which is the dimensionless temperature of the
                         ! electron gas.

!~~~~~~~The cosine of the polar angle of the observer~~~~~~~~~~~~~~~~
        mu_obs = dcos( 85.D0 * dtor )

!~~~~~~~~Set the normalized incident Stokes parameters~~~~~~~~~~~~~~~~~
        !S_in(1) = one
        !S_in(2) = one * zero!one / dsqrt(3.D0)
        !S_in(3) = one * zero!one / dsqrt(3.D0)
        !S_in(4) = one !* zero!one / dsqrt(3.D0)

        theta = 30.D0
        phis = 60.D0
        S_in(1) = one
        S_in(2) = dsin(theta*dtor) * dcos(phis*dtor)!one / dsqrt(3.D0)
        S_in(3) = dsin(theta*dtor) * dsin(phis*dtor)!one / dsqrt(3.D0)
        S_in(4) = dcos(theta*dtor) !* zero!one / dsqrt(3.D0)
        filename = './spectrum/bcs_V_alp=3.txt'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        alp = 3.D0   ! the power law index of electron gas.
        gama1 = 60.D0  ! the Lorentz factor gamma of the power law electron gas is 
        gama2 = 10.D9  ! distributed between gama1 and gama2, i.e. gama1 < gamma < gama2.
        vy1 = 1.D0   ! vy is the index of the observed energy: E_obs, i.e.
        vy2 = 7.D0   ! (E_obs / E_ini) = 10^vy, E_ini is the energy of incident beam.
        E_ini = 2.5D-11 * mec2  ! the energy of incident beam.

!~~~~~~~Set initial paramters for the semi-analytic formulae of BCS 1970~~~~~~~~~~~~~~
        BCS_phot%alp = alp
        BCS_phot%gama1 = gama1
        BCS_phot%gama2 = gama2
        BCS_phot%vy1 = vy1
        BCS_phot%vy2 = vy2
        BCS_phot%E_ini = E_ini
!~~~~~~~ Implement the calculations and the results it saved in file: filename.~~~~~~~~~
        call BCS_phot%BCS_analytical_formula_Powerlaw( Theta_e, S_in, mu_obs, filename )
 
        T_elec = 50.D0 * mec2   ! In unit of MeV
        CrossSec_filename = './data/SigmaArrayT_e=100kev.txt' 
        Spentrum_filename = './spectrum/I_Te=100_V_alp=3.txt'
 

        Total_Phot_Num = 1.1D10
        call mimick_of_ph_Slab_BoundReflc( Total_Phot_Num, T_elec, &
                    CrossSec_filename, Spentrum_filename, S_in, alp, & 
                    gama1, gama2, vy1, vy2, E_ini )
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

    END PROGRAM main








 
