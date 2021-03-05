      module Basic_Variables_And_Methods
      use constants
      use RandUtils
      use SubFunction
      use BLcoordinate 
      use Integrations
      use IntialParameters
      implicit none

      type, abstract, public :: Basic_Variables_And_Methods_Of_Particle
          logical :: my_test = .FALSE. 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The BL coordinates of the Particle
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: r
          real(mcp) :: theta
          real(mcp) :: mucos  ! mucos = cos(theta)
          real(mcp) :: musin  ! musin = sin(theta) 
          real(mcp) :: phi
          real(mcp) :: t
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          real(mcp) :: r_ini
          real(mcp) :: theta_ini
          real(mcp) :: mucos_ini  ! mucos = cos(theta)
          real(mcp) :: musin_ini  ! musin = sin(theta)  
          real(mcp) :: t_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          real(mcp) :: r_p
          real(mcp) :: theta_p
          real(mcp) :: mucos_p  ! mucos = cos(theta)
          real(mcp) :: musin_p  ! musin = sin(theta) 
          real(mcp) :: phi_p
          real(mcp) :: t_p
          real(mcp) :: sign_pr_p
          real(mcp) :: sign_pth_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          real(mcp) :: r_samp
          real(mcp) :: theta_samp
          real(mcp) :: mucos_samp  ! mucos = cos(theta)
          real(mcp) :: musin_samp  ! musin = sin(theta) 
          real(mcp) :: phi_samp
          real(mcp) :: t_samp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! The four-Momentum of the particle(the photon)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1:4) :: Elec4U_BL
          real(mcp), dimension(1:4) :: Phot4k_CtrCF
          real(mcp), dimension(1:4) :: Phot4k_CovCF
          real(mcp), dimension(1:4) :: Phot4k_LNRF
          real(mcp), dimension(1:4) :: Phot4k_BL
          real(mcp), dimension(1:4) :: Phot4k_CovBL
          real(mcp), dimension(1:4) :: Phot4k_CtrBL
          real(mcp), dimension(1:4) :: Phot_Covariant4k_BL
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1:4) :: Phot4k_LNRF_ini
          real(mcp), dimension(1:4) :: Phot4k_CtrBL_ini
          real(mcp), dimension(1:4) :: Phot4k_CovBL_ini
          real(mcp), dimension(1:4) :: Elec4U_CtrBL_At_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1:4) :: Phot4k_CtrCF_ini
          real(mcp), dimension(1:4) :: Phot4k_CovCF_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1: 4) :: Phot4k_CovBL_At_p
          real(mcp), dimension(1: 4) :: Phot4k_CtrBL_At_p
          real(mcp), dimension(1: 4) :: Phot4k_CovCF_At_p
          real(mcp), dimension(1: 4) :: Phot4k_CtrCF_At_p
          real(mcp), dimension(1: 4) :: Elec4U_CtrBL_At_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1: 4) :: Elec4U_CtrBL_At_samp
          real(mcp), dimension(1: 4) :: Phot4k_CovBL_At_samp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1: 4) :: f4_CovBL      !the 4 vector of polarization of Photon
          real(mcp), dimension(1: 4) :: f4_CtrBL      !the 4 vector of polarization of Photon
          real(mcp), dimension(1: 4) :: f4_CovBL_At_p
          real(mcp), dimension(1: 4) :: f4_CtrBL_At_p
          real(mcp), dimension(1: 4) :: f4_CtrCF_At_p    !the 4 vector of polarization of Photon
          real(mcp), dimension(1: 4) :: f4_CovCF_At_p    !the 4 vector of polarization of Photon
          real(mcp), dimension(1: 4) :: f4_scat_CF       !f4 affer scattering
          real(mcp), dimension(1: 4) :: f4_scat_CovCF    !f4 affer scattering
          real(mcp), dimension(1: 4) :: f4_CF
          real(mcp), dimension(1: 4) :: f4_CovCF
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: lambda
          real(mcp) :: q
          real(mcp) :: lambda_ini
          real(mcp) :: q_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: K1_pw_ini
          real(mcp) :: K2_pw_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: sign_pr_samp
          real(mcp) :: sign_pth_samp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: E_ini
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: Q_sp      !sp means Stocks Parameters, and Q_sp = Q / I
          real(mcp) :: U_sp      !sp means Stocks Parameters, and U_sp = U / I
          real(mcp) :: V_sp      !sp means Stocks Parameters, and V_sp = V / I
          real(mcp) :: delta_pd  !pd means Polarization Degree, delta_pd = sqrt[ Q_sp^2 + U_sp^2 ]
          real(mcp) :: Q_sp_scat      !sp means Stocks Parameters, Q_sp = Q / I
          real(mcp) :: U_sp_scat      !sp means Stocks Parameters, U_sp = U / I
          real(mcp) :: V_sp_scat      !sp means Stocks Parameters, V_sp = V / I
          real(mcp) :: delta_pd_scat  !pd means Polarization Degree, delta_pd = sqrt[ Q_sp^2 + U_sp^2 ]
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: w_ini
          real(mcp) :: w_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: p
          real(mcp) :: p_scattering
          real(mcp) :: p_samp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Quantities of Kerr Metric
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: aspin
          real(mcp) :: Delta
          real(mcp) :: Sigma
          real(mcp) :: bigA
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: expnu
          real(mcp) :: exppsi
          real(mcp) :: expmu1
          real(mcp) :: expmu2
          real(mcp) :: somega
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: g_tt
          real(mcp) :: g_rr
          real(mcp) :: g_thth
          real(mcp) :: g_tp
          real(mcp) :: g_pt
          real(mcp) :: g_pp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: g_up_tt
          real(mcp) :: g_up_rr
          real(mcp) :: g_up_thth
          real(mcp) :: g_up_tp
          real(mcp) :: g_up_pt
          real(mcp) :: g_up_pp         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: lambda_p
          real(mcp) :: q_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: Delta_p
          real(mcp) :: Sigma_p
          real(mcp) :: bigA_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: expnu_p
          real(mcp) :: exppsi_p
          real(mcp) :: expmu1_p
          real(mcp) :: expmu2_p
          real(mcp) :: somega_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: g_tt_p
          real(mcp) :: g_rr_p
          real(mcp) :: g_thth_p
          real(mcp) :: g_pp_p
          real(mcp) :: g_tp_p
          real(mcp) :: g_pt_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: g_up_tt_p
          real(mcp) :: g_up_rr_p
          real(mcp) :: g_up_thth_p
          real(mcp) :: g_up_pp_p
          real(mcp) :: g_up_tp_p
          real(mcp) :: g_up_pt_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1:4, 1:4) :: Lorentz_Matrix_p
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_BLCF_a_mu_p
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_BLCF_mu_a_p
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_CF2LNRF_p
          real(mcp) :: absorption_coef !absorption coefficient at Frequency nu
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: Delta_samp
          real(mcp) :: Sigma_samp
          real(mcp) :: bigA_samp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: expnu_samp
          real(mcp) :: exppsi_samp
          real(mcp) :: expmu1_samp
          real(mcp) :: expmu2_samp
          real(mcp) :: somega_samp
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: g_tt_samp
          real(mcp) :: g_rr_samp
          real(mcp) :: g_thth_samp
          real(mcp) :: g_pp_samp
          real(mcp) :: g_tp_samp
          real(mcp) :: g_pt_samp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The Descartes coordinates of the Particle in Flat space
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: x
          real(mcp) :: y
          real(mcp) :: z 
          !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: x_ini
          real(mcp) :: y_ini
          real(mcp) :: z_ini
          real(mcp) :: x_p
          real(mcp) :: y_p
          real(mcp) :: z_p
          real(mcp) :: x_samp
          real(mcp) :: y_samp
          real(mcp) :: z_samp
          !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: rg
          real(mcp) :: BigM
          !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          real(mcp), dimension(1:3) :: Vector_of_position
          real(mcp), dimension(1:3) :: Vector_of_position_ini
          real(mcp), dimension(1:3) :: Vector_of_position_p
          real(mcp), dimension(1:3) :: Vector_of_Momentum
          real(mcp), dimension(1:3) :: Vector_of_Momentum_ini
          real(mcp), dimension(1:3) :: Vector_of_Momentum_p
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^ The BL coordinates of the photon at p, which determines the BL coordinates  ^
!^ by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  ^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_BLCF_a_mu
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_BLCF_mu_a
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_CF2LNRF 
          real(mcp), dimension(1:5) :: KerrMetric_exp
          real(mcp), dimension(1:6) :: KerrMetric_g_ij
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Four Velocity and Four Momentum in Kerr SpaceTime
          real(mcp), dimension(1:4) :: ParticleU4_BL     !BL means Boyer Lindquist Coordinate.
          real(mcp), dimension(1:4) :: ParticleU4_CF     !CF means Comoving Frame
          real(mcp), dimension(1:4) :: ParticleU4_LNRF
! The Comoving Frame's Four Velocity With Respect to BL Coordinate
          real(mcp), dimension(1:4) :: ComFrameU4_BL
          real(mcp), dimension(1:4) :: ComFrameU4_BL_At_p
! The Comoving Frame's Three Physical Velocity With Respect to LNRF
          real(mcp), dimension(1:3) :: ComFrameV3_LNRF
! LNRF's Tetrad Matrix
          real(mcp), dimension(1:4, 1:4) :: eLNRF
          real(mcp), dimension(1:4, 1:4) :: Inverse_eLNRF
 
! Lorentz Transformation Matrix
          real(mcp), dimension(1:4, 1:4) :: Lorentz_Matrix
          real(mcp), dimension(1:4, 1:4) :: Lorentz_Inverse_Matrix

! Transformation Matrix from CF to BL coordinate. 
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_CtrCF2CtrBL
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_CovCF2CovBL
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_CovBL2CF
          real(mcp), dimension(1:4, 1:4) :: Matrix_Of_ContraBL2CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Quantities in terms of a Observer at a distance r_obs, he/her has a virtual screen to 
! Image the picture of the black hole and accretion disks. the screen has a coordiantes 
! alpha, beta, each photon starts from the screen corresponding to a special coordiante.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: alpha
          real(mcp) :: beta
          real(mcp) :: sin_theta_obs
          real(mcp) :: cos_theta_obs
          real(mcp) :: theta_obs
          real(mcp) :: r_obs
          real(mcp) :: E_obs
          real(mcp) :: E_em
          real(mcp) :: phot4k_obs(1: 4)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: Delta_obs
          real(mcp) :: Sigma_obs
          real(mcp) :: bigA_obs
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: expnu_obs
          real(mcp) :: exppsi_obs
          real(mcp) :: expmu1_obs
          real(mcp) :: expmu2_obs
          real(mcp) :: somega_obs
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp), dimension(0: N_H3) :: H3Array = zero
          integer :: PS_Nums
          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          real(mcp) :: Qksi3
          real(mcp) :: Uksi1
          real(mcp) :: Vksi2
          real(mcp) :: delta_PolDeg
          real(mcp) :: r_one_hvlmec2_one_cosE
          real(mcp) :: N_scat
          real(mcp) :: z_tau
          real(mcp) :: tau_max
          real(mcp) :: I_IQ
          real(mcp) :: Q_IQ
          real(mcp) :: Psi_I
          real(mcp) :: Psi_Q
          real(mcp) :: Psi_U
          real(mcp) :: Psi_V 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: phi_ini 
          real(mcp) :: cosphi_ini
          real(mcp) :: sinphi_ini
          real(mcp) :: cos2phi_ini
          real(mcp) :: sin2phi_ini
          real(mcp) :: Scat_phi
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          integer :: medium_case
          real(mcp) :: mu_estimat(1: 4)
          real(mcp) :: mu_estimates(0: Num_mu)
          real(mcp) :: smu_estimates(0: Num_mu)
          real(mcp) :: phi_estimat(1: 4)
          real(mcp) :: phi_estimates(0: Num_Phi)
          real(mcp) :: f_IQUV_estimat(1: 4) 
          real(mcp) :: logE_low
          real(mcp) :: logE_up
          real(mcp) :: dE_esti
          real(mcp) :: dindexE
          real(mcp) :: mu_esti(1: 2)
          real(mcp) :: smu_esti(1: 2)
          real(mcp) :: phi_esti(0: N_esti_phi)
          real(mcp) :: cos_phi_esti(0: Num_Phi)
          real(mcp) :: sin_phi_esti(0: Num_Phi)
          real(mcp) :: E_array_esti(0: vL_sc_up) 
          real(mcp) :: Phot3k_CF_esti(1: 3) 
          real(mcp) :: Phot3k_EF_esti(1: 3) 
          real(mcp) :: Phot3k_ECF_esti(1: 3)
          real(mcp) :: Phot3k_ECF1_esti(1: 3) 
          real(mcp) :: Phot3k_ECF_ini(1: 3) 

          real(mcp) :: Scat_Phot3k_CF(1: 3) 
          !real(mcp) :: nu_vs_nup
          !real(mcp) :: nup_vs_nu

          real(mcp) :: PhotE_ECF_esti
          real(mcp) :: E_esti_vs_mu_phi
          real(mcp) :: y1, y2, dy 
          real(mcp) :: Matrix_Elec_2_CF(1:3, 1:3)
          real(mcp) :: Matrix_CF_2_Elec(1:3, 1:3)
          real(mcp) :: Matrix_ECF_2_ECF1(1:3, 1:3)
          real(mcp) :: Matrix_ECF1_2_ECF(1:3, 1:3)
          real(mcp) :: z_times_SigmaT
          real(mcp) :: ne_times_SigmaT 
          real(mcp) :: ne_times_SigmaKN
          real(mcp) :: ne_times_Sigma_a
          real(mcp) :: z_max
          real(mcp) :: z_max2
          real(mcp) :: tau_max2
          real(mcp) :: n_e1
          real(mcp) :: n_e2
          real(mcp) :: Sigma_KN_E_ini
          real(mcp) :: Sigma_a_E_ini
          real(mcp) :: CROS_absorption
          real(mcp) :: log10_Tbb
          real(mcp) :: P_mu_normal
          real(mcp) :: P_nu_normal
          !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
          real(mcp) :: Phot_Tetrad1_ECF_axisX(1:3)
          real(mcp) :: Phot_Tetrad1_ECF_axisY(1:3)
          real(mcp) :: Phot_Tetrad1_ECF_axisZ(1:3)
          real(mcp) :: Vector_Stokes3_CF(1:3)
          real(mcp) :: Vector_Stokes4_CF(1:4)
          real(mcp) :: Vector_Stokes4_ECF(1:4)
          real(mcp) :: Vector_Stokes4_ECF_scat(1:4)
          real(mcp) :: Vector_Stokes4_CF_scat(1:4)
            
          real(mcp) :: z_max1
          real(mcp) :: p_reflect
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      contains
          procedure, public :: Set_Kerr_Metric               => Set_Kerr_Metric_Sub
          procedure, public :: Set_Inverse_Kerr_Metric       => Set_Inverse_Kerr_Metric_sub 
          procedure, public :: Set_Transformation_Matrices   => Set_Transformation_Matrices_Sub
          procedure, public :: Set_Transformation_Matrices_At_p => &
                                Set_Transformation_Matrices_At_p_Sub
          procedure, public :: get_lambdaq        =>  get_lambdaq_Sub 
          procedure, public :: get_lambdaq_At_p   =>  get_lambdaq_At_p_Sub
          procedure, public :: Get_Phot_BLCoordinates_At_p       => Get_Phot_BLCoordinates_At_p_Sub
          procedure, public :: Get_Phot_BLCoor_At_samp           => Get_Phot_BLCoor_At_samp_Sub
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Set_Kerr_Metric_At_p              => Set_Kerr_Metric_At_p_Sub
          procedure, public :: Set_Inverse_Kerr_Metric_At_p      => Set_Inverse_Kerr_Metric_At_p_Sub
          procedure, public :: Set_Kerr_Metric_At_samp           => Set_Kerr_Metric_At_samp_Sub
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Get_Phot4k_CovBL_At_p             => Get_Phot4k_CovBL_At_p_Sub
          procedure, public :: Get_Phot4k_CtrBL_At_p             => Get_Phot4k_CtrBL_At_p_Sub
          procedure, public :: Get_Phot4k_CovBL_At_samp          => Get_Phot4k_CovBL_At_samp_Sub
          procedure, public :: get_Phot4k_LNRF                   => get_Phot4k_LNRF_Sub
          procedure, public :: get_Phot4k_CtrBL                  => get_Phot4k_CtrBL_Sub
          procedure, public :: get_Phot4k_CovBL                  => get_Phot4k_CovBL_Sub
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: get_K1_K2_pw                =>  get_K1_K2_pw_Sub
          procedure, public :: get_K1_K2_pw_At_ini         =>  get_K1_K2_pw_At_ini_Sub
          procedure, public :: get_f4_Ctr_BL_At_p          =>  get_f4_Ctr_BL_At_p_Sub
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !procedure, public :: get_Tau                           =>  get_Tau_fn
          !procedure, public :: get_Tau2                          =>  get_Tau2_fn
          !procedure, public :: Set_Cross_Section                 =>  Set_Cross_Section_sub
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Set_Elec4U_CtrBL_At_samp          => Set_Elec4U_CtrBL_At_samp_Sub
          procedure, public :: Set_Elec4U_CtrBL_At_p             => Set_Elec4U_CtrBL_At_p_Sub
          !procedure, public :: sigma_fn
          procedure, public :: Get_filename     =>    Get_filename_Sub
          procedure, public :: Get_filename2    =>    Get_filename2_Sub
          procedure, public :: r_p2  =>   r_p2_Sub
      end type Basic_Variables_And_Methods_Of_Particle

      private :: Set_Kerr_Metric_sub
      private :: Set_Inverse_Kerr_Metric_sub 
      private :: Set_Transformation_Matrices_Sub
      private :: Set_Transformation_Matrices_At_p_Sub
      private :: get_lambdaq_Sub
      private :: get_lambdaq_At_p_Sub
      private :: Get_Phot_BLCoordinates_At_p_Sub
      private :: Get_Phot_BLCoor_At_samp_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Set_Kerr_Metric_At_p_Sub
      private :: Set_Inverse_Kerr_Metric_At_p_Sub
      private :: Set_Kerr_Metric_At_samp_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Get_Phot4k_CovBL_At_p_Sub
      private :: Get_Phot4k_CtrBL_At_p_Sub
      private :: Get_Phot4k_CovBL_At_samp_Sub
      private :: get_Phot4k_LNRF_Sub
      private :: get_Phot4k_CovBL_Sub
      private :: get_Phot4k_CtrBL_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: get_K1_K2_pw_At_ini_Sub
      private :: get_f4_Ctr_BL_At_p_Sub
      private :: get_K1_K2_pw_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Set_Elec4U_CtrBL_At_p_Sub
      private :: Set_Elec4U_CtrBL_At_samp_Sub
      private :: Get_filename_Sub
      private :: Get_filename2_Sub
      private :: r_p2_Sub

      contains

!*******************************************************************************************************
      subroutine Set_Kerr_Metric_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this  
      !real(mcp), intent(in), optional :: rp, mucosp, musinp
      !real(mcp) :: r, mucos, musin, a

      associate( a => this%aspin, &
                 r => this%r, &
                 mucos => this%mucos, &
                 musin => this%musin )
      this%Delta = r*r - two*r + a**2
      this%Sigma = r*r + (a * mucos)**2
      this%bigA = (r*r + a**2)**2 - this%Delta*( a * musin )**2

      this%expnu = dsqrt( this%Sigma*this%Delta/this%bigA )
      this%exppsi = musin * dsqrt( this%bigA/this%Sigma )
      this%expmu1 = dsqrt( this%Sigma / this%Delta )
      this%expmu2 = dsqrt( this%Sigma )
      this%somega = two * a * r / this%bigA

      this%g_tt = - this%Sigma * this%Delta / this%bigA  +  &
                    musin**2 * this%bigA / this%Sigma * this%somega**2 
      this%g_rr = this%Sigma / this%Delta
      this%g_thth = this%Sigma 
      this%g_pp = musin**2 * this%bigA / this%Sigma
      this%g_tp = - this%somega * this%g_pp
      this%g_pt = this%g_tp
      end associate  

      !this%is_Kerr_Metric_Setted = .TRUE.

      end subroutine Set_Kerr_Metric_sub

!*******************************************************************************************************
      subroutine Set_Inverse_Kerr_Metric_Sub(this, rp, mucosp, musinp)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this  
      real(mcp), intent(in), optional :: rp, mucosp, musinp
      real(mcp) :: r, mucos, musin, a

      a = this%aspin
      If ( present(rp) ) then
          r = rp
          mucos = mucosp
          musin = musinp
      Else 
          r     = this%r
          mucos = this%mucos
          musin = this%musin
      Endif
      this%Delta = r*r - two*r + a**2
      this%Sigma = r*r + (a * mucos)**2
      this%bigA = (r*r + a**2)**2 - this%Delta*( a * musin )**2

      this%expnu = dsqrt( this%Sigma*this%Delta / this%bigA )
      this%exppsi = musin * dsqrt( this%bigA / this%Sigma )
      this%expmu1 = dsqrt( this%Sigma / this%Delta )
      this%expmu2 = dsqrt( this%Sigma )
      this%somega = two * a * r / this%bigA

      this%g_up_tt = - this%bigA / this%Sigma / this%Delta
      this%g_up_rr = this%Delta / this%Sigma
      this%g_up_thth = one / this%Sigma 
      this%g_up_pp = this%Sigma / musin**2 / this%bigA + this%g_up_tt * ( this%somega )**2
      this%g_up_tp = this%somega * this%g_up_tt
      this%g_up_pt = this%g_up_tp
   
      end subroutine Set_Inverse_Kerr_Metric_sub

!*******************************************************************************************************
      subroutine Set_Kerr_Metric_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      associate( rp => this%r_p, &
                 a => this%aspin, &
                 mucosp => this%mucos_p, &
                 musinp => this%musin_p )   
      this%Delta_p = rp*rp - two*rp + a**2
      this%Sigma_p = rp*rp + (a * mucosp)**2
      this%bigA_p = (rp*rp + a**2)**2 - this%Delta_p*( a * musinp )**2

      this%expnu_p = dsqrt( this%Sigma_p * this%Delta_p / this%bigA_p )
      this%exppsi_p = musinp * dsqrt( this%bigA_p / this%Sigma_p )
      this%expmu1_p = dsqrt( this%Sigma_p / this%Delta_p )
      this%expmu2_p = dsqrt( this%Sigma_p )
      this%somega_p = two * a * rp / this%bigA_p

      this%g_tt_p = - this%Sigma_p * this%Delta_p / this%bigA_p  +  &
                    musinp**2 * this%bigA_p / this%Sigma_p * this%somega_p**2 
      this%g_rr_p = this%Sigma_p / this%Delta_p
      this%g_thth_p = this%Sigma_p
      this%g_pp_p = musinp**2 * this%bigA_p / this%Sigma_p
      this%g_tp_p = - this%somega_p * this%g_pp_p
      this%g_pt_p = this%g_tp_p

      end associate

      end subroutine Set_Kerr_Metric_At_p_Sub

!*******************************************************************************************************
      subroutine Set_Inverse_Kerr_Metric_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      associate( r => this%r_p, &
                 a => this%aspin, &
                 mucos => this%mucos_p, &
                 musin => this%musin_p )
      this%Delta_p = r*r - two*r + a**2
      this%Sigma_p = r*r + (a * mucos)**2
      this%bigA_p = (r*r + a**2)**2 - this%Delta_p*( a * musin )**2

      this%expnu_p = dsqrt( this%Sigma_p*this%Delta_p / this%bigA_p )
      this%exppsi_p = musin * dsqrt( this%bigA_p / this%Sigma_p )
      this%expmu1_p = dsqrt( this%Sigma_p / this%Delta_p )
      this%expmu2_p = dsqrt( this%Sigma_p )
      this%somega_p = two * a * r / this%bigA_p

      this%g_up_tt_p = - this%bigA_p / this%Sigma_p / this%Delta_p
      this%g_up_rr_p = this%Delta_p / this%Sigma_p
      this%g_up_thth_p = one / this%Sigma_p 
      this%g_up_pp_p = this%Sigma_p / musin**2 / this%bigA_p + &
                       this%g_up_tt_p * ( this%somega_p )**2
      this%g_up_tp_p = this%somega_p * this%g_up_tt_p
      this%g_up_pt_p = this%g_up_tp_p
      end associate
   
      end subroutine Set_Inverse_Kerr_Metric_At_p_sub

!*******************************************************************************************************
      subroutine Set_Kerr_Metric_At_samp_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      !*******************************************************************

      associate( rp => this%r_samp, &
                 a => this%aspin, &
                 mucosp => this%mucos_samp, &
                 musinp => this%musin_samp )   
      this%Delta_samp = rp*rp - two*rp + a**2
      this%Sigma_samp = rp*rp + (a * mucosp)**2
      this%bigA_samp = (rp*rp + a**2)**2 - this%Delta_samp*( a * musinp )**2

      this%expnu_samp = dsqrt( this%Sigma_samp * this%Delta_samp / this%bigA_samp )
      this%exppsi_samp = musinp * dsqrt( this%bigA_samp / this%Sigma_samp )
      this%expmu1_samp = dsqrt( this%Sigma_samp / this%Delta_samp )
      this%expmu2_samp = dsqrt( this%Sigma_samp )
      this%somega_samp = two * a * rp / this%bigA_samp

      this%g_tt_samp = - this%Sigma_samp * this%Delta_samp / this%bigA_samp  +  &
                    musinp**2 * this%bigA_samp / this%Sigma_samp * this%somega_samp**2 
      this%g_rr_samp = this%Sigma_samp / this%Delta_samp
      this%g_thth_samp = this%Sigma_samp
      this%g_pp_samp = musinp**2 * this%bigA_samp / this%Sigma_samp
      this%g_tp_samp = - this%somega_samp * this%g_pp_samp
      this%g_pt_samp = this%g_tp_samp

      end associate   

      end subroutine Set_Kerr_Metric_At_samp_Sub


!*******************************************************************************************************
      SUBROUTINE Set_Transformation_Matrices_Sub(this)
!*******************************************************************************************************
      class(Basic_Variables_And_Methods_Of_Particle) :: this  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp), dimension(1:3) :: CFV3_LNRF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp), dimension(1:4, 1:4) :: Lorentz_Matrix
      real(mcp), dimension(1:4, 1:4) :: Lorentz_Inverse_Matrix
      real(mcp), dimension(1:4, 1:4) :: eLNRF
      real(mcp), dimension(1:4, 1:4) :: Inverse_eLNRF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp) :: gama
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  To Calculate the 3 Physical Velocity of the Comoving Frame with respect to the LNRF.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CFV3_LNRF(1) = this%expmu1 * this%ComFrameU4_BL(2) / &
                     this%ComFrameU4_BL(1) / this%expnu
      CFV3_LNRF(2) = this%expmu2 * this%ComFrameU4_BL(3) / &
                     this%ComFrameU4_BL(1) / this%expnu
      CFV3_LNRF(3) = this%exppsi *( this%ComFrameU4_BL(4) / &
                     this%ComFrameU4_BL(1) - this%somega ) / this%expnu
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  To Calculate the Matrix of CF's Tetrad in the LNRF Frame, which actually are 
!  the Lorentz Tramsformation Matrix.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      associate(  Vr => CFV3_LNRF(1), &
                  Vt => CFV3_LNRF(2), &
                  Vp => CFV3_LNRF(3) )

      gama = one / dsqrt(one - Vr**2 - Vt**2 - Vp**2)
      Lorentz_Matrix(1,1) = gama
      Lorentz_Matrix(1,2) = gama*Vr
      Lorentz_Matrix(1,3) = gama*Vt
      Lorentz_Matrix(1,4) = gama*Vp
      Lorentz_Matrix(2,1) = gama*Vr
      Lorentz_Matrix(2,2) = one + (gama*Vr)**2/(one + gama)
      Lorentz_Matrix(2,3) = gama**2*Vr*Vt/(one + gama)
      Lorentz_Matrix(2,4) = gama**2*Vr*Vp/(one + gama)
      Lorentz_Matrix(3,1) = gama*Vt
      Lorentz_Matrix(3,2) = gama**2*Vt*Vr/(one + gama)
      Lorentz_Matrix(3,3) = one + (gama*Vt)**2/(one + gama)
      Lorentz_Matrix(3,4) = gama**2*Vt*Vp/(one + gama)
      Lorentz_Matrix(4,1) = gama*Vp
      Lorentz_Matrix(4,2) = gama**2*Vp*Vr/(one + gama)
      Lorentz_Matrix(4,3) = gama**2*Vp*Vt/(one + gama)
      Lorentz_Matrix(4,4) = one + (gama*Vp)**2/(one + gama)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Lorentz_Inverse_Matrix(1,1) = gama
      Lorentz_Inverse_Matrix(1,2) = - gama*Vr
      Lorentz_Inverse_Matrix(1,3) = - gama*Vt
      Lorentz_Inverse_Matrix(1,4) = - gama*Vp
      Lorentz_Inverse_Matrix(2,1) = - gama*Vr
      Lorentz_Inverse_Matrix(2,2) = one + (gama*Vr)**2/(one + gama)
      Lorentz_Inverse_Matrix(2,3) = gama**2*Vr*Vt/(one + gama)
      Lorentz_Inverse_Matrix(2,4) = gama**2*Vr*Vp/(one + gama)
      Lorentz_Inverse_Matrix(3,1) = - gama*Vt
      Lorentz_Inverse_Matrix(3,2) = gama**2*Vt*Vr/(one + gama)
      Lorentz_Inverse_Matrix(3,3) = one + (gama*Vt)**2/(one + gama)
      Lorentz_Inverse_Matrix(3,4) = gama**2*Vt*Vp/(one + gama)
      Lorentz_Inverse_Matrix(4,1) = - gama*Vp
      Lorentz_Inverse_Matrix(4,2) = gama**2*Vp*Vr/(one + gama)
      Lorentz_Inverse_Matrix(4,3) = gama**2*Vp*Vt/(one + gama)
      Lorentz_Inverse_Matrix(4,4) = one + (gama*Vp)**2/(one + gama)
      END associate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   To obtain the Matrix of Covariant 4 Tetrad of LNRF Frame, i.e. e_(a)^\mu, where subindex
!   (a) is the base vector index of the tetrad.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      eLNRF(1,1) = one / this%expnu
      eLNRF(1,2) = zero
      eLNRF(1,3) = zero
      eLNRF(1,4) = this%somega / this%expnu
      eLNRF(2,1) = zero
      eLNRF(2,2) = one / this%expmu1
      eLNRF(2,3) = zero
      eLNRF(2,4) = zero
      eLNRF(3,1) = zero
      eLNRF(3,2) = zero
      eLNRF(3,3) = one / this%expmu2
      eLNRF(3,4) = zero
      eLNRF(4,1) = zero
      eLNRF(4,2) = zero
      eLNRF(4,3) = zero
      eLNRF(4,4) = one / this%exppsi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   To obtain the Matrix of Contravariant 4 Tetrad of LNRF Frame, i.e. e_\mu^(a), where the
!   super-index (a) is the base vector index of the tetrad.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Inverse_eLNRF(1,1) = this%expnu
      Inverse_eLNRF(1,2) = zero
      Inverse_eLNRF(1,3) = zero
      Inverse_eLNRF(1,4) = - this%somega * this%exppsi
      Inverse_eLNRF(2,1) = zero
      Inverse_eLNRF(2,2) = this%expmu1
      Inverse_eLNRF(2,3) = zero
      Inverse_eLNRF(2,4) = zero
      Inverse_eLNRF(3,1) = zero
      Inverse_eLNRF(3,2) = zero
      Inverse_eLNRF(3,3) = this%expmu2
      Inverse_eLNRF(3,4) = zero
      Inverse_eLNRF(4,1) = zero
      Inverse_eLNRF(4,2) = zero
      Inverse_eLNRF(4,3) = zero
      Inverse_eLNRF(4,4) = this%exppsi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     to obtain the Matrix_Of_BLCF_(a)^\mu, which can complete the transformation
!     of components of a 4-Vector from BL to CF, or CF to BL, according the following
!     indentities:
! (1) CF to BL, P^\mu(BL) = P^(a)(CF) * Matrix_Of_BLCF_(a)^\mu
! (2) BL to CF, P_(a)(CF) = Matrix_Of_BLCF_(a)^\mu * P_\mu(BL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL Matrix_Multiplication44X44_Sub( Lorentz_Matrix, &
                                    eLNRF, this%Matrix_Of_BLCF_a_mu )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     to obtain the Matrix_OF_BLCF_\mu^(a),which can complete the transformation
!     of components of a 4-Vector from BL to CF, or CF to BL, according the following
!     indentities:
! (1) CF to BL,  P^(a)(CF) = P^\mu(BL) * Matrix_OF_BLCF_\mu^(a)
! (2) BL to CF,  P_\mu(BL) = Matrix_OF_BLCF_\mu^(a) * P_(a)(CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call Matrix_Multiplication44X44_Sub( Inverse_eLNRF, &
                  Lorentz_Inverse_Matrix, this%Matrix_OF_BLCF_mu_a )

      this%Matrix_Of_CF2LNRF = Lorentz_Matrix

      END SUBROUTINE Set_Transformation_Matrices_Sub 

!*******************************************************************************************************
      SUBROUTINE Set_Transformation_Matrices_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp), dimension(1:3) :: CFV3_LNRF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp), dimension(1:4, 1:4) :: Lorentz_Matrix
      real(mcp), dimension(1:4, 1:4) :: Lorentz_Inverse_Matrix
      real(mcp), dimension(1:4, 1:4) :: eLNRF
      real(mcp), dimension(1:4, 1:4) :: Inverse_eLNRF
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp) :: gama
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  To Calculate the 3 Physical Velocity of the Comoving Frame with respect to the LNRF.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      CFV3_LNRF(1) = this%expmu1_p * this%Elec4U_CtrBL_At_p(2) / &
                     this%Elec4U_CtrBL_At_p(1) / this%expnu_p
      CFV3_LNRF(2) = this%expmu2_p * this%Elec4U_CtrBL_At_p(3) / &
                     this%Elec4U_CtrBL_At_p(1) / this%expnu_p
      CFV3_LNRF(3) = this%exppsi_p *( this%Elec4U_CtrBL_At_p(4) / &
                     this%Elec4U_CtrBL_At_p(1) - this%somega_p ) / this%expnu_p
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  To Calculate the Matrix of CF's Tetrad in the LNRF Frame, which actually are 
!  the Lorentz Tramsformation Matrix.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      associate(  Vr => CFV3_LNRF(1), &
                  Vt => CFV3_LNRF(2), &
                  Vp => CFV3_LNRF(3) )

      if (Vr**2 + Vt**2 + Vp**2 >= one)then
          write(*,*)'V is great than 1 ===', Vr**2 + Vt**2 + Vp**2 ,this%Elec4U_CtrBL_At_p, &
           this%g_tt_p, this%g_rr_p, this%g_pp_p
      endif 
      gama = one / dsqrt(one - Vr**2 - Vt**2 - Vp**2)
      Lorentz_Matrix(1,1) = gama
      Lorentz_Matrix(1,2) = gama*Vr
      Lorentz_Matrix(1,3) = gama*Vt
      Lorentz_Matrix(1,4) = gama*Vp
      Lorentz_Matrix(2,1) = gama*Vr
      Lorentz_Matrix(2,2) = one + (gama*Vr)**2/(one + gama)
      Lorentz_Matrix(2,3) = gama**2*Vr*Vt/(one + gama)
      Lorentz_Matrix(2,4) = gama**2*Vr*Vp/(one + gama)
      Lorentz_Matrix(3,1) = gama*Vt
      Lorentz_Matrix(3,2) = gama**2*Vt*Vr/(one + gama)
      Lorentz_Matrix(3,3) = one + (gama*Vt)**2/(one + gama)
      Lorentz_Matrix(3,4) = gama**2*Vt*Vp/(one + gama)
      Lorentz_Matrix(4,1) = gama*Vp
      Lorentz_Matrix(4,2) = gama**2*Vp*Vr/(one + gama)
      Lorentz_Matrix(4,3) = gama**2*Vp*Vt/(one + gama)
      Lorentz_Matrix(4,4) = one + (gama*Vp)**2/(one + gama)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Lorentz_Inverse_Matrix(1,1) = gama
      Lorentz_Inverse_Matrix(1,2) = - gama*Vr
      Lorentz_Inverse_Matrix(1,3) = - gama*Vt
      Lorentz_Inverse_Matrix(1,4) = - gama*Vp
      Lorentz_Inverse_Matrix(2,1) = - gama*Vr
      Lorentz_Inverse_Matrix(2,2) = one + (gama*Vr)**2/(one + gama)
      Lorentz_Inverse_Matrix(2,3) = gama**2*Vr*Vt/(one + gama)
      Lorentz_Inverse_Matrix(2,4) = gama**2*Vr*Vp/(one + gama)
      Lorentz_Inverse_Matrix(3,1) = - gama*Vt
      Lorentz_Inverse_Matrix(3,2) = gama**2*Vt*Vr/(one + gama)
      Lorentz_Inverse_Matrix(3,3) = one + (gama*Vt)**2/(one + gama)
      Lorentz_Inverse_Matrix(3,4) = gama**2*Vt*Vp/(one + gama)
      Lorentz_Inverse_Matrix(4,1) = - gama*Vp
      Lorentz_Inverse_Matrix(4,2) = gama**2*Vp*Vr/(one + gama)
      Lorentz_Inverse_Matrix(4,3) = gama**2*Vp*Vt/(one + gama)
      Lorentz_Inverse_Matrix(4,4) = one + (gama*Vp)**2/(one + gama)
      END associate
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   To obtain the Matrix of Covariant 4 Tetrad of LNRF Frame, i.e. e_(a)^\mu, where subindex
!   (a) is the base vector index of the tetrad.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      eLNRF(1,1) = one / this%expnu_p
      eLNRF(1,2) = zero
      eLNRF(1,3) = zero
      eLNRF(1,4) = this%somega_p / this%expnu_p
      eLNRF(2,1) = zero
      eLNRF(2,2) = one / this%expmu1_p
      eLNRF(2,3) = zero
      eLNRF(2,4) = zero
      eLNRF(3,1) = zero
      eLNRF(3,2) = zero
      eLNRF(3,3) = one / this%expmu2_p
      eLNRF(3,4) = zero
      eLNRF(4,1) = zero
      eLNRF(4,2) = zero
      eLNRF(4,3) = zero
      eLNRF(4,4) = one / this%exppsi_p 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   To obtain the Matrix of Contravariant 4 Tetrad of LNRF Frame, i.e. e_\mu^(a), where the
!   super-index (a) is the base vector index of the tetrad.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Inverse_eLNRF(1,1) = this%expnu_p
      Inverse_eLNRF(1,2) = zero
      Inverse_eLNRF(1,3) = zero
      Inverse_eLNRF(1,4) = - this%somega_p * this%exppsi_p
      Inverse_eLNRF(2,1) = zero
      Inverse_eLNRF(2,2) = this%expmu1_p
      Inverse_eLNRF(2,3) = zero
      Inverse_eLNRF(2,4) = zero
      Inverse_eLNRF(3,1) = zero
      Inverse_eLNRF(3,2) = zero
      Inverse_eLNRF(3,3) = this%expmu2_p
      Inverse_eLNRF(3,4) = zero
      Inverse_eLNRF(4,1) = zero
      Inverse_eLNRF(4,2) = zero
      Inverse_eLNRF(4,3) = zero
      Inverse_eLNRF(4,4) = this%exppsi_p
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     to obtain the Matrix_Of_BLCF_(a)^\mu, which can complete the transformation
!     of components of a 4-Vector from BL to CF, or CF to BL, according the following
!     indentities:
! (1) CF to BL, P^\mu(BL) = P^(a)(CF) * Matrix_Of_BLCF_(a)^\mu
! (2) BL to CF, P_(a)(CF) = Matrix_Of_BLCF_(a)^\mu * P_\mu(BL)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL Matrix_Multiplication44X44_Sub( Lorentz_Matrix, &
                                    eLNRF, this%Matrix_Of_BLCF_a_mu_p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     to obtain the Matrix_OF_BLCF_\mu^(a),which can complete the transformation
!     of components of a 4-Vector from BL to CF, or CF to BL, according the following
!     indentities:
! (1) CF to BL,  P^(a)(CF) = P^\mu(BL) * Matrix_OF_BLCF_\mu^(a)
! (2) BL to CF,  P_\mu(BL) = Matrix_OF_BLCF_\mu^(a) * P_(a)(CF)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call Matrix_Multiplication44X44_Sub( Inverse_eLNRF, &
                  Lorentz_Inverse_Matrix, this%Matrix_OF_BLCF_mu_a_p )

      this%Matrix_Of_CF2LNRF_p = Lorentz_Matrix

      END SUBROUTINE Set_Transformation_Matrices_At_p_Sub

!*******************************************************************************************************
      subroutine get_lambdaq_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this  
      real(mcp) :: tempA, x, y

      !tempA = this%Phot4k_LNRF(4) / ( - dsqrt(this%Delta) * this%Sigma / &
      !            this%bigA*this%Phot4k_LNRF(1) + &
      !            this%Phot4k_LNRF(4) * this%somega * this%musin )
      !this%lambda = tempA * this%musin
      !this%q = ( tempA**2 - this%aspin**2 ) * this%mucos**2 + &
      !         ( this%Phot4k_LNRF(3) / this%Phot4k_LNRF(1) * &
      !         (one - this%lambda * this%somega) )**2 * this%bigA / this%Delta
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'cttttttss === ', this%lambda, this%q, this%Phot4k_LNRF
      !write(unit = *, fmt = *)'************************************************************' 
      associate( expnu  => this%expnu, &
                 exppsi => this%exppsi, &
                 expmu1 => this%expmu1, &
                 expmu2 => this%expmu2, &
                 somega => this%somega )
      x = this%Phot4k_LNRF(4) / this%Phot4k_LNRF(1)
      y = this%Phot4k_LNRF(3) / this%Phot4k_LNRF(1)
      this%lambda = x / ( - expnu / exppsi + x * somega )
      this%q = ( y*( one - this%lambda * somega ) * expmu2 / expnu )**2 + &
                 ( ( this%lambda / this%musin )**2 - this%aspin**2 ) * this%mucos**2
      END associate
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'cttttttss2222 === ', this%lambda, this%q
      !write(unit = *, fmt = *)'************************************************************'
      

      end subroutine get_lambdaq_Sub

!*******************************************************************************************************
      subroutine get_lambdaq_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp) :: tempA, x, y

      !tempA = this%Phot4k_LNRF(4) / ( - dsqrt(this%Delta) * this%Sigma / this%bigA*this%Phot4k_LNRF(1) + &
      !            this%Phot4k_LNRF(4) * this%somega * this%musin_p )
      !this%lambda = tempA * this%musin_p
      !this%q = ( tempA**2 - this%aspin**2 ) * this%mucos_p**2 + &
      !         ( this%Phot4k_LNRF(3) / this%Phot4k_LNRF(1) * &
      !         (one - this%lambda * this%somega) )**2 * this%bigA / this%Delta 
      associate( expnu  => this%expnu_p, &
                 exppsi => this%exppsi_p, &
                 expmu1 => this%expmu1_p, &
                 expmu2 => this%expmu2_p, &
                 somega => this%somega_p )
      x = this%Phot4k_LNRF_ini(4)/this%Phot4k_LNRF_ini(1)
      y = this%Phot4k_LNRF_ini(3)/this%Phot4k_LNRF_ini(1)
      this%lambda_p = x / ( - expnu / exppsi + x * somega )
      this%q_p = ( y*( one - this%lambda_p * somega ) * expmu2 / expnu )**2 + &
               ( ( this%lambda_p / this%musin_p )**2 - this%aspin**2 ) * this%mucos_p**2
      END associate
      
      end subroutine get_lambdaq_At_p_Sub

!*******************************************************************************************************
      subroutine Get_Phot_BLCoordinates_At_p_Sub(this, p)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), intent(in) :: p
      !real(mcp) :: sign_pth, sign_pr
 
      this%r_p = radius( p, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
                             this%q_ini, this%aspin, this%r_ini, one, this%sign_pr_p )
      this%mucos_p = mucos( p, this%Phot4k_LNRF_ini(4), &
             this%Phot4k_LNRF_ini(3), this%lambda_ini, &
             this%q_ini, this%musin_ini, this%mucos_ini, this%aspin, one, this%sign_pth_p )
      this%theta_p = DACOS( this%mucos_p )
      this%musin_p = DSQRT( one - this%mucos_p**2 )
      this%phi_p = zero
      this%t_p = zero

      end subroutine Get_Phot_BLCoordinates_At_p_Sub

!*******************************************************************************************************
      subroutine Get_Phot_BLCoor_At_samp_Sub(this, p)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), intent(in) :: p
      real(mcp), dimension(1:4) :: f1234
      real(mcp) :: rdp, mudp, phip, timep, sigmap 
 
      this%r_samp = radius(p,this%Phot4k_LNRF_ini(2), &
                             this%lambda_ini, &
                             this%q_ini, &
                             this%aspin, &
                             this%r_ini, one, this%sign_pr_samp)
      this%mucos_samp = mucos(p, this%Phot4k_LNRF_ini(4), &
                                 this%Phot4k_LNRF_ini(3), &
                                 this%lambda_ini, &
                                 this%q_ini,  &
                                 this%musin_ini, &
                                 this%mucos_ini, &
                                 this%aspin, one, this%sign_pth_samp)
      this%p_samp = p
      !this%r_samp = radius(p, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !               this%q_ini, this%aspin, this%r_ini, one)
      !this%mucos_samp = mucos(p, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
      !                  this%lambda_ini, this%q_ini, this%musin_ini, &
      !                         this%mucos_ini, this%aspin, one, sign_pth)

      !f1234(1) = this%Phot4k_LNRF(2)
      !f1234(2) = this%Phot4k_LNRF(3)
      !f1234(3) = this%Phot4k_LNRF(4)
      !f1234(4) = this%Phot4k_LNRF(1)
      !call YNOGK(p,f1234,this%lambda,this%q,this%musin,this%mucos,this%aspin,this%r,one,&
      !                  rdp,mudp,phip,timep,sigmap) 

      !this%r_samp = this%r_samp
      !this%mucos_samp = mudp 
      this%theta_samp = DACOS( this%mucos_samp )
      this%musin_samp = DSQRT( one - this%mucos_samp**2 )
      this%phi_samp = zero
      this%t_samp = zero 

    !write(unit = *, fmt = *)'************************************************************' 
    !write(unit = *, fmt = *)'**hhhhhhhh === ', p, this%r_samp
    !write(unit = *, fmt = *)'**hhhhhhhh === ',p, this%mucos_samp, mudp
    !write(unit = *, fmt = *)'************************************************************' 
     ! write(unit = *, fmt = *)'************************************************************'
     !write(unit = *, fmt = *)'*====================ddddddd*', this%r_ini - this%r, &
    !this%theta_ini - this%theta, &
    !this%mucos_ini - this%mucos,&
    !this%musin_ini - this%musin,&
    !this%phi_ini - this%phi,&
    !this%t_ini - this%t,&
    !this%Phot4k_LNRF_ini - this%Phot4k_LNRF, &
    !this%lambda_ini - this%lambda,&
    !this%q_ini - this%q 
     ! write(unit = *, fmt = *)'************************************************************'

      end subroutine Get_Phot_BLCoor_At_samp_Sub

!*******************************************************************************************************
      subroutine Get_Phot4k_CovBL_At_p_Sub(this, p)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), intent(in) :: p
      real(mcp) :: R_p, Theta_p, Delta_p, dp, rdp, mudp
      !real(mcp) :: sign_pth, sign_pr

      associate( rp => this%r_p, &
                 ap => this%aspin, &
                 mucosp => this%mucos_p, &
                 musinp => this%musin_p, &
                 lp => this%lambda_ini, &
                 qp => this%q_ini )
      Delta_p = rp*rp - two*rp + ap**2
      R_p = rp**4 - ( qp + lp**2 - ap**2 )*rp**2 + &
            two*( qp + (lp - ap)**2 )*rp - ap**2*qp
      if (lp .NE. zero) then
          Theta_p = qp + ( ap*mucosp )**2 - ( lp * mucosp / musinp )**2
      else
          Theta_p = qp + ( ap*mucosp )**2
      endif
      !dp = p*1.D-5
      !rdp = radius(p+dp,this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                         this%q_ini, this%aspin, this%r_ini, one, sign_pr)
      !mudp = mucos(p+dp, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
      !                   this%lambda_ini, this%q_ini, this%musin_ini, &
      !                   this%mucos_ini, this%aspin, one, sign_pth)
      If (R_p < zero  ) then
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'R_p ff === ', rdp - rp , R_p , DSQRT(dabs(R_p)), rp
          write(unit = *, fmt = *)'************************************************************' 
           R_p = - R_p
      endif
      If (Theta_p < zero )then! .and. dabs(Theta_p)<1.D-9) then
          write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', rdp ,rp ,R_p, Delta_p, Theta_p, musinp, mucosp, lp, qp
          write(*,*)'rdp sssss222 === ', Theta_p, DSQRT(dabs(Theta_p)), p, lp, qp
          write(unit = *, fmt = *)'************************************************************' 
           Theta_p = -Theta_p
      endif
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', rdp ,rp ,R_p, Delta_p, Theta_p, musinp, mucosp, lp, qp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', qp - ( qp + lp**2 -ap**2 ) * mucosp**2 - mucosp**4
      !write(unit = *, fmt = *)'************************************************************'
      this%Phot4k_CovBL_At_p(2) = this%sign_pr_p * DSQRT(R_p) / Delta_p
      !this%Phot4k_CovBL_At_p(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
      this%Phot4k_CovBL_At_p(3) = this%sign_pth_p * DSQRT(Theta_p)
      this%Phot4k_CovBL_At_p(4) = lp
      this%Phot4k_CovBL_At_p(1) = - one
      this%Phot4k_CovBL_At_p = this%Phot4k_CovBL_At_p * this%E_ini
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'Phot_Covariant4k_BL_At_p****', this%Phot_Covariant4k_BL_At_p
  !write(unit = *, fmt = *)'Phot_Covariant4k_BL_At_p****', DSIGN(one, rdp - rp), - DSIGN(one, mudp - mucosp)  
      !write(unit = *, fmt = *)'************************************************************' 
      end associate
      end subroutine Get_Phot4k_CovBL_At_p_Sub
 
!*******************************************************************************************************
      subroutine Get_Phot4k_CtrBL_At_p_Sub(this, p)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), intent(in), optional :: p
      real(mcp) :: R_p, Theta_p, Delta_p, bigT, dp, rdp, mudp
      real(mcp) :: a, l, q, rp, mucosp, musinp 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      associate( rp => this%r_p, &
                 a => this%aspin, &
                 mucosp => this%mucos_p, &
                 musinp => this%musin_p, &
                 l => this%lambda_ini, &
                 q => this%q_ini )
      Delta_p = rp**2 - two*rp + a**2
      bigT = rp**2 + a**2 - l*a
      R_p = rp**4 - ( q + l**2 - a**2 )*rp**2 + &
            two*( q + (l - a)**2 )*rp - a**2*q
      if (l .NE. zero) then
          Theta_p = q + ( a*mucosp )**2 - ( l * mucosp / musinp )**2
          this%Phot4k_CtrBL_At_p(4) = - (a - l / musinp**2) + a * bigT / Delta_p
      else
          Theta_p = q + ( a*mucosp )**2
          this%Phot4k_CtrBL_At_p(4) = - a + a * bigT / Delta_p
      endif
      !dp = p*1.D-5
      !rdp = radius(p+dp, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                   this%q_ini, this%aspin, this%r_ini, one, sign_pr)
      !mudp = mucos(p+dp, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
      !                   this%lambda_ini, this%q_ini, this%musin_ini, &
      !                   this%mucos_ini, this%aspin, one, sign_pth)     
      If (R_p < zero  ) then
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'R_p ff === ', rdp - rp , R_p , DSQRT(dabs(R_p)), rp 
          write(unit = *, fmt = *)'************************************************************'
           R_p = - R_p
      endif
      If (Theta_p < zero )then! .and. dabs(Theta_p)<1.D-9) then
          write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', rdp ,rp ,R_p, Delta_p, Theta_p, musinp, mucosp, lp, qp
          write(*,*)'rdp sssss222 === ', Theta_p, DSQRT(dabs(Theta_p)), p, l, q 
          write(unit = *, fmt = *)'************************************************************' 
           Theta_p = -Theta_p
      endif  
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'Get_Phot4k_CtrBL_At_p_Sub = ', this%sign_pr_p, DSQRT(R_p), this%sign_pth_p, DSQRT(Theta_p)
          !- DSIGN(one, mudp - mucosp), sign_pth
          !write(unit = *, fmt = *)'************************************************************'
!*******************************************************************************************************
      this%Phot4k_CtrBL_At_p(1) = - a*(a*musinp**2 - l) + (rp**2 + a**2) * bigT / Delta_p
      this%Phot4k_CtrBL_At_p(2) = this%sign_pr_p * DSQRT(R_p)
      !this%Phot4k_CtrBL_At_p(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
      this%Phot4k_CtrBL_At_p(3) = this%sign_pth_p * DSQRT(Theta_p)
      this%Phot4k_CtrBL_At_p = this%Phot4k_CtrBL_At_p / ( rp**2 + (a*mucosp)**2 ) * this%E_ini
!*******************************************************************************************************  
      end associate
      end subroutine Get_Phot4k_CtrBL_At_p_Sub

!*******************************************************************************************************
      subroutine Get_Phot4k_CovBL_At_samp_Sub(this, p)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), intent(in) :: p
      real(mcp) :: R_p, Theta_p, Delta_p, dp, rdp, mudp
      real(mcp) :: phip,timep,sigmap
      real(mcp), dimension(1:4) :: f1234
      real(mcp) :: mu_tp1, mu_tp2, f12341, r_tp1, r_tp2 
      logical :: mobseqmtp,robs_eq_rtp,indrhorizon
      integer :: reals, cases
      complex*16 r1,r2,r3,r4
      complex*16 bb(4)
      !********************************************************************* 

      associate( rp => this%r_samp, &
                 a => this%aspin, &
                 mucosp => this%mucos_samp, &
                 musinp => this%musin_samp, &
                 l => this%lambda_ini, &
                 q => this%q_ini )
      Delta_p = rp*rp - two*rp + a**2
      R_p = rp**4 - ( q + l**2 - a**2 )*rp**2 + two*( q + (l - a)**2 )*rp - a**2*q
      if (l .NE. zero) then
          Theta_p = q + ( a*mucosp )**2 - ( l * mucosp / musinp )**2
      else
          Theta_p = q + ( a*mucosp )**2
      endif
      !dp = p*1.D-5
      !rdp = radius(p+dp, this%Phot4k_LNRF_ini(2), this%lambda_ini, &
      !                      this%q_ini, this%aspin, this%r_ini, one, sign_pr)
      !mudp = mucos(p+dp, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
      !                   this%lambda_ini, this%q_ini, this%musin_ini, &
      !                   this%mucos_ini, this%aspin, one, sign_pth)
      ! f12341 = this%Phot4k_LNRF_ini(2)
      !f1234(1) = this%Phot4k_LNRF_ini(2)
      !f1234(2) = this%Phot4k_LNRF_ini(3)
      !f1234(3) = this%Phot4k_LNRF_ini(4)
      !f1234(4) = this%Phot4k_LNRF_ini(1)
      !call YNOGK(p+dp,f1234,this%lambda_ini,this%q_ini,this%musin_ini,&
      !           this%mucos_ini,this%aspin,this%r_ini,one,&
      !                  rdp,mudp,phip,timep,sigmap,sign_pr,sign_pth) 
      !call radiustp( this%Phot4k_LNRF_ini(2),this%aspin,this%r_ini,this%lambda_ini,this%q_ini,r_tp1,&
      !                  r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
      !call mutp(f1234(2),f1234(3),this%musin_ini, this%mucos_ini,this%aspin,this%lambda_ini,&
      !                    this%q_ini,mu_tp1,mu_tp2,reals,mobseqmtp)
      !call root4(zero, ( this%q_ini + this%lambda_ini**2 - this%aspin**2 ), zero, &
      !                                  -this%q_ini ,r1,r2,r3,r4,reals)
    
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt= ', - DSIGN(one, mudp - mucosp), sign_pth
          !If (- DSIGN(one, mudp - mucosp) /= sign_pth) then
              !write(*,*)'sssssstttttt= ', mu_tp1, mu_tp2, mucosp, mudp
              !write(*,*)'mussss===', mucos(p, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
              !           this%lambda_ini, this%q_ini, this%musin_ini, &
              !           this%mucos_ini, this%aspin, -two, sign_pth)
              !stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'
          !write(unit = *, fmt = *)'************************************************************'
          !write(*,*)'sssssssssssssssssstttttt= ',DSIGN(one, rdp - rp), this%sign_pr_samp
          !If (DSIGN(one, rdp - rp) /= this%sign_pr_samp) then
          !    write(*,*)'sssssstttttt= ', r_tp1, r_tp2, rp, rdp, this%r_ini
              !write(*,*)'mussss===', mucos(p, this%Phot4k_LNRF_ini(4), this%Phot4k_LNRF_ini(3), &
              !           this%lambda_ini, this%q_ini, this%musin_ini, &
              !           this%mucos_ini, this%aspin, -two, sign_pth)
              !stop
          !endif
          !write(unit = *, fmt = *)'************************************************************'

      If (R_p < zero  ) then
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'R_p ff === ', rdp - rp , R_p , DSQRT(dabs(R_p)), rp
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'r_ini  P_r === ', this%r_ini, this%Phot4k_LNRF_ini(2), r_tp1, r_tp2
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'r_tp1, r_tp2 === ',&
          r_tp1**4 - ( q + l**2 - a**2 )*r_tp1**2 + two*( q + (l - a)**2 )*r_tp1 - a**2*q, & 
          r_tp2**4 - ( q + l**2 - a**2 )*r_tp2**2 + two*( q + (l - a)**2 )*r_tp2 - a**2*q
           R_p = - R_p
      endif
      If (Theta_p < zero )then! .and. dabs(Theta_p)<1.D-9) then
          write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', rdp ,rp ,R_p, Delta_p, Theta_p, musinp, mucosp, lp, qp
          write(*,*)'rdp sssss222 === ', Theta_p, DSQRT(dabs(Theta_p)), p, l, q
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'rdp sssss222 === ', rdp - rp , R_p  , Delta_p, this%r_ini, r_tp1, r_tp2
          write(unit = *, fmt = *)'************************************************************'
          write(*,*)'rdp sssss222 === ', mudp - mucosp , Theta_p 
          write(unit = *, fmt = *)'************************************************************' 
          write(*,*)'rdp sssss === ', this%mucos_ini, mucosp, mu_tp1, mu_tp2, dabs(mucosp) - dabs(mu_tp1)
          write(unit = *, fmt = *)'************************************************************' 
          !write(*,*)'kkkkkkkkkk === ', this%q_ini - ( this%q_ini + this%lambda_ini**2 - this%aspin**2 ) * &
          !                          mu_tp1**2 - this%aspin**2 * mu_tp1**4, reals
          write(unit = *, fmt = *)'************************************************************' 
          write(*,*)'rdp sssss === ',this%Phot4k_LNRF_ini, &
                                     this%Phot4k_LNRF_ini(3)/this%Phot4k_LNRF_ini(1)
           Theta_p = - Theta_p
      endif
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ', this%mucos, mucosp, mu_tp1, dabs(mucosp) - dabs(mu_tp1)
      !write(unit = *, fmt = *)'************************************************************'
      !write(*,*)'rdp sssss === ', r1,r2,r3,r4, dabs(mu_tp1)
      !write(unit = *, fmt = *)'************************************************************'
      !write(*,*)'rdp sssss === ', this%r_ini, rp,rdp, this%mucos_ini, mucosp
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'kkkkkkkkkk === ', this%q - ( this%q + this%lambda**2 - this%aspin**2 ) * &
      !                             mu_tp1**2 - this%aspin**2 * mu_tp1**4
      !write(unit = *, fmt = *)'************************************************************'
      this%Phot4k_CovBL_At_samp(1) = - one
      this%Phot4k_CovBL_At_samp(2) = this%sign_pr_samp * DSQRT(R_p) / Delta_p
      !this%Phot4k_CovBL_At_samp(3) = - DSIGN(one, mudp - mucosp) * DSQRT(Theta_p)
      this%Phot4k_CovBL_At_samp(3) = this%sign_pth_samp * DSQRT(Theta_p)
      this%Phot4k_CovBL_At_samp(4) = l
      this%Phot4k_CovBL_At_samp = this%Phot4k_CovBL_At_samp * this%E_ini
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(*,*)'rdp sssss === ',this%Phot4k_LNRF, this%Phot4k_LNRF(3)/this%Phot4k_LNRF(1)
      !write(unit = *, fmt = *)'Phot_Covariant4k_BL_At_p****', this%Phot_Covariant4k_BL_At_p
  !write(unit = *, fmt = *)'Phot_Covariant4k_BL_At_p****', DSIGN(one, rdp - rp), - DSIGN(one, mudp - mucosp)  
      !write(unit = *, fmt = *)'************************************************************' 
      end associate
      end subroutine Get_Phot4k_CovBL_At_samp_Sub

 
!*******************************************************************************************************
      subroutine get_K1_K2_pw_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), dimension(1:4) :: detaf, gamaf
!****************************************************************************

      detaf(1) = - this%r * this%Phot4k_BL(2) - &
                   this%aspin**2 * this%mucos * &
                   this%musin * this%Phot4k_BL(3)
      detaf(2) = this%r * this%Phot4k_BL(1) - &
                 this%r * this%aspin * &
                 this%musin**2 * this%Phot4k_BL(4)
      detaf(3) = this%aspin * this%mucos * this%musin * ( &
                 this%aspin * this%Phot4k_BL(1) - &
                (this%r**2 + this%aspin**2) * this%Phot4k_BL(4) )
      detaf(4) = this%aspin * this%r * this%musin**2 * this%Phot4k_BL(2) + &
                (this%r**2 + this%aspin**2) * this%aspin * &
                 this%mucos * this%musin * this%Phot4k_BL(3)
      !***********************************************************************************
      gamaf(1) = - this%aspin * this%mucos * this%Phot4k_BL(2) + &
                   this%aspin * this%r * this%musin * this%Phot4k_BL(3)
      gamaf(2) = this%aspin * this%mucos * this%Phot4k_BL(1) - &
                 this%aspin**2 * this%musin**2 * this%mucos * this%Phot4k_BL(4)
      gamaf(3) = (this%r**2 + this%aspin**2) * this%r * this%musin * this%Phot4k_BL(4) - &
                 this%aspin * this%r * this%musin * this%Phot4k_BL(1)
      gamaf(4) = this%aspin**2 * this%mucos * this%musin**2 * this%Phot4k_BL(2) - &
                (this%r**2 + this%aspin**2) * this%r * this%musin * this%Phot4k_BL(3)
      !***********************************************************************************
      this%K1_pw_ini = this%f4_CtrBL(1) * detaf(1) + this%f4_CtrBL(2) * detaf(2) + &
                   this%f4_CtrBL(3) * detaf(3) + this%f4_CtrBL(4) * detaf(4)
      this%K2_pw_ini = this%f4_CtrBL(1) * gamaf(1) + this%f4_CtrBL(2) * gamaf(2) + &
                   this%f4_CtrBL(3) * gamaf(3) + this%f4_CtrBL(4) * gamaf(4)
      end subroutine get_K1_K2_pw_Sub

!*******************************************************************************************************
      subroutine get_K1_K2_pw_At_ini_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp), dimension(1:4) :: detaf, gamaf
!****************************************************************************

      associate( a  => this%aspin, &
                 r  => this%r_ini, &
                 mu => this%mucos_ini, &
                 su => this%musin_ini, &
                 k0 => this%Phot4k_CtrBL_ini(1), &
                 k1 => this%Phot4k_CtrBL_ini(2), &
                 k2 => this%Phot4k_CtrBL_ini(3), &
                 k3 => this%Phot4k_CtrBL_ini(4) )
      detaf(1) = - r * k1 - a**2 * mu * su * k2
      detaf(2) = r * k0 - r * a * su**2 * k3
      detaf(3) = a * mu * su * ( a * k0 - ( r**2 + a**2 ) * k3 )
      detaf(4) = a * r * su**2 * k1 + ( r**2 + a**2 ) * a * mu * su * k2
      !*********************************************************************************** 
      gamaf(1) = - a * mu * k1 + a * r * su * k2
      gamaf(2) = a * mu * k0 - a**2 * su**2 * mu * k3
      gamaf(3) = ( r**2 + a**2 ) * r * su * k3 - a * r * su * k0
      gamaf(4) = a**2 * mu * su**2 * k1 - ( r**2 + a**2 ) * r * su * k2
      !***********************************************************************************
      this%K1_pw_ini = this%f4_CtrBL(1) * detaf(1) + this%f4_CtrBL(2) * detaf(2) + &
                   this%f4_CtrBL(3) * detaf(3) + this%f4_CtrBL(4) * detaf(4)
      this%K2_pw_ini = this%f4_CtrBL(1) * gamaf(1) + this%f4_CtrBL(2) * gamaf(2) + &
                   this%f4_CtrBL(3) * gamaf(3) + this%f4_CtrBL(4) * gamaf(4)
      END associate

    If ( this%my_test ) then
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******kkkkkkkkkk1 === ', this%Phot4k_CtrBL_ini
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******kkkkkkkkkk2 === ', detaf, gamaf
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******kkkkkkkkkk3 === ', this%Phot4k_CovBL_At_p, &
                     this%r_ini, this%mucos_ini,  this%musin_ini, this%aspin
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******kkkkkkkkkk4 === ',this%f4_CtrBL
      write(unit = *, fmt = *)'************************************************************'
   endif
      end subroutine get_K1_K2_pw_At_ini_Sub

!*******************************************************************************************************
      subroutine get_f4_Ctr_BL_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      real(mcp) :: df1, df2, df3, gf1, gf2, gf3, N_temp, N1, N2, N3, Normal_f4
      !****************************************************************************

      associate( r  => this%r_p, &
                 mu => this%mucos_p, &
                 su => this%musin_p, &
                 a  => this%aspin, &
                 k0 => this%Phot4k_CovBL_At_p(1), &
                 k1 => this%Phot4k_CovBL_At_p(2), &
                 k2 => this%Phot4k_CovBL_At_p(3), &
                 k3 => this%Phot4k_CovBL_At_p(4), &
                 K1pw => this%K1_pw_ini, &
                 K2pw => this%K2_pw_ini )
 
      df1 = r * k0 - r * a * su**2 * k3
      df2 = a * mu * su * (  a * k0 - (r**2 + a**2) * k3 )
      df3 = a * r * su**2 * k1 + (r**2 + a**2) * a * mu * su * k2
      !*********************************************************************************** 
      gf1 = a * mu * k0 - a**2 * su**2 * mu * k3
      gf2 = (r**2 + a**2) * r * su * k3 - a * r * su * k0
      gf3 = a**2 * mu * su**2 * k1 - (r**2 + a**2) * r * su * k2
      !***********************************************************************************
      N_temp = df1 * gf2 * k3 + df2 * gf3 * k1 + df3 * gf1 * k2 - &
               df3 * gf2 * k1 - df2 * gf1 * k3 - df1 * gf3 * k2
      N1 = - k2 * (K1pw * gf3 - df3 * K2pw) + k3 * (K1pw * gf2 - K2pw * df2)
      N2 =   k1 * (K1pw * gf3 - df3 * K2pw) + k3 * (K2pw * df1 - K1pw * gf1)
      N3 =   k1 * (K2pw * df2 - gf2 * K1pw) - k2 * (K2pw * df1 - K1pw * gf1)
      end associate

      this%f4_CtrBL_At_p(1) = zero
      this%f4_CtrBL_At_p(2) = N1 / N_temp
      this%f4_CtrBL_At_p(3) = N2 / N_temp
      this%f4_CtrBL_At_p(4) = N3 / N_temp

      this%f4_CovBL_At_p(1) = this%g_tt_p * this%f4_CtrBL_At_p(1) + &
                              this%g_tp_p * this%f4_CtrBL_At_p(4)
      this%f4_CovBL_At_p(2) = this%g_rr_p * this%f4_CtrBL_At_p(2)
      this%f4_CovBL_At_p(3) = this%g_thth_p * this%f4_CtrBL_At_p(3)
      this%f4_CovBL_At_p(4) = this%g_pt_p * this%f4_CtrBL_At_p(1) + &
                              this%g_pp_p * this%f4_CtrBL_At_p(4)

      Normal_f4 = DSQRT( Vector4D_Inner_Product( this%f4_CovBL_At_p, &
                         this%f4_CtrBL_At_p ) )

      this%f4_CtrBL_At_p = this%f4_CtrBL_At_p / Normal_f4
      this%f4_CovBL_At_p = this%f4_CovBL_At_p / Normal_f4
    If ( this%my_test ) then
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******ssssssssssssdddddd1 === ', this%Phot4k_CovBL_At_p,&
                  this%f4_CtrBL_At_p(1) * this%Phot4k_CovBL_At_p(1) + &
                  this%f4_CtrBL_At_p(2) * this%Phot4k_CovBL_At_p(2) + &
                  this%f4_CtrBL_At_p(3) * this%Phot4k_CovBL_At_p(3) + &
                  this%f4_CtrBL_At_p(4) * this%Phot4k_CovBL_At_p(4)
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******ssssssssssssdddddd2 === ', N1,N2,N3,N_temp, Normal_f4
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******ssssssssssssdddddd3 === ', this%Phot4k_CovBL_At_p, &
                     this%r_p, this%mucos_p,  this%musin_p, this%aspin
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******ssssssssssssdddddd4 === ', this%K1_pw_ini, this%K2_pw_ini
      write(unit = *, fmt = *)'************************************************************' 
      write(unit = *, fmt = *)'******ssssssssssssdddddd5 === ', this%g_tt_p, this%g_rr_p, this%g_thth_p,&
                          this%g_pp_p,this%g_tp_p,this%g_pt_p
      write(unit = *, fmt = *)'************************************************************'
   endif
      
      end subroutine get_f4_Ctr_BL_At_p_Sub



!*******************************************************************************************************
      subroutine Set_Elec4U_CtrBL_At_p_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      real(mcp) :: Omega_K, Omega
  
      associate( r => this%r_p, &
                 theta => this%theta_p )
      Omega_K = one / ( this%aspin + r**1.5D0 ) !Kepler velocity
      Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
                    ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega_p
      this%Elec4U_CtrBL_At_p(1) = one / DSQRT( - ( this%g_tt_p + two * this%g_tp_p * Omega + &
                                                     this%g_pp_p * Omega**2 ) )
      end associate
  
      this%Elec4U_CtrBL_At_p(2) = zero
      this%Elec4U_CtrBL_At_p(3) = zero
      this%Elec4U_CtrBL_At_p(4) = this%Elec4U_CtrBL_At_p(1) * Omega
      end subroutine Set_Elec4U_CtrBL_At_p_Sub

!*******************************************************************************************************
      subroutine Set_Elec4U_CtrBL_At_samp_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this 
      !real(mcp), intent(in), optional :: p
      real(mcp) :: Omega_K, Omega, r, theta
      !**************************************************************
 
      r = this%r_samp
      theta = this%theta_samp 
      Omega_K = one / ( this%aspin + r**1.5D0 ) !Kepler velocity
      Omega = ( theta / (pi / two) )**( one / three ) * Omega_K + &
              ( one - ( theta / (pi / two) )**( one / three ) ) * this%somega_samp
      this%Elec4U_CtrBL_At_samp(1) = one / DSQRT( - ( this%g_tt_samp + two * this%g_tp_samp * Omega + &
                                                     this%g_pp_samp * Omega**2 ) )
      
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**ss Set_Elec4U_BL === ',this%g_tt_samp,this%g_tp_samp ,this%g_pp_samp! ,&
               !this%g_up_tt_samp,this%g_up_tp_samp ,this%g_up_pp_samp ,'***************'
      !write(*,*)'tttttt==', Omega, this%g_tt_samp + two * this%g_tp_samp * Omega + &
      !                                               this%g_pp_samp * Omega**2,r,theta
      !write(unit = *, fmt = *)'************************************************************' 
  
      this%Elec4U_CtrBL_At_samp(2) = zero
      this%Elec4U_CtrBL_At_samp(3) = zero
      this%Elec4U_CtrBL_At_samp(4) = this%Elec4U_CtrBL_At_samp(1) * Omega
      end subroutine Set_Elec4U_CtrBL_At_samp_Sub


!*******************************************************************************************************
      subroutine get_Phot4k_LNRF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      call Matrix_Multiplication14X44_Sub( this%Phot4k_CtrCF, &
                        this%Matrix_Of_CF2LNRF, this%Phot4k_LNRF )

      end subroutine get_Phot4k_LNRF_Sub

!*******************************************************************************************************
      subroutine get_Phot4k_CtrBL_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
  
      call Matrix_Multiplication14X44_Sub( this%Phot4k_CtrCF, &
                        this%Matrix_Of_BLCF_a_mu, this%Phot4k_CtrBL ) 

      end subroutine get_Phot4k_CtrBL_Sub

!*******************************************************************************************************
      subroutine get_Phot4k_CovBL_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this

      call Matrix_Multiplication44X41_Sub( this%Matrix_OF_BLCF_mu_a, &
                      this%Phot4k_CovCF, this%Phot4k_CovBL ) 

      end subroutine get_Phot4k_CovBL_Sub

!*******************************************************************************************
      subroutine Get_filename_Sub(this, dir, midname, id, ext, filename)
!*******************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      !=== variables for int2str ====
      integer :: id1, len
      real(mcp), intent(in) :: id
      character(*), intent(in) :: dir, midname, ext
      character(*), intent(inout) :: filename
      integer :: i
      real(mcp) :: id2
      !- - - - - - - - -- - - - - - -
      id2 = id
      if( id2 < 100.D0 )then
          do
              id2 = id2 * 10.D0
              if( id2 > 100.D0 )exit
          enddo
      endif
      id1 = id2
      write(filename, '(i8)') id1 !midname
      !write(*, *)'ss=', id1, id2
      filename=adjustl(filename)
      !len=len_trim(filename)
      do i=1, 2 !8-len
          !filename=filename//'0'
      enddo
      filename=trim(dir)//trim(midname)//trim(filename)//trim(ext)
      end subroutine Get_filename_Sub

!*******************************************************************************************
      subroutine Get_filename2_Sub(this, filename )
!*******************************************************************************************
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      !=== variables for int2str ==== 
      character(*), intent(inout) :: filename 
      !- - - - - - - - -- - - - - - - 
      filename=trim('./image/images/')//trim(filename)
      end subroutine Get_filename2_Sub


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      subroutine r_p2_Sub( this, r_ini, vector_of_momentum, p, vector_of_position, r_p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      implicit none
      class(Basic_Variables_And_Methods_Of_Particle) :: this
      real(mcp), intent(in) :: r_ini(1: 3), vector_of_momentum(1: 3), p
      real(mcp), intent(out) :: vector_of_position(1: 3), r_p
      real(mcp) :: vector_p(1: 3)
 
      vector_p = r_ini + p * vector_of_momentum
      vector_of_position = vector_p 
      r_p = dsqrt( vector_p(1)**2 + vector_p(2)**2 + vector_p(3)**2 )

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      end subroutine r_p2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Basic_Variables_And_Methods



