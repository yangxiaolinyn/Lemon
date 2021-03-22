    module IntialParameters
    use constants
    implicit none
    integer, parameter :: N_sigma = 2000 
    integer, parameter :: vL_sc_up = 1000
    integer, parameter :: Num_mu = 50
    integer, parameter :: Num_Phi = 30
    integer, parameter :: Num_J_tau_z = 1
    integer, parameter :: Num_J_tau_the = 1
    integer, parameter :: Num_J_tau_phi = 1
    real(mcp), dimension(0: 20) :: psi_array, phi_array, chi_array, &
                      zeta_array, H2_array, H1_array, H_nu_array, H_r_array, H_l_array


    contains 

    Subroutine Set_psi_phi_chi_zera_array()
    implicit none

      psi_array(0) = zero
      psi_array(1) = 0.04075D0
      psi_array(2) = 0.09144D0
      psi_array(3) = 0.15105D0
      psi_array(4) = 0.21916D0
      psi_array(5) = 0.29557D0
      psi_array(6) = 0.38014D0
      psi_array(7) = 0.47276D0
      psi_array(8) = 0.57338D0
      psi_array(9) = 0.68195D0
      psi_array(10) = 0.79843D0
      psi_array(11) = 0.92277D0
      psi_array(12) = 1.05497D0
      psi_array(13) = 1.19501D0
      psi_array(14) = 1.34286D0
      psi_array(15) = 1.49852D0
      psi_array(16) = 1.66198D0
      psi_array(17) = 1.83321D0
      psi_array(18) = 2.01223D0
      psi_array(19) = 2.19902D0
      psi_array(20) = 2.39357D0

      phi_array(0) = one
      phi_array(1) = 1.12988D0
      phi_array(2) = 1.20976D0
      phi_array(3) = 1.26850D0
      phi_array(4) = 1.31108D0
      phi_array(5) = 1.33973D0
      phi_array(6) = 1.35569D0
      phi_array(7) = 1.35971D0
      phi_array(8) = 1.35228D0
      phi_array(9) = 1.33374D0
      phi_array(10) = 1.30437D0
      phi_array(11) = 1.26432D0
      phi_array(12) = 1.21375D0
      phi_array(13) = 1.15279D0
      phi_array(14) = 1.08153D0
      phi_array(15) = 1.00003D0
      phi_array(16) = 0.90836D0
      phi_array(17) = 0.80655D0
      phi_array(18) = 0.69468D0
      phi_array(19) = 0.57276D0
      phi_array(20) = 0.44083D0

      chi_array(0) = one
      chi_array(1) = 1.10352D0
      chi_array(2) = 1.18638D0
      chi_array(3) = 1.26329D0
      chi_array(4) = 1.33687D0
      chi_array(5) = 1.40821D0
      chi_array(6) = 1.47801D0
      chi_array(7) = 1.54664D0
      chi_array(8) = 1.61435D0
      chi_array(9) = 1.68132D0
      chi_array(10) = 1.74772D0
      chi_array(11) = 1.81362D0
      chi_array(12) = 1.87911D0
      chi_array(13) = 1.94425D0
      chi_array(14) = 2.00907D0
      chi_array(15) = 2.07365D0
      chi_array(16) = 2.13799D0
      chi_array(17) = 2.20213D0
      chi_array(18) = 2.26609D0
      chi_array(19) = 2.32990D0
      chi_array(20) = 2.39356D0

      zeta_array(0) = zero
      zeta_array(1) = 0.01824D0
      zeta_array(2) = 0.03764D0
      zeta_array(3) = 0.05780D0
      zeta_array(4) = 0.07852D0
      zeta_array(5) = 0.09969D0
      zeta_array(6) = 0.12121D0
      zeta_array(7) = 0.14303D0
      zeta_array(8) = 0.16510D0
      zeta_array(9) = 0.18738D0
      zeta_array(10) = 0.20984D0
      zeta_array(11) = 0.23247D0
      zeta_array(12) = 0.25523D0
      zeta_array(13) = 0.27812D0
      zeta_array(14) = 0.30112D0
      zeta_array(15) = 0.32421D0
      zeta_array(16) = 0.34739D0
      zeta_array(17) = 0.37065D0
      zeta_array(18) = 0.39398D0
      zeta_array(19) = 0.41738D0
      zeta_array(20) = 0.44038D0


      H2_array(0) = one
      H2_array(1) = 1.04967D0
      H2_array(2) = 1.08621D0
      H2_array(3) = 1.11762D0
      H2_array(4) = 1.14552D0
      H2_array(5) = 1.17075D0
      H2_array(6) = 1.19383D0
      H2_array(7) = 1.21508D0
      H2_array(8) = 1.23476D0
      H2_array(9) = 1.25308D0
      H2_array(10) = 1.27019D0
      H2_array(11) = 1.28624D0
      H2_array(12) = 1.30132D0
      H2_array(13) = 1.31554D0
      H2_array(14) = 1.32895D0
      H2_array(15) = 1.34166D0
      H2_array(16) = 1.35371D0
      H2_array(17) = 1.36515D0
      H2_array(18) = 1.37601D0
      H2_array(19) = 1.38638D0
      H2_array(20) = 1.39625D0
 
      H1_array(0) = one
      H1_array(1) = 1.07301D0
      H1_array(2) = 1.12164D0
      H1_array(3) = 1.16151D0
      H1_array(4) = 1.19571D0
      H1_array(5) = 1.22577D0
      H1_array(6) = 1.25256D0
      H1_array(7) = 1.27674D0
      H1_array(8) = 1.29872D0
      H1_array(9) = 1.31882D0
      H1_array(10) = 1.33733D0
      H1_array(11) = 1.35444D0
      H1_array(12) = 1.37030D0
      H1_array(13) = 1.38507D0
      H1_array(14) = 1.39886D0
      H1_array(15) = 1.41179D0
      H1_array(16) = 1.42392D0
      H1_array(17) = 1.43533D0
      H1_array(18) = 1.44608D0
      H1_array(19) = 1.45625D0
      H1_array(20) = 1.46586D0

      H_nu_array(0) = 1.D0
      H_nu_array(1) = 1.02031D0
      H_nu_array(2) = 1.03834D0
      H_nu_array(3) = 1.05458D0
      H_nu_array(4) = 1.06934D0
      H_nu_array(5) = 1.08284D0
      H_nu_array(6) = 1.09525D0
      H_nu_array(7) = 1.10671D0
      H_nu_array(8) = 1.11735D0
      H_nu_array(9) = 1.12724D0
      H_nu_array(10) = 1.13647D0
      H_nu_array(11) = 1.14510D0
      H_nu_array(12) = 1.15320D0
      H_nu_array(13) = 1.16082D0
      H_nu_array(14) = 1.16799D0
      H_nu_array(15) = 1.17476D0
      H_nu_array(16) = 1.18116D0
      H_nu_array(17) = 1.18722D0
      H_nu_array(18) = 1.19297D0
      H_nu_array(19) = 1.19843D0
      H_nu_array(20) = 1.20362D0

      H_l_array(0) = 1.D0
      H_l_array(1) = 1.1814D0
      H_l_array(2) = 1.3255D0
      H_l_array(3) = 1.4596D0
      H_l_array(4) = 1.5884D0
      H_l_array(5) = 1.7137D0
      H_l_array(6) = 1.8367D0
      H_l_array(7) = 1.9579D0
      H_l_array(8) = 2.0778D0
      H_l_array(9) = 2.1966D0
      H_l_array(10) = 2.3146D0
      H_l_array(11) = 2.4319D0
      H_l_array(12) = 2.5486D0
      H_l_array(13) = 2.6649D0
      H_l_array(14) = 2.7807D0
      H_l_array(15) = 2.8962D0
      H_l_array(16) = 3.0113D0
      H_l_array(17) = 3.1262D0
      H_l_array(18) = 3.2408D0
      H_l_array(19) = 3.3552D0
      H_l_array(20) = 3.4695D0


      H_r_array(0) = 1.D0
      H_r_array(1) = 1.05737D0
      H_r_array(2) = 1.09113D0
      H_r_array(3) = 1.11703D0
      H_r_array(4) = 1.13816D0
      H_r_array(5) = 1.15594D0
      H_r_array(6) = 1.17128D0
      H_r_array(7) = 1.18468D0
      H_r_array(8) = 1.19654D0
      H_r_array(9) = 1.20713D0
      H_r_array(10) = 1.21668D0
      H_r_array(11) = 1.22532D0
      H_r_array(12) = 1.23320D0
      H_r_array(13) = 1.24042D0
      H_r_array(14) = 1.24705D0
      H_r_array(15) = 1.25318D0
      H_r_array(16) = 1.25886D0
      H_r_array(17) = 1.26414D0
      H_r_array(18) = 1.26906D0
      H_r_array(19) = 1.27366D0
      H_r_array(20) = 1.27797D0

    end Subroutine Set_psi_phi_chi_zera_array

    end module IntialParameters









