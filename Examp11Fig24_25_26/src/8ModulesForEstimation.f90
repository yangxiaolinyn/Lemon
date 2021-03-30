      module ModuleForCompScatEsti 
      use ScatDistance_FlatSP 
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photon_With_ScatDistance_FlatSP) :: Photon_ForEstimation
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
          real(mcp) :: PolArrIQUVpmu11(1: 5, 0: vL_sc_up) = zero
          real(mcp) :: PolArrIQUVpmu50(1: 5, 0: vL_sc_up) = zero
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: PolArrImu11(0: 6, 0: vL_sc_up) = zero 
          real(mcp) :: PolArrImu50(0: 6, 0: vL_sc_up) = zero  
          real(mcp) :: PolArrIQUV(1: 4, 0: 6, 1: Num_mu, 0: vL_sc_up) = zero   
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      contains 
!******************************************************************************************************* 
          procedure, public :: get_J_emissivity_for_estimation_Phiarr   =>   &
                               get_J_emissivity_for_estimation_Phiarr_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          procedure, public :: Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF   =>    &
                               Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF_Sub
          procedure, public :: Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withPol   =>   &
                               Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withPol_Sub
          procedure, public :: Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withpol_medium2  =>   &
                               Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withpol_medium2_Sub
          procedure, public :: Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation    =>   &
                               Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_Sub 
          procedure, public :: Get_PolaVector_And_IQU_InCF_for_Estimation    =>   &
                               Get_PolaVector_And_IQU_InCF_for_Estimation_Sub 
          procedure, public :: Get_Phot4k_phi_As_Egiven_in_CF   =>    &
                               Get_Phot4k_phi_As_Egiven_in_CF_Sub
          procedure, public :: Set_Phot_f4_Tetrad_In_Phot_CFrame    =>   &
                               Set_Phot_f4_Tetrad_In_Phot_CFrame_Sub 
          procedure, public :: Get_Observed_Stokes_Parameters   =>   &
                               Get_Observed_Stokes_Parameters_Sub
          procedure, public :: Get_Observed_Energy_Bin_index_i  =>  &
                               Get_Observed_Energy_Bin_index_i_Sub
          procedure, public :: Making_An_Estimation_One_MC_Component_1   =>   &
                               Making_An_Estimation_One_MC_Component_1_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon_ForEstimation
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      private :: get_J_emissivity_for_estimation_Phiarr_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF_Sub 
      private :: Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_Sub
      private :: Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withPol_Sub
      private :: Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withpol_medium2_Sub
      private :: Get_PolaVector_And_IQU_InCF_for_Estimation_Sub 
      private :: Get_Phot4k_phi_As_Egiven_in_CF_Sub
      private :: Set_Phot_f4_Tetrad_In_Phot_CFrame_Sub 
      private :: Get_Observed_Stokes_Parameters_Sub
      private :: Get_Observed_Energy_Bin_index_i_Sub
      private :: Making_An_Estimation_One_MC_Component_1_Sub
 
      contains   

!*******************************************************************************************************
      subroutine Making_An_Estimation_One_MC_Component_1_Sub( this, sphot )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this 
      TYPE(ScatPhoton_KN), INTENT(INOUT) :: sphot 
      integer :: i_1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Vector_Stokes4_ECF = sphot%Vector_Stokes4_ECF
      this%Phot4k_In_Elec_CF = sphot%Phot4k_In_Elec_CF
      this%Elec_Phot_mu = sphot%Elec_Phot_mu
      this%Elec_Phot_sin = sphot%Elec_Phot_sin
      this%Elec_Phot_phi = sphot%Elec_Phot_phi
      this%Elec_gama = sphot%Elec_gama
      this%Elec_V = sphot%Elec_V

      this%Elec_Phot_sin_In_Elec_CF = sphot%Elec_Phot_sin_In_Elec_CF 
      this%Elec_Phot_mu_In_Elec_CF = sphot%Elec_Phot_mu_In_Elec_CF

      this%Matrix_Of_Tetrad_Of_ElecAxis = sphot%Matrix_Of_Tetrad_Of_ElecAxis
      this%Matrix_Of_Tetrad_Of_PhotAxis = sphot%Matrix_Of_Tetrad_Of_PhotAxis  
      this%Matrix_Of_Tetrad1_Of_photAxis = sphot%Matrix_Of_Tetrad1_Of_photAxis
      this%Matrix_ECF_2_ECF1 = sphot%Matrix_ECF_2_ECF1
      this%Matrix_ECF1_2_ECF = sphot%Matrix_ECF1_2_ECF

      CALL this%Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation() 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      return
      end subroutine Making_An_Estimation_One_MC_Component_1_Sub
 

!*******************************************************************************************************
      subroutine Get_Observed_Energy_Bin_index_i_Sub( this, v_esti, E_i )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      real(mcp), intent(in) :: v_esti
      integer, intent(out) :: E_i
      integer :: i_1
        
    
      i_1 = floor( ( dlog10(v_esti) - this%log10_Tbb - this%y1 ) / (this%dy/two) )
      if(mod(i_1, 2)==0)then
          E_i = i_1 / 2
      else
          E_i = (i_1 + 1) / 2
      endif 
 
      return
      end subroutine Get_Observed_Energy_Bin_index_i_Sub


!*******************************************************************************************************
      subroutine get_J_emissivity_for_estimation_Phiarr_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      !integer, intent(in) :: mu_i 
      integer :: i_E, mu_i, i_phi
      real(mcp) :: vLv, v_esti, J_esti, Sigam_a, nu
       
      do mu_i = 1, Num_mu
        do i_phi = 0, Num_phi
              !v_esti = this%E_array_esti(i_E)
              nu = this%ln_nu1 + ( this%ln_nu2 - this%ln_nu1 ) * ranmar()
              v_esti = 10.D0**( nu ) * h_ev * 1.D-6
              Sigam_a = this%Sigma_fn(v_esti) * this%n_e1

              if(this%mu_estimates(mu_i) > zero)then
                  vLv = dexp( - this%z_tau * Sigam_a / this%mu_estimates(mu_i) ) / &
                           this%mu_estimates(mu_i)  * Sigam_a
              else
                  vLv = - dexp( (this%z_max - this%z_tau) * Sigam_a / this%mu_estimates(mu_i) ) / &
                           this%mu_estimates(mu_i)  * Sigam_a
              endif

              !write(*, *)'ffs=', mu_i, this%mu_estimates(mu_i), this%z_tau * Sigam_a / this%mu_estimates(mu_i)
              J_esti = ( v_esti / this%T_s )**3 / ( dexp( v_esti / this%T_s ) - one ) * &
                   ( this%mu_estimates(mu_i) + two * this%mu_estimates(mu_i)**2 ) * &
                     this%P_mu_normal! * this%P_nu_normal
    
              call this%Get_Observed_Energy_Bin_index_i( v_esti, i_E )
 
              if( i_E > vL_sc_up .or. i_E < 0 )cycle

              this%PolArrIQUV(1, 0, mu_i, i_E) = this%PolArrIQUV(1, 0, mu_i, i_E) + J_esti * vLv! * v_esti
              this%PolArrIQUV(1, 6, mu_i, i_E) = this%PolArrIQUV(1, 6, mu_i, i_E) + J_esti * vLv! * v_esti 
          enddo 
      enddo 
      return
      end subroutine get_J_emissivity_for_estimation_Phiarr_Sub
  

!*******************************************************************************************************
      subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this 
      integer :: i, i_phi, i_E, i_1, i_2, stimes, mu_i
      real(mcp) :: vLv, Lv, v_esti, J_esti, Qpsi, Upsi
      real(mcp) :: mu_ep, smu_ep, sin_phi_ep, mu_psi, epsip, epsi, KN_CrossSection, &
                   nup_vs_nu, chi, factor1, factor2, factor, Cos_Theta, Sin_Theta, &
                   KN_CrossSection_ECF, chi1, E_temp, Sin_Theta2, SigmaKN, Sigma_E1, &
                   temp_P1(1: 3), temp_value
      real(mcp) :: cos_phi, sin_phi, cos2phi, sin2phi 
      
      do mu_i = 1, Num_mu
        do i_phi = 0, Num_phi
          this%Phot3k_CF_esti(1) = this%smu_estimates(mu_i) * this%cos_phi_esti(i_phi)
          this%Phot3k_CF_esti(2) = this%smu_estimates(mu_i) * this%sin_phi_esti(i_phi)
          this%Phot3k_CF_esti(3) = this%mu_estimates(mu_i) 

          CALL this%Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          this%Phot3k_ECF_ini(1: 3) = this%Phot4k_In_Elec_CF(2: 4) / dabs( this%Phot4k_In_Elec_CF(1) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          factor = one - this%Phot3k_EF_esti(3) * this%Elec_v 
          !Cos_Theta = Vector3D_Inner_Product( this%Phot3k_ECF_ini, this%Phot3k_ECF_esti )
          Cos_Theta = this%Phot3k_ECF1_esti(3) 
           
          !write(*, *)'s1=', Cos_Theta, temp_value, Cos_Theta - temp_value

          Sin_Theta2 = one - Cos_Theta**2
          Sin_Theta = dsqrt( Sin_Theta2 )
           
          nup_vs_nu = one / ( one + dabs( this%Phot4k_In_Elec_CF(1) ) / mec2 * (one - Cos_Theta) )
 
          cos_phi = this%Phot3k_ECF1_esti(1) / Sin_Theta
          sin_phi = this%Phot3k_ECF1_esti(2) / Sin_Theta
          cos2phi = cos_phi**2 - sin_phi**2
          sin2phi = two * sin_phi * cos_phi
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          chi = nup_vs_nu + one / nup_vs_nu - Sin_Theta2
          this%Q_sp_scat = - Sin_Theta2 / chi
          this%U_sp_scat = zero 
          this%delta_pd_scat = dabs( this%Q_sp_scat ) 

          this%Vector_Stokes4_ECF_scat(1) = chi
          this%Vector_Stokes4_ECF_scat(2) = Sin_Theta2
          this%Vector_Stokes4_ECF_scat(3) = zero
          this%Vector_Stokes4_ECF_scat(4) = zero
  
          this%Scat_Phot3k_CF = this%Phot3k_CF_esti
          call this%Get_Scattered_Stokes_Vector_And_f_CF( this%Phot3k_ECF1_esti )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !write(*, *)'s1=', KN_CrossSection_ECF 
          epsi = ( two * this%E_ini / mec2 ) * this%Elec_gama * &
                         ( one - this%Elec_V * this%Elec_Phot_mu )
          if ( epsi > 8.19D-4 ) then
              SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
          else 
              SigmaKN = one - epsi
          end if 
          KN_CrossSection_ECF = half_re2 * chi * (nup_vs_nu / this%Elec_gama / factor)**2 / SigmaKN
          !write(*, *)'s1=', epsi, this%Phot4k_In_Elec_CF(1) * two / mec2, KN_CrossSection_ECF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

          E_temp = dabs( this%Phot4k_In_Elec_CF(1) ) * nup_vs_nu
          this%E_esti_vs_mu_phi = this%Elec_gama * E_temp * ( one + this%Elec_v * this%Phot3k_ECF_esti(3) )
          !write(*, *)'ffs=', this%E_esti_vs_mu_phi, this%E_array_esti(0), this%E_array_esti(vL_sc_up)    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Sigma_E1 = this%sigma_fn( this%E_esti_vs_mu_phi ) * this%n_e1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(this%mu_estimates(mu_i) > zero)then
              Lv = this%w_ini * dexp( - this%z_tau * Sigma_E1 / this%mu_estimates(mu_i) ) / &
                          this%mu_estimates(mu_i) * KN_CrossSection_ECF * Sigma_E1 
          else
              Lv = - this%w_ini * dexp( (this%z_max - this%z_tau) * Sigma_E1 / this%mu_estimates(mu_i) ) / &
                          this%mu_estimates(mu_i) * KN_CrossSection_ECF * Sigma_E1 
          endif
          vLv = Lv! * dabs( this%E_esti_vs_mu_phi )

          call this%Get_Observed_Energy_Bin_index_i( this%E_esti_vs_mu_phi, i_E )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
          if( i_E > vL_sc_up .or. i_E < 0 )cycle
          stimes = this%scatter_times + 1  

          !write(*, *)'s2=',  i_E, this%E_esti_vs_mu_phi, this%cos_phi_esti(i_phi)
          if(stimes == 0)stop
          this%PolArrIQUV(1, 6, mu_i, i_E) = this%PolArrIQUV(1, 6, mu_i, i_E) + vLv * &
                                           this%Vector_Stokes4_ECF_scat(1)

          if(stimes <= 4)then
              this%PolArrIQUV(1, stimes, mu_i, i_E) = this%PolArrIQUV(1, stimes, mu_i, i_E) + &
                                              vLv * this%Vector_Stokes4_ECF_scat(1)
          else
              this%PolArrIQUV(1, 5, mu_i, i_E) = this%PolArrIQUV(1, 5, mu_i, i_E) + &
                                              vLv * this%Vector_Stokes4_ECF_scat(1)
          endif 

          !if( this%delta_pd_scat /= zero )then 
              call this%Get_Observed_Stokes_Parameters( )
              this%PolArrIQUV(2, 6, mu_i, i_E) = this%PolArrIQUV(2, 6, mu_i, i_E) + &
                                  Lv * this%Vector_Stokes4_ECF_scat(2)
              this%PolArrIQUV(3, 6, mu_i, i_E) = this%PolArrIQUV(3, 6, mu_i, i_E) + &
                                  Lv * this%Vector_Stokes4_ECF_scat(3) 
          !endif 

        enddo
        !write(*, *)'#######################################' 
      enddo  
      return
      end subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_Sub
!*******************************************************************************************************



!*******************************************************************************************************
      Subroutine Get_Observed_Stokes_Parameters_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this   
      real(mcp), dimension(1: 3) :: Temp_Vec, Scat_Phot_Tetrad_AxisX, Scat_Phot_Tetrad_AxisY, &
                 Scat_Phot_Tetrad_AxisZ, f3_x, f3_y
      real(mcp), dimension(1: 4) :: f4_scat_e
      real(mcp) :: mu_0, normal_L, phi_obs, mu_signs, mu_ya_yb


      Scat_Phot_Tetrad_AxisZ = this%Phot3k_CF_esti 
      Temp_Vec(1) = zero
      Temp_Vec(2) = zero
      Temp_Vec(3) = one

      mu_0 = Vector3D_Inner_Product( Scat_Phot_Tetrad_AxisZ, Temp_Vec )
      Scat_Phot_Tetrad_AxisY = Temp_Vec - mu_0 * Scat_Phot_Tetrad_AxisZ
      normal_L = Vector3D_Length( Scat_Phot_Tetrad_AxisY )
      Scat_Phot_Tetrad_AxisY = Scat_Phot_Tetrad_AxisY / normal_L
      call this%Vector_Cross_Product( Scat_Phot_Tetrad_AxisY, Scat_Phot_Tetrad_AxisZ, &
                                      Scat_Phot_Tetrad_AxisX ) 
  
      f3_x = this%f4_scat_CF(2: 4)
      call this%Vector_Cross_Product( Scat_Phot_Tetrad_AxisZ, f3_x, f3_y ) 

      mu_ya_yb = Vector3D_Inner_Product( f3_y, Scat_Phot_Tetrad_AxisY )
      mu_signs = Vector3D_Inner_Product( f3_x, Scat_Phot_Tetrad_AxisY )
      phi_obs = - sign(one, mu_signs) * dacos( mu_ya_yb )
      
      CALL this%StokesPara_Rotation_Matrix( phi_obs, this%Vector_Stokes4_ECF_scat )

      return
      end subroutine Get_Observed_Stokes_Parameters_Sub
 
 

!*******************************************************************************************************
      subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withPol_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      !integer, intent(in) :: mu_i 
      integer :: i, i_phi, i_E, i_1, i_2, scat_times, mu_i
      real(mcp) :: vLv, Lv, v_esti, J_esti, Qpsi, Upsi
      real(mcp) :: mu_ep, smu_ep, sin_phi_ep, mu_psi, epsip, epsi, KN_CrossSection, &
                   nup_vs_nu, chi, factor, Cos_Theta, Sin_Theta, SigmaKN, &
                   KN_CrossSection_ECF, chi1, Sin_Theta_Cos_Phi, E_temp, Sigma_E1
      real(mcp) :: N_temp, Sin_Theta2, x, y, cosphi, sinphi, cos2phi, sin2phi
      real(mcp) :: F0, F11, F1, F22, F33, tp_vales, Q_xi, U_xi, I_xi, V_xi, gam
       
      do mu_i = 1, Num_mu
        do i_phi = 0, Num_phi - 1
          this%Phot3k_CF_esti(1) = this%smu_estimates(mu_i) * this%cos_phi_esti(i_phi)
          this%Phot3k_CF_esti(2) = this%smu_estimates(mu_i) * this%sin_phi_esti(i_phi)
          this%Phot3k_CF_esti(3) = this%mu_estimates(mu_i) 

          CALL this%Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          this%Phot3k_ECF_ini(1: 3) = this%Phot4k_In_Elec_CF(2: 4) / dabs( this%Phot4k_In_Elec_CF(1) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          factor = one - this%Phot3k_EF_esti(3) * this%Elec_v  
          !Cos_Theta = Vector3D_Inner_Product( this%Phot3k_ECF_ini, this%Phot3k_ECF_esti )
          Cos_Theta = this%Phot3k_ECF1_esti(3) 
           
          !write(*, *)'s1=', Cos_Theta, temp_value, Cos_Theta - temp_value

          Sin_Theta2 = one - Cos_Theta**2
          Sin_Theta = dsqrt( Sin_Theta2 )
           
          nup_vs_nu = one / ( one + dabs( this%Phot4k_In_Elec_CF(1) ) / mec2 * (one - Cos_Theta) ) 
          !write(*, fmt="(' ', A5, 3ES17.6)")'s1=', Cos_Theta, Cos_Theta-one, &
          !        dabs( this%Phot4k_In_Elec_CF(1) ) / mec2 
 
          cosphi = this%Phot3k_ECF1_esti(1) / Sin_Theta
          sinphi = this%Phot3k_ECF1_esti(2) / Sin_Theta
          cos2phi = cosphi**2 - sinphi**2
          sin2phi = two * sinphi * cosphi

          !write(*, fmt=*)'s1=', cosphi, sinphi
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          F0 = nup_vs_nu + one / nup_vs_nu - Sin_Theta2
          F1 = Sin_Theta2
          F11 = one + Cos_Theta**2
          F22 = two * Cos_Theta
          F33 = ( nup_vs_nu + one / nup_vs_nu ) * Cos_Theta
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          I_xi = this%Vector_Stokes4_ECF(1) * F0 + F1 * ( this%Vector_Stokes4_ECF(2) * cos2phi + &
                                                          this%Vector_Stokes4_ECF(3) * sin2phi )
          Q_xi = this%Vector_Stokes4_ECF(1) * F1 + F11 * (this%Vector_Stokes4_ECF(2) * cos2phi + &
                                                          this%Vector_Stokes4_ECF(3) * sin2phi )
          U_xi = F22 * (-this%Vector_Stokes4_ECF(2) * sin2phi + this%Vector_Stokes4_ECF(3) * cos2phi)  
          V_xi = F33 * this%Vector_Stokes4_ECF(4)
 

          this%Vector_Stokes4_ECF_scat(1) = I_xi
          this%Vector_Stokes4_ECF_scat(2) = Q_xi
          this%Vector_Stokes4_ECF_scat(3) = U_xi
          this%Vector_Stokes4_ECF_scat(4) = V_xi
          !this%Vector_Stokes4_ECF_scat(4) = zero
  
          this%Scat_Phot3k_CF = this%Phot3k_CF_esti
          call this%Get_Scattered_Stokes_Vector_And_f_CF( this%Phot3k_ECF1_esti )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Q_xi = this%Vector_Stokes4_ECF(2) / this%Vector_Stokes4_ECF(1)
          U_xi = this%Vector_Stokes4_ECF(3) / this%Vector_Stokes4_ECF(1)
          N_temp = I_xi / this%Vector_Stokes4_ECF(1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          !write(*, *)'s1=', KN_CrossSection_ECF 
          epsi = ( two * this%E_ini / mec2 ) * this%Elec_gama * &
                         ( one - this%Elec_V * this%Elec_Phot_mu )
          if ( epsi > 8.19D-4 ) then
              SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
          else 
              SigmaKN = one - epsi
          end if  
          !KN_CrossSection = half_re2 * N_temp * (nup_vs_nu / this%Elec_gama / factor)**2 / SigmaKN
          KN_CrossSection = half_re2 * (nup_vs_nu / this%Elec_gama / factor)**2 / SigmaKN
          !write(*, *)'s1=', epsi, this%Phot4k_In_Elec_CF(1) * two / mec2, KN_CrossSection_ECF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

          E_temp = dabs( this%Phot4k_In_Elec_CF(1) ) * nup_vs_nu
          this%E_esti_vs_mu_phi = this%Elec_gama * E_temp * ( one + this%Elec_v * this%Phot3k_ECF_esti(3) )
          !write(*, fmt="(' ', A8, 5ES20.10)")'ffs = ', this%E_esti_vs_mu_phi, nup_vs_nu, &
          !              E_temp, this%Elec_v, this%Phot3k_ECF_esti(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !Sigma_E1 = this%sigma_fn( this%E_esti_vs_mu_phi )! * this%n_e1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !if(this%mu_estimates(mu_i) > zero)then
          !    Lv = this%w_ini * dexp( - this%z_tau * Sigma_E1 / this%mu_estimates(mu_i) ) / &
          !                this%mu_estimates(mu_i) * KN_CrossSection * Sigma_E1 
          !else
          !    Lv = - this%w_ini * dexp( (this%z_max - this%z_tau) * Sigma_E1 / this%mu_estimates(mu_i) ) / &
          !                this%mu_estimates(mu_i) * KN_CrossSection * Sigma_E1 
          !endif
          vLv = KN_CrossSection !* Sigma_E1

          call this%Get_Observed_Energy_Bin_index_i( this%E_esti_vs_mu_phi, i_E )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
           !write(*, fmt=300)'ffs=', this%w_ini, dexp( - this%z_tau / this%mu_estimates(mu_i) ), &
          !      one / this%mu_estimates(mu_i),  KN_CrossSection, this%scatter_times
          !300 FORMAT (' ', A5, 5ES15.6)

          !gam = dsqrt( this%E_esti_vs_mu_phi / (2.5D-11 * mec2) ) / &
          !      dsqrt( two * ( one - dcos( pi * 85.D0 / 180.D0 ) ) )

          !write(*, *)'ffs= ', gam, this%Elec_gama
          if( i_E > vL_sc_up .or. i_E < 0 )cycle
          scat_times = this%scatter_times + 1 

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          this%PolArrIQUV(1, 6, mu_i, i_E) = this%PolArrIQUV(1, 6, mu_i, i_E) + &
                                              vLv * this%Vector_Stokes4_ECF_scat(1)

          if(scat_times <= 4)then
              this%PolArrIQUV(1, scat_times, mu_i, i_E) = this%PolArrIQUV(1, scat_times, mu_i, i_E) + &
                                              vLv * this%Vector_Stokes4_ECF_scat(1)
          else
              this%PolArrIQUV(1, 5, mu_i, i_E) = this%PolArrIQUV(1, 5, mu_i, i_E) + &
                                              vLv * this%Vector_Stokes4_ECF_scat(1)
          endif   

          !if( this%delta_pd_scat /= zero )then 
              !call this%Get_PolaVector_And_IQU_InCF_for_Estimation(Cos_Theta, Sin_Theta, &
              !                            sinphi, cosphi, Qpsi, Upsi)
              call this%Get_Observed_Stokes_Parameters( ) 

              this%PolArrIQUV(2, 6, mu_i, i_E) = this%PolArrIQUV(2, 6, mu_i, i_E) + &
                                  vLv * this%Vector_Stokes4_ECF_scat(2)
              this%PolArrIQUV(3, 6, mu_i, i_E) = this%PolArrIQUV(3, 6, mu_i, i_E) + &
                                  vLv * this%Vector_Stokes4_ECF_scat(3) 
              this%PolArrIQUV(4, 6, mu_i, i_E) = this%PolArrIQUV(4, 6, mu_i, i_E) + &
                                  vLv * this%Vector_Stokes4_ECF_scat(4) 

!write(*, fmt=*)'ffs=',this%Vector_Stokes4_ECF_scat(2),this%Vector_Stokes4_ECF_scat(3), vLv, this%Vector_Stokes4_ECF_scat(1)
          !endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        enddo  
      enddo  
      return
      end subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Estimation_withPol_Sub
!*******************************************************************************************************



!*******************************************************************************************************
      Subroutine Get_PolaVector_And_IQU_InCF_for_Estimation_Sub(this, Cos_Theta, Sin_Theta, &
                                          sinphi, cosphi, Qpsi, Upsi)
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this 
      real(mcp), intent(in) :: Cos_Theta, Sin_Theta, sinphi, cosphi
      real(mcp), intent(out) :: Qpsi, Upsi
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, InvMatrix
      real(mcp) :: beta, Psi_pp, Vector_Length, cos_psi, fx, fy, sin_psi
      real(mcp), dimension(1: 3) :: Scat_Phot_Tetrad_AxisX, Scat_Phot_Tetrad_AxisY, &
                   Scat_Phot_Tetrad_AxisZ, e_p, e_f, e_y, e_xp, e_yp, f3_scat_pp, &
                   f3_scat_e_tilde, f4, Axis_x, Axis_y, Axis_z
      real(mcp), dimension(1: 4) :: f4_scat_e


      Scat_Phot_Tetrad_AxisZ(1) = Sin_Theta * cosphi !this%Phot3k_ECF_esti
      Scat_Phot_Tetrad_AxisZ(2) = Sin_Theta * sinphi !this%Phot3k_ECF_esti
      Scat_Phot_Tetrad_AxisZ(3) = Cos_Theta !this%Phot3k_ECF_esti
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one

      call this%Vector_Cross_Product( Scat_Phot_Tetrad_AxisZ, e_p, Scat_Phot_Tetrad_AxisX )
      Vector_Length = Vector3D_Length( Scat_Phot_Tetrad_AxisX )
      If ( Vector_Length /= zero ) then
            Scat_Phot_Tetrad_AxisX = &
          - Scat_Phot_Tetrad_AxisX / Vector_Length
          call this%Vector_Cross_Product( Scat_Phot_Tetrad_AxisZ, Scat_Phot_Tetrad_AxisX, &
                                      Scat_Phot_Tetrad_AxisY )
      Else
          e_f(1) = one
          e_f(2) = zero
          e_f(3) = zero
          e_y(1) = zero
          e_y(2) = one
          e_y(3) = zero
          e_yp = cosphi * e_f + sinphi * e_y
          CALL this%Vector_Cross_Product( e_yp, e_p, e_xp )
          Scat_Phot_Tetrad_AxisX = e_xp
          Scat_Phot_Tetrad_AxisY = DSIGN( one, Cos_Theta ) * e_yp
      Endif

      Psi_pp = DATAN( DABS( this%U_sp_scat / this%Q_sp_scat ) ) / two
      If ( this%U_sp_scat < zero ) then
          Psi_pp = pi - Psi_pp
      Endif 
      beta = DATAN( DABS( this%U_sp_scat / this%Q_sp_scat ) ) / two
      If ( this%Q_sp_scat > zero .and. this%U_sp_scat > zero ) then
          Psi_pp = beta
      else If ( this%Q_sp_scat > zero .and. this%U_sp_scat < zero ) then
          Psi_pp = pi - beta
      else If ( this%Q_sp_scat < zero .and. this%U_sp_scat > zero ) then
          Psi_pp = pi / two - beta
      else If ( this%Q_sp_scat < zero .and. this%U_sp_scat < zero ) then
          Psi_pp = beta + pi / two
      Endif

      f3_scat_pp = DCOS(psi_pp) * Scat_Phot_Tetrad_AxisX + &
                   DSIN(psi_pp) * Scat_Phot_Tetrad_AxisY

      CALL Matrix_Multiplication13X33_Sub( f3_scat_pp, &
                  this%Matrix_Of_Tetrad_Of_Phot_f4_Axis, f3_scat_e_tilde )

      f4_scat_e(1) = this%Elec_gama * this%Elec_V * f3_scat_e_tilde(3)
      f4_scat_e(2) = f3_scat_e_tilde(1)
      f4_scat_e(3) = f3_scat_e_tilde(2)
      f4_scat_e(4) = this%Elec_gama * f3_scat_e_tilde(3)

      !CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
      !                   this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_3X3 )
      CALL Matrix_Multiplication13X33_Sub( f4_scat_e(2:4), this%Matrix_Elec_2_CF, &
                                                         this%f4_scat_CF(2:4) )
      this%f4_scat_CF(1) = f4_scat_e(1) 
      !this%f4_scat_CovCF = this%f4_scat_CF 
      !this%f4_scat_CovCF(1) = - this%f4_scat_CF(1) 

      Axis_z = this%Phot3k_CF_esti

      if( .true. )then
          Axis_y(1) = zero
          Axis_y(2) = zero
          Axis_y(3) = one 
          call this%Vector_Cross_Product( Axis_y, Axis_z, Axis_x )
          Axis_x = Axis_x / Vector3D_Length( Axis_x )
          call this%Vector_Cross_Product( Axis_z, Axis_x, Axis_y ) 
      else
          Axis_x(1) = zero
          Axis_x(2) = zero
          Axis_x(3) = one 
          call this%Vector_Cross_Product( Axis_z, Axis_x, Axis_y )
          Axis_y = Axis_y / Vector3D_Length( Axis_y )
          call this%Vector_Cross_Product( Axis_y, Axis_z, Axis_x )
      endif
 
      f4 = this%f4_scat_CF(2: 4) - this%f4_scat_CF(1) * this%Phot3k_CF_esti 

      cos_psi = Vector3D_Inner_Product( f4, Axis_x )
      sin_psi = dsqrt( dabs( one - cos_psi**2 ) )
      fx = Vector3D_Inner_Product( f4, Axis_x )
      fy = Vector3D_Inner_Product( f4, Axis_y )
      if( fx > zero .and. fy > zero  )then
          Qpsi = this%delta_pd_scat * ( two * cos_psi ** 2 - one )
          Upsi = this%delta_pd_scat * two * sin_psi * cos_psi 
      else if( fx < zero .and. fy > zero  )then
          Qpsi = this%delta_pd_scat * ( two * cos_psi ** 2 - one )
          Upsi = this%delta_pd_scat * two * sin_psi * cos_psi 
      else if( fx < zero .and. fy < zero  )then
          Qpsi = this%delta_pd_scat * ( two * cos_psi ** 2 - one )
          Upsi = this%delta_pd_scat * two * sin_psi * dabs( cos_psi ) 
      else if( fx > zero .and. fy < zero  )then
          Qpsi = this%delta_pd_scat * ( two * cos_psi ** 2 - one )
          Upsi = - this%delta_pd_scat * two * sin_psi * cos_psi 
      endif
      return
      end subroutine Get_PolaVector_And_IQU_InCF_for_Estimation_Sub



!*******************************************************************************************************
      Subroutine Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      !real(mcp), intent(in) :: Phot3k_CF(1: 3) 
      !real(mcp), intent(inout) :: Phot3k_EF(1: 3) 
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, InvMatrix
      real(mcp) :: factor

      CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
                         this%Matrix_Of_Tetrad_Of_PhotAxis, this%Matrix_Elec_2_CF ) 

      CALL Inverse_Matrix_Of_Matrix_3X3_Sub(this%Matrix_Elec_2_CF, this%Matrix_CF_2_Elec)

      !CALL Matrix_Multiplication33X33_Sub( Temp_Matrix_3X3, InvMatrix, Temp_Matrix_3X3_1 ) 
 
      CALL Matrix_Multiplication13X33_Sub( this%Phot3k_CF_esti, this%Matrix_CF_2_Elec, this%Phot3k_EF_esti )
      factor = one - this%Elec_v * this%Phot3k_EF_esti(3)
      this%Phot3k_ECF_esti(1) = this%Phot3k_EF_esti(1) / this%Elec_gama / factor
      this%Phot3k_ECF_esti(2) = this%Phot3k_EF_esti(2) / this%Elec_gama / factor
      this%Phot3k_ECF_esti(3) = (this%Phot3k_EF_esti(3) - this%Elec_v) / factor

      CALL Matrix_Multiplication13X33_Sub( this%Phot3k_ECF_esti, this%Matrix_ECF_2_ECF1, &
                                           this%Phot3k_ECF1_esti )

      !this%PhotE_ECF_esti = this%Elec_gama * dabs( this%Phot4k_CtrCF(1) ) * &
      !                            ( one - this%Elec_V * this%Elec_Phot_mu )
      
      !write(*, fmt="(' ', A5, 4ES18.7)")'s0=', this%Phot3k_CF_esti, Vector3D_Length( this%Phot3k_CF_esti )
      !write(*, fmt="(' ', A5, 4ES18.7)")'s1=', this%Phot3k_EF_esti, Vector3D_Length( this%Phot3k_EF_esti )
      !write(*, fmt="(' ', A5, 5ES18.7)")'s2=', this%Phot3k_ECF_esti, factor, this%Elec_v
  
      end subroutine Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF_Sub
!*******************************************************************************************************
  

!*******************************************************************************************************
      subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withPol_medium2_Sub(this, mu_i)
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      integer, intent(in) :: mu_i 
      integer :: i, i_phi, i_E, i_1, i_2, scat_times
      real(mcp) :: vLv, Lv, v_esti, J_esti, Qpsi, Upsi
      real(mcp) :: mu_ep, smu_ep, sin_phi_ep, mu_psi, epsip, epsi, KN_CrossSection, &
                   nup_vs_nu, chi, factor, Cos_Theta, Sin_Theta, SigmaKN, &
                   KN_CrossSection_ECF, chi1, Sin_Theta_Cos_Phi, E_temp, Sigma_E1, Sigma_E2
      real(mcp) :: N_temp, Sin_Theta2, x, y, cosphi, sinphi, cos2phi, sin2phi, &
                   tau1, tau2
       
      call this%Set_Phot_f4_Tetrad_In_Phot_CF()
      do i_phi = 0, N_esti_phi
          this%Phot3k_CF_esti(1) = this%smu_esti(mu_i) * this%cos_phi_esti(i_phi)
          this%Phot3k_CF_esti(2) = this%smu_esti(mu_i) * this%sin_phi_esti(i_phi)
          this%Phot3k_CF_esti(3) = this%mu_esti(mu_i)

          CALL Matrix_Multiplication13X33_Sub( this%Phot3k_CF_esti, this%Matrix_CF_2_Elec, &
                         this%Phot3k_EF_esti )
          !CALL this%Get_Phot4k_in_ElecFrame_From_Phot4k_in_CF()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !this%Phot3k_ECF_ini(1: 3) = this%Phot4k_In_Elec_CF(2: 4) / dabs( this%Phot4k_In_Elec_CF(1) )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !Cos_Theta = this%Phot3k_ECF_ini(1) * this%Phot3k_ECF_esti(1) + this%Phot3k_ECF_ini(2) * &
          !            this%Phot3k_ECF_esti(2) + this%Phot3k_ECF_ini(3) * this%Phot3k_ECF_esti(3)
          !Cos_Theta = Vector3D_Inner_Product( this%Phot3k_ECF_ini, this%Phot3k_ECF_esti )
          Cos_Theta = this%Phot3k_EF_esti(3)
          Sin_Theta2 = one - Cos_Theta**2
          Sin_Theta = dsqrt( Sin_Theta2 )
          x = Vector3D_Inner_Product( this%Phot_f4_AxisX_PCF, this%Phot3k_EF_esti )
          y = Vector3D_Inner_Product( this%Phot_f4_AxisY_PCF, this%Phot3k_EF_esti )

          cosphi = x / Sin_Theta
          sinphi = y / Sin_Theta
          cos2phi = ( x**2 - y**2 ) / Sin_Theta2
          sin2phi = two * x * y / Sin_Theta2

          nup_vs_nu = one / ( one + this%E_ini / mec2 * (one - Cos_Theta) ) 
          this%E_esti_vs_mu_phi = this%E_ini * nup_vs_nu 
          !write(*, *)'s1=', this%E_esti_vs_mu_phi, this%E_array_esti(0), this%E_array_esti(vL_sc_up)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Sigma_E1 = this%sigma_fn( this%E_esti_vs_mu_phi ) * this%n_e1
          Sigma_E2 = this%sigma_KNs( this%E_esti_vs_mu_phi ) * this%n_e2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          this%Q_sp = this%delta_pd * cos2phi
          this%U_sp = - this%delta_pd * sin2phi
          !**********************************************************************************
          N_temp = nup_vs_nu + one / nup_vs_nu - Sin_Theta2 * ( one + this%Q_sp )
          this%Q_sp_scat = - ( Sin_Theta2 - ( one + Cos_Theta**2) * this%Q_sp ) / N_temp
          this%U_sp_scat =  two * Cos_Theta * this%U_sp / N_temp
          this%delta_pd_scat = DSQRT( this%Q_sp_scat**2 + this%U_sp_scat**2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          factor = one - this%Phot3k_EF_esti(3) * this%Elec_v
          !KN_CrossSection = half_re2 * N_temp * (nup_vs_nu / this%Elec_gama / factor)**2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !write(*, *)'s1=', KN_CrossSection_ECF 
          epsi = two * this%E_ini / mec2
          if ( epsi > 8.19D-4 ) then
              SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
          else 
              SigmaKN = one - epsi
          end if 
          KN_CrossSection = half_re2 * N_temp * nup_vs_nu**2 / SigmaKN
          !write(*, *)'s1=', epsi, this%Phot4k_In_Elec_CF(1) * two / mec2, KN_CrossSection
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          tau1 = this%z_max * Sigma_E1 / this%mu_esti(mu_i)
          tau2 = ( this%z_tau - this%z_max ) / this%mu_esti(mu_i) * Sigma_E2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Lv = this%w_ini * dexp( - (tau1 + tau2) ) /  &
                               this%mu_esti(mu_i) * KN_CrossSection
          vLv = Lv * DABS( this%E_esti_vs_mu_phi )
          i_1 = floor( ( dlog10(this%E_esti_vs_mu_phi) - log10_mec2 - this%y1 ) / (this%dy/two) )
          if(mod(i_1, 2)==0)then
              i_E = i_1 / 2
          else
              i_E = (i_1 + 1) / 2
          endif
          write(*, fmt=300)'ffs=', this%w_ini, Lv, vLv
           !write(*, fmt=300)'ffs=', this%w_ini, dexp( - this%z_tau / this%mu_esti(mu_i) ), &
           !!     one / this%mu_esti(mu_i),  KN_CrossSection, this%scatter_times
          300 FORMAT (' ', A5, 5ES15.6)
          if( i_E > vL_sc_up .or. i_E < 0 )cycle
          scat_times = this%scatter_times + 1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if(mu_i == 1)then 
              this%PolArrIQUVpmu11(1, i_E) = this%PolArrIQUVpmu11(1, i_E) + Lv 
              this%PolArrImu11(6, i_E) = this%PolArrImu11(6, i_E) + vLv
              !write(*, *)'fs2==',this%w_ini, dexp( - this%z_tau / this%mu_esti(mu_i) ) / this%mu_esti(mu_i), &
              !             KN_CrossSection, this%E_esti_vs_mu_phi, Lv, vLv

              if( scat_times <= 4 )then
                  this%PolArrImu11(scat_times, i_E ) = this%PolArrImu11(scat_times, i_E) + vLv
              else 
                  this%PolArrImu11(5, i_E) = this%PolArrImu11(5, i_E) + vLv
              endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(this%delta_pd_scat /= zero)then
                  call this%Get_PolaVector_And_IQU_InCF_for_Estimation(Cos_Theta, Sin_Theta, &
                                          sinphi, cosphi, Qpsi, Upsi)
                  this%PolArrIQUVpmu11(2, i_E) = this%PolArrIQUVpmu11(2, i_E) + Lv * Qpsi !2 represents Q
                  this%PolArrIQUVpmu11(3, i_E) = this%PolArrIQUVpmu11(3, i_E) + Lv * Upsi !3 represents U
              endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if(mu_i == 2)then    
              this%PolArrIQUVpmu50(1, i_E) = this%PolArrIQUVpmu50(1, i_E) + Lv
              this%PolArrImu50(6, i_E) = this%PolArrImu50(6, i_E) + vLv

              if( scat_times <= 4 )then
                  this%PolArrImu50(scat_times, i_E ) = this%PolArrImu50(scat_times, i_E) + vLv
              else 
                  this%PolArrImu50(5, i_E) = this%PolArrImu50(5, i_E) + vLv
              endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(this%delta_pd_scat /= zero)then
                  call this%Get_PolaVector_And_IQU_InCF_for_Estimation(Cos_Theta, Sin_Theta, &
                                          sinphi, cosphi, Qpsi, Upsi)
                  this%PolArrIQUVpmu50(2, i_E) = this%PolArrIQUVpmu50(2, i_E) + Lv * Qpsi !2 represents Q
                  this%PolArrIQUVpmu50(3, i_E) = this%PolArrIQUVpmu50(3, i_E) + Lv * Upsi !3 represents U
              endif 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
      enddo  
      return
      end subroutine Get_K_P1_P2_Scat_Kernel_InECF_for_Esti_withPol_medium2_Sub
!*******************************************************************************************************



!*******************************************************************************************************
      Subroutine Get_Phot4k_phi_As_Egiven_in_CF_Sub(this, mu_i, Enger_i, x1, x2)
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this
      integer, intent(in) :: mu_i
      real(mcp), intent(in) :: Enger_i
      real(mcp), intent(out) :: x1, x2
      !real(mcp), intent(inout) :: Phot3k_EF(1: 3) 
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, InvMatrix
      real(mcp) :: factor
      real(mcp) :: a1, b1, c1, a2, b2, c2, a3, b3, c3, cos_Psi, &
                   B, C, d1, d2, d3, delta1, delta2, delta3, phi_t, phi01, phi02
 
      CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
                         this%Matrix_Of_Tetrad_Of_PhotAxis, this%Matrix_Elec_2_CF ) 

      CALL Inverse_Matrix_Of_Matrix_3X3_Sub(this%Matrix_Elec_2_CF, this%Matrix_CF_2_Elec)
 
      a1 = this%Matrix_CF_2_Elec(1, 1) * this%smu_esti(mu_i)
      b1 = this%Matrix_CF_2_Elec(2, 1) * this%smu_esti(mu_i)
      c1 = this%Matrix_CF_2_Elec(3, 1) * this%mu_esti(mu_i)
 
      a2 = this%Matrix_CF_2_Elec(1, 2) * this%smu_esti(mu_i)
      b2 = this%Matrix_CF_2_Elec(2, 2) * this%smu_esti(mu_i)
      c2 = this%Matrix_CF_2_Elec(3, 2) * this%mu_esti(mu_i)
 
      a3 = this%Matrix_CF_2_Elec(1, 3) * this%smu_esti(mu_i)
      b3 = this%Matrix_CF_2_Elec(2, 3) * this%smu_esti(mu_i)
      c3 = this%Matrix_CF_2_Elec(3, 3) * this%mu_esti(mu_i)

      C = (one + mec2 / dabs(this%Phot4k_In_Elec_CF(1)) ) * this%Elec_gama * this%Elec_v
      B = (one + mec2 / dabs( this%Phot4k_In_Elec_CF(1) ) ) * this%Elec_gama - mec2 / Enger_i
      delta1 = C * a3 + this%Phot3k_ECF_ini(1) * a1 + this%Phot3k_ECF_ini(2) * a2 + &
                        this%Phot3k_ECF_ini(3) * a3
      delta2 = C * b3 + this%Phot3k_ECF_ini(1) * b1 + this%Phot3k_ECF_ini(2) * b2 + &
                        this%Phot3k_ECF_ini(3) * b3
      delta3 = C * c3 + this%Phot3k_ECF_ini(1) * c1 + this%Phot3k_ECF_ini(2) * c2 + &
                        this%Phot3k_ECF_ini(3) * c3
 
      d1 = delta1**2 + delta2**2
      d2 = - two * delta2 * (B - delta3)
      d3 = (B - delta3)**2 - delta1**2

      x1 = ( - d2 + dsqrt( d2**2 - four * d1 * d3 ) ) / two / d1
      x2 = ( - d2 - dsqrt( d2**2 - four * d1 * d3 ) ) / two / d1
      phi01 = dasin( dabs( delta1 ) / dsqrt(delta1**2+delta2**2) )
      !phi02 = dacos( delta2 / dsqrt(delta1**2+delta2**2) )
      if(delta1 > zero .and. delta2 > zero)then
          phi02 = phi01
      else if(delta1 > zero .and. delta2 < zero)then
          phi02 = pi - phi01
      else if(delta1 < zero .and. delta2 > zero)then
          phi02 = - phi01
      else if(delta1 < zero .and. delta2 < zero)then
          phi02 = - pi + phi01
      endif
      phi_t = dasin( ( B-delta3 ) / dsqrt(delta1**2+delta2**2) )
      if( B-delta3 > zero )then
          !write(*, fmt="(' ', A10, 4ES20.10)")'ss1=', phi_t - phi02, pi - phi_t - phi02, twopi + phi_t - phi02
      else 
         !write(*, fmt="(' ', A10, 4ES20.10)")'ss2=', twopi + phi_t - phi02, pi - phi_t - phi02
      endif
      if(dabs( ( B-delta3 )  ) / dsqrt(delta1**2+delta2**2) <= one + 1.D-2) &
              !write(*, fmt="(' ', A10, 5ES20.10)")'sdd=', x1, x2, dsin( phi_t - phi02 ), &
              !dsin(pi - phi_t - phi02 ), dsin( twopi + phi_t - phi02 )
 
      CALL Matrix_Multiplication13X33_Sub( this%Phot3k_CF_esti, this%Matrix_CF_2_Elec, this%Phot3k_EF_esti )
      factor = one - this%Elec_v * this%Phot3k_EF_esti(3)
      this%Phot3k_ECF_esti(1) = this%Phot3k_EF_esti(1) / this%Elec_gama / factor
      this%Phot3k_ECF_esti(2) = this%Phot3k_EF_esti(2) / this%Elec_gama / factor
      this%Phot3k_ECF_esti(3) = (this%Phot3k_EF_esti(3) - this%Elec_v) / factor

      !this%PhotE_ECF_esti = this%Elec_gama * dabs( this%Phot4k_CtrCF(1) ) * &
      !                            ( one - this%Elec_V * this%Elec_Phot_mu )
      
      !write(*, *)'s4=', this%Phot3k_ECF_esti(1)**2 + this%Phot3k_ECF_esti(2)**2 + this%Phot3k_ECF_esti(3)**2

      end subroutine Get_Phot4k_phi_As_Egiven_in_CF_Sub
!******************************************************************************************************* 

 

!*******************************************************************************************************
      subroutine Set_Phot_f4_Tetrad_In_Phot_CFrame_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photon_ForEstimation) :: this 
      real(mcp) :: vector1(1: 3)
      real(mcp), dimension(1:3, 1:3) :: InvMatrix_Phot, temM
      real(mcp), dimension(1:4) :: f4_scat_e, f4_In_MF, phot_In_MF
 
      !*******************************************************************************************
 
      this%Phot_f4_AxisZ_PCF = this%Phot4k_CtrCF_At_p(2: 4) / Vector3D_Length( this%Phot4k_CtrCF_At_p(2: 4) ) 

      !vector1 = this%f4_CF(2: 4) !f4_In_MF(2: 4) 
      vector1 = this%f4_CF(2: 4) - this%f4_CF(1) / this%Phot4k_CtrCF_At_p(1) * &
                this%Phot4k_CtrCF_At_p(2: 4)
      this%Phot_f4_AxisX_PCF = vector1 / Vector3D_Length( vector1 )
      !this%Phot_f4_AxisX(1) = this%f4_In_Elec_CF(2)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(2) = this%f4_In_Elec_CF(3)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(3) = this%f4_In_Elec_CF(4)! / this%f4_In_Elec_CF(1)
       
      call this%Vector_Cross_Product( this%Phot_f4_AxisZ_PCF, &
                   this%Phot_f4_AxisX_PCF, this%Phot_f4_AxisY_PCF ) 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**f5555 === ', &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisX_PCF, this%Phot_f4_AxisY_PCF), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisX_PCF, this%Phot_f4_AxisZ_PCF), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisY_PCF, this%Phot_f4_AxisX_PCF), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisZ_PCF, this%Phot_f4_AxisX_PCF), &
      !      Vector4D_Inner_Product_Mski( this%Phot4k_In_Elec_CF, this%f4_In_Elec_CF ) 
      !write(unit = *, fmt = *)'************************************************************'
      
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(1, 1) = this%Phot_f4_AxisX_PCF(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(1, 2) = this%Phot_f4_AxisX_PCF(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(1, 3) = this%Phot_f4_AxisX_PCF(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(2, 1) = this%Phot_f4_AxisY_PCF(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(2, 2) = this%Phot_f4_AxisY_PCF(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(2, 3) = this%Phot_f4_AxisY_PCF(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(3, 1) = this%Phot_f4_AxisZ_PCF(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(3, 2) = this%Phot_f4_AxisZ_PCF(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF(3, 3) = this%Phot_f4_AxisZ_PCF(3)

      end subroutine Set_Phot_f4_Tetrad_In_Phot_CFrame_Sub

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module ModuleForCompScatEsti




 
