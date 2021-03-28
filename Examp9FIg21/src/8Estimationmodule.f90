      module EstimationsModule
      use ScatDistance_FlatSP 
      use PhotonEmitterBB 
      implicit none 


      type, public, extends(Photon_With_ScatDistance_FlatSP) :: PhotonEstimate

          real(mcp) :: PolArrImu(1: 8, 0: 6, 0: vL_sc_up) = zero   
      contains
          procedure, public :: Get_K_P1_P2_Scat_Kernel_for_Estimation   =>   &
                               Get_K_P1_P2_Scat_Kernel_for_Estimation_Sub
      end type PhotonEstimate

      private :: Get_K_P1_P2_Scat_Kernel_for_Estimation_Sub
   
      contains

!*******************************************************************************************************
      subroutine Get_K_P1_P2_Scat_Kernel_for_Estimation_Sub(this, mu_i)
!*******************************************************************************************************
      implicit none
      class(PhotonEstimate) :: this
      integer, intent(in) :: mu_i 
      integer :: i, i_phi, i_E, i_1, i_2, scat_times
      real(mcp) :: vLv, Lv, v_esti, J_esti
      real(mcp) :: mu_ep, smu_ep, mu_psi, epsip, epsi, KN_CrossSection, &
                   nup_vs_nu, chi, factor1, factor2, SigmaKN, sigma_temp
        
      do i_phi = 0,  N_esti_phi
          this%Phot3k_CF_esti(1) = this%smu_esti(mu_i) * this%cos_phi_esti(i_phi)
          this%Phot3k_CF_esti(2) = this%smu_esti(mu_i) * this%sin_phi_esti(i_phi)
          this%Phot3k_CF_esti(3) = this%mu_esti(mu_i) 
 
   
          CALL Matrix_Multiplication13X33_Sub( this%Phot3k_CF_esti, this%Matrix_CF_2_Elec, &
                         this%Phot3k_EF_esti )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          mu_ep = this%Phot3k_EF_esti(3)
          smu_ep = dsqrt( one - this%Phot3k_EF_esti(3)**2 )  
          mu_psi = this%Phot3k_EF_esti(3)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          nup_vs_nu = one / (one + this%E_ini / mec2 * (one - mu_psi) )
          this%E_esti_vs_mu_phi = this%E_ini * nup_vs_nu
          sigma_temp = this%sigma_KN( this%E_esti_vs_mu_phi ) * this%n_e
 
          chi = nup_vs_nu + one / nup_vs_nu - ( one - this%Phot3k_EF_esti(3)**2 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          epsi = two * this%E_ini / mec2 
          if ( epsi > 8.19D-4 ) then
              SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
          else 
              SigmaKN = one - epsi
          end if 
          KN_CrossSection = half_re2 * chi / SigmaKN * (nup_vs_nu )**2
          !write(*, *)'s1=', epsi, this%Phot4k_In_Elec_CF(1) * two / mec2, KN_CrossSection_ECF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          Lv = this%w_ini * dexp( - this%z_tau * sigma_temp / this%mu_esti(mu_i) ) / &
                        this%mu_esti(mu_i) * KN_CrossSection  

          i_1 = floor( ( this%E_esti_vs_mu_phi - this%y1 ) / (this%dy/two) )
          if(mod(i_1, 2)==0)then
              i_E = i_1 / 2
          else
              i_E = (i_1 + 1) / 2
          endif
          if( i_E > vL_sc_up .or. i_E < 0 )cycle
          scat_times = this%scatter_times + 1
 
          this%PolArrImu(mu_i, 6, i_E) = this%PolArrImu(mu_i, 6, i_E) + Lv
          if(scat_times <= 4)then 
              this%PolArrImu(mu_i, scat_times, i_E) = this%PolArrImu(mu_i, scat_times, i_E) + Lv
          else 
              this%PolArrImu(mu_i, 5, i_E) = this%PolArrImu(mu_i, 5, i_E) + Lv
          endif 
      enddo  
      return
      end subroutine Get_K_P1_P2_Scat_Kernel_for_Estimation_Sub
!*******************************************************************************************************

 
       end module EstimationsModule




