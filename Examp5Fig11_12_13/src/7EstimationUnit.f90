      module PhotonsEstimation
      use ScatDistance_FlatSP 
      implicit none 

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      type, public, extends(Photon_With_ScatDistance_FlatSP) :: Photons_Esti 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(1: vL_sc_up)  
          real(mcp) :: nu_obs 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !real(mcp) :: PolarArrayd0(0: Num_PolDeg)
          real(mcp) :: PolarArrayI0(0: Num_PolDeg)
          real(mcp) :: PolarArrayQ0(0: Num_PolDeg)
          real(mcp) :: PolarArrayU0(0: Num_PolDeg)
          real(mcp) :: PolarArrayV0(0: Num_PolDeg)
          real(mcp) :: PolarArrayIr0(0: Num_PolDeg)
          real(mcp) :: PolarArrayIl0(0: Num_PolDeg)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !real(mcp) :: PolarArrayd90(0: Num_PolDeg)
          real(mcp) :: PolarArrayI90(0: Num_PolDeg)
          real(mcp) :: PolarArrayQ90(0: Num_PolDeg)
          real(mcp) :: PolarArrayU90(0: Num_PolDeg)
          real(mcp) :: PolarArrayV90(0: Num_PolDeg)
          real(mcp) :: PolarArrayIr90(0: Num_PolDeg)
          real(mcp) :: PolarArrayIl90(0: Num_PolDeg)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !real(mcp) :: PolarArrayd180(0: Num_PolDeg)
          real(mcp) :: PolarArrayI180(0: Num_PolDeg)
          real(mcp) :: PolarArrayQ180(0: Num_PolDeg)
          real(mcp) :: PolarArrayU180(0: Num_PolDeg)
          real(mcp) :: PolarArrayV180(0: Num_PolDeg)
          real(mcp) :: PolarArrayIr180(0: Num_PolDeg)
          real(mcp) :: PolarArrayIl180(0: Num_PolDeg)
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          real(mcp) :: PolarArrayIQUV10(1: 4, 0: Num_phi) = zero
          real(mcp) :: PolarArrayIQUV30(1: 4, 0: Num_Phi) = zero
          real(mcp) :: PolarArrayIQUV60(1: 4, 0: Num_Phi) = zero
          real(mcp) :: PolarArrayIQUV80(1: 4, 0: Num_Phi) = zero
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: PolarArrayIQUV0(1: 4, 0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayIQUV90(1: 4, 0: Num_PolDeg) = zero
          real(mcp) :: PolarArrayIQUV180(1: 4, 0: Num_PolDeg) = zero
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          integer :: times_counter(1: 4) = 0
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          real(mcp) :: d_theta, d_phi, d_tau
          logical :: first_time_recording

      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
          procedure, public :: Calc_Phot_Informations_At_Observer_Flat_SPTM_Reflec   =>   &
                               Calc_Phot_Informations_At_Observer_Flat_SPTM_Reflec_Sub 
          procedure, public :: Calc_Phot_Informations_At_Observer_Diffuse_Reflec   =>   &
                               Calc_Phot_Informations_At_Observer_Diffuse_Reflec_Sub 
          procedure, public :: Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi   =>   &
                               Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub
          procedure, public :: Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat   =>   &
                               Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat_Sub
          procedure, public :: Get_Psi_IQUV_for_Estimat_With_muphi_Geiven    =>    &
                               Get_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
          procedure, public :: Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven    =>   &
                               Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
          procedure, public :: Calc_Phot_Inform_At_Observer_with_mu_phi_Given    =>   &
                               Calc_Phot_Inform_At_Observer_with_mu_phi_Given_Sub
          procedure, public :: Calc_Phot_Inform_At_Observer_with_mu_phi_Given2    =>   &
                               Calc_Phot_Inform_At_Observer_with_mu_phi_Given2_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photons_Esti
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Calc_Phot_Informations_At_Observer_Flat_SPTM_Reflec_Sub
      private :: Calc_Phot_Informations_At_Observer_Diffuse_Reflec_Sub
      private :: Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub
      private :: Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat_Sub
      private :: Get_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
      private :: Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
      private :: Calc_Phot_Inform_At_Observer_with_mu_phi_Given_Sub
      private :: Calc_Phot_Inform_At_Observer_with_mu_phi_Given2_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains   
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observer_Flat_SPTM_Reflec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer :: i, j, h, k, mu_i, N_low, N_low1, i_obs, i_phi, case_phi
      real(mcp) :: temp_P_X_f, vector1(1: 3), f4(1: 3), Axis_x(1: 3), &
                    Axis_y(1: 3), Axis_z(1: 3), cos_psi, Qpsi, Upsi, fx, fy, &
                    Axis_z1(1: 3), phi_p, vLv  

      this%d_theta = one / Num_PolDeg
      this%d_phi =  twopi / Num_Phi 

      if( this%InterSection_Cases == -1 .or. this%InterSection_Cases == -3 )then
          return
      else

          phi_p = datan( dabs( this%Vector_of_Momentum_ini(2) / this%Vector_of_Momentum_ini(1) ) )
          if( this%Vector_of_Momentum_ini(1) > zero .and. this%Vector_of_Momentum_ini(2) >= zero )then
          else if( this%Vector_of_Momentum_ini(1) < zero .and. this%Vector_of_Momentum_ini(2) > zero )then
              phi_p = pi - phi_p
          else if( this%Vector_of_Momentum_ini(1) < zero .and. this%Vector_of_Momentum_ini(2) < zero )then
              phi_p = pi + phi_p
          else if( this%Vector_of_Momentum_ini(1) > zero .and. this%Vector_of_Momentum_ini(2) < zero )then
              phi_p = twopi - phi_p
          else if( this%Vector_of_Momentum_ini(1) == zero .and. this%Vector_of_Momentum_ini(2) > zero )then
              phi_p = pi / two
          else if( this%Vector_of_Momentum_ini(1) == zero .and. this%Vector_of_Momentum_ini(2) < zero )then
              phi_p = pi * three / two
          else if( this%Vector_of_Momentum_ini(1) > zero .and. this%Vector_of_Momentum_ini(2) == zero )then
              phi_p = zero
          else if( this%Vector_of_Momentum_ini(1) < zero .and. this%Vector_of_Momentum_ini(2) == zero )then
              phi_p = pi 
          endif
 
          i_phi = floor( phi_p / this%d_phi ) 
   
      !if( i_phi == 0 .or. i_phi == 1 .or. i_phi == Num_phi-1 )then
      !    case_phi = 1
      !    this%effect_number = this%effect_number + 1 
      !else if( i_phi == Num_phi / 4 .or. abs( i_phi - Num_phi / 4 ) == 1 )then
      !    case_phi = 2
      !else if( i_phi == Num_phi * 3 / 4 .or. abs( i_phi - Num_phi * 3 / 4 ) == 1 )then
      !    case_phi = 3
      !else if( i_phi == Num_phi / 2 .or. abs( i_phi - Num_phi / 2 ) == 1 )then
      !    case_phi = 4
      !else
      !    return
      !endif 

      if( i_phi == 0 )then
          case_phi = 1
          this%effect_number = this%effect_number + 1 
      else if( i_phi == Num_phi / 4 )then
          case_phi = 2 
      else if( i_phi == Num_phi / 2 )then
          case_phi = 4
      else
          return
      endif 

          i_obs = floor( dabs(this%Vector_of_Momentum_ini(3)) / this%d_theta )
  
          if(i_obs < 0)write(*, *)'ss===', i_obs, dabs(this%Vector_of_Momentum_ini(3)), &
                     Dacos( dabs(this%Vector_of_Momentum_ini(3)) ), this%d_theta, &
                      this%Vector_of_Momentum_ini(3), &
                 Vector4D_Inner_Product_Mski( this%Phot4k_CtrCF_ini, this%f4_CF ), &
                 Vector4D_Inner_Product_Mski( this%Phot4k_CtrCF_ini, this%Phot4k_CtrCF_ini ), &
                 Vector4D_Inner_Product_Mski( this%f4_CF, this%f4_CF ), this%scatter_times

          if( this%r_one_hvlmec2_one_cosE /= zero )this%w_ini = this%w_ini * this%r_one_hvlmec2_one_cosE 
          vLv = this%w_ini * dexp( - this%Optical_Depth_scatter ) / dabs(this%Vector_of_Momentum_ini(3))

          if( case_phi == 1 )then
              this%PolarArrayI0( i_obs ) = this%PolarArrayI0( i_obs ) + vLv
          else if( case_phi == 2 .or. case_phi == 3 )then 
              this%PolarArrayI90( i_obs ) = this%PolarArrayI90( i_obs ) + vLv
          else if( case_phi == 4 )then 
              this%PolarArrayI180( i_obs ) = this%PolarArrayI180( i_obs ) + vLv
          endif   
  
          if( this%delta_pd /= zero )then
              !Axis_z1 = this%Phot4k_CtrCF_ini(2: 4) / dabs( this%Phot4k_CtrCF_ini(1) ) 
              Axis_z = this%Phot4k_CtrCF_ini(2: 4) / Vector3D_Length( this%Phot4k_CtrCF_ini(2: 4) )
              !write(*, *)'ff', Axis_z1 - Axis_z

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

              f4 = this%f4_CF(2: 4) - this%f4_CF(1) / this%Phot4k_CtrCF_ini(1) * &
                      this%Phot4k_CtrCF_ini(2: 4)  !! | f4 | = 1 

              cos_psi = Vector3D_Inner_Product( f4, Axis_x )
              fx = Vector3D_Inner_Product( f4, Axis_x )
              fy = Vector3D_Inner_Product( f4, Axis_y )
              if( fx > zero .and. fy > zero  )then
                  Qpsi = this%delta_pd * ( two * cos_psi ** 2 - one )
                  Upsi = this%delta_pd * two * dsqrt( dabs( one - cos_psi**2 ) ) * cos_psi
                  h = 1
              else if( fx < zero .and. fy > zero  )then
                  Qpsi = this%delta_pd * ( two * cos_psi ** 2 - one )
                  Upsi = this%delta_pd * two * dsqrt( dabs( one - cos_psi**2 ) ) * cos_psi
                  h = 2
              else if( fx < zero .and. fy < zero  )then
                  Qpsi = this%delta_pd * ( two * cos_psi ** 2 - one )
                  Upsi = this%delta_pd * two * dsqrt( dabs( one - cos_psi**2 ) ) * dabs( cos_psi )
                  h = 3
              else if( fx > zero .and. fy < zero  )then
                  Qpsi = this%delta_pd * ( two * cos_psi ** 2 - one )
                  Upsi = - this%delta_pd * two * dsqrt( dabs( one - cos_psi**2 ) ) * cos_psi
                  h = 4
              endif

              if( case_phi == 1 )then
                  this%PolarArrayQ0( i_obs ) = this%PolarArrayQ0( i_obs ) + Qpsi * vLv
                  this%PolarArrayU0( i_obs ) = this%PolarArrayU0( i_obs ) + Upsi * vLv
                  !this%PolarArrayV0( i_obs ) = this%PolarArrayV0( i_obs ) + this%V_sp * vLv
              else if( case_phi == 2 )then
                  this%PolarArrayQ90( i_obs ) = this%PolarArrayQ90( i_obs ) + Qpsi * vLv
                  this%PolarArrayU90( i_obs ) = this%PolarArrayU90( i_obs ) + Upsi * vLv
                  !this%PolarArrayV90( i_obs ) = this%PolarArrayV90( i_obs ) + this%V_sp * vLv
              else if( case_phi == 3 )then
                  this%PolarArrayQ90( i_obs ) = this%PolarArrayQ90( i_obs ) + Qpsi * vLv
                  this%PolarArrayU90( i_obs ) = this%PolarArrayU90( i_obs ) - Upsi * vLv
                  !this%PolarArrayV90( i_obs ) = this%PolarArrayV90( i_obs ) + this%V_sp * vLv
              else if( case_phi == 4 )then
                  this%PolarArrayQ180( i_obs ) = this%PolarArrayQ180( i_obs ) + Qpsi * vLv
                  this%PolarArrayU180( i_obs ) = this%PolarArrayU180( i_obs ) + Upsi * vLv
                  !this%PolarArrayV180( i_obs ) = this%PolarArrayV180( i_obs ) + this%V_sp * vLv
              endif  
 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          !this%effect_number = this%effect_number + 1 

          endif 
      endif
      return
      end subroutine Calc_Phot_Informations_At_Observer_Flat_SPTM_Reflec_Sub

!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observer_Diffuse_Reflec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer :: i, j, h, k, i_obs, i_phi, case_phi  
      real(mcp) :: vLv

      !this%d_theta = one / Num_PolDeg
      !this%d_phi =  twopi / Num_Phi 
      !write(*, *)'ffs=', one / Num_PolDeg, twopi / Num_Phi, this%d_theta, this%d_phi

      if( this%InterSection_Cases == -1 .or. this%InterSection_Cases == -3 )then
          return
      else
  
          i_phi = floor( this%phi_ini / this%d_phi ) 
    

          if( i_phi == 0 .or. i_phi == Num_phi - 1 )then
              case_phi = 1
              this%effect_number = this%effect_number + 1 
          else if( i_phi == Num_phi / 4 .or. i_phi == Num_phi / 4 - 1 )then
              case_phi = 2 
          else if( i_phi == 3 * Num_phi / 2 .or. i_phi == 3 * Num_phi / 2 - 1 )then
              case_phi = 2 
          else if( i_phi == Num_phi / 2 .or. i_phi == Num_phi / 2 - 1  )then
              case_phi = 4
          else
              return
          endif 

          i_obs = floor( dabs(this%Vector_of_Momentum_ini(3)) / this%d_theta )
  
          if(i_obs < 0)write(*, *)'ss===', i_obs, dabs(this%Vector_of_Momentum_ini(3)), &
                     Dacos( dabs(this%Vector_of_Momentum_ini(3)) ), this%d_theta, &
                      this%Vector_of_Momentum_ini(3), &
                 Vector4D_Inner_Product_Mski( this%Phot4k_CtrCF_ini, this%f4_CF ), &
                 Vector4D_Inner_Product_Mski( this%Phot4k_CtrCF_ini, this%Phot4k_CtrCF_ini ), &
                 Vector4D_Inner_Product_Mski( this%f4_CF, this%f4_CF ), this%scatter_times

          !if( this%r_one_hvlmec2_one_cosE /= zero )this%w_ini = this%w_ini * this%r_one_hvlmec2_one_cosE 
          vLv = dexp( - this%Optical_Depth_scatter ) / dabs(this%Vector_of_Momentum_ini(3))

          !write(*, *)'ffs=', vLv, this%Psi_I, this%Psi_Q, this%Psi_U, this%Psi_V, &
           !    dexp( - this%Optical_Depth_scatter ), dabs(this%Vector_of_Momentum_ini(3))
          if( case_phi == 1 )then
              this%PolarArrayI0( i_obs ) = this%PolarArrayI0( i_obs ) + this%Psi_I * vLv 
              this%PolarArrayQ0( i_obs ) = this%PolarArrayQ0( i_obs ) + this%Psi_Q * vLv
              this%PolarArrayU0( i_obs ) = this%PolarArrayU0( i_obs ) + this%Psi_U * vLv
              this%PolarArrayV0( i_obs ) = this%PolarArrayV0( i_obs ) + this%Psi_V * vLv
          else if( case_phi == 2 .or. case_phi == 3 )then 
              this%PolarArrayI90( i_obs ) = this%PolarArrayI90( i_obs ) + this%Psi_I * vLv
              this%PolarArrayQ90( i_obs ) = this%PolarArrayQ90( i_obs ) + this%Psi_Q * vLv
              this%PolarArrayU90( i_obs ) = this%PolarArrayU90( i_obs ) + this%Psi_U * vLv
              this%PolarArrayV90( i_obs ) = this%PolarArrayV90( i_obs ) + this%Psi_V * vLv
          else if( case_phi == 4 )then 
              this%PolarArrayI180( i_obs ) = this%PolarArrayI180( i_obs ) + this%Psi_I * vLv
              this%PolarArrayQ180( i_obs ) = this%PolarArrayQ180( i_obs ) + this%Psi_Q * vLv
              this%PolarArrayU180( i_obs ) = this%PolarArrayU180( i_obs ) + this%Psi_U * vLv
              this%PolarArrayV180( i_obs ) = this%PolarArrayV180( i_obs ) + this%Psi_V * vLv
          endif    
          !this%effect_number = this%effect_number + 1  
      endif
      return
      end subroutine Calc_Phot_Informations_At_Observer_Diffuse_Reflec_Sub

  
!*******************************************************************************************************
      subroutine Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub(this)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer :: i, j, k, i_obs, i_phi, case_mu
      real(mcp) :: i_mu  
      real(mcp) :: vLv

      !this%d_theta = one / Num_PolDeg
      !this%d_phi =  twopi / Num_Phi 
      !write(*, *)'ffs=', one / Num_PolDeg, twopi / Num_Phi, this%d_theta, this%d_phi

      if( this%InterSection_Cases == -1 .or. this%InterSection_Cases == -3 )then
          return
      else

          i_mu = floor( dabs(this%Vector_of_Momentum_ini(3)) / this%d_theta ) 

          if( i_mu == 5 )then
              case_mu = 1 
              !this%times_counter(1) = this%times_counter(1) + 1
          else if( i_mu == 30 )then
              case_mu = 2 
              !this%times_counter(2) = this%times_counter(2) + 1
          else if( i_mu == 55 )then
              case_mu = 3 
              !this%times_counter(3) = this%times_counter(3) + 1
          else if( i_mu == 85 )then
              case_mu = 4
              !this%times_counter(4) = this%times_counter(4) + 1
          else
              return
          endif 

          i_phi = floor( this%phi_ini / this%d_phi ) 
    
          vLv = dexp( - this%Optical_Depth_scatter ) / dabs(this%Vector_of_Momentum_ini(3))

          write(*, *)'ffs=', vLv, this%Psi_I, this%Psi_Q, this%Psi_U, this%Psi_V, &
              dexp( - this%Optical_Depth_scatter ), dabs(this%Vector_of_Momentum_ini(3))
          if( case_mu == 1 )then
              this%PolarArrayIQUV10(1, i_phi ) = this%PolarArrayIQUV10(1, i_phi ) + this%Psi_I * vLv 
              this%PolarArrayIQUV10(2, i_phi ) = this%PolarArrayIQUV10(2, i_phi ) + this%Psi_Q * vLv
              this%PolarArrayIQUV10(3, i_phi ) = this%PolarArrayIQUV10(3, i_phi ) + this%Psi_U * vLv
              this%PolarArrayIQUV10(4, i_phi ) = this%PolarArrayIQUV10(4, i_phi ) + this%Psi_V * vLv
          else if( case_mu == 2 )then 
              this%PolarArrayIQUV30(1, i_phi ) = this%PolarArrayIQUV30(1, i_phi ) + this%Psi_I * vLv 
              this%PolarArrayIQUV30(2, i_phi ) = this%PolarArrayIQUV30(2, i_phi ) + this%Psi_Q * vLv
              this%PolarArrayIQUV30(3, i_phi ) = this%PolarArrayIQUV30(3, i_phi ) + this%Psi_U * vLv
              this%PolarArrayIQUV30(4, i_phi ) = this%PolarArrayIQUV30(4, i_phi ) + this%Psi_V * vLv
          else if( case_mu == 3 )then 
              this%PolarArrayIQUV60(1, i_phi ) = this%PolarArrayIQUV60(1, i_phi ) + this%Psi_I * vLv 
              this%PolarArrayIQUV60(2, i_phi ) = this%PolarArrayIQUV60(2, i_phi ) + this%Psi_Q * vLv
              this%PolarArrayIQUV60(3, i_phi ) = this%PolarArrayIQUV60(3, i_phi ) + this%Psi_U * vLv
              this%PolarArrayIQUV60(4, i_phi ) = this%PolarArrayIQUV60(4, i_phi ) + this%Psi_V * vLv
          else if( case_mu == 4 )then 
              this%PolarArrayIQUV80(1, i_phi ) = this%PolarArrayIQUV80(1, i_phi ) + this%Psi_I * vLv 
              this%PolarArrayIQUV80(2, i_phi ) = this%PolarArrayIQUV80(2, i_phi ) + this%Psi_Q * vLv
              this%PolarArrayIQUV80(3, i_phi ) = this%PolarArrayIQUV80(3, i_phi ) + this%Psi_U * vLv
              this%PolarArrayIQUV80(4, i_phi ) = this%PolarArrayIQUV80(4, i_phi ) + this%Psi_V * vLv
          endif    
          !this%effect_number = this%effect_number + 1  
      endif
      return
      end subroutine Calc_Phot_Inform_At_Observer_Diffuse_Reflec_phi_Sub

!*******************************************************************************************************

  
!*******************************************************************************************************
      subroutine Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat_Sub(this, mu_i)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer, intent(in) :: mu_i
      integer :: i, j, k, i_obs, i_phi 
      real(mcp) :: i_mu  
      real(mcp) :: vLv
 
          !i_mu = floor( dabs(this%Vector_of_Momentum_ini(3)) / this%d_theta ) 
 

      i_phi = floor( this%phi_estimat(mu_i) / this%d_phi ) 
    
      this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(mu_i)
      vLv = dexp( - this%Optical_Depth_scatter ) / this%mu_estimat(mu_i)
      !this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(1)!p_out1 

      !write(*, *)'tt3==', this%f_IQUV_estimat, vLv, this%Optical_Depth_scatter 
      if(mu_i == 1)then
          this%PolarArrayIQUV10(1: 4, i_phi ) = this%PolarArrayIQUV10(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv  
      else if(mu_i == 2)then
          this%PolarArrayIQUV30(1: 4, i_phi ) = this%PolarArrayIQUV30(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv  
      else if(mu_i == 3)then
          this%PolarArrayIQUV60(1: 4, i_phi ) = this%PolarArrayIQUV60(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv   
      else if(mu_i == 4)then
          this%PolarArrayIQUV80(1: 4, i_phi ) = this%PolarArrayIQUV80(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv 
      endif 
          !write(*, *)'tt3==', this%PolarArrayIQUV10(1: 4, i_phi ), i_phi
      return
      end subroutine Calc_Phot_Inform_At_Observer_DiffRefl_phi_Estimat_Sub

!*******************************************************************************************************

  
!*******************************************************************************************************
      subroutine Calc_Phot_Inform_At_Observer_with_mu_phi_Given_Sub(this, mu_i)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer, intent(in) :: mu_i
      integer :: i, j, k, i_phi 
      real(mcp) :: i_mu  
      real(mcp) :: vLv
  
      vLv = dexp( - this%z_tau / this%mu_estimat(mu_i) ) / this%mu_estimat(mu_i)
      do i_phi = 0, Num_Phi 
          call this%Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven( this%mu_estimat(mu_i), &
                                                             this%phi_estimates( i_phi ) )
          !this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(1)!p_out1 

          !write(*, *)'tt3==', this%f_IQUV_estimat, vLv, this%Optical_Depth_scatter 
          if(mu_i == 1)then
              this%PolarArrayIQUV10(1: 4, i_phi ) = this%PolarArrayIQUV10(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv  
          else if(mu_i == 2)then
              this%PolarArrayIQUV30(1: 4, i_phi ) = this%PolarArrayIQUV30(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv  
          else if(mu_i == 3)then
              this%PolarArrayIQUV60(1: 4, i_phi ) = this%PolarArrayIQUV60(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv   
          else if(mu_i == 4)then
              this%PolarArrayIQUV80(1: 4, i_phi ) = this%PolarArrayIQUV80(1: 4, i_phi ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv 
          endif  
      enddo
      return
      end subroutine Calc_Phot_Inform_At_Observer_with_mu_phi_Given_Sub

  
!*******************************************************************************************************
      subroutine Calc_Phot_Inform_At_Observer_with_mu_phi_Given2_Sub(this, phi_i)
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this   
      integer, intent(in) :: phi_i
      integer :: i, j, k, i_mu  
      real(mcp) ::  vLv, phi_esti
  
      phi_esti = this%phi_estimat( phi_i )
      do i_mu = 0, Num_PolDeg
          call this%Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven( this%mu_estimates( i_mu ), &
                                                              phi_esti ) 
          vLv = dexp( - this%z_tau / this%mu_estimates(i_mu) ) / this%mu_estimates(i_mu)

          !write(*, *)'tt3==', this%f_IQUV_estimat, vLv, this%Optical_Depth_scatter 
          if(phi_i == 1)then
              this%PolarArrayIQUV0(1: 4, i_mu ) = this%PolarArrayIQUV0(1: 4, i_mu ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv  
          else if(phi_i == 2)then
              this%PolarArrayIQUV90(1: 4, i_mu ) = this%PolarArrayIQUV90(1: 4, i_mu ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv 
          !write(*, *)'tt3==', this%f_IQUV_estimat, vLv!, this%Optical_Depth_scatter  
          !else if(phi_i == 3)then
          !    this%PolarArrayIQUV60(1: 4, i_mu ) = this%PolarArrayIQUV60(1: 4, i_mu ) + &
          !                                       this%f_IQUV_estimat(1: 4) * vLv   
          else if(phi_i == 3)then
              this%PolarArrayIQUV180(1: 4, i_mu ) = this%PolarArrayIQUV180(1: 4, i_mu ) + &
                                                 this%f_IQUV_estimat(1: 4) * vLv 
          endif  
      enddo
      return
      end subroutine Calc_Phot_Inform_At_Observer_with_mu_phi_Given2_Sub


!*******************************************************************************************************
      Subroutine Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub( this, mu_estimat, scat_estimat )
!*******************************************************************************************************
      implicit none
      class(Photons_Esti) :: this
      !integer, intent(in) :: mu_i
      real(mcp), intent(in) :: mu_estimat, scat_estimat
      real(mcp) :: r, phip, smu, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp, beta, &
                   A_const, B_const, C_const, phi_ini, N_normal, N_temp2, &
                   cosphi_ini, sinphi_ini, cos2phi_ini, sin2phi_ini, smup, &
                   mu2, mup2, mup, mu

      real(mcp) :: A1, B1, t1, t2, Norm_C, f1, f2, f3, f4, f5, F0, g1, g2, g3, g4, g5
      real(mcp) :: UlI, QlI, sinDphi, cosDphi, sin2Dphi, cos2Dphi 
      real(mcp) :: P0(1: 3, 1: 3) = zero, P1(1: 3, 1: 3) = zero, P2(1: 3, 1: 3) = zero, &
                   PT(1: 3, 1: 3) = zero
      real(mcp) :: f_Q, f_U, f_V
      !****************************************************************************
      !real(mcp) :: a, b, c, d
      !Complex*16 :: roots3(1: 3), rts(1: 3)
      !integer :: del, cases_f
      !**************************************************************************** 
   
      mup = this%Vector_of_Momentum_ini(3)
      mup2 = mup**2
      smup = dsqrt(one - mup2)
      QlI = this%Psi_Q / this%Psi_I
      UlI = this%Psi_U / this%Psi_I
      !t1 = one + mup2 - (one - mup2) * QlI
      !t2 = two * (one - mup2) * (one + QlI)
      !A1 = t1 - t2
      !B1 = t1 + t2  
      mu = mu_estimat !this%mu_estimat( mu_i )
      mu2 = mu**2 
      smu = dsqrt( one - mu2 ) 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !f1 = (mu2 * A1 + B1) !t1 * (one + mu2) + t2 * (one - mu2)
      !f2 = four * mup * smup * mu * smu * (one + QlI)
      !f3 = four * smup * mu * smu * UlI
      !f4 = (one - mu2) * (one - mup2) - (one - mu2) * (one + mup2) * QlI
      !f5 = - two * (one - mu2) * mup * UlI 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !g2 = f2 * this%cosphi_ini + f3 * this%sinphi_ini
      !g3 = f2 * this%sinphi_ini - f3 * this%cosphi_ini
      !g4 = f4 * this%cos2phi_ini + f5 * this%sin2phi_ini
      !g5 = f4 * this%sin2phi_ini - f5 * this%cos2phi_ini 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      sinDphi =  dsin(this%phi_ini - scat_estimat)
      cosDphi =  dcos(this%phi_ini - scat_estimat)
      sin2Dphi =  dsin(two*(this%phi_ini - scat_estimat))
      cos2Dphi =  dcos(two*(this%phi_ini - scat_estimat))   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      P0(1, 1) = two * (one - mu2)*(one - mup2) + mu2*mup2
      P0(1, 2) = mu2
      P0(2, 1) = mup2
      P0(2, 2) = one
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      P1(1, 1) = four * mu * mup * CosDphi
      P1(1, 3) = two * mu * SinDphi
      P1(3, 1) = - two * mup * SinDphi
      P1(3, 3) = CosDphi 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      P2(1, 1) = mu2 * mup2 * Cos2Dphi
      P2(1, 2) = - mu2 * Cos2Dphi
      P2(1, 3) = mu2 * mup * Sin2Dphi
      P2(2, 1) = - mup2 * Cos2Dphi
      P2(2, 2) = Cos2Dphi
      P2(2, 3) = - mup * Sin2Dphi
      P2(3, 1) = - mu * mup2 * Sin2Dphi
      P2(3, 2) = mu * Sin2Dphi
      P2(3, 3) = mu * mup * Cos2Dphi
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      PT = P0 + P1 * smu * smup + P2
      PT(3, 1:3) = PT(3, 1:3) * two
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      !N_temp = f1 + g2 * dcos(scat_estimat) + g3 * dsin(scat_estimat) + &
      !              g4 * dcos(two * scat_estimat) + g5 * dsin(two * scat_estimat) 
      !Norm_phi = N_temp / twopi / f1 
      !N_temp2 = f1 + f2 * CosDphi + f3 * SinDphi + f4 * Cos2Dphi + f5 * Sin2Dphi
      N_temp = PT(1, 1) + PT(1, 2) + PT(2, 1) + PT(2, 2) + &
                            ( PT(1, 1) - PT(1, 2) + PT(2, 1) - PT(2, 2) ) * QlI + &
                            two * ( PT(1, 3) + PT(2, 3) ) * UlI 
    
      !write(*, fmt="(' ', F22.15)") N_temp2- N_temp!  p2 / Norm_phi * 3.D0 / 32.D0 /pi
      this%f_IQUV_estimat(1) = N_temp * this%Psi_I !N_temp2 * this%Psi_I / Norm_phi
      if( N_temp < zero )then
           write(*, *)'mms55==',  UlI, QlI !, mu, mup, scat_estimat
           !write(*, *)'mmsa77==', pb1, pb2, pb3, pb4, pb5, pb1 + pb2 + pb3 + pb4 + pb5, f4, f5
           !write(*, *)'mmsa88==', f1, f2, f3, f4, f5, UlI, QlI
           stop
      endif

      !temp3 = ( ( one -mu2 )*( one - mup2 )*(two + 3.D0*QlI) - &
      !                ( one + mup2 )*( one - mu2 ) + &
      !                four * mup * mu * smu * smup * (one + QlI) * CosDphi + &
      !                ( one + mu2 )*( mup2 - one + (mup2 + one)*QlI ) * Cos2Dphi + &
      !                four * mu * smu * smup * UlI * SinDphi + &
      !                two*( one + mu2 ) * mup * UlI * Sin2Dphi &
      !              ) !/ N_temp * this%Psi_I 

      f_Q =(         PT(1, 1) + PT(1, 2) - PT(2, 1) - PT(2, 2) + &
             QlI * ( PT(1, 1) - PT(1, 2) - PT(2, 1) + PT(2, 2) )  + &
             UlI * ( PT(1, 3) - PT(2, 3) ) * two ) ! / N_temp 
     ! write(*, *)'tt3==', temp3 , f_Q
      this%f_IQUV_estimat(2) = f_Q * this%Psi_I ! / Norm_phi

      !this%Psi_U = ( - two * mup * smu * smup * SinDphi + &
      !               UlI * smu * smup * CosDphi + &
      !               (QlI - mup**2) * mu * Sin2Dphi + &
      !               UlI * mup * mu * Cos2Dphi  ) * &
      !               four * this%Psi_I / N_temp 
      !write(*, *)'tt3==', this%Psi_U / this%Psi_I , this%Psi_U

      f_U =  ( PT(3, 1) + PT(3, 2) + QlI * ( PT(3, 1) - PT(3, 2) ) + UlI * PT(3, 3)*two ) ! / N_temp  
      this%f_IQUV_estimat(3) = f_U * this%Psi_I  ! / Norm_phi

      f_V = ( mup * mu + smu * smup * cosDphi )! * three / 8.D0 / pi ! / N_temp * four
      this%f_IQUV_estimat(4) = f_V * this%Psi_V !* 3.D0 * f1 / 16.D0 / N_temp !/ Norm_phi 
      !**********************************************************************************    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      !Call this%Get_Psi_IQUV_for_Estimat_With_muphi_Geiven( mu, mu2, smu, mup, mup2, smup, &
      !            scat_estimat, f1, g2, g3, g4, g5, QlI, UlI )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
      !*************************************************************************************************  
      end Subroutine Calcu_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
!*******************************************************************************************************


!*******************************************************************************************************
      Subroutine Get_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub( this, mu, mu2, smu, mup, mup2, smup, &
                  scat_phi, f1, g2, g3, g4, g5, QlI, UlI )
!*******************************************************************************************************
      class(Photons_Esti) :: this 
      real(mcp), intent(in) ::  mu, mu2, smu, mup, mup2, smup, scat_phi, f1, g2, g3, g4, g5, QlI, UlI 

      end Subroutine Get_Psi_IQUV_for_Estimat_With_muphi_Geiven_Sub
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module PhotonsEstimation




 
