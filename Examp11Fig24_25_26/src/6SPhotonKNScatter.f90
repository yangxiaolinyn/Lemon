      module ScatterPhoton_KN 
      use PhotonEmitter
      implicit none 

      type, public, extends(Photon_Emitter) :: ScatPhoton_KN
          real(mcp), dimension(1:3) :: PhotAxisX
          real(mcp), dimension(1:3) :: PhotAxisY
          real(mcp), dimension(1:3) :: PhotAxisZ
          real(mcp), dimension(1:4) :: Phot4k_In_Elec_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_CovCF
          real(mcp), dimension(1:3, 1:3) :: Matrix_Of_Tetrad_Of_PhotAxis
          real(mcp), dimension(1:4) :: f4_In_Elec_CF
          !*****************************************************************
          real(mcp), dimension(1:3) :: temp_px = -one
          real(mcp), dimension(1:3) :: temp_py = -one
          real(mcp), dimension(1:3) :: temp_pz = -one
          real(mcp), dimension(1:3) :: temp_f4 = -one
          !*****************************************************************
          real(mcp), dimension(1:3) :: Scat_Phot_AxisX
          real(mcp), dimension(1:3) :: Scat_Phot_AxisY
          real(mcp), dimension(1:3) :: Scat_Phot_AxisZ
          real(mcp), dimension(1:3) :: Scat_Phot_Pe_AxisX
          real(mcp), dimension(1:3) :: Scat_Phot_Pe_AxisY
          real(mcp), dimension(1:3) :: Scat_Phot_Pe_AxisZ
          real(mcp), dimension(1:3) :: Phot_f4_AxisX
          real(mcp), dimension(1:3) :: Phot_f4_AxisY
          real(mcp), dimension(1:3) :: Phot_f4_AxisZ
          real(mcp), dimension(1:3) :: Phot_f4_AxisX_PCF
          real(mcp), dimension(1:3) :: Phot_f4_AxisY_PCF
          real(mcp), dimension(1:3) :: Phot_f4_AxisZ_PCF
          real(mcp), dimension(1:3, 1:3) :: Matrix_Of_Tetrad_Of_Phot_f4_Axis
          real(mcp), dimension(1:3, 1:3) :: Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF
          !*********************************************************************
          real(mcp) :: Cos_Theta_Scat
          real(mcp) :: Sin_Theta_Scat
          real(mcp) :: Phi_Scat
          real(mcp) :: Cos_Phi_Scat
          real(mcp) :: Sin_Phi_Scat 
!******************************************************************************************************* 
          real(mcp) :: Elec_gama
          real(mcp) :: Elec_V
          real(mcp) :: Elec_Phot_mu
          real(mcp) :: Elec_Phot_mu_In_Elec_CF
          real(mcp) :: Elec_Phot_sin
          real(mcp) :: Elec_Phot_sin_In_Elec_CF
          real(mcp) :: Elec_Phot_Phi
          real(mcp) :: w_before
          real(mcp) :: w_after 
          real(mcp), dimension(1:3) :: Elec4P_CF
          real(mcp), dimension(1:3) :: ElecAxisX
          real(mcp), dimension(1:3) :: ElecAxisY
          real(mcp), dimension(1:3) :: ElecAxisZ
          real(mcp), dimension(1:3, 1:3) :: Matrix_Of_Tetrad_Of_ElecAxis
          real(mcp), dimension(1:3, 1:3) :: Matrix_Of_Tetrad1_Of_PhotAxis
!*******************************************************************************************************
          real(mcp) :: mutilde 
          real(mcp) :: sinmutil 
          real(mcp) :: mupsi 
          real(mcp) :: sinpsi 
          real(mcp) :: mutildep 
          real(mcp) :: sintilphip 
          real(mcp) :: epsi 
          logical :: test = .False.
          real(mcp) :: r11
          real(mcp), dimension(1:10) :: test_quantity
          real(mcp) :: random_number_phi 
      contains 
      procedure, public :: Set_Photon_Tetrad_In_CF      => Set_Photon_Tetrad_In_CF_Sub
      procedure, public :: Set_Photon_f3_Tetrad_In_CF   => Set_Photon_f3_Tetrad_In_CF_Sub
      procedure, public :: Set_Elec_Tetrad_In_CF        => Set_Elec_Tetrad_In_CF_Sub
      procedure, public :: Vector_Cross_Product         => Vector_Cross_Product_Sub
      procedure, public :: Get_gama_mu_phi_Of_Scatter_Electron => &
                           Get_gama_mu_phi_Of_Scatter_Electron_Sub
      procedure, public :: Get_gama_mu_phi_Of_Scatter_Electron_Power    =>  &
                           Get_gama_mu_phi_Of_Scatter_Electron_Power_Sub
      procedure, public :: Set_Phot4k_In_Elec_CF        => Set_Phot4k_In_Elec_CF_Sub 


      procedure, public :: Static_Compton_Scattering_With_Polarization => &
                           Static_Compton_Scattering_With_Polarization_Sub
      procedure, public :: Set_f4_In_Elec_CF            => Set_f4_In_Elec_CF_Sub
      procedure, public :: Set_Phot_f4_Tetrad_In_Elec_CF   => &
                           Set_Phot_f4_Tetrad_In_Elec_CF_Sub
      procedure, public :: Compton_Scattering_With_Polarization  =>  &
                           Compton_Scattering_With_Polarization_Sub 
      procedure, public :: Tompson_Scattering_WithOut_Polarization   =>  &
                           Tompson_Scattering_WithOut_Polarization_Sub 
      procedure, public :: Set_Phot_f4_Tetrad_In_Phot_CF   =>    &
                           Set_Phot_f4_Tetrad_In_Phot_CF_Sub 
      procedure, public :: ComptonScaPhoSampling
      !~~~~~~~~~~~~~~~~~~~~~Subroutines to set H3 array~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      procedure, public :: Calc_H3Array  =>  Calc_H3Array_Sub
      procedure, public :: do_dgm  =>  do_dgm_fn
      procedure, public :: Get_gama_mu_phi_Of_Scatter_Electron_HXM  =>  &
                           Get_gama_mu_phi_Of_Scatter_Electron_HXM_Sub
      procedure, public :: Compton_Scat_Sampling_Phi_temp
      procedure, public :: Compton_Scat_Sampling_Phi_Kepler
      procedure, public :: Compton_Scat_Sampling_Phi_StokesVec
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      procedure, public :: Set_Phot_Tetrad1_In_Elec_CF  =>  Set_Phot_Tetrad1_In_Elec_CF_Sub
      procedure, public :: StokesPara_Rotation_Matrix   =>  StokesPara_Rotation_Matrix_Sub
      procedure, public :: Compton_Scattering_With_Zero_QU_StokesVec  =>   &
                           Compton_Scattering_With_Zero_QU_StokesVec_Sub
      procedure, public :: Get_Scattered_Phot4k_CF   =>   Get_Scattered_Phot4k_CF_Sub
      procedure, public :: Get_Scattered_Stokes_Vector_And_f_CF   =>   &
                           Get_Scattered_Stokes_Vector_And_f_CF_Sub
      procedure, public :: Compton_Scattering_With_Polar_StokesVec   =>   &
                           Compton_Scattering_With_Polar_StokesVec_Sub
      end type ScatPhoton_KN

      private :: Set_Photon_Tetrad_In_CF_Sub
      private :: Vector_Cross_Product_Sub
      private :: Get_gama_mu_phi_Of_Scatter_Electron_Sub
      private :: Set_Elec_Tetrad_In_CF_Sub
      private :: Set_Phot4k_In_Elec_CF_Sub
      !private :: Compton_Scattering_With_Zero_QU_Sub
      !private :: Compton_Scattering_WithOut_Polarizations_Sub 
      private :: Set_f4_In_Elec_CF_Sub
      private :: Set_Phot_f4_Tetrad_In_Elec_CF_Sub
      private :: Compton_Scattering_With_Polarization_Sub
      !private :: Diffuse_Reflection_With_Polarization_Sub 
      private :: Tompson_Scattering_WithOut_Polarization_Sub 
      private :: Set_Phot_f4_Tetrad_In_Phot_CF_Sub
      private :: Static_Compton_Scattering_With_Polarization_Sub 
      !private :: Static_Compton_Scattering_WithOut_Polarizations_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      private :: Calc_H3Array_Sub
      private :: do_dgm_fn
      private :: Get_gama_mu_phi_Of_Scatter_Electron_HXM_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private :: Set_Phot_Tetrad1_In_Elec_CF_Sub
      private :: StokesPara_Rotation_Matrix_Sub
      private :: Compton_Scattering_With_Zero_QU_StokesVec_Sub
      private :: Get_Scattered_Phot4k_CF_Sub
      private :: Get_Scattered_Stokes_Vector_And_f_CF_Sub
      private :: Set_Photon_f3_Tetrad_In_CF_Sub
      private :: Compton_Scattering_With_Polar_StokesVec_Sub

      contains
!********************************************************************************
      subroutine Calc_H3Array_Sub( this )
!********************************************************************************
      Implicit None
      class(ScatPhoton_KN) :: this
      real(mcp) :: g1, g2, gc, f1, f2, fc, &
                   Sigma_e, v, e, Ep
      integer :: i, j
!**********************************************************************

      !write(*, *)'mn=', E_dw_H3 * 1.D6 / h_ev
      !stop
      do i = 0, N_H3
          Ep = 10.D0 ** ( dlog_nu_H3 * i + log_nu_dw_H3 ) * h_ev * 1.D-6
          g1 = two
          j = 0
          do
              f1 = this%do_dgm( g1, Ep )
              if( f1 > zero )exit
              g1 = one + one * 0.9D0**j
              j = j + 1
              if(g1 < one + 0.1D-10)stop
              !write(*, *)'H3_Calc01=', j, g1, f1
          enddo
          g2 = one + 1.D-3
          !write(*, *)'H3_Calc1=', g1, f1, Ep / mec2
          do
              f2 = this%do_dgm( g2, Ep )
              if( f2 < zero )exit
              g2 = g2 + 1.D-3
              !write(*, *)'H3_Calc02=', j, g2, f2
          enddo
          !write(*, *)'H3_Calc2=', g2, f2, Ep / mec2
          do
              gc = half * ( g1 + g2 )
              fc = this%do_dgm( gc, Ep )
              !write(*, *)'H3_Calc3=', g1, gc, g2, i
              !write(*, *)'H3_Calc4=', f1, fc, f2, Ep / mec2
              if(fc * f1 > zero)then
                  g1 = gc
                  f1 = fc
              else if(fc * f2 > zero)then
                  g2 = gc
                  f2 = fc
              else if( fc == zero )then
                  exit
              endif
              if( dabs(g1-g2) <= 1.D-13 )exit
          enddo
          v = dsqrt( one - one / gc**2 )
          e = two * Ep / mec2 * gc * ( one + v )
          if( e > 8.19D-5 )then
              Sigma_e = ( ( one - four / e - eight / e**2 ) * dlog(one + e) + half + &
                      eight / e - one / two / (one + e)**2 ) / e
          else
              Sigma_e = one - e
          endif
          this%H3Array( i ) = ( one + v ) * Sigma_e * three * Sigma_T / four
          !write(*, *)'H3_Calc=', i, this%H3Array( i ), v, e, Sigma_e
      enddo
      end subroutine Calc_H3Array_Sub
 
      
!********************************************************************************
      real(mcp) function do_dgm_fn( this, gm, Ep )
!********************************************************************************
      Implicit None
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: gm, Ep
      real(mcp) :: Sigma_e, dSigmaeff_dgm, e, v, dhde
!**********************************************************************

      v = dsqrt( one - one / gm**2 )
      e = two * Ep / mec2 * gm * ( one + v )
      if( e > 8.19D-5 )then
          Sigma_e = ( ( one - four / e - eight / e**2 )*dlog(one + e) + half + &
                  eight / e - one / two / (one + e)**2 ) / e
          dhde = ( four / e**2 + 16.D0 / e**3 ) * dlog(one + e) + &
                 ( one - four / e - eight / e**2 ) / ( one + e ) - &
                 eight / e**2 + one / (one + e)**3
      else
          Sigma_e = one - e 
          !dhde = four / three - eight / three * e + 26.D0 / 5.D0 * e**2 - &
          !      ( 28.D0 / 15.D0 - 7.D0 ) * e**3
          dhde = ( four / e**2 + 16.D0 / e**3 ) * dlog(one + e) + &
                 ( one - four / e - eight / e**2 ) / ( one + e ) - &
                 eight / e**2 + one / (one + e)**3
      endif
      do_dgm_fn = - v * Sigma_e / gm + dhde / gm / v
      !write(*, *)'HHHH=', e, v, Sigma_e, dhde
      end function do_dgm_fn

  

!*******************************************************************************************************
      subroutine Get_gama_mu_phi_Of_Scatter_Electron_HXM_Sub( this, T_e )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: T_e 
      real(mcp) :: p11, p22, p33, p44, S3, K1, K2, K3, Y1
      real(mcp) :: AcceptProbability, epsi, SigmaKN, t, r1, r2
      real(mcp) :: Ep, H2, a2, b2, H1, a1, b1, prob1
      integer(kind=8) :: i, j, N_H3
      !**********************************************************************

      Ep = dabs( this%Phot4k_CtrCF(1) )
      
      write(*, *)'HXM1 = ', Ep, E_critical
      write(*, *)'f1= ', this%H3Array(1:100) 
      if( this%H3Array(100) == zero )then
          call this%Calc_H3Array()
      endif
      if( Ep >= E_critical )then
          t = T_e / mec2

          b1 = dsqrt( four * t**2 + one / three )
          a1 = two * t + two * b1 * dcos( dacos( t * ( 16.D0 * t**2 - one ) / two / b1**3 ) / three )
          H1 = a1 * dsqrt( a1**2 - one ) * dexp( - ( a1 - one ) / three / t )

          b2 = dsqrt( t**2 + one / three )
          a2 = t + two * b2 * dcos( dacos( t * ( four * t**2 - one ) / four / b2**3 ) / three )
          H2 = a2 * dsqrt( a2**2 - one ) * dexp( - two * ( a2 - one ) / three / t )
          
          i = 0
          j = 0
          do
              i = i + 1
              do
                  r1 = ranmar()
                  r2 = ranmar()
                  j = j + 1
                  if( t <= 0.01D0 )then
                      if( r2**2 <= - ln_e * r1 * dlog( r1 ) )then
                          this%Elec_V = DSQRT( - three * t * dlog( r1 ) )
                          this%Elec_gama = one / DSQRT( one - this%Elec_V ** 2 )
                          exit
                      endif
                  else if( t <= 0.25D0 )then
                      this%Elec_gama = one - 1.5D0 * t * dlog( r1 )
                      if( r2 * H1 <= this%Elec_gama * dsqrt( r1 * ( this%Elec_gama ** 2 - one ) ) )then
                          this%Elec_V = DSQRT( one - one / this%Elec_gama**2 )
                          exit
                      endif
                  else
                      this%Elec_gama = one - three * t * dlog( r1 )
                      this%Elec_V = DSQRT( one - one / this%Elec_gama**2 )
                      !write(*, *)'HXM2 = ', j, t, r2*H2, ( this%Elec_gama * r1 )**2 * this%Elec_V
                      if( r2 * H2 <= ( this%Elec_gama * r1 )**2 * this%Elec_V )exit
                  endif
                  !if( j > 1000 )then
                      !write(*, *)'HXM3 = ', t, Ep, epsi, prob1
                  !endif
              enddo
              !do
              this%Elec_Phot_mu = - one + two * ranmar()
                !If ( dabs( this%Elec_Phot_mu ) /= one ) exit
              !Enddo 
              epsi = ( two * Ep / mec2 ) * this%Elec_gama * ( one - this%Elec_V * this%Elec_Phot_mu )
              if ( epsi > 8.19D-4 ) then
                  SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                        + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
              else 
                  SigmaKN = one - epsi
              end if

              N_H3 = floor( ( dlog10( Ep * 1.D6 / h_ev ) - log_nu_dw_H3 ) / dlog_nu_H3 ) 
              prob1 = ( one - this%Elec_V * this%Elec_Phot_mu ) * SigmaKN * Sigma_T / this%H3Array(N_H3)
              !if( ranmar() * this%H3Array(N_H3) <= ( one - this%Elec_V * &
              !this%Elec_Phot_mu ) * SigmaKN )exit
              !write(*, *)'f1= ', prob1, this%H3Array(N_H3), N_H3
              !stop 
              if( ranmar() <= prob1 )exit
              if( i > 1000 )then
                  write(*, *)'HXM = ', t, Ep, epsi, prob1
              endif
          enddo
          this%Elec_Phot_Phi = twopi * ranmar()
          If ( dabs(this%Elec_Phot_mu) < one ) then
              this%Elec_Phot_sin = DSQRT( one - this%Elec_Phot_mu**2 )
          else
              this%Elec_Phot_sin = zero
          endif
          this%Elec4P_CF(1) = this%Elec_Phot_sin * DCOS( this%Elec_Phot_Phi )
          this%Elec4P_CF(2) = this%Elec_Phot_sin * DSIN( this%Elec_Phot_Phi )
          this%Elec4P_CF(3) = this%Elec_Phot_mu

      else 

          call this%Get_gama_mu_phi_Of_Scatter_Electron( T_e ) 

      endif

      end subroutine Get_gama_mu_phi_Of_Scatter_Electron_HXM_Sub
  


!*******************************************************************************************************
      subroutine Set_Photon_Tetrad_In_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), dimension(1:3) :: e_p
      real(mcp) :: Vector_Length
!**********************************************************************
 
      !this%PhotAxisZ(1) = this%Phot4k_CtrCF(2) / dabs( this%Phot4k_CtrCF(1) )
      !this%PhotAxisZ(2) = this%Phot4k_CtrCF(3) / dabs( this%Phot4k_CtrCF(1) )
      !this%PhotAxisZ(3) = this%Phot4k_CtrCF(4) / dabs( this%Phot4k_CtrCF(1) )

      this%PhotAxisZ = this%Phot4k_CtrCF(2: 4) / Vector3D_Length( this%Phot4k_CtrCF(2: 4) )
 
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one
      call this%Vector_Cross_Product( this%PhotAxisZ, e_p, this%PhotAxisX )
      Vector_Length = Vector3D_Length( this%PhotAxisX )
      If ( Vector_Length /= zero ) then
          this%PhotAxisX = this%PhotAxisX / Vector_Length
          call this%Vector_Cross_Product( this%PhotAxisZ, &
                           this%PhotAxisX, this%PhotAxisY )
      Else
          this%PhotAxisX(1) = one
          this%PhotAxisX(2) = zero
          this%PhotAxisX(3) = zero
          this%PhotAxisY(1) = zero
          this%PhotAxisY(2) = DSIGN( one, this%PhotAxisZ(3) ) * one
          this%PhotAxisY(3) = zero
      Endif
    
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 1) = this%PhotAxisX(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 2) = this%PhotAxisX(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 3) = this%PhotAxisX(3)

      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 1) = this%PhotAxisY(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 2) = this%PhotAxisY(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 3) = this%PhotAxisY(3)

      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 1) = this%PhotAxisZ(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 2) = this%PhotAxisZ(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 3) = this%PhotAxisZ(3)
      !write(*, *)'kk=' ,this%Phot4k_CtrCF, this%PhotAxisZ, Vector3D_Length( this%Phot4k_CtrCF(2: 4) ), &
      !          this%Matrix_Of_Tetrad_Of_PhotAxis
      end subroutine Set_Photon_Tetrad_In_CF_Sub



!*******************************************************************************************************
      subroutine Set_Photon_f3_Tetrad_In_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), dimension(1:3) :: e_p
      real(mcp) :: Vector_Length
!********************************************************************** 

      this%PhotAxisZ = this%Phot4k_CtrCF(2: 4) / Vector3D_Length( this%Phot4k_CtrCF(2: 4) )
      this%PhotAxisX = this%f4_CF(2: 4)

      !write(*, fmt="(' ', A5, 3ES18.7)")'s0=', Vector3D_Length( this%PhotAxisZ ), &
      !    Vector3D_Inner_Product( this%PhotAxisX, this%PhotAxisZ ), Vector3D_Length( this%PhotAxisX )

      call this%Vector_Cross_Product( this%PhotAxisZ, this%PhotAxisX, this%PhotAxisY )
      Vector_Length = Vector3D_Length( this%PhotAxisX ) 
    
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 1) = this%PhotAxisX(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 2) = this%PhotAxisX(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(1, 3) = this%PhotAxisX(3)

      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 1) = this%PhotAxisY(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 2) = this%PhotAxisY(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(2, 3) = this%PhotAxisY(3)

      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 1) = this%PhotAxisZ(1)
      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 2) = this%PhotAxisZ(2)
      this%Matrix_Of_Tetrad_Of_PhotAxis(3, 3) = this%PhotAxisZ(3) 

      end subroutine Set_Photon_f3_Tetrad_In_CF_Sub


!********************************************************************************
      subroutine Vector_Cross_Product_Sub(this, e1, e2, e3)
!********************************************************************************
      Implicit None
      class(ScatPhoton_KN) :: this
      real(mcp),dimension(1:3), intent(in) :: e1, e2
      real(mcp),dimension(1:3), intent(out) :: e3
!**********************************************************************

      e3(1) = e1(2) * e2(3) - e1(3) * e2(2)    
      e3(2) = e1(3) * e2(1) - e1(1) * e2(3)    
      e3(3) = e1(1) * e2(2) - e1(2) * e2(1)  
      end subroutine Vector_Cross_Product_Sub
 
!*******************************************************************************************************
      real(mcp) function ComptonScaPhoSampling(this, epsi)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: epsi
      real(mcp) :: r, r1 
      !******************************************************************** 

      do
          if( ranmar() <= 27.D0 / ( 2.D0*epsi + 29.D0 ) )then
              r1 = ranmar()
              r = ( epsi + one ) / ( epsi*r1 + one )
              !if(epsi <= 1.D-10)then
              !    r = one + epsi * ( one - r1 )
              !    this%r11 = r1
              !endif
              if( ranmar() <= ( ( (epsi + two - two*r)/epsi )**2 + one) / two )exit 
          else
              r = epsi*ranmar() + one
              if( ranmar() <= 6.75D0*(r - one)**2 / r**3 ) exit
          endif
      end do 
      !if(epsi<1.D-15)r=one
      ComptonScaPhoSampling = r
      return
      end function ComptonScaPhoSampling
 
!*******************************************************************************************************
      real(mcp) function Compton_Scat_Sampling_Phi1( B )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: B
      real(mcp) :: r, A, m, Phi1
      integer :: m1

      A = one / twopi
      Phi1 = ranmar() * twopi
      If ( Phi1 <= pi / two ) then
          m1 = 0!zero
      Else If ( Phi1 <= pi ) then
          m1 = 1!one
      Else If ( Phi1 <= pi * three / two ) then
          m1 = 2!two
      Else If ( Phi1 <= twopi ) then
          m1 = 3!three
      Endif
      If(m1/=floor(Phi1/pi*two))then
          write(*,*)'tttt==', m1, floor(Phi1/pi*two),m1-floor(Phi1/pi*two)
          stop
      endif
      If ( ranmar() * two * A <= A - B * DCOS( two*Phi1 ) ) then
          Compton_Scat_Sampling_Phi1 = Phi1
      Else
          Compton_Scat_Sampling_Phi1 = (two*m + one) * pi / two - Phi1
      Endif
      return
      end function Compton_Scat_Sampling_Phi1

!*******************************************************************************************************
      real(mcp) function Compton_Scat_Sampling_Phi( B )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: B
      real(mcp) :: r, A, m, Phi1 

      A = one / twopi
      !this%random_number_phi = ranmar()
      Phi1 = ranmar() * twopi 
      If ( ranmar() * two * A <= A - B * DCOS( two*Phi1 ) ) then
          Compton_Scat_Sampling_Phi = Phi1
      Else
          m = floor( Phi1 / ( pi/two ) ) 
          Compton_Scat_Sampling_Phi = (two * m + one) * pi / two - Phi1
          !write(*, *)'kk=' , m, A - B * DCOS( two*Phi1 ), A * two
      Endif
      return
      end function Compton_Scat_Sampling_Phi

!*******************************************************************************************************
      real(mcp) function Compton_Scat_Sampling_Phi_temp( this, B )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: B
      real(mcp) :: r, A, m, Phi1 

      A = one / twopi
      Phi1 = ranmar() * twopi 
      If ( ranmar() * two * A <= A - B * DCOS( two*Phi1 ) ) then
          Compton_Scat_Sampling_Phi_temp = Phi1
      Else
          m = floor( Phi1 / ( pi/two ) ) 
          Compton_Scat_Sampling_Phi_temp = (two * m + one) * pi / two - Phi1
          !write(*, *)'kk=' , m, A - B * DCOS( two*Phi1 ), A * two
      Endif
      return
      end function Compton_Scat_Sampling_Phi_temp


!*******************************************************************************************************
      real(mcp) function Compton_Scat_Sampling_Phi_Kepler( this, B )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: B
      real(mcp) :: r, A, m, Phi1, e, q, x0, x, x1
      integer :: i

      e = B * twopi
      q = twopi * ranmar() * two
      x0 = twopi
      i = 0
      do
          x = e * dsin( x0 ) + q
          ! write(*, *)'kks=', x, x0, dabs( x - x0 ), dabs( x - x0 ) <= 1.D-10, i
          if( dabs( x - x0 ) <= 1.D-13 )then
              !write(*, *)'kl11=', i, x, x0, dabs( x - x0 ), dabs( x - x0 ) <= 1.D-10
              exit
          endif
          x1 = x0
          x0 = x
          i = i + 1
          if(i >= 1000)then
              !write(*, *)'kl22=', i, x, x1, dabs( x - x1 ), dabs( x - x1 ) <= 1.D-13
              exit
          endif
      enddo
      Compton_Scat_Sampling_Phi_Kepler = x / two
      return
      end function Compton_Scat_Sampling_Phi_Kepler


!*******************************************************************************************************
      real(mcp) function Compton_Scat_Sampling_Phi_StokesVec( this, B, C )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: B, C
      real(mcp) :: r, A, m, Phi1, temp_phi, F_phi, H_phi
      integer :: m2

      !A = one / twopi
      !Phi1 = ranmar() * twopi 
      !If ( ranmar() * two * A <= A + B * DCOS( two*Phi1 ) ) then
      !    Compton_Scat_Sampling_Phi_StokesVec = Phi1
      !Else
      !    m = floor( Phi1 / ( pi/two ) ) 
      !    Compton_Scat_Sampling_Phi_StokesVec = (two * m + one) * pi / two - Phi1
      !    !write(*, *)'kk=' , m, A - B * DCOS( two*Phi1 ), A * two
      !Endif 
      If(B == zero .and. C == zero)then
          Compton_Scat_Sampling_Phi_StokesVec = twopi * ranmar()
          return
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            A = one / twopi
            Phi1 = ranmar() * twopi
            If ( ranmar() * two * A <= A + B * DCOS( two*Phi1 ) ) then
                temp_phi = Phi1
            Else
                m2 = floor( Phi1 / ( pi/two ) ) 
                temp_phi = (two * m2 + one) * pi / two - Phi1
                !write(*, *)'kk=' , m, A - B * DCOS( two*Phi1 ), A * two
            Endif
            F_phi = A + B * dcos(two * temp_phi) 
            H_phi = C * dsin(two * temp_phi)
            If( ranmar() - one <= H_phi / F_phi )then
                Compton_Scat_Sampling_Phi_StokesVec = temp_phi
            else
                Compton_Scat_Sampling_Phi_StokesVec = twopi - temp_phi
            endif  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      return
      end function Compton_Scat_Sampling_Phi_StokesVec


!*******************************************************************************************************
      real(mcp) function chisquare(n)
!*******************************************************************************************************
      implicit none
      integer, intent(in) :: n
      integer :: i
     !***************************************************** 

      chisquare = zero
      i = 0
      do 
          chisquare = chisquare + GAUSSIAN1()**2
          i = i + 1
          if ( i == n) exit
      end do
      chisquare = dsqrt( chisquare / two )
      return
      end function chisquare

!*******************************************************************************************************
      subroutine Get_gama_mu_phi_Of_Scatter_Electron_Sub( this, T_e )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), intent(in) :: T_e 
      real(mcp) :: p11, p22, p33, p44, S3, K1, K2, K3, Y1
      real(mcp) :: AcceptProbability, epsi, SigmaKN, t, r1, r2
      integer(kind=8) :: i
      !**********************************************************************

      t = T_e / mec2
      p11 = sqrtpi * 0.25D0
      p22 = dsqrt( 0.5D0 * t )*0.5D0
      p33 = 0.375D0 * sqrtpi * t
      p44 = t * dsqrt( 0.5D0 * t )
      S3 = p11 + p22 + p33 + p44
      i = 0
 
      gamamuphi : do
        i = i + 1
        do 
          K1 = ranmar() * S3 - p11
          K2 = K1 - p22
          K3 = K2 - p33
          if (K1 <= zero) then
              Y1 = chisquare(3)
          else if (K2 <= zero) then  
              Y1 = chisquare(4)  
          else if (K3 <= zero) then
              Y1 = chisquare(5)
          else
              Y1 = chisquare(6)
          end if 
          AcceptProbability = DSQRT( one + half * t * Y1**2 ) /&
                              ( one + DSQRT( half * t ) * Y1 )
          if ( ranmar() <= AcceptProbability ) exit
        end do      

        this%Elec_gama = one + t * Y1**2
        this%Elec_V = DSQRT( one - one / this%Elec_gama**2 )

        !Do
            this%Elec_Phot_mu = - one + two * ranmar()
            !If ( dabs( this%Elec_Phot_mu ) /= one ) exit 
        !Enddo
        if ( ranmar() >= half * ( one - this%Elec_V * this%Elec_Phot_mu ) ) then  
          this%Elec_Phot_mu = - this%Elec_Phot_mu
        end if

        epsi = ( two * dabs( this%Phot4k_CtrCF(1) ) / mec2 ) * this%Elec_gama * &
               ( one - this%Elec_V * this%Elec_Phot_mu )
        if ( epsi > 8.19D-4 ) then
          SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                    + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
        else 
          SigmaKN = one - epsi
        end if
 
        r1 = ranmar() 
        if( isnan(epsi) )then
            write(*, *)'Get_gama_mu_phi==epsi', epsi, dabs( this%Phot4k_CtrCF(1) ), &
                  this%Elec_V, this%Elec_Phot_mu
            stop
        endif
        if( i > 200000 )write(*, *)'Get_gama_mu_phi_Times ==', r1, SigmaKn, AcceptProbability, epsi, &
                       this%Phot4k_CtrCF(1)
        if ( r1 <= SigmaKN ) exit
      end do gamamuphi
      this%Elec_Phot_Phi = twopi * ranmar()
      If ( dabs(this%Elec_Phot_mu) < one ) then
          this%Elec_Phot_sin = DSQRT( one - this%Elec_Phot_mu**2 )
      else
          this%Elec_Phot_sin = zero
      endif

      this%Elec4P_CF(1) = this%Elec_Phot_sin * DCOS( this%Elec_Phot_Phi )
      this%Elec4P_CF(2) = this%Elec_Phot_sin * DSIN( this%Elec_Phot_Phi )
      this%Elec4P_CF(3) = this%Elec_Phot_mu
      end subroutine Get_gama_mu_phi_Of_Scatter_Electron_Sub



!*******************************************************************************************************
      subroutine Get_gama_mu_phi_Of_Scatter_Electron_Power_Sub( this )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp) :: K1, K2, K3, Y1
      real(mcp) :: AcceptProbability, epsi, SigmaKN, r1, r2
      integer(kind=8) :: i
      !**********************************************************************
 
      i = 0
 
      gamamuphi : do
        i = i + 1    

        r1 = ranmar() 
        this%Elec_gama = ( this%gama1**(one - this%alp) - r1*&
                           (this%alp-one)/this%N_coef )**( one / (one - this%alp) )
        !write(*, *)'sdf=', this%Elec_gama, r1, this%gama1**(one - this%alp), r1*&
        !                   (this%alp-one)/this%N_coef
        this%Elec_V = DSQRT( one - one / this%Elec_gama**2 )

        !Do
            this%Elec_Phot_mu = - one + two * ranmar()
            !If ( dabs( this%Elec_Phot_mu ) /= one ) exit
        !Enddo
        if ( ranmar() >= half * ( one - this%Elec_V * this%Elec_Phot_mu ) ) then  
          this%Elec_Phot_mu = - this%Elec_Phot_mu
        end if

        epsi = ( two * dabs( this%Phot4k_CtrCF(1) ) / mec2 ) * this%Elec_gama * &
               ( one - this%Elec_V * this%Elec_Phot_mu )
        if ( epsi > 8.19D-4 ) then
          SigmaKN = ( ( one - four / epsi - eight / epsi / epsi ) * DLOG( one + epsi )&
                    + half + eight / epsi - half / ( one + epsi )**2 ) / epsi * (three/four)
        else 
          SigmaKN = one - epsi
        end if
 
        r1 = ranmar()   
        if ( r1 <= SigmaKN ) exit
      end do gamamuphi
      this%Elec_Phot_Phi = twopi * ranmar()
      If ( dabs(this%Elec_Phot_mu) < one ) then
          this%Elec_Phot_sin = DSQRT( one - this%Elec_Phot_mu**2 )
      else
          this%Elec_Phot_sin = zero
      endif

      this%Elec4P_CF(1) = this%Elec_Phot_sin * DCOS( this%Elec_Phot_Phi )
      this%Elec4P_CF(2) = this%Elec_Phot_sin * DSIN( this%Elec_Phot_Phi )
      this%Elec4P_CF(3) = this%Elec_Phot_mu
      end subroutine Get_gama_mu_phi_Of_Scatter_Electron_Power_Sub


!*******************************************************************************************************
      subroutine Set_Elec_Tetrad_In_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), dimension(1:3) :: e_p
      real(mcp) :: Vector_Length
      !**********************************************************************

      this%ElecAxisZ(1) = this%Elec4P_CF(1)
      this%ElecAxisZ(2) = this%Elec4P_CF(2)
      this%ElecAxisZ(3) = this%Elec4P_CF(3)
      
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one
 
      call this%Vector_Cross_Product( e_p, this%ElecAxisZ, this%ElecAxisY ) 
      Vector_Length = Vector3D_Length( this%ElecAxisY )
      If ( Vector_Length /= zero ) then
          this%ElecAxisY = this%ElecAxisY / Vector_Length
          call this%Vector_Cross_Product( this%ElecAxisY, &
                           this%ElecAxisZ, this%ElecAxisX )
      Else
          stop
          this%ElecAxisX(1) = one
          this%ElecAxisX(2) = zero
          this%ElecAxisX(3) = zero
          this%ElecAxisY(1) = zero
          this%ElecAxisY(2) = one !DSIGN( one, this%ElecAxisZ(3) ) * one
          this%ElecAxisY(3) = zero
      Endif
      
      this%Matrix_Of_Tetrad_Of_ElecAxis(1, 1) = this%ElecAxisX(1)
      this%Matrix_Of_Tetrad_Of_ElecAxis(1, 2) = this%ElecAxisX(2)
      this%Matrix_Of_Tetrad_Of_ElecAxis(1, 3) = this%ElecAxisX(3)

      this%Matrix_Of_Tetrad_Of_ElecAxis(2, 1) = this%ElecAxisY(1)
      this%Matrix_Of_Tetrad_Of_ElecAxis(2, 2) = this%ElecAxisY(2)
      this%Matrix_Of_Tetrad_Of_ElecAxis(2, 3) = this%ElecAxisY(3)

      this%Matrix_Of_Tetrad_Of_ElecAxis(3, 1) = this%ElecAxisZ(1)
      this%Matrix_Of_Tetrad_Of_ElecAxis(3, 2) = this%ElecAxisZ(2)
      this%Matrix_Of_Tetrad_Of_ElecAxis(3, 3) = this%ElecAxisZ(3)
      end subroutine Set_Elec_Tetrad_In_CF_Sub

!**************************************************************************** 
      subroutine Set_Phot4k_In_Elec_CF_Sub(this)
!**************************************************************************** 
      implicit none
      class(ScatPhoton_KN) :: this 
      !**********************************************************************

      this%Elec_Phot_mu_In_Elec_CF = ( this%Elec_Phot_mu - this%Elec_V ) / &
                                     ( one - this%Elec_V * this%Elec_Phot_mu )
      this%Elec_Phot_sin_In_Elec_CF = DSQRT( one - this%Elec_Phot_mu_In_Elec_CF**2 )

      this%Phot4k_In_Elec_CF(1) = this%Elec_gama * dabs( this%Phot4k_CtrCF(1) ) * &
                                  ( one - this%Elec_V * this%Elec_Phot_mu )
      !write(*, *)'Elec_cf=', this%Phot4k_In_Elec_CF(1), this%Elec_gama, this%Phot4k_CtrCF(1), &
      !                            ( one - this%Elec_V * this%Elec_Phot_mu )
      this%Phot4k_In_Elec_CF(2) = - dabs( this%Phot4k_CtrCF(1) ) * this%Elec_Phot_sin
      this%Phot4k_In_Elec_CF(3) = zero !- dabs( this%Phot4k_CtrCF(1) ) * this%Elec_Phot_sin
      this%Phot4k_In_Elec_CF(4) = this%Elec_gama * dabs( this%Phot4k_CtrCF(1) ) * &
                                  ( this%Elec_Phot_mu - this%Elec_V ) 
      !this%Scal_Factor_before = this%Phot4k_In_Elec_CF(1) / this%Phot4k_CtrCF(1)
      
      end subroutine Set_Phot4k_In_Elec_CF_Sub



!**************************************************************************** 
      subroutine Set_Phot_Tetrad1_In_Elec_CF_Sub(this)
!**************************************************************************** 
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp), dimension(1:3) :: e_p
      real(mcp) :: Vector_Length
      !**********************************************************************

      this%Phot_Tetrad1_ECF_axisZ(1:3) = this%Phot4k_In_Elec_CF(2: 4) &
                                 / dabs( this%Phot4k_In_Elec_CF(1) )
  
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one
 
      call this%Vector_Cross_Product( this%Phot_Tetrad1_ECF_axisZ, e_p, this%Phot_Tetrad1_ECF_axisY ) 
      Vector_Length = Vector3D_Length( this%Phot_Tetrad1_ECF_axisY )
      If ( Vector_Length /= zero ) then
          this%Phot_Tetrad1_ECF_axisY = this%Phot_Tetrad1_ECF_axisY / Vector_Length
          call this%Vector_Cross_Product( this%Phot_Tetrad1_ECF_axisY, &
                           this%Phot_Tetrad1_ECF_axisZ, this%Phot_Tetrad1_ECF_axisX )
      Else
          stop
          this%Phot_Tetrad1_ECF_axisX(1) = one
          this%Phot_Tetrad1_ECF_axisX(2) = zero
          this%Phot_Tetrad1_ECF_axisX(3) = zero
          this%Phot_Tetrad1_ECF_axisY(1) = zero
          this%Phot_Tetrad1_ECF_axisY(2) = one  
          this%Phot_Tetrad1_ECF_axisY(3) = zero
      Endif

      this%Matrix_Of_Tetrad1_Of_photAxis(1, 1:3) = this%Phot_Tetrad1_ECF_axisX(1:3)
      this%Matrix_Of_Tetrad1_Of_photAxis(2, 1:3) = this%Phot_Tetrad1_ECF_axisY(1:3)
      this%Matrix_Of_Tetrad1_Of_photAxis(3, 1:3) = this%Phot_Tetrad1_ECF_axisZ(1:3)
      this%Matrix_ECF1_2_ECF = this%Matrix_Of_Tetrad1_Of_photAxis
     
      CALL Inverse_Matrix_Of_Matrix_3X3_Sub(this%Matrix_ECF1_2_ECF, this%Matrix_ECF_2_ECF1)
   
      end subroutine Set_Phot_Tetrad1_In_Elec_CF_Sub

 
!**************************************************************************** 
      subroutine StokesPara_Rotation_Matrix_Sub(this, angle_phi, Vec_Stokes)
!**************************************************************************** 
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp), dimension(1:3), intent(inout) :: Vec_Stokes
      !real(mcp), dimension(1:3), intent(out) :: Rot_Vec_Stokes
      real(mcp), intent(in) :: angle_phi
      real(mcp), dimension(2: 3) :: Rot_Vec_Stokes
      !**********************************************************************

      Rot_Vec_Stokes(2) =  dcos(two*angle_phi) * Vec_Stokes(2) + &
                           dsin(two*angle_phi) * Vec_Stokes(3)
      Rot_Vec_Stokes(3) = -dsin(two*angle_phi) * Vec_Stokes(2) + &
                           dcos(two*angle_phi) * Vec_Stokes(3)
      Vec_Stokes(2) = Rot_Vec_Stokes(2)
      Vec_Stokes(3) = Rot_Vec_Stokes(3)
   
      end subroutine StokesPara_Rotation_Matrix_Sub



!*************************************************************************** 
      subroutine Set_f4_In_Elec_CF_Sub(this)
!***************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this
      real(mcp), dimension(1:3, 1:3) :: InvMatrix_Elec, InvMatrix_Phot
      real(mcp), dimension(1:4) :: f4_In_Elec_MF   !MF means in Electron Moving Frame
      real(mcp), dimension(1:3, 1:3) :: Matrix_CF2Elec_MF
      !*********************************************************************************

      CALL Inverse_Matrix_Of_Matrix_3X3_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, InvMatrix_Elec )
      CALL Inverse_Matrix_Of_Matrix_3X3_Sub( this%Matrix_Of_Tetrad_Of_PhotAxis, InvMatrix_Phot )
      !******************************************************************************************* 
      CALL Matrix_Multiplication33X33_Sub( InvMatrix_Phot, InvMatrix_Elec, Matrix_CF2Elec_MF )
      CALL Matrix_Multiplication13X33_Sub( this%f4_CF(2:4), Matrix_CF2Elec_MF, f4_In_Elec_MF(2:4) )
      f4_In_Elec_MF(1) = this%f4_CF(1)
      !*******************************************************************************************
      this%f4_In_Elec_CF(1) = this%Elec_gama * ( f4_In_Elec_MF(1) - &
                                   this%Elec_V * f4_In_Elec_MF(4) )
      this%f4_In_Elec_CF(2) = f4_In_Elec_MF(2)
      this%f4_In_Elec_CF(3) = f4_In_Elec_MF(3)
      this%f4_In_Elec_CF(4) = this%Elec_gama * ( f4_In_Elec_MF(4) - &
                                   this%Elec_V * f4_In_Elec_MF(1) )
      
      end subroutine Set_f4_In_Elec_CF_Sub
 

!*******************************************************************************************************
      subroutine Compton_Scattering_With_Zero_QU_StokesVec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp) :: epsi, r, phip, mupsi, mu_tilde_p, sinmu_tilde, sinpsi, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp
      real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      real(mcp), dimension(1:3) :: Temp_Vector, f3_scat_e_tilde, P_e, k_p
      real(mcp), dimension(1:4) :: f4_scat_e
      real(mcp) :: psi_p, beta, temp_cos, temp_sin
      real(mcp) :: mu_kp_pe_ECF, mu_kp_pe_CF, sin_kp_pe_CF, cos_phi, sin_phi, &
                   mu_ya_yb, mu_signs, Rot_angle_k_pe
      !****************************************************************************
      real(mcp) :: a, b, c, mutilde, Etilde
      Complex*16 :: roots3(1: 3)
      integer :: del, cases_f
      !****************************************************************************

      mutilde = this%Elec_Phot_mu_In_Elec_CF
      Etilde  = this%Phot4k_In_Elec_CF(1) 
      epsi = two * dabs( Etilde ) / mec2
      !write(*, *)'mmmk=', epsi, Etilde
      if( epsi >= 1.D-6 )then
          r = this%ComptonScaPhoSampling(epsi)
 
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r  
          mupsi = one - ( r - one ) * two / epsi 
      else 
          a = zero
          b = three
          c = four - ranmar() * eight
          call root3(a, b, c, roots3(1), roots3(2), roots3(3), del)
          if( del == 1 )then
              mupsi = real( roots3(1) )
          else
              stop
          endif
          r = one + epsi * ( one - mupsi ) / two
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r 
      endif
      this%r_one_hvlmec2_one_cosE = r

  
      phip = twopi * ranmar()
   
      !**********************************************************************
      if ( dabs(mutilde) < one ) then
          sinmu_tilde = dsqrt( one - mutilde**2 )
      else
          sinmu_tilde = zero
      endif
      if ( dabs(mupsi) < one ) then
          sinpsi = dsqrt( one - mupsi**2 )
      else
          sinpsi = zero
      endif 
   
      !********************************************************************************** 
      if(this%delta_pd /= zero)then
          write(*, *)'mmmkss=',  this%delta_pd
          stop
      else
          N_temp = r + one / r - sinpsi**2
          this%Q_sp_scat = sinpsi**2 / N_temp
          this%U_sp_scat = zero 
          this%delta_pd_scat = DSQRT( this%Q_sp_scat**2 + this%U_sp_scat**2 )
          this%Vector_Stokes4_ECF_scat(1) = this%Vector_Stokes4_CF(1)
          this%Vector_Stokes4_ECF_scat(2) = sinpsi**2 / N_temp
          this%Vector_Stokes4_ECF_scat(3) = zero
          this%Vector_Stokes4_ECF_scat(4) = zero
      endif
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      k_p(1) = sinpsi * dcos( phip )
      k_p(2) = sinpsi * dsin( phip )
      k_p(3) = mupsi
   
      call this%Get_Scattered_Phot4k_CF(k_p)
      !write(*, *)'ffs2=',this%Scattered_Phot4k_CF
      !write(*, *)'ffs=', dacos( - 0.5D0 ), dacos( - 1.D0 ), dacos( - 0.D0 )
      ! stop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      call this%Get_Scattered_Stokes_Vector_And_f_CF( k_p )
      !write(*, *)'ffs2=',this%f4_scat_CF(2:4)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      end subroutine Compton_Scattering_With_Zero_QU_StokesVec_Sub



!*******************************************************************************************************
      subroutine Get_Scattered_Phot4k_CF_Sub(this, k_p)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp), intent(in) :: k_p(1: 3)
      real(mcp) :: Temp_Matrix_1X3(1: 3), temp_sin, temp_cos, cos_phi, sin_phi, mu_kp_pe_CF, &
          sin_kp_pe_CF, Scattered_Phot4k_In_Elec(1: 4), Temp_Matrix_3X3(1:3, 1:3), Temp_vec(1: 3)
      !********************************************************************** 

      this%Scattered_Phot4k_In_Elec_CF(2: 4) = k_p(1: 3) * this%Scattered_Phot4k_In_Elec_CF(1)

      CALL Matrix_Multiplication13X33_Sub( k_p, this%Matrix_Of_Tetrad1_Of_photAxis, Temp_Matrix_1X3 ) 
      temp_sin = dsqrt( one - Temp_Matrix_1X3(3)**2 )
      temp_cos = Temp_Matrix_1X3(3)
       
      cos_phi = Temp_Matrix_1X3(1) / temp_sin
      sin_phi = Temp_Matrix_1X3(2) / temp_sin
  
      mu_kp_pe_CF = ( Temp_Matrix_1X3(3) + this%Elec_V ) / ( one + Temp_Matrix_1X3(3) * this%Elec_V )
      sin_kp_pe_CF = dsqrt( one - mu_kp_pe_CF**2 )

      Scattered_Phot4k_In_Elec(1) = this%Elec_gama * this%Scattered_Phot4k_In_Elec_CF(1) * ( one + &
                                         this%Elec_V * Temp_Matrix_1X3(3) )
      Temp_vec(1) = sin_kp_pe_CF * cos_phi
      Temp_vec(2) = sin_kp_pe_CF * sin_phi
      Temp_vec(3) = mu_kp_pe_CF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
                         this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_3X3 )
      CALL Matrix_Multiplication13X33_Sub( Temp_vec, Temp_Matrix_3X3, Temp_Matrix_1X3 )
      
      !write(*, fmt="(' ', A5, 2ES18.7)")'s5=', Vector3D_Length( Temp_vec ), Vector3D_Length( Temp_Matrix_1X3 )
      this%Scat_Phot3k_CF = Temp_Matrix_1X3
      this%Scattered_Phot4k_CF(1) = Scattered_Phot4k_In_Elec(1)
      this%Scattered_Phot4k_CF(2:4) = Temp_Matrix_1X3 * Scattered_Phot4k_In_Elec(1)
  
      end subroutine Get_Scattered_Phot4k_CF_Sub
 

!*******************************************************************************************************
      subroutine Get_Scattered_Stokes_Vector_And_f_CF_Sub(this, k_p )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp), intent(in) :: k_p(1: 3)
      real(mcp) :: Temp_Matrix_1X3(1: 3), temp_sin, temp_cos, cos_phi, sin_phi, mu_kp_pe_CF, &
               sin_kp_pe_CF, Scattered_Phot4k_In_Elec(1: 4), Temp_Matrix_3X3(1:3, 1:3), P_e(1: 3), &
               Temp_Vector(1: 3), mu_ya_yb, mu_signs, Rot_angle_k_Pe, k_p2(1: 3)
      !********************************************************************** 
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      this%Scat_Phot_AxisZ = k_p 
 
      Temp_Vector(1) = zero 
      Temp_Vector(2) = zero 
      Temp_Vector(3) = one

      call this%Vector_Cross_Product( this%Scat_Phot_AxisZ, Temp_Vector, &
                                      this%Scat_Phot_AxisY )
      this%Scat_Phot_AxisY = this%Scat_Phot_AxisY / Vector3D_Length( this%Scat_Phot_AxisY )
      call this%Vector_Cross_Product( this%Scat_Phot_AxisY, &
                                      this%Scat_Phot_AxisZ, this%Scat_Phot_AxisX )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      this%Scat_Phot_Pe_AxisZ = k_p 

      P_e(1) = this%Elec_Phot_sin_In_Elec_CF
      P_e(2) = zero
      P_e(3) = this%Elec_Phot_mu_In_Elec_CF

      !mu_kp_pe_ECF = Vector3D_Inner_Product( P_e, k_p )

      call this%Vector_Cross_Product( this%Scat_Phot_Pe_AxisZ, P_e, this%Scat_Phot_Pe_AxisY )
      this%Scat_Phot_Pe_AxisY = this%Scat_Phot_Pe_AxisY / Vector3D_Length( this%Scat_Phot_Pe_AxisY )
      call this%Vector_Cross_Product( this%Scat_Phot_Pe_AxisY, &
                                      this%Scat_Phot_Pe_AxisZ, this%Scat_Phot_Pe_AxisX )

      mu_ya_yb = Vector3D_Inner_Product( this%Scat_Phot_AxisY, this%Scat_Phot_Pe_AxisY )
      mu_signs = Vector3D_Inner_Product( this%Scat_Phot_AxisX, this%Scat_Phot_Pe_AxisY )
      Rot_angle_k_Pe = - sign(one, mu_signs) * dacos( mu_ya_yb )

      !write(*, *)'tts=', ( P_e(3) - temp_cos * mupsi ) / temp_sin / sinpsi, mu_ya_yb 
      CALL this%StokesPara_Rotation_Matrix( Rot_angle_k_Pe, this%Vector_Stokes4_ECF_scat )
  
      P_e(1) = this%Elec_Phot_sin * dcos( this%Elec_Phot_phi )
      P_e(2) = this%Elec_Phot_sin * dsin( this%Elec_Phot_phi )
      P_e(3) = this%Elec_Phot_mu

      !CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
      !                   this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_3X3 )

      CALL Matrix_Multiplication13X33_Sub( P_e, &
                         this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_1X3 )
      mu_signs = Vector3D_Inner_Product( this%Scat_Phot3k_CF, Temp_Matrix_1X3 )
      k_p2 = Temp_Matrix_1X3 - mu_signs * this%Scat_Phot3k_CF
      this%f4_scat_CF(2: 4) = k_p2 / Vector3D_Length( k_p2 )
      this%f4_scat_CF(1) = zero

      !this%f4_scat_CovCF = this%f4_scat_CF 
      !this%f4_scat_CovCF(1) = - this%f4_scat_CF(1)
      !write(*, fmt="(' ', A5, 2ES18.7)")'s3=', Vector3D_Inner_Product( &
          !this%f4_scat_CF(2:4), this%Scat_Phot3k_CF), Vector3D_Length( this%Scat_Phot3k_CF )

      end subroutine Get_Scattered_Stokes_Vector_And_f_CF_Sub


!*******************************************************************************************************
      subroutine Set_Phot_f4_Tetrad_In_Elec_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp) :: vector1(1: 3)
      !**********************************************************************

      !this%Phot_f4_AxisZ(1) = this%Phot4k_In_Elec_CF(2) / this%Phot4k_In_Elec_CF(1)
      !this%Phot_f4_AxisZ(2) = this%Phot4k_In_Elec_CF(3) / this%Phot4k_In_Elec_CF(1)
      !this%Phot_f4_AxisZ(3) = this%Phot4k_In_Elec_CF(4) / this%Phot4k_In_Elec_CF(1)
      this%Phot_f4_AxisZ = this%Phot4k_In_Elec_CF(2: 4) / Vector3D_Length( this%Phot4k_In_Elec_CF(2: 4) )

      vector1 = this%f4_In_Elec_CF(2: 4) - this%f4_In_Elec_CF(1) / this%Phot4k_In_Elec_CF(1) * &
                this%Phot4k_In_Elec_CF(2: 4)
      this%Phot_f4_AxisX = vector1 / Vector3D_Length( vector1 )
      !this%Phot_f4_AxisX(1) = this%f4_In_Elec_CF(2)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(2) = this%f4_In_Elec_CF(3)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(3) = this%f4_In_Elec_CF(4)! / this%f4_In_Elec_CF(1)
       
      call this%Vector_Cross_Product( this%Phot_f4_AxisZ, &
                   this%Phot_f4_AxisX, this%Phot_f4_AxisY ) 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**f5555 === ', &
            !this%Phot4k_In_Elec_CF, this%f4_In_Elec_CF
            !Vector3D_Inner_Product(this%Phot_f4_AxisY, this%Phot_f4_AxisY), &
            !Vector3D_Inner_Product(this%Phot_f4_AxisZ, this%Phot_f4_AxisX), &
            !Vector3D_Inner_Product(this%Phot_f4_Axisz, this%Phot_f4_Axisz), &
            !Vector4D_Inner_Product_Mski( this%Phot4k_In_Elec_CF, this%f4_In_Elec_CF ) 
      !write(unit = *, fmt = *)'************************************************************'
      
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 1) = this%Phot_f4_AxisX(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 2) = this%Phot_f4_AxisX(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 3) = this%Phot_f4_AxisX(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 1) = this%Phot_f4_AxisY(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 2) = this%Phot_f4_AxisY(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 3) = this%Phot_f4_AxisY(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 1) = this%Phot_f4_AxisZ(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 2) = this%Phot_f4_AxisZ(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 3) = this%Phot_f4_AxisZ(3)

      end subroutine Set_Phot_f4_Tetrad_In_Elec_CF_Sub
 

!*******************************************************************************************************
      subroutine Compton_Scattering_With_Polar_StokesVec_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      !logical, intent(in) :: test
      real(mcp) :: epsi, r, phip, Cos_Theta1, Sin_Theta1, Sin_Theta12, Phi1, B_coef, C_coef
      real(mcp) :: N_temp, Vector_Length
      real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ
      real(mcp), dimension(1:3) :: e_p, e_f, e_y, e_yp, e_xp
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      real(mcp), dimension(1:3) :: Temp_Vector
      real(mcp), dimension(1:3) :: f3_scat_e_tilde
      real(mcp), dimension(1:3) :: f3_scat_pp
      real(mcp), dimension(1:3) :: Scat_Phot4k_In_Phot_Frame, k_p
      real(mcp), dimension(1:4) :: f4_scat_e, tempV
      real(mcp) :: psi_pp,  beta
      real(mcp) :: a, b, c
      real(mcp) :: F0, F11, F1, F22, F33, tp_vales, Q_xi, U_xi
      Complex*16 :: roots3(1: 3)
      integer :: del
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epsi = two * dabs( this%Phot4k_In_Elec_CF(1) ) / mec2 
      if( epsi >= 1.D-6 )then
          r = this%ComptonScaPhoSampling(epsi)
 
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r  
          Cos_Theta1 = one - ( r - one ) * two / epsi 
      else 
          a = zero
          b = three
          c = four - ranmar() * eight
          call root3(a, b, c, roots3(1), roots3(2), roots3(3), del)
          if( del == 1 )then
              Cos_Theta1 = real( roots3(1) )
              !write(*,*)'sss=', Cos_Theta1 , roots3
          else
              stop
          endif
          r = one + epsi * ( one - Cos_Theta1 ) / two
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r  
          !write(*,*)'sss=', Cos_Theta1 , roots3
      endif 
  
      if ( dabs( Cos_Theta1 ) < one ) then
          Sin_Theta12 = one - Cos_Theta1**2
          Sin_Theta1 = dsqrt( Sin_Theta12 )
      else
          Sin_Theta1 = zero
          Sin_Theta12 = zero
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      F0 = r + one /r - Sin_Theta12
      F1 = Sin_Theta12
      F11 = one + Cos_Theta1**2
      F22 = two * Cos_Theta1
      F33 = ( r + one / r ) * Cos_Theta1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      tp_vales = Sin_Theta12 / F0 / twopi
      Q_xi = this%Vector_Stokes4_ECF(2) / this%Vector_Stokes4_ECF(1)
      U_xi = this%Vector_Stokes4_ECF(3) / this%Vector_Stokes4_ECF(1)
      B_coef = tp_vales * Q_xi
      C_coef = tp_vales * U_xi
      Phi1 = this%Compton_Scat_Sampling_Phi_StokesVec( B_coef, C_coef ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      N_temp = F0 + F1 * ( Q_xi * dcos( two*phi1 ) + U_xi * dsin( two*phi1 ) )
      this%Vector_Stokes4_ECF_scat(1) = this%Vector_Stokes4_ECF(2)
      this%Vector_Stokes4_ECF_scat(2) = ( F1 + F11 * ( Q_xi * dcos( two*phi1 ) + &
                     U_xi * dsin( two*phi1 ) ) ) / N_temp
      this%Vector_Stokes4_ECF_scat(3) =  F22 * ( - Q_xi * dsin( two*phi1 ) + &
                     U_xi * dcos( two*phi1 ) ) / N_temp 
      this%Vector_Stokes4_ECF_scat(4) =  F33 * this%Vector_Stokes4_ECF(4) / N_temp 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      !**********************************************************************
      !*** To obtain the Scattered 4 momentum of photon in the CF 
      !********************************************************************** 
      k_p(1) = Sin_Theta1 * DCOS( Phi1 )
      k_p(2) = Sin_Theta1 * DSIN( Phi1 )
      k_p(3) = Cos_Theta1
       
      call this%Get_Scattered_Phot4k_CF( k_p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      call this%Get_Scattered_Stokes_Vector_And_f_CF( k_p )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
      end subroutine Compton_Scattering_With_Polar_StokesVec_Sub


!*******************************************************************************************************
      subroutine Compton_Scattering_With_Polarization_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      !logical, intent(in) :: test
      real(mcp) :: epsi, r, phip, Cos_Theta1, Sin_Theta1, Phi1, B_coef
      real(mcp) :: N_temp, Vector_Length
      real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ
      real(mcp), dimension(1:3) :: e_p, e_f, e_y, e_yp, e_xp
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      real(mcp), dimension(1:3) :: Temp_Vector
      real(mcp), dimension(1:3) :: f3_scat_e_tilde
      real(mcp), dimension(1:3) :: f3_scat_pp
      real(mcp), dimension(1:3) :: Scat_Phot4k_In_Phot_Frame
      real(mcp), dimension(1:4) :: f4_scat_e, tempV
      real(mcp) :: psi_pp,  beta
      real(mcp) :: a, b, c
      Complex*16 :: roots3(1: 3)
      integer :: del
      !****************************************************************************

      epsi = two * dabs( this%Phot4k_In_Elec_CF(1) ) / mec2 
      if( epsi >= 1.D-6 )then
          r = this%ComptonScaPhoSampling(epsi)
 
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r  
          Cos_Theta1 = one - ( r - one ) * two / epsi 
      else 
          a = zero
          b = three
          c = four - ranmar() * eight
          call root3(a, b, c, roots3(1), roots3(2), roots3(3), del)
          if( del == 1 )then
              Cos_Theta1 = real( roots3(1) )
              !write(*,*)'sss=', Cos_Theta1 , roots3
          else
              stop
          endif
          r = one + epsi * ( one - Cos_Theta1 ) / two
          this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_In_Elec_CF(1) / r  
          !write(*,*)'sss=', Cos_Theta1 , roots3
      endif
      this%r_one_hvlmec2_one_cosE = r
      !write(*,*)'sss=',  r, epsi, this%Phot4k_In_Elec_CF(1), this%Scattered_Phot4k_In_Elec_CF(1), &
      !       this%Elec_gama, this%Elec_V
      !Cos_Theta1 = one - ( r - one ) * two / epsi
      !**********************************************************************************
      if ( dabs( Cos_Theta1 ) < one ) then
          Sin_Theta1 = dsqrt( one - Cos_Theta1**2 )
      else
          Sin_Theta1 = zero
      endif
      !**********************************************************************************
      B_coef = Sin_Theta1**2 * this%delta_pd / ( r + one / r - Sin_Theta1**2 ) / twopi
      Phi1 = Compton_Scat_Sampling_Phi( B_coef )
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'** Cos_Theta1, Sin_Theta1, Phi1 === ',Cos_Theta1, Sin_Theta1, Phi1
      !write(unit = *, fmt = *)'************************************************************' 
      
      this%Cos_Theta_Scat = Cos_Theta1
      this%Sin_Theta_Scat = Sin_Theta1

      this%Q_sp = this%delta_pd * DCOS( two*Phi1 )
      this%U_sp = - this%delta_pd * DSIN( two*Phi1 )
      !**********************************************************************************
      N_temp = r + one / r - Sin_Theta1**2 * ( one + this%Q_sp )
      this%Q_sp_scat = - ( Sin_Theta1**2 - ( one + Cos_Theta1**2) * this%Q_sp ) / N_temp
      this%U_sp_scat =  two * Cos_Theta1 * this%U_sp / N_temp
      this%delta_pd_scat = DSQRT( this%Q_sp_scat**2 + this%U_sp_scat**2 )
      this%N_scat = N_temp
      !**********************************************************************************
   !write(unit = *, fmt = *)'************************************************************' 
   !write(unit = *, fmt = *)'** Delta_pd === ', this%delta_pd_scat,N_temp, this%delta_pd_scat , this%Q_sp_scat, &
    !this%U_sp_scat, r, Sin_Theta1, this%delta_pd, Phi1
  ! write(unit = *, fmt = *)'************************************************************' 
      
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(1) = Sin_Theta1 * DCOS( Phi1 )
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(2) = Sin_Theta1 * DSIN( Phi1 )
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(3) = Cos_Theta1
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one

      call this%Vector_Cross_Product( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ, e_p, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX )
      Vector_Length = Vector3D_Length( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX )
      If ( Vector_Length /= zero ) then
            Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX = &
          - Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX / Vector_Length
          call this%Vector_Cross_Product( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY )
      Else
          e_f(1) = one
          e_f(2) = zero
          e_f(3) = zero
          e_y(1) = zero
          e_y(2) = one
          e_y(3) = zero
          e_yp = DCOS( Phi1 ) * e_f + DSIN( Phi1 ) * e_y
          CALL this%Vector_Cross_Product( e_yp, e_p, e_xp )
          Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX = e_xp
          Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY = DSIGN( one, Cos_Theta1 ) * e_yp
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

      f3_scat_pp = DCOS(psi_pp) * Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX + &
                   DSIN(psi_pp) * Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY

      CALL Matrix_Multiplication13X33_Sub( f3_scat_pp, &
                  this%Matrix_Of_Tetrad_Of_Phot_f4_Axis, f3_scat_e_tilde )

      f4_scat_e(1) = this%Elec_gama * this%Elec_V * f3_scat_e_tilde(3)
      f4_scat_e(2) = f3_scat_e_tilde(1)
      f4_scat_e(3) = f3_scat_e_tilde(2)
      f4_scat_e(4) = this%Elec_gama * f3_scat_e_tilde(3)

      CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
                         this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_3X3 )
      CALL Matrix_Multiplication13X33_Sub( f4_scat_e(2:4), Temp_Matrix_3X3, &
                                                         Temp_Matrix_1X3 )
      this%f4_scat_CF(1) = f4_scat_e(1)
      this%f4_scat_CF(2:4) = Temp_Matrix_1X3
      this%f4_scat_CovCF = this%f4_scat_CF 
      this%f4_scat_CovCF(1) = - this%f4_scat_CF(1)
       
      If ( Vector_Length == 80.d0 ) then
          write(unit = *, fmt = *)'************************************************************' 
          write(unit = *, fmt = *)'************************************************************' 
         this%test_quantity(8) = Vector_Length
         this%temp_px = Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX
         this%temp_py = Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY
         this%temp_pz = Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ
         this%temp_f4 = f3_scat_pp
          write(unit = *, fmt = *)'sdfsfsfsfsfsfsfsfsf1=  === ', Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX,&
                Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY,  Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ
          write(unit = *, fmt = *)'sdfsfsfsfsfsfsfsfsf2=  === ', f3_scat_pp, this%f4_scat_CF
          write(unit = *, fmt = *)'sdfsfsfsfsfsfsfsfsf3=  === ', DSIN(psi_pp), DCOS(psi_pp), psi_pp
          write(unit = *, fmt = *)'sdfsfsfsfsfsfsfsfsf4=  === ', this%U_sp_scat, this%Q_sp_scat, Sin_Theta1, &
             Cos_Theta1, phi1, this%delta_pd
          write(unit = *, fmt = *)'************************************************************'   
          write(unit = *, fmt = *)'************************************************************' 
          !stop
      endif
      !**********************************************************************
      !*** To obtain the Scattered 4 momentum of photon in the CF 
      !**********************************************************************
      !*************************************************************************************************
      !************************************************************************************************* 
      Scat_Phot4k_In_Phot_Frame(1) = Sin_Theta1 * DCOS( Phi1 )
      Scat_Phot4k_In_Phot_Frame(2) = Sin_Theta1 * DSIN( Phi1 )
      Scat_Phot4k_In_Phot_Frame(3) = Cos_Theta1
      
      Scat_Phot4k_In_Phot_Frame = Scat_Phot4k_In_Phot_Frame * &
                             this%Scattered_Phot4k_In_Elec_CF(1)

      CALL Matrix_Multiplication13X33_Sub( Scat_Phot4k_In_Phot_Frame, &
                  this%Matrix_Of_Tetrad_Of_Phot_f4_Axis, this%Scattered_Phot4k_In_Elec_CF(2:4) )
      If ( Vector_Length == 880.d0 ) then
        write(unit = *, fmt = *)'************************************************************' 
        write(unit = *, fmt = *)'**gggggggg  ssss === ', &
            Vector3D_Inner_Product( this%Scattered_Phot4k_In_Elec_CF(2:4), &
          this%Scattered_Phot4k_In_Elec_CF(2:4) ) /  this%Scattered_Phot4k_In_Elec_CF(1)**2, &
          Vector3D_Inner_Product( Scat_Phot4k_In_Phot_Frame, Scat_Phot4k_In_Phot_Frame) / &
              this%Scattered_Phot4k_In_Elec_CF(1)**2
        write(unit = *, fmt = *)'************************************************************'
      endif

      Scattered_Phot4k_In_Elec(1) = this%Elec_gama * ( this%Scattered_Phot4k_In_Elec_CF(1) + &
                                         this%Elec_V * this%Scattered_Phot4k_In_Elec_CF(4) )
      Scattered_Phot4k_In_Elec(4) = this%Elec_gama * ( this%Scattered_Phot4k_In_Elec_CF(4) + &
                                         this%Elec_V * this%Scattered_Phot4k_In_Elec_CF(1) )
      Scattered_Phot4k_In_Elec(2) = this%Scattered_Phot4k_In_Elec_CF(2) ! x component
      Scattered_Phot4k_In_Elec(3) = this%Scattered_Phot4k_In_Elec_CF(3) ! y component
      !*************************************************************************************************
      CALL Matrix_Multiplication13X33_Sub( Scattered_Phot4k_In_Elec(2:4), &
                                       Temp_Matrix_3X3, Temp_Matrix_1X3 ) 
      !*************************************************************************************************
      this%Scattered_Phot4k_CF(1) = Scattered_Phot4k_In_Elec(1)
      this%Scattered_Phot4k_CF(2:4) = Temp_Matrix_1X3
      this%Scattered_Phot4k_CovCF = this%Scattered_Phot4k_CF
      this%Scattered_Phot4k_CovCF(1) = - this%Scattered_Phot4k_CF(1)
      tempV = Scattered_Phot4k_In_Elec
      if( this%Scattered_Phot4k_CF(4) / dabs(this%Scattered_Phot4k_CF(1)) > 1.D0 )then
          write(*, *)'*****************************************************************************'
          write(*, *)'mms1111==', this%Scattered_Phot4k_CF / dabs(this%Scattered_Phot4k_CF(1))
          write(*, *)'mms2222==', Scattered_Phot4k_In_Elec / dabs(Scattered_Phot4k_In_Elec(1))
          write(*, *)'mms3333==', this%Scattered_Phot4k_In_Elec_CF / &
                                 dabs(this%Scattered_Phot4k_In_Elec_CF(1)) 
          write(*, *)'mms4444==', Scat_Phot4k_In_Phot_Frame / this%Scattered_Phot4k_In_Elec_CF(1), &
                Sin_Theta1, Cos_Theta1, Sin_Theta1**2 + Cos_Theta1**2
          write(*, *)'mms5555==', Vector4D_Inner_Product_Mski( this%Scattered_Phot4k_CF, &
                              this%Scattered_Phot4k_CF ) / this%Scattered_Phot4k_CF(1)**2, &
             Vector3D_Inner_Product(this%Phot_f4_AxisX, this%Phot_f4_AxisX), &
             Vector3D_Inner_Product(this%Phot_f4_AxisY, this%Phot_f4_AxisY), &
             Vector3D_Inner_Product(this%Phot_f4_AxisZ, this%Phot_f4_AxisX), &
             Vector3D_Inner_Product(this%Phot_f4_Axisz, this%Phot_f4_Axisz), &
             Vector4D_Inner_Product_Mski( this%Phot4k_In_Elec_CF, this%f4_In_Elec_CF )
          write(*, *)'mms6666=', Vector3D_Inner_Product( this%PhotAxisX, this%PhotAxisX ), &
                                 Vector3D_Inner_Product( this%PhotAxisy, this%PhotAxisy ), &
                                 Vector3D_Inner_Product( this%PhotAxisz, this%PhotAxisz ), &
                                 Vector3D_Inner_Product( this%ELecAxisX, this%ELecAxisX ), &
                                 Vector3D_Inner_Product( this%ELecAxisy, this%ELecAxisy ), &
                                 Vector3D_Inner_Product( this%ELecAxisz, this%ELecAxisz ), &
                  ( r - one ) * two / epsi, r, epsi
          write(*, *)'*****************************************************************************'
      endif
      !this%Scal_Factor_after = this%Scattered_Phot4k_CF(1) / this%Scattered_Phot4k_In_Elec_CF(1)
      !*************************************************************************************************
      !************************************************************************************************* 
    !write(unit = *, fmt = *)'************************************************************' 
      If ( Vector_Length == 7770.d0 ) then
    write(unit = *, fmt = *)'** f44 * k44   === ', -this%Scattered_Phot4k_CF(1)*this%f4_scat_CF(1) + &
                                                   this%Scattered_Phot4k_CF(2)*this%f4_scat_CF(2) + &
                                                   this%Scattered_Phot4k_CF(3)*this%f4_scat_CF(3) + &
                                                   this%Scattered_Phot4k_CF(4)*this%f4_scat_CF(4), &
                            -this%Scattered_Phot4k_CF(1)*this%Scattered_Phot4k_CF(1) + &
                            this%Scattered_Phot4k_CF(2)*this%Scattered_Phot4k_CF(2) + &
                            this%Scattered_Phot4k_CF(3)*this%Scattered_Phot4k_CF(3) + &
                            this%Scattered_Phot4k_CF(4)*this%Scattered_Phot4k_CF(4), &
                           -tempV(1)*tempV(1) +  tempV(2)*tempV(2) +  tempV(3)*tempV(3) +  tempV(4)*tempV(4)  
     write(unit = *, fmt = *)'************************************************************'
      this%test_quantity(1) = r
      this%test_quantity(2) = this%Phot4k_In_Elec_CF(1)
      this%test_quantity(3) = this%Scattered_Phot4k_In_Elec_CF(1)
      this%test_quantity(4) = Cos_Theta1
      this%test_quantity(5) = Sin_Theta1
      this%test_quantity(6) = B_coef
      this%test_quantity(7) = Phi1
   endif
      !****************************************************************************

      end subroutine Compton_Scattering_With_Polarization_Sub



!*******************************************************************************************************
      subroutine Static_Compton_Scattering_With_Polarization_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      !logical, intent(in) :: test
      real(mcp) :: epsi, r, phip, Cos_Theta1, Sin_Theta1, Phi1, B_coef
      real(mcp) :: N_temp, Vector_Length, E_ini
      real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY
      real(mcp), dimension(1:3) :: Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ
      real(mcp), dimension(1:3) :: e_p, e_f, e_y, e_yp, e_xp
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      real(mcp), dimension(1:3) :: Temp_Vector
      real(mcp), dimension(1:3) :: f3_scat_e_tilde
      real(mcp), dimension(1:3) :: f3_scat_pp
      real(mcp), dimension(1:3) :: Scat_Phot4k_In_Phot_Frame
      real(mcp), dimension(1:4) :: f4_scat_e, tempV
      real(mcp) :: psi_pp,  beta
      real(mcp) :: a, b, c
      Complex*16 :: roots3(1: 3)
      integer :: del
      !****************************************************************************
 
      E_ini = dabs( this%Phot4k_CtrCF(1) )
      epsi = two * E_ini / mec2 
      if( epsi >= 1.D-6 )then
          r = this%ComptonScaPhoSampling(epsi)
 
          this%Scattered_Phot4k_In_Elec_CF(1) = E_ini / r  
          Cos_Theta1 = one - ( r - one ) * two / epsi 
      else 
          a = zero
          b = three
          c = four - ranmar() * eight
          call root3(a, b, c, roots3(1), roots3(2), roots3(3), del)
          if( del == 1 )then
              Cos_Theta1 = real( roots3(1) )
              !write(*,*)'sss=', Cos_Theta1 , roots3
          else
              stop
          endif
          r = one + epsi * ( one - Cos_Theta1 ) / two
          this%Scattered_Phot4k_In_Elec_CF(1) = E_ini / r  
          !write(*,*)'sss=', Cos_Theta1 , roots3
      endif
      this%r_one_hvlmec2_one_cosE = r
      !write(*,*)'sss=',  r, epsi, this%Phot4k_In_Elec_CF(1), this%Scattered_Phot4k_In_Elec_CF(1), &
      !       this%Elec_gama, this%Elec_V 
      !**********************************************************************************
      if ( dabs( Cos_Theta1 ) < one ) then
          Sin_Theta1 = dsqrt( one - Cos_Theta1**2 )
      else
          Sin_Theta1 = zero
      endif
      !**********************************************************************************
      B_coef = Sin_Theta1**2 * this%delta_pd / ( r + one / r - Sin_Theta1**2 ) / twopi
      Phi1 = Compton_Scat_Sampling_Phi( B_coef )
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'** Cos_Theta1, Sin_Theta1, Phi1 === ',Cos_Theta1, Sin_Theta1, Phi1
      !write(unit = *, fmt = *)'************************************************************' 
      
      this%Cos_Theta_Scat = Cos_Theta1
      this%Sin_Theta_Scat = Sin_Theta1

      this%Q_sp = this%delta_pd * DCOS( two*Phi1 )
      this%U_sp = - this%delta_pd * DSIN( two*Phi1 )
      !**********************************************************************************
      N_temp = r + one / r - Sin_Theta1**2 * ( one + this%Q_sp )
      this%Q_sp_scat = - ( Sin_Theta1**2 - ( one + Cos_Theta1**2) * this%Q_sp ) / N_temp
      this%U_sp_scat =  two * Cos_Theta1 * this%U_sp / N_temp
      this%delta_pd_scat = DSQRT( this%Q_sp_scat**2 + this%U_sp_scat**2 )
      this%N_scat = N_temp
      !**********************************************************************************
   !write(unit = *, fmt = *)'************************************************************' 
   !write(unit = *, fmt = *)'** Delta_pd === ', this%delta_pd_scat,N_temp, this%delta_pd_scat , this%Q_sp_scat, &
    !this%U_sp_scat, r, Sin_Theta1, this%delta_pd, Phi1
  ! write(unit = *, fmt = *)'************************************************************' 
      
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(1) = Sin_Theta1 * DCOS( Phi1 )
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(2) = Sin_Theta1 * DSIN( Phi1 )
      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ(3) = Cos_Theta1
      e_p(1) = zero
      e_p(2) = zero
      e_p(3) = one

      call this%Vector_Cross_Product( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ, e_p, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX )
      Vector_Length = Vector3D_Length( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX )
      If ( Vector_Length /= zero ) then
            Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX = &
          - Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX / Vector_Length
          call this%Vector_Cross_Product( Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisZ, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX, &
                                      Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY )
      Else
          e_f(1) = one
          e_f(2) = zero
          e_f(3) = zero
          e_y(1) = zero
          e_y(2) = one
          e_y(3) = zero
          e_yp = DCOS( Phi1 ) * e_f + DSIN( Phi1 ) * e_y
          CALL this%Vector_Cross_Product( e_yp, e_p, e_xp )
          Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX = e_xp
          Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY = DSIGN( one, Cos_Theta1 ) * e_yp
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

      f3_scat_pp = DCOS(psi_pp) * Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisX + &
                   DSIN(psi_pp) * Scat_Phot_Tetrad_In_Phot_f4_Tetrad_AxisY

      CALL Matrix_Multiplication13X33_Sub( f3_scat_pp, &
                  this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF, f3_scat_e_tilde )
  
      !CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_ElecAxis, &
      !                   this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_3X3 )
      CALL Matrix_Multiplication13X33_Sub( f3_scat_e_tilde, this%Matrix_Of_Tetrad_Of_PhotAxis, &
                                                         Temp_Matrix_1X3 )
      this%f4_scat_CF(1) = zero
      this%f4_scat_CF(2:4) = Temp_Matrix_1X3
      this%f4_scat_CovCF = this%f4_scat_CF 
      this%f4_scat_CovCF(1) = - this%f4_scat_CF(1)
 
      !**********************************************************************
      !*** To obtain the Scattered 4 momentum of photon in the CF 
      !**********************************************************************   
      Scat_Phot4k_In_Phot_Frame(1) = Sin_Theta1 * DCOS( Phi1 )
      Scat_Phot4k_In_Phot_Frame(2) = Sin_Theta1 * DSIN( Phi1 )
      Scat_Phot4k_In_Phot_Frame(3) = Cos_Theta1
      
      Scat_Phot4k_In_Phot_Frame = Scat_Phot4k_In_Phot_Frame * &
                             this%Scattered_Phot4k_In_Elec_CF(1)

      CALL Matrix_Multiplication13X33_Sub( Scat_Phot4k_In_Phot_Frame, &
                  this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF, Scattered_Phot4k_In_Elec(2:4) ) 
 
      !*************************************************************************************************
      CALL Matrix_Multiplication13X33_Sub( Scattered_Phot4k_In_Elec(2:4), &
                                       this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_1X3 ) 
      !*************************************************************************************************
      !write(*, *)'sc1==', Scat_Phot4k_In_Phot_Frame, Scattered_Phot4k_In_Elec, Temp_Matrix_1X3
      !write(*, *)'sc2==', this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF, this%Matrix_Of_Tetrad_Of_PhotAxis

      this%Scattered_Phot4k_CF(1) = this%Scattered_Phot4k_In_Elec_CF(1)
      this%Scattered_Phot4k_CF(2:4) = Temp_Matrix_1X3 
      this%Scattered_Phot4k_CovCF = this%Scattered_Phot4k_CF
      this%Scattered_Phot4k_CovCF(1) = - this%Scattered_Phot4k_CF(1) 
      !************************************************************************************************* 
      !write(*, *)'sc==', Scat_Phot4k_In_Phot_Frame, Scattered_Phot4k_In_Elec(2:4), Temp_Matrix_1X3
      !write(*, *)'sc==', this%Matrix_Of_Tetrad_Of_Phot_f4_Axis_PCF, this%Matrix_Of_Tetrad_Of_PhotAxis, &
      !                 this%Scattered_Phot4k_In_Elec_CF(1)
 
      end subroutine Static_Compton_Scattering_With_Polarization_Sub
 

!*******************************************************************************************************
      Subroutine Tompson_Scattering_WithOut_Polarization_Sub( this )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp) :: epsi, r, phip, mupsi, mu_tilde_p, sinmu_tilde, sinpsi, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp, beta
      real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      real(mcp), dimension(1:3) :: Temp_Vector, f3_scat_e_tilde
      real(mcp), dimension(1:4) :: f4_scat_e
      real(mcp) :: psi_p
      !****************************************************************************

      !****************************************************************************
      real(mcp) :: a, b, c, Etilde
      Complex*16 :: roots3(1: 3)
      integer :: del, cases_f
      !****************************************************************************

      !mutilde = this%Elec_Phot_mu_In_Elec_CF
      Etilde  = this%Phot4k_CtrCF(1) 
      epsi = two * dabs( Etilde ) / mec2
      if( epsi >= 1.D-4 )then
          r = this%ComptonScaPhoSampling(epsi)
 
          !this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_CtrCF(1) / r  
          mupsi = one - ( r - one ) * two / epsi 
      else 
          a = zero
          b = three
          c = four - ranmar() * eight
          call root3(a, b, c, roots3(1), roots3(2), roots3(3), del)
          if( del == 1 )then
              mupsi = real( roots3(1) )
          else
              stop
          endif
          r = one! + epsi * ( one - mupsi ) / two
      endif
      phip = twopi * ranmar()

      !this%Scattered_Phot4k_In_Elec_CF(1) = this%Phot4k_CtrCF(1) / r

      this%r_one_hvlmec2_one_cosE = r
 
      !if( epsi > 1.D-10 )then
      !    mupsi = one - ( r - one ) * two / epsi
      !else
      !    mupsi = one - (one - this%r11) * two
      !endif
      !**********************************************************************
 
      if ( dabs(mupsi) < one ) then
          sinpsi = dsqrt( one - mupsi**2 )
      else
          sinpsi = zero
      endif
      this%Cos_Theta_Scat = mupsi
      this%Sin_Theta_Scat = sinpsi
      
      !**********************************************************************************
      N_temp = r + one / r - this%Sin_Theta_Scat**2 *( one + this%delta_pd *  DCOS( two*phip ) )
      this%Q_sp_scat = - ( this%Sin_Theta_Scat**2 - (one + this%Cos_Theta_Scat**2) &
                                         * this%delta_pd * DCOS( two*phip ) ) / N_temp
      this%U_sp_scat = two * this%Cos_Theta_Scat * this%delta_pd * DSIN( two*phip ) / N_temp
      this%delta_pd_scat = DSQRT( this%Q_sp_scat**2 + this%U_sp_scat**2 )  
      this%N_scat = N_temp! * three / 16.D0 / pi
      !**********************************************************************
      !*** To obtain the Scattered 4 momentum of photon in the CF 
      !**********************************************************************
      Scattered_Phot4k_In_Elec(2) = sinpsi * dcos( phip )
      Scattered_Phot4k_In_Elec(3) = sinpsi * dsin( phip )
      Scattered_Phot4k_In_Elec(4) = mupsi
      CALL Matrix_Multiplication13X33_Sub( Scattered_Phot4k_In_Elec(2:4), &
                            this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_1X3 )
      
      this%Scattered_Phot4k_CF(1) = this%Phot4k_CtrCF(1) / r
      this%Scattered_Phot4k_CF(2:4) = Temp_Matrix_1X3 * this%Scattered_Phot4k_CF(1)
      this%Scattered_Phot4k_CovCF = this%Scattered_Phot4k_CF
      this%Scattered_Phot4k_CovCF(1) = - this%Scattered_Phot4k_CF(1)

      if( this%Scattered_Phot4k_CF(4) / dabs(this%Scattered_Phot4k_CF(1)) > 1.D0 )then
          write(*, *)'mms0000==', this%Scattered_Phot4k_CF / dabs(this%Scattered_Phot4k_CF(1))
      endif
      !*************************************************************************************************
      !*************************************************************************************************
      this%Scat_Phot_AxisZ(1:3) = Scattered_Phot4k_In_Elec(2:4)
 
      Temp_Vector(1) = zero !this%Phot4k_In_Elec_CF(2:4) / this%Phot4k_In_Elec_CF(1)
      Temp_Vector(2) = zero
      Temp_Vector(3) = one

      call this%Vector_Cross_Product( Temp_Vector, this%Scat_Phot_AxisZ, &
                                      this%Scat_Phot_AxisY )
      this%Scat_Phot_AxisY = this%Scat_Phot_AxisY / Vector3D_Length( this%Scat_Phot_AxisY )
      call this%Vector_Cross_Product( this%Scat_Phot_AxisY, &
                                      this%Scat_Phot_AxisZ, this%Scat_Phot_AxisX )
      beta = DATAN( DABS( this%U_sp_scat / this%Q_sp_scat ) ) / two
      If ( this%Q_sp_scat > zero .and. this%U_sp_scat > zero ) then
          Psi_p = beta
          cases_f = 1
      else If ( this%Q_sp_scat > zero .and. this%U_sp_scat < zero ) then
          Psi_p = pi - beta
          cases_f = 2
      else If ( this%Q_sp_scat < zero .and. this%U_sp_scat > zero ) then
          Psi_p = pi / two - beta
          cases_f = 3
      else If ( this%Q_sp_scat < zero .and. this%U_sp_scat < zero ) then
          Psi_p = beta + pi / two
          cases_f = 4
      else If ( this%Q_sp_scat > zero .and. this%U_sp_scat == zero ) then
          Psi_p = zero
          cases_f = 5
      else If ( this%Q_sp_scat < zero .and. this%U_sp_scat == zero ) then
          Psi_p = Pi / two
          cases_f = 6
      else If ( this%Q_sp_scat == zero .and. this%U_sp_scat > zero ) then
          Psi_p = Pi / four
          cases_f = 7
      else If ( this%Q_sp_scat == zero .and. this%U_sp_scat < zero ) then
          Psi_p = - Pi / four
          cases_f = 8
      Endif

      !psi_p = DATAN( this%U_sp_scat / this%Q_sp_scat ) / two

      f3_scat_e_tilde = DCOS(psi_p) * this%Scat_Phot_AxisX + &
                        DSIN(psi_p) * this%Scat_Phot_AxisY
      f4_scat_e(1) = zero !this%Elec_gama * this%Elec_V * f3_scat_e_tilde(3)
      f4_scat_e(2) = f3_scat_e_tilde(1)
      f4_scat_e(3) = f3_scat_e_tilde(2)
      f4_scat_e(4) = f3_scat_e_tilde(3)
      CALL Matrix_Multiplication13X33_Sub( f4_scat_e(2:4), &
                      this%Matrix_Of_Tetrad_Of_PhotAxis, Temp_Matrix_1X3 )

      this%f4_scat_CF(1) = f4_scat_e(1)
      this%f4_scat_CF(2:4) = Temp_Matrix_1X3
      this%f4_scat_CovCF = this%f4_scat_CF 
      this%f4_scat_CovCF(1) = - this%f4_scat_CF(1)
      !write(unit = *, fmt = *)'**f111 === ', &
      !      Vector4D_Inner_Product( this%f4_scat_CF, this%Scattered_Phot4k_CF)

      end Subroutine Tompson_Scattering_WithOut_Polarization_Sub


!*******************************************************************************************************
      subroutine Set_Phot_f4_Tetrad_In_Phot_CF_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton_KN) :: this 
      real(mcp) :: vector1(1: 3)
      real(mcp), dimension(1:3, 1:3) :: InvMatrix_Phot, temM
      real(mcp), dimension(1:4) :: f4_scat_e, f4_In_MF, phot_In_MF

      !CALL Inverse_Matrix_Of_Matrix_3X3_Sub( this%Matrix_Of_Tetrad_Of_PhotAxis, InvMatrix_Phot )
      !CALL Matrix_Multiplication13X33_Sub( this%f4_CF(2:4), InvMatrix_Phot, f4_In_MF(2:4) )
      !CALL Matrix_Multiplication13X33_Sub( this%Phot4k_CtrCF(2: 4) / &
      !        Vector3D_Length( this%Phot4k_CtrCF(2: 4) ), InvMatrix_Phot, phot_In_MF(2:4) )
      !*******************************************************************************************

      !CALL Matrix_Multiplication33X33_Sub( this%Matrix_Of_Tetrad_Of_PhotAxis, &
      !            InvMatrix_Phot, temM )
      !write(unit = *, fmt = *)'ss==', phot_In_MF
 
      this%Phot_f4_AxisZ = this%Phot4k_CtrCF(2: 4) / Vector3D_Length( this%Phot4k_CtrCF(2: 4) )
      !this%Phot_f4_AxisZ(1) = zero!this%Phot4k_CtrCF(2: 4) / Vector3D_Length( this%Phot4k_CtrCF(2: 4) )
      !this%Phot_f4_AxisZ(2) = zero
      !this%Phot_f4_AxisZ(3) = one

      vector1 = this%f4_CF(2:4) !f4_In_MF(2: 4) 
      this%Phot_f4_AxisX = vector1 / Vector3D_Length( vector1 )
      !this%Phot_f4_AxisX(1) = this%f4_In_Elec_CF(2)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(2) = this%f4_In_Elec_CF(3)! / this%f4_In_Elec_CF(1)
      !this%Phot_f4_AxisX(3) = this%f4_In_Elec_CF(4)! / this%f4_In_Elec_CF(1)
       
      call this%Vector_Cross_Product( this%Phot_f4_AxisZ, &
                   this%Phot_f4_AxisX, this%Phot_f4_AxisY ) 
      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'**f5555 === ', &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisX, this%Phot_f4_AxisY), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisX, this%Phot_f4_AxisZ), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisY, this%Phot_f4_AxisX), &
      !      Vector3D_Inner_Product(this%Phot_f4_AxisZ, this%Phot_f4_AxisX), &
      !      Vector4D_Inner_Product_Mski( this%Phot4k_In_Elec_CF, this%f4_In_Elec_CF ) 
      !write(unit = *, fmt = *)'************************************************************'
      
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 1) = this%Phot_f4_AxisX(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 2) = this%Phot_f4_AxisX(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(1, 3) = this%Phot_f4_AxisX(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 1) = this%Phot_f4_AxisY(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 2) = this%Phot_f4_AxisY(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(2, 3) = this%Phot_f4_AxisY(3)

      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 1) = this%Phot_f4_AxisZ(1)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 2) = this%Phot_f4_AxisZ(2)
      this%Matrix_Of_Tetrad_Of_Phot_f4_Axis(3, 3) = this%Phot_f4_AxisZ(3)

      end subroutine Set_Phot_f4_Tetrad_In_Phot_CF_Sub
 

!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************

      end module ScatterPhoton_KN






