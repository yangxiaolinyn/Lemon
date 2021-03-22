      module ScatterPhoton
      !use Basic_Variables_And_Methods
      use ChandraDiffuseReflection
      implicit none 

!*******************************************************************************************************
      type, public, extends(ChandDR) :: ScatPhoton 
!*******************************************************************************************************
          real(mcp), dimension(1:4) :: Phot4k_In_Elec_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_CF
          real(mcp), dimension(1:4) :: Scattered_Phot4k_CovCF 
          !*********************************************************************
          real(mcp) :: Cos_Theta_Scat
          real(mcp) :: Sin_Theta_Scat
          real(mcp) :: Phi_Scat
          real(mcp) :: Cos_Phi_Scat
          real(mcp) :: Sin_Phi_Scat  
          real(mcp) :: ds_phi 
          real(mcp) :: PolarArrayIQUVs10(1: 4)
          real(mcp) :: PolarArrayIQUVs30(1: 4)
          real(mcp) :: PolarArrayIQUVs60(1: 4)
          real(mcp) :: PolarArrayIQUVs80(1: 4)
      contains 
!*******************************************************************************************************
      procedure, public :: Tompson_Scat_With_Polarized_Diffuse_Reflection   =>   &
                           Tompson_Scat_With_Polarized_Diffuse_Reflection_Sub
      procedure, public :: Vector_Cross_Product   =>   Vector_Cross_Product_Sub
      procedure, public :: Get_Psi_IQUV_for_Estimation    =>   Get_Psi_IQUV_for_Estimation_Sub
      procedure, public :: Tompson_Scat_With_Polarized_Diffuse_Reflection2    =>   &
                           Tompson_Scat_With_Polarized_Diffuse_Reflection2_Sub
      procedure, public :: Calc_Phot_Inform_At_Observor_DiffRefl_phi_Estimats   =>   &
                           Calc_Phot_Inform_At_Observor_DiffRefl_phi_Estimats_Sub
      end type ScatPhoton
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      private :: Tompson_Scat_With_Polarized_Diffuse_Reflection_Sub
      private :: Vector_Cross_Product_Sub
      private :: Get_Psi_IQUV_for_Estimation_Sub
      private :: Tompson_Scat_With_Polarized_Diffuse_Reflection2_Sub
      private :: Calc_Phot_Inform_At_Observor_DiffRefl_phi_Estimats_Sub
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      contains 
!********************************************************************************
      subroutine Vector_Cross_Product_Sub(this, e1, e2, e3)
!********************************************************************************
      Implicit None
      class(ScatPhoton) :: this
      real(mcp),dimension(1:3), intent(in) :: e1, e2
      real(mcp),dimension(1:3), intent(out) :: e3
!**********************************************************************

      e3(1) = e1(2) * e2(3) - e1(3) * e2(2)    
      e3(2) = e1(3) * e2(1) - e1(1) * e2(3)    
      e3(3) = e1(1) * e2(2) - e1(2) * e2(1)  
      end subroutine Vector_Cross_Product_Sub

!*******************************************************************************************************
      Subroutine Tompson_Scat_With_Polarized_Diffuse_Reflection_Sub( this )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton) :: this 
      real(mcp) :: r, phip, smu, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp, beta, &
                   A_const, B_const, C_const, phi_ini, N_normal, N_temp2, &
                   cosphi_ini, sinphi_ini, cos2phi_ini, sin2phi_ini, smup, &
                   mu2, mup2, mu, mup

      real(mcp) :: A1, B1, t1, t2, Norm_C, f1, f2, f3, f4, f5, F0, g1, g2, g3, g4, g5
      real(mcp) :: UlI, QlI, scat_phi, sinDphi, cosDphi, sin2Dphi, cos2Dphi
      real(mcp) :: pb1, pb2, pb3, pb4, pb5, xi1
      real(mcp) :: P0(1: 3, 1: 3) = zero, P1(1: 3, 1: 3) = zero, P2(1: 3, 1: 3) = zero, &
                   PT(1: 3, 1: 3) = zero
      real(mcp) :: f_Q, f_U, f_V
      !****************************************************************************
      real(mcp) :: a, b, c, d
      Complex*16 :: roots3(1: 3), rts(1: 3)
      integer :: del, cases_f
      !****************************************************************************
 
      mup = this%Vector_of_Momentum_ini(3)
      mup2 = mup**2
      smup = dsqrt(one - mup2)
      QlI = this%Psi_Q / this%Psi_I
      UlI = this%Psi_U / this%Psi_I
      t1 = one + mup2 - (one - mup2) * QlI
      t2 = two * (one - mup2) * (one + QlI)
      A1 = t1 - t2
      B1 = t1 + t2
      a = A1
      b = zero
      c = 3.D0 * B1
      d = a + c - 16.D0 * ranmar()
    
      if( dabs( a ) > 1.D-7 )then
          call root3new(a, b, c, d, roots3, del)  
          if( del == 3 )then
              mu = real( roots3(3) )
              if( dabs( mu ) > one )then
                  mu = real( roots3(2) )
                  if( dabs( mu ) > one )then
                      mu = real( roots3(1) )
                      if( dabs( mu ) > one )then
                          write(*, *)'mms11==', roots3, A1, B1, a, b, c, d 
                          !call root3new(B_const, a, b, c, rts, del)  
                          write(*, *)'mms11==', - c / b
                          stop
                      endif
                  endif
              endif
          else
              mu = real( roots3(1) )
              if( dabs( mu ) > one .or. isnan( mu ) )then
                  write(*, *)'mms22==', A1, B1, mup, a, b, c
                  write(*, *)'mms33==', del, roots3  
                  write(*, *)'mms11==', rts, B_const, a, b, c,  QlI, this%Q_IQ, this%I_IQ
                  write(*, *)'mms11==', - c / b
                  stop
              endif 
          endif
      else 
          mu = - d / c
          if( dabs( mu ) > one  )then
              write(*, *)'mms41==', A1, B1, mup, a, b, c, d, mu
              stop 
          endif 
      endif  
      mu2 = mu**2 
 
      if ( dabs(mu) < one ) then
          smu = dsqrt( one - mu2 )
      else
          smu = zero
      endif  

      Norm_C = twopi * (mu2 * A1 + B1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      f1 = (mu2 * A1 + B1) !t1 * (one + mu2) + t2 * (one - mu2)
      f2 = four * mup * smup * mu * smu * (one + QlI)
      f3 = four * smup * mu * smu * UlI
      f4 = (one - mu2) * (one - mup2) - (one - mu2) * (one + mup2) * QlI
      f5 = - two * (one - mu2) * mup * UlI
      !F0 = f1 - dabs(f2) - dabs(f3) - dabs(f4) - dabs(f5)
      !p1 = F0 / Norm_C
      !p2 = dabs(f2) / Norm_C
      !p3 = dabs(f3) / Norm_C
      !p4 = dabs(f4) / Norm_C
      !p5 = dabs(f5) / Norm_C
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      g2 = f2 * this%cosphi_ini + f3 * this%sinphi_ini
      g3 = f2 * this%sinphi_ini - f3 * this%cosphi_ini
      g4 = f4 * this%cos2phi_ini + f5 * this%sin2phi_ini
      g5 = f4 * this%sin2phi_ini - f5 * this%cos2phi_ini
      !pb1 = ( f1 - dabs(g2) - dabs(g3)  - dabs(g4)  - dabs(g5) ) / f1
      !pb2 = dabs(g2) / f1
      !pb3 = dabs(g3) / f1
      !pb4 = dabs(g4) / f1
      !pb5 = dabs(g5) / f1
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !if(pb1 * pb2 * pb3 * pb4 * pb5 < zero)then
           !write(*, *)'mms55==', pb1, pb2, pb3, pb4, pb5
           !pb1 = ( f1 - dabs(g2) - dabs(g4) ) / f1
           !pb2 = dabs(g2) / f1
           !pb3 = dabs(g4) / f1
           !write(*, *)'mms55==', pb1, pb2, pb3 
           if(pb1 * pb2 * pb3 < zero)then 
               write(*, *)'mms55==', pb1, pb2, pb3, pb1 + pb2 + pb3 
               stop
           endif
      !endif
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~Now we sample the phi~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      scat_phi = Sampling_Phi( f1, g2, g3, g4, g5 )
      !********************************************************************************** 
      !write(*, *)'mms66==',  scat_phi, xi1
      !**********************************************************************************
      sinDphi =  dsin(this%phi_ini - scat_phi)
      cosDphi =  dcos(this%phi_ini - scat_phi)
      sin2Dphi =  dsin(two*(this%phi_ini - scat_phi))
      cos2Dphi =  dcos(two*(this%phi_ini - scat_phi))  
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
      N_temp = f1 + g2 * dcos(scat_phi) + g3 * dsin(scat_phi) + &
                    g4 * dcos(two * scat_phi) + g5 * dsin(two * scat_phi) 
      !N_temp2 = f1 + f2 * CosDphi + f3 * SinDphi + f4 * Cos2Dphi + f5 * Sin2Dphi

      !write(*, *)'ppp==', N_temp, N_temp2, PT(1, 1) + PT(1, 2) + PT(2, 1) + PT(2, 2) + &
      !                      ( PT(1, 1) - PT(1, 2) + PT(2, 1) - PT(2, 2) ) * QlI + &
      !                      two * ( PT(1, 3) + PT(2, 3) ) * UlI 

      if( N_temp < zero )then
           write(*, *)'mms55==',  UlI, QlI !, mu, mup, scat_phi
           write(*, *)'mmsa77==', pb1, pb2, pb3, pb4, pb5, pb1 + pb2 + pb3 + pb4 + pb5, f4, f5
           write(*, *)'mmsa88==', f1, f2, f3, f4, f5, UlI, QlI
           stop
      endif
  
      !this%Psi_Q = ( ( one -mu2 )*( one - mup2 )*(two + 3.D0*QlI) - &
      !                ( one + mup2 )*( one - mu2 ) + &
      !                four * mup * mu * smu * smup * (one + QlI) * CosDphi + &
      !                ( one + mu2 )*( mup2 - one + (mup2 + one)*QlI ) * Cos2Dphi + &
      !                four * mu * smu * smup * UlI * SinDphi + &
      !                two*( one + mu2 ) * mup * UlI * Sin2Dphi &
      !              ) / N_temp * this%Psi_I 
           !write(*, *)'tt1==',  this%Psi_Q, N_temp * this%Psi_I !  
      f_Q =(         PT(1, 1) + PT(1, 2) - PT(2, 1) - PT(2, 2) + &
             QlI * ( PT(1, 1) - PT(1, 2) - PT(2, 1) + PT(2, 2) )  + &
             UlI * ( PT(1, 3) - PT(2, 3) ) * two ) / N_temp 
       !write(*, *)'tt1==',  f_Q, this%Psi_Q / this%Psi_I
      this%Psi_Q = f_Q * this%Psi_I
           !write(*, *)'tt2==', this%Psi_Q,N_temp * this%Psi_I !  
      
      !this%Psi_U = ( - two * mup * smu * smup * SinDphi + &
      !               UlI * smu * smup * CosDphi + &
      !               (QlI - mup**2) * mu * Sin2Dphi + &
      !               UlI * mup * mu * Cos2Dphi  ) * &
      !               four * this%Psi_I / N_temp 
      !write(*, *)'tt3==', this%Psi_U / this%Psi_I , this%Psi_U
      !f_U =  ( PT(3, 1) + QlI * PT(3, 2) + UlI * PT(3, 3) )  / N_temp 
      f_U =  ( PT(3, 1) + PT(3, 2) + QlI * ( PT(3, 1) - PT(3, 2) ) + UlI * PT(3, 3)*two )  / N_temp 
      !write(*, *)'tt2==',  f_U, this%Psi_U / this%Psi_I
      this%Psi_U = f_U * this%Psi_I 
       !write(*, *)'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
      !write(*, *)'tt4==', PT(3, 1) + PT(3, 2)!this%Psi_U / this%Psi_I , this%Psi_U

      f_V = ( mup * mu + smu * smup * cosDphi ) / N_temp * four
      this%Psi_V = f_V * this%Psi_V
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !if( dabs(f_Q) > one .or. dabs(f_Q) > one)then
          !write(*, *)'mms11==', QlI, UlI, this%Psi_V, this%Psi_I, dsqrt(QlI**2 + UlI**2 + &
          !     ( this%Psi_V / this%Psi_I )**2 )
          !write(*, *)'mms22==',  f_Q, f_U, f_V
         ! write(*, *)'mms33==',  dsqrt(f_Q**2 + f_U**2 + ( this%Psi_V / this%Psi_I )**2)
      !endif
      this%phi_ini = scat_phi
      this%cosphi_ini = dcos(scat_phi)
      this%sinphi_ini = dsin(scat_phi)
      this%cos2phi_ini = dcos(two*scat_phi)
      this%sin2phi_ini = dsin(two*scat_phi)
      this%Vector_of_Momentum_ini(3) = mu
 
      !*************************************************************************************************  
      end Subroutine Tompson_Scat_With_Polarized_Diffuse_Reflection_Sub
!*******************************************************************************************************


!*******************************************************************************************************
      Subroutine Tompson_Scat_With_Polarized_Diffuse_Reflection2_Sub( this, mu_i )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton) :: this
      integer, intent(in) :: mu_i
      !real(mcp), intent(in) :: mu
      real(mcp) :: r, phip, smu, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp, beta, &
                   A_const, B_const, C_const, phi_ini, N_normal, N_temp2, &
                   cosphi_ini, sinphi_ini, cos2phi_ini, sin2phi_ini, smup, &
                   mu2, mup2, mup, mu

      real(mcp) :: A1, B1, t1, t2, Norm_C, f1, f2, f3, f4, f5, F0, g1, g2, g3, g4, g5
      real(mcp) :: UlI, QlI, scat_phi, sinDphi, cosDphi, sin2Dphi, cos2Dphi
      real(mcp) :: pb1, pb2, pb3, pb4, pb5, xi1
      real(mcp) :: P0(1: 3, 1: 3) = zero, P1(1: 3, 1: 3) = zero, P2(1: 3, 1: 3) = zero, &
                   PT(1: 3, 1: 3) = zero
      real(mcp) :: f_Q, f_U, f_V
      !****************************************************************************
      real(mcp) :: a, b, c, d
      Complex*16 :: roots3(1: 3), rts(1: 3)
      integer :: del, cases_f
      !**************************************************************************** 
   
      mu = this%mu_estimat( mu_i )
      mup = this%Vector_of_Momentum_ini(3)
      mup2 = mup**2
      smup = dsqrt(one - mup2)
      QlI = this%Psi_Q / this%Psi_I
      UlI = this%Psi_U / this%Psi_I
      t1 = one + mup2 - (one - mup2) * QlI
      t2 = two * (one - mup2) * (one + QlI)
      A1 = t1 - t2
      B1 = t1 + t2
      a = A1
      b = zero
      c = 3.D0 * B1
      d = a + c - 16.D0 * ranmar()
      
      mu2 = mu**2 
 
      if ( dabs(mu) < one ) then
          smu = dsqrt( one - mu2 )
      else
          smu = zero
      endif  

      !Norm_C = twopi * (mu2 * A1 + B1)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      f1 = (mu2 * A1 + B1) !t1 * (one + mu2) + t2 * (one - mu2)
      f2 = four * mup * smup * mu * smu * (one + QlI)
      f3 = four * smup * mu * smu * UlI
      f4 = (one - mu2) * (one - mup2) - (one - mu2) * (one + mup2) * QlI
      f5 = - two * (one - mu2) * mup * UlI
      !F0 = f1 - dabs(f2) - dabs(f3) - dabs(f4) - dabs(f5)
      !p1 = F0 / Norm_C
      !p2 = dabs(f2) / Norm_C
      !p3 = dabs(f3) / Norm_C
      !p4 = dabs(f4) / Norm_C
      !p5 = dabs(f5) / Norm_C
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      g2 = f2 * this%cosphi_ini + f3 * this%sinphi_ini
      g3 = f2 * this%sinphi_ini - f3 * this%cosphi_ini
      g4 = f4 * this%cos2phi_ini + f5 * this%sin2phi_ini
      g5 = f4 * this%sin2phi_ini - f5 * this%cos2phi_ini
      !pb1 = ( f1 - dabs(g2) - dabs(g3)  - dabs(g4)  - dabs(g5) ) / f1
      !pb2 = dabs(g2) / f1
      !pb3 = dabs(g3) / f1
      !pb4 = dabs(g4) / f1
      !pb5 = dabs(g5) / f1
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !if(pb1 * pb2 * pb3 * pb4 * pb5 < zero)then
           !write(*, *)'mms55==', pb1, pb2, pb3, pb4, pb5
           !pb1 = ( f1 - dabs(g2) - dabs(g4) ) / f1
           !pb2 = dabs(g2) / f1
           !pb3 = dabs(g4) / f1
           !write(*, *)'mms55==', pb1, pb2, pb3 
           if(pb1 * pb2 * pb3 < zero)then 
               write(*, *)'mms55==', pb1, pb2, pb3, pb1 + pb2 + pb3 
               stop
           endif
      !endif
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !~~~~~~~~~~~~~~~~~~~~~Now we sample the phi~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      scat_phi = Sampling_Phi( f1, g2, g3, g4, g5 )
      !**********************************************************************************   
      this%phi_estimat( mu_i ) = scat_phi
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      !N_temp = f1 + g2 * dcos(scat_phi) + g3 * dsin(scat_phi) + &
      !              g4 * dcos(two * scat_phi) + g5 * dsin(two * scat_phi)   
      Call this%Get_Psi_IQUV_for_Estimation( mu, mu2, smu, mup, mup2, smup, &
                  scat_phi, f1, g2, g3, g4, g5, QlI, UlI )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
 
      !*************************************************************************************************  
      end Subroutine Tompson_Scat_With_Polarized_Diffuse_Reflection2_Sub
!*******************************************************************************************************


!*******************************************************************************************************
      Subroutine Get_Psi_IQUV_for_Estimation_Sub( this, mu, mu2, smu, mup, mup2, smup, &
                  scat_phi, f1, g2, g3, g4, g5, QlI, UlI )
!*******************************************************************************************************
      class(ScatPhoton) :: this 
      real(mcp), intent(in) ::  mu, mu2, smu, mup, mup2, smup, scat_phi, f1, g2, g3, g4, g5, QlI, UlI
      real(mcp) :: sinDphi, cosDphi, sin2Dphi, cos2Dphi
      real(mcp) :: P0(1: 3, 1: 3) = zero, P1(1: 3, 1: 3) = zero, P2(1: 3, 1: 3) = zero, &
                   PT(1: 3, 1: 3) = zero
      real(mcp) :: f_Q, f_U, f_V, N_temp, N_temp2, temp3, Norm_phi

      sinDphi =  dsin(this%phi_ini - scat_phi)
      cosDphi =  dcos(this%phi_ini - scat_phi)
      sin2Dphi =  dsin(two*(this%phi_ini - scat_phi))
      cos2Dphi =  dcos(two*(this%phi_ini - scat_phi)) 
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
      !N_temp = f1 + g2 * dcos(scat_phi) + g3 * dsin(scat_phi) + &
      !              g4 * dcos(two * scat_phi) + g5 * dsin(two * scat_phi) 
      !Norm_phi = N_temp / twopi / f1 
      !N_temp2 = f1 + f2 * CosDphi + f3 * SinDphi + f4 * Cos2Dphi + f5 * Sin2Dphi
      N_temp = PT(1, 1) + PT(1, 2) + PT(2, 1) + PT(2, 2) + &
                            ( PT(1, 1) - PT(1, 2) + PT(2, 1) - PT(2, 2) ) * QlI + &
                            two * ( PT(1, 3) + PT(2, 3) ) * UlI 
    
      !write(*, fmt="(' ', F22.15)") N_temp2- N_temp!  p2 / Norm_phi * 3.D0 / 32.D0 /pi
      this%f_IQUV_estimat(1) = 3.D0 * f1 / 16.D0 * this%Psi_I !N_temp2 * this%Psi_I / Norm_phi
      if( N_temp < zero )then
           write(*, *)'mms55==',  UlI, QlI !, mu, mup, scat_phi
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
      this%f_IQUV_estimat(2) = f_Q * this%Psi_I * 3.D0 * f1 / 16.D0 / N_temp ! / Norm_phi

      !this%Psi_U = ( - two * mup * smu * smup * SinDphi + &
      !               UlI * smu * smup * CosDphi + &
      !               (QlI - mup**2) * mu * Sin2Dphi + &
      !               UlI * mup * mu * Cos2Dphi  ) * &
      !               four * this%Psi_I / N_temp 
      !write(*, *)'tt3==', this%Psi_U / this%Psi_I , this%Psi_U

      f_U =  ( PT(3, 1) + PT(3, 2) + QlI * ( PT(3, 1) - PT(3, 2) ) + UlI * PT(3, 3)*two ) ! / N_temp  
      this%f_IQUV_estimat(3) = f_U * this%Psi_I * 3.D0 * f1 / 16.D0 / N_temp ! / Norm_phi

      f_V = ( mup * mu + smu * smup * cosDphi )! * three / 8.D0 / pi ! / N_temp * four
      this%f_IQUV_estimat(4) = f_V * this%Psi_V * 3.D0 * f1 / 16.D0 / N_temp !/ Norm_phi
 
      !write(*, *)'tt3==', this%f_IQUV_estimat

      end Subroutine Get_Psi_IQUV_for_Estimation_Sub
!*******************************************************************************************************

!*******************************************************************************************************
      real(mcp) function  Sampling_Phi( f1, g2, g3, g4, g5 )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: f1, g2, g3, g4, g5
      real(mcp) :: xi1, scat_phi, pb1, pb2, pb4, H_s, F_s

      pb1 = ( f1 - dabs(g2) - dabs(g4) ) / f1
      pb2 = dabs(g2) / f1
      pb4 = dabs(g4) / f1

      xi1 = ranmar()
      if( xi1 <= pb1 )then
          scat_phi = ranmar() * twopi
      else if( xi1 <= pb1 + pb2 )then 
          scat_phi = Cos_phi_2phi_Sampling( dsign(one, g2), one ) 
      else if( xi1 <= pb1 + pb2 + pb4 )then 
          scat_phi = Cos_phi_2phi_Sampling( dsign(one, g4), two )  
      endif 
 
      H_s =  g3 * dsin(scat_phi) + g5 * dsin(two*scat_phi) 
      F_s = f1 + g2 * dcos(scat_phi) + g4 * dcos(two*scat_phi)
      if( ranmar() - one <= H_s / F_s )then
          Sampling_Phi = scat_phi
      else
          Sampling_Phi = twopi - scat_phi
      endif 
      return
      end function Sampling_Phi


!*******************************************************************************************************
      real(mcp) function Cos_phi_2phi_Sampling( B, rank )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: B, rank
      real(mcp) :: r, A, m, Phi1 

      A = one / twopi
      !this%random_number_phi = ranmar()
      Phi1 = ranmar() * twopi 
      If ( ranmar() * two <= one + B * DCOS( rank*Phi1 ) ) then
          Cos_phi_2phi_Sampling = Phi1
      Else
          m = floor( Phi1 / ( pi/rank ) ) 
          Cos_phi_2phi_Sampling = (two * m + one) * pi / rank - Phi1
          !write(*, *)'kk=' , m, A - B * DCOS( two*Phi1 ), A * two
      Endif
      return
      end function Cos_phi_2phi_Sampling

!*******************************************************************************************************
      real(mcp) function Sin_phi_2phi_Sampling( B, rank )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: B, rank
      real(mcp) :: r, A, m, Phi1, Del

      A = one / twopi
      Del = pi/rank
      !this%random_number_phi = ranmar()
      Phi1 = ranmar() * twopi 
      If ( ranmar() * two <= one + B * DSIN( rank*Phi1 ) ) then
          Sin_phi_2phi_Sampling = Phi1
      Else
          If( Phi1 <= Del/two .or. Phi1 >= twopi - Del/two )then
              Sin_phi_2phi_Sampling = twopi - Phi1
          else
              m = floor( (Phi1 - Del/two) / ( Del ) ) 
              Sin_phi_2phi_Sampling = (two * m + one) * two * Del - Phi1 
          endif
      Endif
      return
      end function Sin_phi_2phi_Sampling


!*******************************************************************************************************
      subroutine Calc_Phot_Inform_At_Observor_DiffRefl_phi_Estimats_Sub(this)
!*******************************************************************************************************
      implicit none
      class(ScatPhoton) :: this   
      integer :: i, j, k, i_obs, i_phi, case_mu
      real(mcp) :: i_mu, Optical_Depth_scatter
      integer, parameter :: n = 128
      real(mcp) :: xx(0: n-1), ww(0: n-1), temp_m, vLv

      !this%d_theta = one / Num_PolDeg
      !this%d_phi =  twopi / Num_Phi 
      !write(*, *)'ffs=', one / Num_PolDeg, twopi / Num_Phi, this%d_theta, this%d_phi
 

          !i_mu = floor( dabs(this%Vector_of_Momentum_ini(3)) / this%d_theta ) 
 

          i_phi = floor( this%phi_estimat(1) / this%ds_phi ) 
    
          vLv = dexp( - Optical_Depth_scatter ) / this%mu_estimat(1)
          !this%Optical_Depth_scatter = this%z_tau / this%mu_estimat(1)!p_out1 

          !write(*, *)'tt3==', this%f_IQUV_estimat, vLv, this%Optical_Depth_scatter 
          this%PolarArrayIQUVs10(1: 4) = this%PolarArrayIQUVs10(1: 4) + &
                                        this%f_IQUV_estimat(1: 4) * vLv  
          this%PolarArrayIQUVs30(1: 4) = this%PolarArrayIQUVs30(1: 4) + &
                                        this%f_IQUV_estimat(1: 4) * vLv  
          this%PolarArrayIQUVs60(1: 4) = this%PolarArrayIQUVs60(1: 4) + &
                                        this%f_IQUV_estimat(1: 4) * vLv   
          this%PolarArrayIQUVs80(1: 4) = this%PolarArrayIQUVs80(1: 4) + &
                                        this%f_IQUV_estimat(1: 4) * vLv  
          !write(*, *)'tt3==', this%PolarArrayIQUV10(1: 4, i_phi ), i_phi
      return
      end subroutine Calc_Phot_Inform_At_Observor_DiffRefl_phi_Estimats_Sub

      end module ScatterPhoton






