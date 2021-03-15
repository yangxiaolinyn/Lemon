      module ScatterPhoton
      use PhotonEmitter
      implicit none 

!*******************************************************************************************************
      type, public, extends(Photon_Emitter) :: ScatPhoton 
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
      contains 
!*******************************************************************************************************
      procedure, public :: Tompson_Scattering_WithOut_Polarization_IQ  =>   &
                            Tompson_Scattering_WithOut_Polarization_IQ_Sub 
      end type ScatPhoton
 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
      private :: Tompson_Scattering_WithOut_Polarization_IQ_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      contains 

!*******************************************************************************************************
      Subroutine Tompson_Scattering_WithOut_Polarization_IQ_Sub( this )
!*******************************************************************************************************
      implicit none
      class(ScatPhoton) :: this 
      real(mcp) :: epsi, r, phip, mupsi, mu_tilde_p, sinmu_tilde, sinpsi, &
                   sinmu_tilde_p, sin_tilphi_p, cos_tilphi_p, N_temp, beta, &
                   A_const, B_const, C_const, QlI, mu_ini, N_normal
      !real(mcp), dimension(1:4) :: Scattered_Phot4k_In_Elec
      !real(mcp), dimension(1:3, 1:3) :: Temp_Matrix_3X3, TTemp_Matrix_3X3, PTemp_Matrix_3X3
      !real(mcp), dimension(1:3) :: Temp_Matrix_1X3
      !real(mcp), dimension(1:3) :: Temp_Vector, f3_scat_e_tilde
      !real(mcp), dimension(1:4) :: f4_scat_e
      real(mcp) :: psi_p
      !****************************************************************************

      !****************************************************************************
      real(mcp) :: a, b, c 
      Complex*16 :: roots3(1: 3), rts(1: 3)
      integer :: del, cases_f
      !****************************************************************************

      !mutilde = this%Elec_Phot_mu_In_Elec_CF
      !Etilde  = this%Phot4k_CtrCF(1) 
      !epsi = two * dabs( Etilde ) / mec2
      mu_ini = this%Vector_of_Momentum_ini(3)
      QlI = this%Q_IQ / this%I_IQ
      A_const = three - mu_ini**2 + ( one - mu_ini**2 ) * QlI ! - mu_ini
      B_const = three * mu_ini**2 - ( one + three * ( one - mu_ini**2 ) * QlI )
      C_const = zero !- mu_ini  ! * dexp( - this%z_tau )
      N_normal = A_const * two + two / three * B_const
      
      !The cdf function is given by: mu^3 + a * mu^2 + b * mu + c = 0, and 
          a = C_const * three / two / B_const
          b = three * A_const / B_const
          c = b + one - a - ranmar() * three * N_normal / B_const
          if( dabs( B_const ) > 1.D-7 )then
            !call root3(a, b, c, roots3(1), roots3(2), roots3(3), del) 
            call root3new(B_const, a*B_const, b*B_const, c*B_const, roots3, del) 
            if( del == 3 )then
              mupsi = real( roots3(3) )
              if( dabs( mupsi ) > one )then
                  mupsi = real( roots3(2) )
                  if( dabs( mupsi ) > one )then
                      mupsi = real( roots3(1) )
                      if( dabs( mupsi ) > one )then
                          write(*, *)'mms11==', roots3, A_const, B_const, mu_ini, a, b, c
                          a = a * B_const
                          b = b * B_const
                          c = c * B_const
                          call root3new(B_const, a, b, c, rts, del) 
                          write(*, *)'mms11==', rts, B_const, a, b, c
                          write(*, *)'mms11==', (- b + dsqrt(b**2 - four*a*c)) / two / a, &
                                       (- b - dsqrt(b**2 - four*a*c)) / two / a
                          stop
                      endif
                  endif
              endif
            else
              mupsi = real( roots3(1) )
              if( dabs( mupsi ) > one .or. isnan( mupsi ) )then
                  write(*, *)'mms22==', A_const, B_const, mu_ini, a, b, c
                  write(*, *)'mms33==', del, roots3 
                          a = a * B_const
                          b = b * B_const
                          c = c * B_const
                          call root3new(B_const, a, b, c, rts, del) 
                          write(*, *)'mms11==', rts, B_const, a, b, c,  QlI, this%Q_IQ, this%I_IQ
                          write(*, *)'mms11==', (- b + dsqrt(b**2 - four*a*c)) / two / a, &
                                       (- b - dsqrt(b**2 - four*a*c)) / two / a
                  stop
              endif 
            endif
          else
              ! At this case B = 0, 
              !mupsi = ( - b + dsqrt(b**2 - four*a*c) ) / two / a
              mupsi = - c / b
              if( dabs( mupsi ) > one  )then
                  !write(*, *)'mms41==', A_const, B_const, mu_ini, a, b, c, mupsi
                  mupsi = ( - b - dsqrt(b**2 - four*a*c) ) / two / a
                  if( dabs( mupsi ) > one  )then
                      mupsi = - c / b
                  else
                      !write(*, *)'mms42==', A_const, B_const, mu_ini, a, b, c 
                      !write(*, *)'mms43==', (- b + dsqrt(b**2 - four*a*c)) / two / a, &
                      !                 (- b - dsqrt(b**2 - four*a*c)) / two / a, - c / b, this%Q_IQ, this%I_IQ
                      !stop
                  endif
              endif 
          endif
          !r = one + epsi * ( one - mupsi ) / two
   
      !this%r_one_hvlmec2_one_cosE = r
  
      !**********************************************************************
 
      if ( dabs(mupsi) < one ) then
          sinpsi = dsqrt( one - mupsi**2 )
      else
          sinpsi = zero
      endif  
      !********************************************************************************** 
      N_temp = ( A_const + B_const * mupsi**2 + C_const * mupsi ) / N_normal
      if( N_temp < zero )then
          write(*, *)'mms55==',  this%I_IQ, this%Q_IQ, N_temp, mupsi, mu_ini, QlI
          stop
      endif
      this%Q_IQ = this%I_IQ * ( ( one - mupsi**2 ) * ( one - three * mu_ini**2 ) + QlI * three * &
                    ( one - mupsi**2 ) * ( one - mu_ini**2 ) ) / N_temp * three / 16.D0
      !this%I_IQ = this%I_IQ * N_normal * three / 16.D0

      !write(*, *)'mms44==',this%I_IQ,   N_normal * three / 16.D0
      !this%I_IQ = this%I_IQ * ( 3.D0 - mupsi**2 - mu_ini**2 + QlI * &
      !              ( one - three * mupsi**2 ) * ( one - mu_ini**2 ) ) / N_temp * three / 16.D0 
      !this%I_IQ = this%I_IQ * N_normal * three / 16.D0
      !this%N_scat = N_normal
      ! write(*, *)'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
      ! write(*, *)'mms44==',  this%I_IQ, this%Q_IQ
      !stop
      !if( QlI > one )then
          !write(*, *)'mms44==',  this%I_IQ, this%Q_IQ, N_temp, mupsi, mu_ini, QlI
         ! stop
      !endif 
      !**********************************************************************
      !*** To obtain the Scattered 4 momentum of photon in the CF 
      !********************************************************************** 
      
      this%Scattered_Phot4k_CF(1) = this%Phot4k_CtrCF(1) 
      this%Scattered_Phot4k_CF(2:3) = zero
      this%Scattered_Phot4k_CF(4) = mupsi
      !this%Scattered_Phot4k_CovCF = this%Scattered_Phot4k_CF
      !this%Scattered_Phot4k_CovCF(1) = - this%Scattered_Phot4k_CF(1)

      !if( this%Scattered_Phot4k_CF(4) / dabs(this%Scattered_Phot4k_CF(1)) > 1.D0 )then
         !write(*, *)'mms0000==', this%Scattered_Phot4k_CF / dabs(this%Scattered_Phot4k_CF(1))
      !endif 
      !*************************************************************************************************  
      end Subroutine Tompson_Scattering_WithOut_Polarization_IQ_Sub 
!*******************************************************************************************************

      end module ScatterPhoton






