!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE SemiAnalyMethod
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      USE Photons 
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Semi_Analytical_Calculations1( phot1 )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      TYPE(Photon), INTENT(INOUT) :: phot1
      REAL(mcp) :: j_nu, D, L, R_sphere, r_xy, r_max, dr, dtheta, theta_xy
      REAL(mcp) :: dnu, nu_low, nu_up, v_Lv(0: 500)
      REAL(mcp) :: alpha, nu, costheta, sintheta
      integer(kind=8) :: i, j, k, Nr, Nt, Nn

      v_Lv = zero
      R_sphere = one
      Nr = 300
      Nt = 300
      Nn = 500
      D = one
      L = 1.D4
      r_max = D * R_sphere / dsqrt( L**2 - R_sphere**2 )
      dr = r_max / Nr
      dtheta = twopi / Nt
      nu_low = 8.D0
      nu_up = 15.D0
      dnu = ( nu_up - nu_low ) / Nn
      DO i = 0, Nr
          r_xy = dr * i
          write(*, *)'sss====', i, r_xy
          DO j = 0, Nt
              Do k = 0, Nn
                  nu = 10**( k * dnu + nu_low )
                  theta_xy = dtheta * j
                  costheta = r_xy * dsin( theta_xy ) / dsqrt( r_xy**2 + D**2 )
                  sintheta = dsqrt( one - costheta**2 ) 
                  v_Lv( k ) = v_Lv( k ) + nu * phot1%j_theta_nu_emissity&
                              ( nu, costheta, sintheta ) * two * &
                              dsqrt( R_sphere**2 - ( L*r_xy )**2 / ( D**2 + r_xy**2 ) )
                  !write(*, *)'sss====', i, j, j_nu, nu, costheta, sintheta
              enddo
          ENDDO
      ENDDO
      open(unit=16, file='./image/v_Lv1.txt', status="replace") 
      do k = 0,Nn
          write(unit = 16, fmt = *)v_Lv( k )
      enddo 
      RETURN
      END SUBROUTINE Semi_Analytical_Calculations1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Semi_Analytical_Calculations2( phot )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
      !TYPE(Photon_Emitter), INTENT(INOUT) :: Emitter
      TYPE(Photon), INTENT(INOUT) :: Phot 
      REAL(mcp) :: j_nu, D, L, R_sphere, r_xy, r_max, dr, dtheta, theta_xy
      REAL(mcp) :: dnu, nu_low, nu_up, v_Lv(0: 500), Length_0
      REAL(mcp) :: alpha, nu, costhetaB, sinthetaB, alpha_nu 
      integer :: i, j, k, Nr, Nt, Nn

      v_Lv = zero
      R_sphere = phot%R_out
      Nr = 300
      Nt = 300
      Nn = 500
      D = one
      L = 1.D10 
      r_max = D * R_sphere / dsqrt( L**2 - R_sphere**2 )
      dr = r_max / Nr
      dtheta = twopi / Nt 
      dnu = ( phot%ln_nu2 - phot%ln_nu1 ) / Nn
      DO i = 0, Nr
          r_xy = dr * i
          write(*, *)'sss====', i, r_xy
          DO j = 0, Nt
              Do k = 0, Nn
                  nu = 10**( k * dnu + phot%ln_nu1 )
                  theta_xy = dtheta * j
                  costhetaB = ( phot%cos_theta_obs + r_xy * dsin( theta_xy ) * &
                                phot%sin_theta_obs ) / dsqrt( r_xy**2 + D**2 )
                  sinthetaB = dsqrt( one - costhetaB**2 ) 
                  Length_0 = two * dsqrt( R_sphere**2 - ( L*r_xy )**2 / ( D**2 + r_xy**2 ) ) 
                  alpha_nu = phot%j_theta_nu_emissity( nu, costhetaB, sinthetaB ) / &
                             ( two*planck_h*nu**3/Cv**2 / ( dexp( h_ev*nu*1.D-6/phot%T_e ) - one ) ) 
                  if( Length_0 * alpha_nu <= 1.D-10 )then
                      v_Lv( k ) = v_Lv( k ) + nu * phot%j_theta_nu_emissity &
                              ( nu, costhetaB, sinthetaB ) * Length_0  
                  else
                      v_Lv( k ) = v_Lv( k ) + nu * phot%j_theta_nu_emissity &
                              ( nu, costhetaB, sinthetaB ) / alpha_nu * &
                              ( one - dexp( - Length_0 * alpha_nu ) )  
                  endif 
                  !write(*, *)'sss====', phot%n_e * phot%sigma_fn( phot%T_e, nu*h_ev*1.D-6 ),alpha_nu, &
                  !       alpha_nu - phot%n_e * phot%sigma_fn( phot%T_e, nu*h_ev*1.D-6 )
                  !       phot%j_theta_nu_emissity( one, nu, costheta, sintheta )
              enddo
          ENDDO
      ENDDO
      open(unit=16, file='./spectrum/vLv_Integral3.txt', status="replace") 
      do k = 0,Nn
          write(unit = 16, fmt = *)v_Lv( k )
      enddo 
      RETURN
      END SUBROUTINE Semi_Analytical_Calculations2
 
!**************************************************************************************
      END MODULE SemiAnalyMethod
!**************************************************************************************



