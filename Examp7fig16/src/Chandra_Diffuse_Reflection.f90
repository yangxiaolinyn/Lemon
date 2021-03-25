      module ChandraDiffuseReflection
      use PhotonEmitter
      implicit none
      
      type, public, extends(Photon_Emitter) :: ChandDR
      contains
      procedure, public :: Diffuse_Reflection_Ir_Il_U_Setting    =>   &
                           Diffuse_Reflection_Ir_Il_U_Setting_Sub
      procedure, public :: Diffuse_Reflection_QS_Matrix_Setting  =>    &
                           Diffuse_Reflection_QS_Matrix_Setting_Sub
      procedure, public :: Diffuse_Reflection_QS_Matrix_of_IQUV_Setting  =>    &
                           Diffuse_Reflection_QS_Matrix_of_IQUV_Setting_Sub

      end type ChandDR
!*******************************************************************************************************

      private :: Diffuse_Reflection_Ir_Il_U_Setting_Sub
      private :: Diffuse_Reflection_QS_Matrix_Setting_Sub
      private :: Diffuse_Reflection_QS_Matrix_of_IQUV_Setting_Sub
!*******************************************************************************************************

      contains
!*******************************************************************************************************
      subroutine Diffuse_Reflection_Ir_Il_U_Setting_Sub( this, mu, mu0, phi, phi0, I_r, I_l, U, V )
!*******************************************************************************************************
      implicit none
      class(ChandDR) :: this
      real(mcp), intent(in) :: mu, mu0, phi, phi0
      real(mcp), intent(out) :: I_r, I_l, U, V
      real(mcp) :: psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu, H_nu_mu, H_r_mu, H_l_mu
      real(mcp) :: psi_mu0, phi_mu0, chi_mu0, zeta_mu0, H1_mu0, H2_mu0, H_nu_mu0, H_r_mu0, H_l_mu0
      real(mcp) :: smu, smu0, cosDphi, sinDphi, cos2Dphi, sin2Dphi

      smu  = dsqrt( one - mu**2 )
      smu0 = dsqrt( one - mu0**2 )
      cosDphi = dcos( phi0 - phi )
      sinDphi = dsin( phi0 - phi )
      cos2Dphi = dcos( two*(phi0 - phi) )
      sin2Dphi = dsin( two*(phi0 - phi) )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, chi_mu, zeta_mu, &
                            H1_mu, H2_mu, H_nu_mu, H_r_mu, H_l_mu )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu0, psi_mu0, phi_mu0, chi_mu0, zeta_mu0, &
                            H1_mu0, H2_mu0, H_nu_mu0, H_r_mu0, H_l_mu0 )
  
      I_l = three / ( mu + mu0 ) / 32.D0 * ( psi_mu * ( psi_mu0 + chi_mu0 ) + two * phi_mu * &
         ( phi_mu0 + zeta_mu0 ) - mu**2 * (one - mu0**2) * H2_mu * H2_mu0 * cos2Dphi - &
      four * mu * mu0 * smu * smu0 * H1_mu * H1_mu0 * cosDphi )

      I_r = three / ( mu + mu0 ) / 32.D0 * ( chi_mu * ( psi_mu0 + chi_mu0 ) + two * zeta_mu * &
               ( phi_mu0 + zeta_mu0 ) + (one - mu0**2) * H2_mu * H2_mu0 * cos2Dphi )
 
      U = three / ( mu + mu0 ) / 16.D0 * ( two * smu *smu0 * &
                    mu0 * H1_mu * H1_mu0 * sinDphi + mu * ( one - mu0**2 ) * &
                    H2_mu * H2_mu0 * sin2Dphi )
    
      V = three / ( mu + mu0 ) / 8.D0 * ( - mu * mu0 * H_nu_mu * H_nu_mu0 + &
                   smu * smu0 * H_r_mu * H_r_mu0 * cosDphi ) 

      end subroutine Diffuse_Reflection_Ir_Il_U_Setting_Sub



!*******************************************************************************************************
      real(mcp) Function psi_mu( mu )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: mu
      integer :: i

      i = floor( mu / 0.05D0 )

      psi_mu = ( psi_array(i + 1) - psi_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + psi_array( i )
      end Function psi_mu


!*******************************************************************************************************
      Subroutine Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, &
                                      chi_mu, zeta_mu, H1_mu, H2_mu, H_nu_mu, H_r_mu, H_l_mu )
!*******************************************************************************************************
      implicit none
      real(mcp), intent(in) :: mu
      real(mcp), intent(out) :: psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu, &
                         H_nu_mu, H_r_mu, H_l_mu
      integer :: i

      if( mu == one )then
          psi_mu = psi_array( 20 )
          phi_mu = phi_array( 20 )
          chi_mu = chi_array( 20 )
          zeta_mu = zeta_array( 20 )
          H1_mu = H1_array( 20 )
          H2_mu = H2_array( 20 )
          H_nu_mu = H_nu_array(20)
          return
      endif

      i = floor( mu / 0.05D0 )

      psi_mu = ( psi_array(i + 1) - psi_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + psi_array( i )

      phi_mu = ( phi_array(i + 1) - phi_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + phi_array( i )

      chi_mu = ( chi_array(i + 1) - chi_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + chi_array( i )

      zeta_mu = ( zeta_array(i + 1) - zeta_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + zeta_array( i )

      H1_mu = ( H1_array(i + 1) - H1_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + H1_array( i )

      H2_mu = ( H2_array(i + 1) - H2_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + H2_array( i )

      H_nu_mu = ( H_nu_array(i + 1) - H_nu_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + H_nu_array( i )

      H_r_mu = ( H_r_array(i + 1) - H_r_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + H_r_array( i )

      H_l_mu = ( H_l_array(i + 1) - H_l_array(i) ) / 0.05D0 * &
               ( mu - i * 0.05D0 ) + H_l_array( i )

      end Subroutine Get_psi_Phi_chi_zeta_H1_H2_mu
      
!*******************************************************************************************************
      subroutine Diffuse_Reflection_QS_Matrix_Setting_Sub( this, mu, mu0, phi, phi0, QS_Matrix )
!*******************************************************************************************************
      implicit none
      class(ChandDR) :: this
      real(mcp), intent(in) :: mu, mu0, phi, phi0
      real(mcp), dimension(1:3, 1:3), intent(out) :: QS_Matrix
      real(mcp) :: psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu
      real(mcp) :: psi_mu0, phi_mu0, chi_mu0, zeta_mu0, H1_mu0, H2_mu0
      real(mcp) :: H_nu_mu, H_r_mu, H_l_mu, H_nu_mu0, H_r_mu0, H_l_mu0
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part1 = zero
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part2 = zero
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part3 = zero
      real(mcp), dimension(1:3, 1:3) :: Q_Matrix = zero
      real(mcp), dimension(1:3, 1:3) :: QS_Matrix1 = zero
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix1 = zero
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix2 = zero
      real(mcp) :: cosDphi, sinDphi, cos2Dphi, sin2Dphi


      !CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu )
      !CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu0, psi_mu0, phi_mu0, chi_mu0, zeta_mu0, H1_mu0, H2_mu0 )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, chi_mu, zeta_mu, &
                            H1_mu, H2_mu, H_nu_mu, H_r_mu, H_l_mu )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu0, psi_mu0, phi_mu0, chi_mu0, zeta_mu0, &
                            H1_mu0, H2_mu0, H_nu_mu0, H_r_mu0, H_l_mu0 )
  

      Temp_Matrix1(1, 1) = psi_mu
      Temp_Matrix1(1, 2) = phi_mu * sqrt2
      Temp_Matrix1(2, 1) = chi_mu
      Temp_Matrix1(2, 2) = zeta_mu * sqrt2

      Temp_Matrix2(1, 1) = psi_mu0
      Temp_Matrix2(1, 2) = chi_mu0
      Temp_Matrix2(2, 1) = phi_mu0 * sqrt2
      Temp_Matrix2(2, 2) = zeta_mu0 * sqrt2

      CALL Matrix_Multiplication33X33_Sub( Temp_Matrix1, Temp_Matrix2, S_Matrix_part1 )
      S_Matrix_part1 = S_Matrix_part1 * three / four

      cosDphi = dcos(phi0 - phi)
      sinDphi = dsin(phi0 - phi)
      cos2Dphi = dcos( two * (phi0 - phi) )
      sin2Dphi = dsin( two * (phi0 - phi) )

      S_Matrix_part2(1, 1) = - four * mu0 * mu * cosDphi
      S_Matrix_part2(1, 3) = two * mu * sinDphi
      S_Matrix_part2(3, 1) = two * mu0 * sinDphi
      S_Matrix_part2(3, 3) = cosDphi
      S_Matrix_part2 = S_Matrix_part2 * dsqrt( one - mu**2 ) * dsqrt( one - mu0**2 ) *&
                       H1_mu * H1_mu0 * three / four
 
      S_Matrix_part3(1, 1) = ( mu0 * mu )**2 * cos2Dphi
      S_Matrix_part3(1, 2) = - mu**2 * cos2Dphi
      S_Matrix_part3(1, 3) = - mu**2 * mu0 * sin2Dphi

      S_Matrix_part3(2, 1) = - mu0**2 * cos2Dphi
      S_Matrix_part3(2, 2) = cos2Dphi
      S_Matrix_part3(2, 3) = mu0 * sin2Dphi

      S_Matrix_part3(3, 1) = - mu * mu0**2 * sin2Dphi
      S_Matrix_part3(3, 2) = mu * sin2Dphi
      S_Matrix_part3(3, 3) = - mu * mu0 * cos2Dphi
      S_Matrix_part3 = S_Matrix_part3 * H2_mu * H2_mu0 * three / four

      QS_Matrix1 = ( S_Matrix_part1 + S_Matrix_part2 + S_Matrix_part3 ) / ( mu + mu0 )
      Q_Matrix(1, 1) = one
      Q_Matrix(2, 2) = one
      Q_Matrix(3, 3) = two
      CALL Matrix_Multiplication33X33_Sub( Q_Matrix, QS_Matrix1, QS_Matrix ) 
      QS_Matrix = QS_Matrix / four

      end subroutine Diffuse_Reflection_QS_Matrix_Setting_Sub
 
!*******************************************************************************************************
      subroutine Diffuse_Reflection_QS_Matrix_of_IQUV_Setting_Sub( this, mu, mu0, phi, phi0, QS_Matrix )
!*******************************************************************************************************
      implicit none
      class(ChandDR) :: this
      real(mcp), intent(in) :: mu, mu0, phi, phi0
      real(mcp), dimension(1:3, 1:3), intent(out) :: QS_Matrix
      real(mcp) :: psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu
      real(mcp) :: psi_mu0, phi_mu0, chi_mu0, zeta_mu0, H1_mu0, H2_mu0
      real(mcp) :: H_nu_mu, H_r_mu, H_l_mu, H_nu_mu0, H_r_mu0, H_l_mu0
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part1 = zero
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part2 = zero
      real(mcp), dimension(1:3, 1:3) :: S_Matrix_part3 = zero
      real(mcp), dimension(1:3, 1:3) :: Q_Matrix = zero
      real(mcp), dimension(1:3, 1:3) :: QS_Matrix1 = zero
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix1 = zero
      real(mcp), dimension(1:3, 1:3) :: Temp_Matrix2 = zero
      real(mcp) :: cosDphi, sinDphi, cos2Dphi, sin2Dphi


      !CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, chi_mu, zeta_mu, H1_mu, H2_mu )
      !CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu0, psi_mu0, phi_mu0, chi_mu0, zeta_mu0, H1_mu0, H2_mu0 )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu, psi_mu, phi_mu, chi_mu, zeta_mu, &
                            H1_mu, H2_mu, H_nu_mu, H_r_mu, H_l_mu )
      CALL Get_psi_Phi_chi_zeta_H1_H2_mu( mu0, psi_mu0, phi_mu0, chi_mu0, zeta_mu0, &
                            H1_mu0, H2_mu0, H_nu_mu0, H_r_mu0, H_l_mu0 )
  

      Temp_Matrix1(1, 1) = psi_mu
      Temp_Matrix1(1, 2) = phi_mu * sqrt2
      Temp_Matrix1(2, 1) = chi_mu
      Temp_Matrix1(2, 2) = zeta_mu * sqrt2

      Temp_Matrix2(1, 1) = psi_mu0
      Temp_Matrix2(1, 2) = chi_mu0
      Temp_Matrix2(2, 1) = phi_mu0 * sqrt2
      Temp_Matrix2(2, 2) = zeta_mu0 * sqrt2

      CALL Matrix_Multiplication33X33_Sub( Temp_Matrix1, Temp_Matrix2, S_Matrix_part1 )
      S_Matrix_part1 = S_Matrix_part1 * three / four

      cosDphi = dcos(phi0 - phi)
      sinDphi = dsin(phi0 - phi)
      cos2Dphi = dcos( two * (phi0 - phi) )
      sin2Dphi = dsin( two * (phi0 - phi) )
      S_Matrix_part2(1, 1) = - four * mu0 * mu * cosDphi
      S_Matrix_part2(1, 3) = two * mu * sinDphi
      S_Matrix_part2(3, 1) = two * mu0 * sinDphi
      S_Matrix_part2(3, 3) = cosDphi
      S_Matrix_part2 = S_Matrix_part2 * dsqrt( one - mu**2 ) * dsqrt( one - mu0**2 ) *&
                       H1_mu * H1_mu0 * three / four
 
      S_Matrix_part3(1, 1) = ( mu0 * mu )**2 * cos2Dphi
      S_Matrix_part3(1, 2) = - mu**2 * cos2Dphi
      S_Matrix_part3(1, 3) = - mu**2 * mu0 * sin2Dphi

      S_Matrix_part3(2, 1) = - mu0**2 * cos2Dphi
      S_Matrix_part3(2, 2) = cos2Dphi
      S_Matrix_part3(2, 3) = mu0 * sin2Dphi

      S_Matrix_part3(3, 1) = - mu * mu0**2 * sin2Dphi
      S_Matrix_part3(3, 2) = mu * sin2Dphi
      S_Matrix_part3(3, 3) = - mu * mu0 * cos2Dphi
      S_Matrix_part3 = S_Matrix_part3 * H2_mu * H2_mu0 * three / four

      QS_Matrix1 = ( S_Matrix_part1 + S_Matrix_part2 + S_Matrix_part3 ) * mu0 / ( mu + mu0 )
      QS_Matrix1(3, 1: 3) = QS_Matrix1(3, 1: 3) * two
      !Q_Matrix(1, 1) = one
      !Q_Matrix(2, 2) = one
      !Q_Matrix(3, 3) = two
      !CALL Matrix_Multiplication33X33_Sub( Q_Matrix, QS_Matrix1, QS_Matrix ) 
      QS_Matrix(1, 1) = ( QS_Matrix1(1, 1) + QS_Matrix1(1, 2) + QS_Matrix1(2, 1) + QS_Matrix1(2, 2) ) / two
      QS_Matrix(1, 2) = ( QS_Matrix1(1, 1) - QS_Matrix1(1, 2) + QS_Matrix1(2, 1) - QS_Matrix1(2, 2) ) / two
      QS_Matrix(1, 3) =   QS_Matrix1(1, 3) + QS_Matrix1(2, 3)

      QS_Matrix(2, 1) = ( QS_Matrix1(1, 1) + QS_Matrix1(1, 2) - QS_Matrix1(2, 1) - QS_Matrix1(2, 2) ) / two
      QS_Matrix(2, 2) = ( QS_Matrix1(1, 1) - QS_Matrix1(1, 2) - QS_Matrix1(2, 1) + QS_Matrix1(2, 2) ) / two
      QS_Matrix(2, 3) =   QS_Matrix1(1, 3) - QS_Matrix1(2, 3)

      QS_Matrix(3, 1) = ( QS_Matrix1(3, 1) + QS_Matrix1(3, 2) ) / two
      QS_Matrix(3, 2) = ( QS_Matrix1(3, 1) - QS_Matrix1(3, 2) ) / two
      QS_Matrix(3, 3) =   QS_Matrix1(3, 3)
      QS_Matrix = QS_Matrix / four 

      end subroutine Diffuse_Reflection_QS_Matrix_of_IQUV_Setting_Sub
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************

      end module ChandraDiffuseReflection







