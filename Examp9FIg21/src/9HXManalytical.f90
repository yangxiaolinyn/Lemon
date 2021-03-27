      module HXManalyticalresults
      !use Basic_Variables_And_Methods
      use PhotonEmitterBB
      implicit none
      
      type, public, extends(Photon_Emitter_BB) :: Hxmformulae
      contains 
      procedure, public :: Hxm_alaytical_results  =>  Hxm_alaytical_results_Sub
      procedure, public :: Sigma_alpha => Sigma_alpha_Sub

      end type Hxmformulae
!******************************************************************************************************* 
      private :: Hxm_alaytical_results_Sub
      private :: Sigma_alpha_Sub
!*******************************************************************************************************

      contains
!*******************************************************************************************************
      subroutine Hxm_alaytical_results_Sub( this, E0, E, mu2, g1Emu2 )
!*******************************************************************************************************
      implicit none
      class(Hxmformulae) :: this
      real(mcp), intent(in) :: E0, E, mu2
      real(mcp), intent(out) :: g1Emu2
      real(mcp) :: f_r, Imumu2, b, c, d
      real(mcp) :: mu, sig1, sig2, alpha, r, f1chi0, f2chi0, coschi0, sinchi0, &
                   sinmu, sinmu2, alpha0, cass
  
      call this%Sigma_alpha( E0, sig1 )
      call this%Sigma_alpha( E, sig2 )

      r = E / E0
      alpha = E / mec2
      alpha0 = E0 / mec2
      mu = ( alpha0 * r + r - one ) / alpha0 / r
      b = one - sig2 / sig1 * mu
      sinmu = dsqrt(one - mu**2)
      sinmu2 = dsqrt(one - mu2**2)
      c = - sig2 / sig1 / mu2 * sinmu * sinmu2

      if( one / ( two * alpha0 + one ) <= r .and. r <= 1 )then
          f_r = three * sigma_T / eight / sig1 / alpha0 * ( one / r + r + mu**2 - one ) 
      else
          f_r = zero
      endif
      !write(*, *)'fss= ', three, sigma_T, eight, sig1, alpha0, mu, r, (one / r + r + mu**2 - one)
 
      if( mu > - sinmu2 )then
          coschi0 = - mu * mu2 / sinmu / sinmu2
          sinchi0 = dsqrt( one - coschi0**2 )
      else
          coschi0 = one
          sinchi0 = zero
      endif

      if( c + one < b .and. b < -c )then
          d = dsqrt( c**2 - b**2 )
          f1chi0 = d * sinchi0 / ( c + b * coschi0 )
          Imumu2 = - one / d * dlog( (one + f1chi0) / (one-f1chi0) )
          cass = 1
      else if( b == -c )then 
          Imumu2 = two * sinchi0 / b / (one - coschi0)
          cass = 2
      else if( b > -c )then
          d = dsqrt( b**2 - c**2 )
          if( coschi0 == one )then
              Imumu2 = twopi / d
              cass = 3
          else if( coschi0 < one .and. coschi0 >= - c / b )then
              f2chi0 = d * sinchi0 / ( c + b * coschi0 )
              Imumu2 = two / d * ( pi - datan(f2chi0) )
              cass = 4
          else
              f2chi0 = d * sinchi0 / ( c + b * coschi0 )
              Imumu2 = - two / d * ( datan(f2chi0) )
              cass = 5
          endif
      endif
      !write(*, *)'fss= ', cass, coschi0, sinchi0, b, sig2, sig1, mu2 , sinmu , sinmu2
      g1Emu2 = f_r / E0 * Imumu2
      return
      end subroutine Hxm_alaytical_results_Sub


!*******************************************************************************************************
      subroutine Sigma_alpha_Sub( this, E, sig )
!*******************************************************************************************************
      implicit none
      class(Hxmformulae) :: this
      real(mcp), intent(in) :: E 
      real(mcp), intent(out) :: sig 
      real(mcp) :: alpha
 
      alpha = E / mec2
      sig = three * sigma_T / eight / alpha * ( ( one - two *(one + alpha)/alpha**2 )*&
             dlog(two*alpha + one) + 0.5D0 + four / alpha - one /two/(two*alpha + one)**2 ) 
      return
      end subroutine Sigma_alpha_Sub

  
!*******************************************************************************************************

      end module HXManalyticalresults







