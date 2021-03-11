      module Photons 
      use PhotonEmitter
      !use ScatDistance
      implicit none 
      integer, parameter :: N_sigma = 1000 

      type, public, extends(Basic_Variables_And_Methods_Of_Particle) :: Photon
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!* The BL coordinates of the photon at p, which determines the BL coordinates  *
!* by YNOGK functions: r(p), mucos(p), phi(p), t(p), sigma(p)                  *
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
          real(mcp) :: v_L_v_i(0:500) = zero
          real(mcp) :: v_L_v_i_ET(1:500, 1:1000)
          real(mcp) :: frequency_v
          real(mcp) :: Optical_Depth_scatter
          real(mcp) :: Optical_Depth_absorption
          real(mcp) :: j_Enu_theta
          real(mcp) :: R_out
          real(mcp) :: logE_low
          real(mcp) :: logE_up
          real(mcp) :: dindexE 
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_85
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_230
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_50
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE_100
          real(mcp), dimension(0:N_sigma) :: sigmaaTeE
          real(mcp) :: NormalA
          real(mcp) :: n_e 
          integer(kind=8) :: effect_number
          integer(kind=8) :: scatter_times

      contains 
!******************************************************************************************************* 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      Subroutines for non-Kerr space-time
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
       procedure, public :: Calc_Phot_Informations_At_Observor_2zones => &
                            Calc_Phot_Informations_At_Observor_2zones_Sub 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end type Photon
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      private ::  Calc_Phot_Informations_At_Observor_2zones_Sub 
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
      contains   
 
 
!*******************************************************************************************************
      subroutine Calc_Phot_Informations_At_Observor_2zones_Sub( this )
!*******************************************************************************************************
      implicit none
      class(Photon) :: this
      !integer, intent(in) :: cases   
      real(mcp) :: index_i1, del_x1, r_times_p, p_out1
      integer :: k, N_low1
  
      del_x1 = (this%ln_nu2 - this%ln_nu1) / 500.D0
      !N_low1 = - floor( this%ln_nu1 / del_x1 )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%frequency_v = DABS( this%Phot4k_CovCF_ini(1) ) * 1.D6 / h_ev
      index_i1 = dlog10( this%frequency_v / 1.D8 ) 
      k = floor( index_i1 / del_x1 ) + 1 ! + N_low1 + 1 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      r_times_p = Vector3D_Inner_Product( this%Vector_of_Momentum_ini, &
                                          this%Vector_of_position_ini )
      p_out1 = - r_times_p + dsqrt( r_times_p**2 - ( this%r_ini**2 - this%R_out**2 ) )  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
      this%v_L_v_i(k) = this%v_L_v_i(k) + this%w_ini * dexp( - this%Optical_Depth_absorption * p_out1 ) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
     
      return
      end subroutine Calc_Phot_Informations_At_Observor_2zones_Sub
!*******************************************************************************************************
!******************************************************************************************************* 
!*******************************************************************************************************
!*******************************************************************************************************
!*******************************************************************************************************
 
      end module Photons





