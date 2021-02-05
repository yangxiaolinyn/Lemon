    module PhotoElectron 
    use Basic_Variables_And_Methods 
    implicit none 
    real(mcp), dimension(0: 13) :: photo_elect_arr_c0, photo_elect_arr_c1, photo_elect_arr_c2
    real(mcp), dimension(0: 48) :: E_1H_2He, COH_CS_1H, INCOH_CS_1H, COH_CS_2He, INCOH_CS_2He, &
                                   log10_E
    real(mcp) :: half_dy, COH_y1, COH_y2, COH_dy
 
    contains 

!*****************************************************************************************
    Subroutine Set_PhotoElect_CS_array()
!*****************************************************************************************
    implicit none
!*****************************************************************************************

      photo_elect_arr_c0(0) = 17.3D0
      photo_elect_arr_c0(1) = 34.6D0
      photo_elect_arr_c0(2) = 78.1D0
      photo_elect_arr_c0(3) = 71.4D0
      photo_elect_arr_c0(4) = 95.5D0
      photo_elect_arr_c0(5) = 308.9D0
      photo_elect_arr_c0(6) = 120.6D0
      photo_elect_arr_c0(7) = 141.3D0
      photo_elect_arr_c0(8) = 202.7D0
      photo_elect_arr_c0(9) = 342.7D0
      photo_elect_arr_c0(10) = 352.2D0
      photo_elect_arr_c0(11) = 433.9D0
      photo_elect_arr_c0(12) = 629.0D0
      photo_elect_arr_c0(13) = 701.2D0 

      photo_elect_arr_c1(0) = 608.1D0
      photo_elect_arr_c1(1) = 267.9D0
      photo_elect_arr_c1(2) = 18.8D0
      photo_elect_arr_c1(3) = 66.8D0
      photo_elect_arr_c1(4) = 145.8D0
      photo_elect_arr_c1(5) = -380.6D0
      photo_elect_arr_c1(6) =  169.3D0
      photo_elect_arr_c1(7) =  146.8D0
      photo_elect_arr_c1(8) =  104.7D0
      photo_elect_arr_c1(9) =  18.7D0
      photo_elect_arr_c1(10) =  18.7D0
      photo_elect_arr_c1(11) =  -2.4D0
      photo_elect_arr_c1(12) =  30.9D0
      photo_elect_arr_c1(13) =  25.2D0 

      photo_elect_arr_c2(0) = -2150.D0
      photo_elect_arr_c2(1) = -476.1D0
      photo_elect_arr_c2(2) = 4.3D0
      photo_elect_arr_c2(3) = -51.4D0
      photo_elect_arr_c2(4) = -61.1D0
      photo_elect_arr_c2(5) = 294.0D0
      photo_elect_arr_c2(6) = -47.7D0
      photo_elect_arr_c2(7) =  -31.5D0
      photo_elect_arr_c2(8) =  -17.0D0
      photo_elect_arr_c2(9) =  0.0D0
      photo_elect_arr_c2(10) = 0.0D0
      photo_elect_arr_c2(11) = 0.75D0
      photo_elect_arr_c2(12) = 0.0D0
      photo_elect_arr_c2(13) = 0.0D0 
 
 
    end Subroutine Set_PhotoElect_CS_array


!*****************************************************************************************
    real(mcp) function PhotoElect_CSect( E_Mev )
!*****************************************************************************************
    implicit none
    real(mcp), intent(in) :: E_Mev
    real(mcp) :: E
    E = E_Mev * 1.D3 ! trans it into keV
!*****************************************************************************************
 
    if(0.030D0 <= E .and. E < 0.100D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(0) + E * photo_elect_arr_c1(0) &
                            + E**2 * photo_elect_arr_c2(0) ) / E**3

    else if(0.100D0 <= E .and. E < 0.284D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(1) + E * photo_elect_arr_c1(1) &
                            + E**2 * photo_elect_arr_c2(1) ) / E**3

    else if(0.284D0 <= E .and. E < 0.400D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(2) + E * photo_elect_arr_c1(2) &
                            + E**2 * photo_elect_arr_c2(2) ) / E**3

    else if(0.400D0 <= E .and. E < 0.532D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(3) + E * photo_elect_arr_c1(3) &
                            + E**2 * photo_elect_arr_c2(3) ) / E**3

    else if(0.532D0 <= E .and. E < 0.707D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(4) + E * photo_elect_arr_c1(4) &
                            + E**2 * photo_elect_arr_c2(4) ) / E**3

    else if(0.707D0 <= E .and. E < 0.867D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(5) + E * photo_elect_arr_c1(5) &
                            + E**2 * photo_elect_arr_c2(5) ) / E**3

    else if(0.867D0 <= E .and. E < 1.303D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(6) + E * photo_elect_arr_c1(6) &
                            + E**2 * photo_elect_arr_c2(6) ) / E**3

    else if(1.303D0 <= E .and. E < 1.840D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(7) + E * photo_elect_arr_c1(7) &
                            + E**2 * photo_elect_arr_c2(7) ) / E**3

    else if(1.840D0 <= E .and. E < 2.471D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(8) + E * photo_elect_arr_c1(8) &
                            + E**2 * photo_elect_arr_c2(8) ) / E**3

    else if(2.471D0 <= E .and. E < 3.210D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(9) + E * photo_elect_arr_c1(9) &
                            + E**2 * photo_elect_arr_c2(9) ) / E**3

    else if(3.210D0 <= E .and. E < 4.038D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(10) + E * photo_elect_arr_c1(10) &
                            + E**2 * photo_elect_arr_c2(10) ) / E**3

    else if(4.038D0 <= E .and. E < 7.111D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(11) + E * photo_elect_arr_c1(11) &
                            + E**2 * photo_elect_arr_c2(11) ) / E**3

    else if(7.111D0 <= E .and. E < 8.331D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(12) + E * photo_elect_arr_c1(12) &
                            + E**2 * photo_elect_arr_c2(12) ) / E**3

    else if(8.331D0 <= E .and. E <= 10.000D0)then

         PhotoElect_CSect = ( photo_elect_arr_c0(13) + E * photo_elect_arr_c1(13) &
                            + E**2 * photo_elect_arr_c2(13) ) / E**3

    endif
    PhotoElect_CSect = PhotoElect_CSect * barn
!*****************************************************************************************
    end function PhotoElect_CSect
!*****************************************************************************************


!*****************************************************************************************
    Subroutine PhotoElect_1H_COH_INCOH_CSect( )
!*****************************************************************************************
    implicit none
!*****************************************************************************************
    E_1H_2He(0) = 1.0D+02
    E_1H_2He(1) = 1.5D+02
    E_1H_2He(2) = 2.0D+02
    E_1H_2He(3) = 3.0D+02
    E_1H_2He(4) = 4.0D+02
    E_1H_2He(5) = 5.0D+02
    E_1H_2He(6) = 6.0D+02
    E_1H_2He(7) = 8.0D+02

    E_1H_2He(8) = 1.0D+03
    E_1H_2He(9) = 1.5D+03
    E_1H_2He(10) = 2.0D+03
    E_1H_2He(11) = 3.0D+03
    E_1H_2He(12) = 4.0D+03
    E_1H_2He(13) = 5.0D+03
    E_1H_2He(14) = 6.0D+03
    E_1H_2He(15) = 8.0D+03

    E_1H_2He(16) = 1.0D+04
    E_1H_2He(17) = 1.5D+04
    E_1H_2He(18) = 2.0D+04
    E_1H_2He(19) = 3.0D+04
    E_1H_2He(20) = 4.0D+04
    E_1H_2He(21) = 5.0D+04
    E_1H_2He(22) = 6.0D+04
    E_1H_2He(23) = 8.0D+04
 
    E_1H_2He(24) = 1.0D+05
    E_1H_2He(25) = 1.5D+05
    E_1H_2He(26) = 2.0D+05
    E_1H_2He(27) = 3.0D+05
    E_1H_2He(28) = 4.0D+05
    E_1H_2He(29) = 5.0D+05
    E_1H_2He(30) = 6.0D+05
    E_1H_2He(31) = 8.0D+05

    E_1H_2He(32) = 1.0D+06
    E_1H_2He(33) = 1.5D+06
    E_1H_2He(34) = 2.0D+06
    E_1H_2He(35) = 3.0D+06
    E_1H_2He(36) = 4.0D+06
    E_1H_2He(37) = 5.0D+06
    E_1H_2He(38) = 6.0D+06
    E_1H_2He(39) = 8.0D+06

    E_1H_2He(40) = 1.0D+07
    E_1H_2He(41) = 1.5D+07
    E_1H_2He(42) = 2.0D+07
    E_1H_2He(43) = 3.0D+07
    E_1H_2He(44) = 4.0D+07
    E_1H_2He(45) = 5.0D+07
    E_1H_2He(46) = 6.0D+07
    E_1H_2He(47) = 8.0D+07

    E_1H_2He(48) = 1.0D+08
 

    COH_CS_1H(0) = 6.65D-1
    COH_CS_1H(1) = 6.635D-1
    COH_CS_1H(2) = 6.617D-1
    COH_CS_1H(3) = 6.569D-1
    COH_CS_1H(4) = 6.503D-1
    COH_CS_1H(5) = 6.421D-1
    COH_CS_1H(6) = 6.323D-1
    COH_CS_1H(7) = 6.087D-1

    COH_CS_1H(8) = 5.806D-1
    COH_CS_1H(9) = 4.984D-1
    COH_CS_1H(10) = 4.142D-1
    COH_CS_1H(11) = 2.764D-1
    COH_CS_1H(12) = 1.881D-1
    COH_CS_1H(13) = 1.341D-1
    COH_CS_1H(14) = 9.987D-2
    COH_CS_1H(15) = 6.126D-2

    COH_CS_1H(16) = 4.121D-2
    COH_CS_1H(17) = 1.943D-2
    COH_CS_1H(18) = 1.119D-2
    COH_CS_1H(19) = 5.062D-3
    COH_CS_1H(20) = 2.866D-3
    COH_CS_1H(21) = 1.840D-3
    COH_CS_1H(22) = 1.280D-3
    COH_CS_1H(23) = 7.211D-4
 
    COH_CS_1H(24) = 4.619D-04
    COH_CS_1H(25) = 2.054D-04
    COH_CS_1H(26) = 1.156D-04
    COH_CS_1H(27) = 5.138D-05
    COH_CS_1H(28) = 2.890D-05
    COH_CS_1H(29) = 1.850D-05
    COH_CS_1H(30) = 1.285D-05
    COH_CS_1H(31) = 7.226D-06
 

    COH_CS_1H(32) = 4.625D-06
    COH_CS_1H(33) = 2.056D-06
    COH_CS_1H(34) = 1.156D-06
    COH_CS_1H(35) = 5.139D-07
    COH_CS_1H(36) = 2.891D-07
    COH_CS_1H(37) = 1.850D-07
    COH_CS_1H(38) = 1.285D-07
    COH_CS_1H(39) = 7.227D-08
 
    COH_CS_1H(40) = 4.625D-08
    COH_CS_1H(41) = 2.056D-08
    COH_CS_1H(42) = 1.156D-08
    COH_CS_1H(43) = 5.139D-09
    COH_CS_1H(44) = 2.890D-09
    COH_CS_1H(45) = 1.850D-09
    COH_CS_1H(46) = 1.284D-09
    COH_CS_1H(47) = 7.222D-10

    COH_CS_1H(48) = 4.620D-10 
 
 

    INCOH_CS_1H(0) = 9.552D-04
    INCOH_CS_1H(1) = 2.144D-03
    INCOH_CS_1H(2) = 3.802D-03
    INCOH_CS_1H(3) = 8.494D-03
    INCOH_CS_1H(4) = 1.496D-02
    INCOH_CS_1H(5) = 2.310D-02
    INCOH_CS_1H(6) = 3.279D-02
    INCOH_CS_1H(7) = 5.629D-02 

    INCOH_CS_1H(8) = 8.424D-02
    INCOH_CS_1H(9) = 1.650D-01
    INCOH_CS_1H(10) = 2.478D-01
    INCOH_CS_1H(11) = 3.822D-01
    INCOH_CS_1H(12) = 4.675D-01
    INCOH_CS_1H(13) = 5.187D-1
    INCOH_CS_1H(14) = 5.503D-1
    INCOH_CS_1H(15) = 5.840D-1
 
    INCOH_CS_1H(16) = 5.993D-01
    INCOH_CS_1H(17) = 6.095D-01
    INCOH_CS_1H(18) = 6.068D-01
    INCOH_CS_1H(19) = 5.924D-01
    INCOH_CS_1H(20) = 5.759D-01
    INCOH_CS_1H(21) = 5.597D-01
    INCOH_CS_1H(22) = 5.444D-01
    INCOH_CS_1H(23) = 5.166D-01
  
    INCOH_CS_1H(24) = 4.923D-01
    INCOH_CS_1H(25) = 4.435D-01
    INCOH_CS_1H(26) = 4.064D-01
    INCOH_CS_1H(27) = 3.535D-01
    INCOH_CS_1H(28) = 3.168D-01
    INCOH_CS_1H(29) = 2.893D-01
    INCOH_CS_1H(30) = 2.676D-01
    INCOH_CS_1H(31) = 2.351D-01
 
    INCOH_CS_1H(32) = 2.114D-1
    INCOH_CS_1H(33) = 1.718D-1
    INCOH_CS_1H(34) = 1.466D-1
    INCOH_CS_1H(35) = 1.153D-1
    INCOH_CS_1H(36) = 9.620D-02
    INCOH_CS_1H(37) = 8.308D-02
    INCOH_CS_1H(38) = 7.343D-02
    INCOH_CS_1H(39) = 6.007D-02
 
    INCOH_CS_1H(40) = 5.116D-2
    INCOH_CS_1H(41) = 3.786D-2
    INCOH_CS_1H(42) = 3.039D-2
    INCOH_CS_1H(43) = 2.212D-2
    INCOH_CS_1H(44) = 1.758D-02
    INCOH_CS_1H(45) = 1.467D-02
    INCOH_CS_1H(46) = 1.264D-02
    INCOH_CS_1H(47) = 9.972D-03

    INCOH_CS_1H(48) = 8.276D-3 



    COH_CS_2He(0) = 2.660D0
    COH_CS_2He(1) = 2.658D0
    COH_CS_2He(2) = 2.655D0
    COH_CS_2He(3) = 2.648D0
    COH_CS_2He(4) = 2.637D0
    COH_CS_2He(5) = 2.624D0
    COH_CS_2He(6) = 2.608D0
    COH_CS_2He(7) = 2.567D0

    COH_CS_2He(8) = 2.517D0
    COH_CS_2He(9) = 2.356D0
    COH_CS_2He(10) = 2.162D0
    COH_CS_2He(11) = 1.743D0
    COH_CS_2He(12) = 1.369D0
    COH_CS_2He(13) = 1.072D0
    COH_CS_2He(14) = 8.492D-1
    COH_CS_2He(15) = 5.592D-1

    COH_CS_2He(16) = 3.921D-1
    COH_CS_2He(17) = 1.962D-1
    COH_CS_2He(18) = 1.166D-1
    COH_CS_2He(19) = 5.432D-2
    COH_CS_2He(20) = 3.113D-2
    COH_CS_2He(21) = 2.010D-2
    COH_CS_2He(22) = 1.403D-2
    COH_CS_2He(23) = 7.933D-3
 
    COH_CS_2He(24) = 5.089D-03
    COH_CS_2He(25) = 2.267D-03
    COH_CS_2He(26) = 1.276D-03
    COH_CS_2He(27) = 5.676D-04
    COH_CS_2He(28) = 3.194D-04
    COH_CS_2He(29) = 2.044D-04
    COH_CS_2He(30) = 1.420D-04
    COH_CS_2He(31) = 7.986D-05
 

    COH_CS_2He(32) = 5.111D-05
    COH_CS_2He(33) = 2.272D-05
    COH_CS_2He(34) = 1.278D-05
    COH_CS_2He(35) = 5.679D-06
    COH_CS_2He(36) = 3.194D-06
    COH_CS_2He(37) = 2.044D-06
    COH_CS_2He(38) = 1.420D-06
    COH_CS_2He(39) = 7.986D-07
 
    COH_CS_2He(40) = 5.111D-07
    COH_CS_2He(41) = 2.272D-07
    COH_CS_2He(42) = 1.278D-07
    COH_CS_2He(43) = 5.679D-8
    COH_CS_2He(44) = 3.194D-08
    COH_CS_2He(45) = 2.044D-08
    COH_CS_2He(46) = 1.420D-08
    COH_CS_2He(47) = 7.984D-9

    COH_CS_2He(48) = 5.109D-9
 


    INCOH_CS_2He(0) = 1.524D-03
    INCOH_CS_2He(1) = 2.559D-03
    INCOH_CS_2He(2) = 3.793D-03
    INCOH_CS_2He(3) = 7.081D-03
    INCOH_CS_2He(4) = 1.183D-02
    INCOH_CS_2He(5) = 1.796D-02
    INCOH_CS_2He(6) = 2.543D-02
    INCOH_CS_2He(7) = 4.426D-02

    INCOH_CS_2He(8) = 6.767D-02
    INCOH_CS_2He(9) = 1.419D-01
    INCOH_CS_2He(10) = 2.317D-01
    INCOH_CS_2He(11) = 4.230D-01
    INCOH_CS_2He(12) = 5.948D-01
    INCOH_CS_2He(13) = 7.328D-01
    INCOH_CS_2He(14) = 8.385D-01
    INCOH_CS_2He(15) = 9.785D-01
 

    INCOH_CS_2He(16) = 1.059D0
    INCOH_CS_2He(17) = 1.145D+00
    INCOH_CS_2He(18) = 1.168D+00
    INCOH_CS_2He(19) = 1.163D+00
    INCOH_CS_2He(20) = 1.139D+0
    INCOH_CS_2He(21) = 1.111D+00
    INCOH_CS_2He(22) = 1.083D+00
    INCOH_CS_2He(23) = 1.030D+00
   
    INCOH_CS_2He(24) = 9.825D-01
    INCOH_CS_2He(25) = 8.860D-01
    INCOH_CS_2He(26) = 8.124D-01
    INCOH_CS_2He(27) = 7.068D-01
    INCOH_CS_2He(28) = 6.334D-01
    INCOH_CS_2He(29) = 5.785D-01
    INCOH_CS_2He(30) = 5.352D-01
    INCOH_CS_2He(31) = 4.702D-01
 


    INCOH_CS_2He(32) = 4.228D-01
    INCOH_CS_2He(33) = 3.436D-1
    INCOH_CS_2He(34) = 2.932D-1
    INCOH_CS_2He(35) = 2.307D-1
    INCOH_CS_2He(36) = 1.924D-01
    INCOH_CS_2He(37) = 1.662D-01
    INCOH_CS_2He(38) = 1.469D-01
    INCOH_CS_2He(39) = 1.201D-01
  
    INCOH_CS_2He(40) = 1.023D-01
    INCOH_CS_2He(41) = 7.573D-02
    INCOH_CS_2He(42) = 6.078D-02
    INCOH_CS_2He(43) = 4.425D-02
    INCOH_CS_2He(44) = 3.516D-02
    INCOH_CS_2He(45) = 2.935D-02
    INCOH_CS_2He(46) = 2.529D-02
    INCOH_CS_2He(47) = 1.994D-02

    INCOH_CS_2He(48) = 1.655D-2 

    COH_y1 = dlog10(E_1H_2He(0))
    COH_y2 = dlog10(E_1H_2He(48))
    half_dy = ( COH_y2 - COH_y1 ) / 48.D0 / two
    COH_dy = ( COH_y2 - COH_y1 ) / 48.D0

    log10_E = dlog10( E_1H_2He )
!*****************************************************************************************
    end Subroutine PhotoElect_1H_COH_INCOH_CSect
!*****************************************************************************************


!*****************************************************************************************
    Subroutine COH_INCOH_CS( E_Mev, COH_CS, INCOH_CS )
!*****************************************************************************************
    implicit none
    real(mcp), intent(in) :: E_Mev
    real(mcp), intent(out) :: COH_CS, INCOH_CS
    real(mcp) :: E_ev, Ei, Ei1, COH_CS_H, COH_CS_He, INCOH_CS_H, INCOH_CS_He
    integer :: i
!*****************************************************************************************

    E_ev = E_Mev * 1.D6 ! trans it into eV

    i = floor( ( dlog10(E_ev) - COH_y1 ) / COH_dy ) 
    !E_1H_2He
    !Ei = E_1H_2He(i)
    !Ei1 = E_1H_2He(i+1) 
    write(*, *)'COH_INCOH_CS=', E_ev, dlog10(E_ev), COH_y1, COH_dy, i
    COH_CS_H = ( COH_CS_1H(i+1) - COH_CS_1H(i) ) / ( log10_E(i+1) - log10_E(i) ) *&
             ( dlog10(E_ev) - log10_E(i) ) + COH_CS_1H(i)

    COH_CS_He = ( COH_CS_2He(i+1) - COH_CS_2He(i) ) / ( log10_E(i+1) - log10_E(i) ) *&
             ( dlog10(E_ev) - log10_E(i) ) + COH_CS_2He(i)

    COH_CS = COH_CS_H + COH_CS_He * 1.D-1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    INCOH_CS_H = ( INCOH_CS_1H(i+1) - INCOH_CS_1H(i) ) / ( log10_E(i+1) - log10_E(i) ) *&
             ( dlog10(E_ev) - log10_E(i) ) + INCOH_CS_1H(i)

    INCOH_CS_He = ( INCOH_CS_2He(i+1) - INCOH_CS_2He(i) ) / ( log10_E(i+1) - log10_E(i) ) *&
             ( dlog10(E_ev) - log10_E(i) ) + INCOH_CS_2He(i)

    INCOH_CS = INCOH_CS_H + INCOH_CS_He * 1.D-1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!*****************************************************************************************
    end Subroutine COH_INCOH_CS
!*****************************************************************************************

    end module PhotoElectron









