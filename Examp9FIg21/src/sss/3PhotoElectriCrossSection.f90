    module PhotoElectron 
    use Basic_Variables_And_Methods 
    implicit none 
    real(mcp), dimension(0: 13) :: photo_elect_arr_c0, photo_elect_arr_c1, photo_elect_arr_c2
 
    contains 

    Subroutine Set_PhotoElect_CS_array()
    implicit none

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


    real(mcp) function PhotoElect_CSect( E )
    implicit none
    real(mcp), intent(in) :: E
 
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
    PhotoElect_CSect = PhotoElect_CSect / Sigma_T 
    end function PhotoElect_CSect

    end module PhotoElectron









