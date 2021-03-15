!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MODULE SubFunction
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      USE constants
      USE RandUtils 
      IMPLICIT NONE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer, parameter :: N_H3 = 5000
      real(mcp), parameter :: E_up_H3 = 1.D25 * h_ev * 1.D-6
      real(mcp), parameter :: E_dw_H3 = 1.D4 * h_ev * 1.D-6
      real(mcp), parameter :: log_E_up_H3 = dlog10( E_up_H3 )
      real(mcp), parameter :: log_E_dw_H3 = dlog10( E_dw_H3 )
      real(mcp), parameter :: dlog_E_H3 = ( log_E_up_H3 - log_E_dw_H3 ) / N_H3
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(mcp), parameter :: E_critical = 1.D-5   ! Its unit is MeV. 
      real(mcp), parameter :: nu_critical = E_critical * 1.D6 / h_ev   ! Its unit is MeV. 
      real(mcp), parameter :: nu_dw_H3 = nu_critical
      real(mcp), parameter :: nu_up_H3 = 1.D38
      real(mcp), parameter :: log_nu_up_H3 = dlog10( nu_up_H3 )
      real(mcp), parameter :: log_nu_dw_H3 = dlog10( nu_dw_H3 )
      real(mcp), parameter :: dlog_nu_H3 = ( log_nu_up_H3 - log_nu_dw_H3 ) / N_H3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CONTAINS
 
!******************************************************************* 
      subroutine Matrix_Multiplication44X44_Sub(a, b, c)
!*******************************************************************
      implicit none
      real(mcp), dimension(1:4, 1:4), intent(in) :: a, b 
      real(mcp), dimension(1:4, 1:4), intent(out) :: c
      integer :: i, j, k
!*******************************************************************
 
      c = c * zero
      Do i = 1, 4
          Do j = 1, 4
              Do k = 1, 4 
                  c(i, j) = c(i, j) + a(i, k) * b(k, j)
              Enddo
          Enddo
      Enddo  
 
      end subroutine Matrix_Multiplication44X44_Sub

!******************************************************************** 
      subroutine Matrix_Multiplication33X33_Sub(a, b, c)
!******************************************************************** 
      implicit none
      real(mcp), dimension(1:3, 1:3), intent(in) :: a, b 
      real(mcp), dimension(1:3, 1:3), intent(out) :: c
      integer :: i, j, k
!***********************************************************************
 
      c = c * zero
      Do i = 1, 3
          Do j = 1, 3
              Do k = 1, 3 
                  c(i, j) = c(i, j) + a(i, k) * b(k, j)
              Enddo
          Enddo
      Enddo  
 
      end subroutine Matrix_Multiplication33X33_Sub

!*******************************************************************************************************
      subroutine Matrix_Multiplication14X44_Sub(a, b, c)
!*******************************************************************************************************
      implicit none
      real(mcp), dimension(1:4), intent(in) :: a
      real(mcp), dimension(1:4, 1:4), intent(in) :: b 
      real(mcp), dimension(1:4), intent(out) :: c
      integer :: i, j 
!***********************************************************************

      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'c1 === ', c
      !write(unit = *, fmt = *)'************************************************************' 

      c = c * zero
      Do i = 1, 4
          Do j = 1, 4
              c(i) = c(i) + a(j) * b(j, i)
          Enddo
      Enddo 

      end subroutine Matrix_Multiplication14X44_Sub

!*******************************************************************************************************
      subroutine Matrix_Multiplication44X41_Sub(a, b, c)
!*******************************************************************************************************
      implicit none
      real(mcp), dimension(1:4, 1:4), intent(in) :: a 
      real(mcp), dimension(1:4), intent(in) :: b
      real(mcp), dimension(1:4), intent(out) :: c
      integer :: i, j 

      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'c1 === ', c
      !write(unit = *, fmt = *)'************************************************************' 

      c = c * zero
      Do i = 1, 4
          Do j = 1, 4
              c(i) = c(i) + a(i, j) * b(j)
          Enddo
      Enddo 

      end subroutine Matrix_Multiplication44X41_Sub

!*******************************************************************************************************
      subroutine Matrix_Multiplication13X33_Sub(a, b, c)
!*******************************************************************************************************
      implicit none
      real(mcp), dimension(1:3), intent(in) :: a
      real(mcp), dimension(1:3, 1:3), intent(in) :: b 
      real(mcp), dimension(1:3), intent(out) :: c
      integer :: i, j 
!***********************************************************************

      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'c1 === ', c
      !write(unit = *, fmt = *)'************************************************************' 

      c = c * zero
      Do i = 1, 3
          Do j = 1, 3
              c(i) = c(i) + a(j) * b(j, i)
          Enddo
      Enddo 

      end subroutine Matrix_Multiplication13X33_Sub
 

!*******************************************************************************************************
      subroutine Matrix_Multiplication33X31_Sub(a, b, c)
!*******************************************************************************************************
      implicit none
      real(mcp), dimension(1:3), intent(in) :: b
      real(mcp), dimension(1:3, 1:3), intent(in) :: a 
      real(mcp), dimension(1:3), intent(out) :: c
      integer :: i, j 

      !write(unit = *, fmt = *)'************************************************************' 
      !write(unit = *, fmt = *)'c1 === ', c
      !write(unit = *, fmt = *)'************************************************************' 

      c = c * zero
      Do i = 1, 3
          Do j = 1, 3
              c(i) = c(i) + a(i, j) * b(j)
          Enddo
      Enddo 

      end subroutine Matrix_Multiplication33X31_Sub
  
!*********************************************************************** 
      real(mcp) function Vector3D_Length( Vector )
!*********************************************************************** 
      real(mcp), dimension(1:3), intent(in) :: Vector
      !********************************************************* 

      Vector3D_Length = DSQRT( Vector(1)*Vector(1) + Vector(2)*&
                               Vector(2) + Vector(3)*Vector(3) )
      end function Vector3D_Length

!*********************************************************************** 
      real(mcp) function Vector4D_Inner_Product( Vec1, Vec2 )
!*********************************************************************** 
      real(mcp), dimension(1:4), intent(in) :: Vec1, Vec2
      !********************************************************* 

      Vector4D_Inner_Product = Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2) + &
                               Vec1(3)*Vec2(3) + Vec1(4)*Vec2(4)
      end function Vector4D_Inner_Product

!*********************************************************************** 
      real(mcp) function Vector3D_Inner_Product( Vec1, Vec2 )
!*********************************************************************** 
      real(mcp), dimension(1:3), intent(in) :: Vec1, Vec2
      !********************************************************* 

      Vector3D_Inner_Product = Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2) + &
                               Vec1(3)*Vec2(3)
      end function Vector3D_Inner_Product
      
!*********************************************************************** 
      real(mcp) function Vector4D_Inner_Product_Mski( Vec1, Vec2 )
!*********************************************************************** 
      real(mcp), dimension(1: 4), intent(in) :: Vec1, Vec2
      !********************************************************* 

      Vector4D_Inner_Product_Mski = - Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2) + &
                               Vec1(3)*Vec2(3) + Vec1(4)*Vec2(4)
      end function Vector4D_Inner_Product_Mski

!*******************************************************************************************************
      subroutine Inverse_Matrix_Of_Matrix_3X3_Sub(Matrix, InvMatrix)
!*******************************************************************************************************
      implicit none 
      real(mcp), dimension(1:3, 1:3), intent(in) :: Matrix
      real(mcp), dimension(1:3, 1:3), intent(out) :: InvMatrix
      integer :: i, j
       
      Do i = 1, 3
          Do j = 1, 3
             InvMatrix(j, i) = Matrix(i, j)
          Enddo
      Enddo
      end subroutine Inverse_Matrix_Of_Matrix_3X3_Sub

!*******************************************************************************************************
      subroutine Matrix3X3_Trans_Sub(a, b)
!*******************************************************************************************************
      implicit none
      real(mcp), dimension(1:3, 1:3), intent(in) :: a
      real(mcp), dimension(1:3, 1:3), intent(out) :: b 
      integer :: i, j 
  
      Do i = 1, 3
          Do j = 1, 3
              b(i, j) = a(j, i)
          Enddo
      Enddo 

      end subroutine Matrix3X3_Trans_Sub
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!|       ||
!|       ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END MODULE SubFunction
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|       ||
!|       ||
!|        ||
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~














 
