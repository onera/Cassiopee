! ================================================================
  module RotateBackM
! ================================================================

  implicit none

  private
  public :: RotateBack

  contains


! ===================================================================
  subroutine RotateBack(x,y,g,xo,xt)
! ===================================================================
! Rotate a point back into its original panel plane.

  use IntrType,only: rd

  implicit none

  real(kind=rd),intent(in)                :: x,y
  real(kind=rd),dimension(3,3),intent(in) :: g
  real(kind=rd),dimension(3),intent(in)   :: xo
  real(kind=rd),dimension(3),intent(out)  :: xt


  continue


  xt(1) = g(1,1)*x + g(2,1)*y + xo(1)
  xt(2) = g(1,2)*x + g(2,2)*y + xo(2)
  xt(3) = g(1,3)*x + g(2,3)*y + xo(3)

  return
  end subroutine RotateBack

  end module RotateBackM
