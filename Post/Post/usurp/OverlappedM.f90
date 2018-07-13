! ================================================================
  module OverlappedM
! ================================================================

  implicit none

  private
  public :: Overlapped

  contains

! ========================================================
  pure function overlapped(bound1,bound2) result (overlap)
! ========================================================
! Function to test whether two axis-aligned rectangular 
! prisms overlap by checking for a separation plane 
! perpendicular to one of the Cartesian axes.

  use IntrType,only: rd

  implicit none

  real(kind=rd),dimension(3,2),intent(in) :: bound1,bound2
  integer :: idir
  logical :: overlap


  continue

  
  overlap = .false.

  do idir = 1,3
    if (bound2(idir,1) > bound1(idir,2)) return
    if (bound1(idir,1) > bound2(idir,2)) return
  end do

  overlap = .true.

  return
  end function overlapped

  end module OverlappedM
