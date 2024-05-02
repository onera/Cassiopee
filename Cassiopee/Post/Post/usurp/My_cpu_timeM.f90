! ================================================================
  module my_cpu_timeM
! ================================================================

  implicit none

  private
  public :: my_cpu_time

  contains

! ===========================================
  subroutine my_cpu_time(time)
! ===========================================
! Return cpu time in seconds using intrinsics
! available under different compilers.

  use IntrType,only: rs

  implicit none

  real(kind=rs),intent(out) :: time

  continue


  time = 0.0_rs

  return
  end subroutine my_cpu_time

  end module my_cpu_timeM
