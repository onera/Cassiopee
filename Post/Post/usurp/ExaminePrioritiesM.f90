! ================================================================
  module ExaminePrioritiesM
! ================================================================

  implicit none

  private
  public :: ExaminePriorities

  contains

! ========================================
  subroutine ExaminePriorities()
! ========================================
! check to see whether the patch ranking system
! is already consistent with any priority pair
! specifications

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: rankp,num_patches
  use PriPairs    ,only: ipri
  use UserInput   ,only: cpufmt,showcpu

  implicit none

  real(kind=rs) :: time1,time2
  integer       :: i,j


  continue


  call my_cpu_time(time1)


  do j = 1,num_patches
  do i = 1,num_patches
    if ((ipri(i,j) > 0).and.(rankp(i) < rankp(j))) then
      write(*,"(1x,a,i4,a,i4)")"WARNING! set priority reverses "// &
                            "ranking for patches ",i," and ",j
    end if
  end do
  end do


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in ExaminePriorities:',time2-time1
  end if


  return
  end subroutine ExaminePriorities

  end module ExaminePrioritiesM
