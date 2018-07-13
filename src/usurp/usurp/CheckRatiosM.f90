! ================================================================
  module CheckRatiosM
! ================================================================

  implicit none

  private
  public :: CheckRatios

  contains

! ======================================================================
  subroutine CheckRatios(p)
! ======================================================================
! This procedure checks to make sure the area ratios
! are between 0 and 1.  Ratios below zero are considered
! fatal errors.  Ratios significantly above 1.0 are reported
! and corrected.  Note that because p is a pointer array,
! this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,showcpu

  implicit none

  type(panel),pointer,dimension (:) :: p
  real(kind=rs)                     :: time1,time2
  integer                           :: ipan,ic,iv,ip1,ip2


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)


  do ipan = ip1,ip2

    if (p(ipan)%ratio < 0.0_rd) then
      write(0,*)'ERROR!  polygon ratio < 0, ratio = ',p(ipan)%ratio
      write(0,*)'num_cont = ',p(ipan)%poly%num_cont
      do ic = 1,p(ipan)%poly%num_cont
        write(0,*)'hole = ',p(ipan)%poly%hole(ic)
        write(0,*)'num_vert = ',p(ipan)%poly%cont(ic)%num_vert
        do iv = 1,p(ipan)%poly%cont(ic)%num_vert
          write(0,*)p(ipan)%poly%cont(ic)%x(iv),p(ipan)%poly%cont(ic)%y(iv)
        end do
      end do
      stop
    end if

    if (p(ipan)%ratio > 1.0001_rd) then
      write(*,*)'WARNING! panel, area ratio = ',ipan,p(ipan)%ratio
      p(ipan)%ratio = 1.0_rd
    end if

  end do

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CheckRatios:',time2-time1
  end if


  return
  end subroutine CheckRatios

  end module CheckRatiosM
