! ================================================================
  module WritePanelWeightsM
! ================================================================

  implicit none

  private
  public :: WritePanelWeights

  contains

! ====================================================================
  subroutine WritePanelWeights(p)
! ====================================================================
! This procedure writes the panel weights to a data file, using one
! record for each surface patch.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: num_patches
  use Types       ,only: panel,surf
  use UserInput   ,only: pathname,nphase,showcpu,cpufmt

  implicit none

  type(panel),pointer,dimension (:) :: p
  real(kind=rs)                     :: time1,time2
  integer                           :: ip,ipan


  continue


  call my_cpu_time(time1)

  open( unit   = 3,                                   &
        file   = trim(pathname)//"panel_weights.dat", &
        action = "write")


  if (nphase) then

!   write all panels in a single record
    write(unit=3,fmt=*)(p(ip)%ratio,ip=lbound(p,1),ubound(p,1))

  else

!   use one record for each surface patch
    do ip = 1,num_patches
      write(unit=3,fmt=*)(p(ipan)%ratio, &
           ipan=surf(ip)%first_panel,surf(ip)%last_panel)
    end do

  end if

  close(unit=3)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in WritePanelWeights:',time2-time1
  end if


  return
  end subroutine WritePanelWeights

  end module WritePanelWeightsM
