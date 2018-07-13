! ================================================================
  module ReadPanelWeightsM
! ================================================================

  implicit none

  private
  public :: ReadPanelWeights

  contains

! ====================================================================
  subroutine ReadPanelWeights(p)
! ====================================================================
! This procedure reads the panel weights from a data file, using one
! record for each surface patch.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use PatchInfo   ,only: num_patches
  use Types       ,only: panel,surf
  use UserInput   ,only: cpufmt,nphase,showcpu,workdir

  implicit none

  type(panel),pointer,dimension (:) :: p
  real(kind=rs)                     :: time1,time2
  integer                           :: ip,ipan
  logical                           :: weightfile


  continue


  call my_cpu_time(time1)

  inquire(file=trim(workdir)//"panel_weights.dat",exist=weightfile)

  if (.not.weightfile) then
    write(0,*)
    write(0,*)"ERROR!  option --keep-weights has been specified but"
    write(0,*)"panel_weights.dat file was not found in ",trim(workdir)
    stop
  else


    open( unit   = 3,                                   &
          file   = trim(workdir)//"panel_weights.dat", &
          action = "read")


    if (nphase) then

!     read all panels from a single record
      read(unit=3,fmt=*)(p(ip)%ratio,ip=lbound(p,1),ubound(p,1))

    else

!     use one record for each surface patch
      do ip = 1,num_patches
        read(unit=3,fmt=*)(p(ipan)%ratio, &
             ipan=surf(ip)%first_panel,surf(ip)%last_panel)
      end do

    end if

    close(unit=3)

  end if


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in ReadPanelWeights:',time2-time1
  end if


  return
  end subroutine ReadPanelWeights

  end module ReadPanelWeightsM
