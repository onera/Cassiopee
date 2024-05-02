! ================================================================
  module CalculateProjectedAreasM
! ================================================================

  implicit none

  private
  public :: CalculateProjectedAreas

  contains

! ======================================================================
  subroutine CalculateProjectedAreas(p)
! ======================================================================
! procedure to check the projected areas of the solid surface
! (which should normally sum to zero) as well as the wetted area.
! note that the dummy argument p is a pointer, which requires that
! we have an explicit interface for this procedure.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,showcpu

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rd)     :: areax,areay
  real(kind=rd)     :: areaz,areat
  real(kind=rd)     :: area0,areab,areap
  real(kind=rd)     :: ratio0,ratiob
  real(kind=rs)     :: time1,time2
  integer           :: ipan,ip1,ip2
  character(len=51) :: afmt


  continue


! initialization
  call my_cpu_time(time1)
  areax = 0.0_rd
  areay = 0.0_rd
  areaz = 0.0_rd
  areat = 0.0_rd
  areab = 0.0_rd
  area0 = 0.0_rd
  ip1   = lbound(p,1)
  ip2   = ubound(p,1)
  
  do ipan = ip1,ip2
    areap = p(ipan)%area_ref
    if (p(ipan)%iblank == 1) areab = areab + areap
    area0 = area0 + areap
    areax = areax + areap*p(ipan)%xn(1) * p(ipan)%ratio
    areay = areay + areap*p(ipan)%xn(2) * p(ipan)%ratio
    areaz = areaz + areap*p(ipan)%xn(3) * p(ipan)%ratio
    areat = areat + areap               * p(ipan)%ratio
  end do

  write(0,*)
  write(0,*)"Integrated Areas"
  write(0,"(t8,a1,t22,a1,t36,a1,t49,a5)")"X","Y","Z","Total"
  write(0,*)'----------------------------------------------------------'
  write(0,'(3e14.6,e15.6,/)')areax,areay,areaz,areat
  write(0,*)

  if ((area0 == 0.0_rd).or.(areab == 0.0_rd)) then
    write(0,'(" The original area of all of the panels was ",e15.6)')area0
  else
    ratio0 = (area0-areat)/area0*100.0
    ratiob = (areab-areat)/areab*100.0

    afmt = '(" The original area was reduced by ",f6.2,a)'

    if (nint(ratiob*100) /= nint(ratio0*100)) then
      write(0,afmt)ratio0,"% (all panels)"
      write(0,afmt)ratiob,"% (panels with iblank=1)"
    else
      write(0,afmt)ratio0,"%"
    end if
  end if


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CalculateProjectedAreas:',time2-time1
  end if


  return
  end subroutine CalculateProjectedAreas

  end module CalculateProjectedAreasM
