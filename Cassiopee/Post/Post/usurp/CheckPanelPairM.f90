! ================================================================
  module CheckPanelPairM
! ================================================================

  implicit none

  private
  public :: CheckPanelPair

  contains

! =========================================================
  subroutine CheckPanelPair(p)
! =========================================================
! check whether two panels overlap and tell why or why not

  use IntrType   ,only: rd
  use OverlappedM,only: Overlapped
  use Types      ,only: panel
  use UserInput  ,only: disjoin,icpanel1,icpanel2,dotlimit

  implicit none

  type(panel),pointer,dimension(:) :: p

  type(panel),pointer :: q1,q2
  real(kind=rd)       :: dot
  integer             :: ip1,ip2
  logical             :: ov


  continue


! initialization
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)

  if (min(icpanel1,icpanel2) < ip1) then
    write(0,*)"ERROR! panel in --check-panels is out of range."
    write(0,*)"requested panel ",min(icpanel1,icpanel2)
    write(0,*)"lower bound on panel array is ",ip1
    stop
  end if
  if (max(icpanel1,icpanel2) > ip2) then
    write(0,*)"ERROR! panel in --check-panels is out of range."
    write(0,*)"requested panel ",max(icpanel1,icpanel2)
    write(0,*)"upper bound on panel array is ",ip2
    stop
  end if

  q1 => p(icpanel1)
  q2 => p(icpanel2)

  ov  = overlapped(q1%bound,q2%bound)
  dot = q1%xn(1)*q2%xn(1) + q1%xn(2)*q2%xn(2) + q1%xn(3)*q2%xn(3)

  write(0,*)
  write(0,*)"checking panels "
  write(0,*)
  write(0,10)q1%panel_number,q1%isurf,q1%icomp,q1%iblank
  10 format(1x,"panel ",i9,", surface ",i4,", component ",i3,", i-blank = ",i3)
  write(0,*)"bounding box"
  write(0,*)q1%bound(1,1),q1%bound(2,1),q1%bound(3,1)
  write(0,*)q1%bound(1,2),q1%bound(2,2),q1%bound(3,2)
  write(0,*)
  write(0,10)q2%panel_number,q2%isurf,q2%icomp,q2%iblank
  write(0,*)"bounding box"
  write(0,*)q2%bound(1,1),q2%bound(2,1),q2%bound(3,1)
  write(0,*)q2%bound(1,2),q2%bound(2,2),q2%bound(3,2)
  write(0,*)
  write(0,*)"dot product = ",dot
  write(0,*)
  write(0,*)"area"
  write(0,*)q1%panel_number,q1%area_ref
  write(0,*)q2%panel_number,q2%area_ref


! reasons why the pair is valid
  if (ov) write(0,*)"the panel bounding boxes overlap."
  if (dot > dotlimit) write(0,*)"the dot product is greater than ",dotlimit
  if (q1%isurf /= q2%isurf) write(0,*)"the panels are from different surfaces."
  if (q1%icomp == q2%icomp) then
    write(0,*)"the panels are from the same component."
  else if (.not.disjoin) then
    write(0,*)"the panels are from different components but are allowed ", &
              "to interact because the disjoin option is set to false."
  end if


! reasons why the pair is invalid:

! bounding boxes must overlap for pair to be valid
  if (.not.ov) write(0,*)"the panel pair is invalid because the ", &
                         "bounding boxes do not overlap."

! dot product must be positive for pair to be valid
  if (dot <= dotlimit) write(0,*)"the panel pair is invalid because the ", &
                            "dot product is less than ",dotlimit

! panels must not be from the same surface
  if (q1%isurf == q2%isurf) &
     write(0,*)"the panel pair is invalid because the panels are ", &
               "from the same surface."


! panels must be from the same component is disjoin = .true.
  if ((q1%icomp /= q2%icomp).and.disjoin) then
    write(0,*)"the panel pair is invalid because the panels are from ",&
              "different components and the disjoin option is set to true."
  end if

  nullify(q1,q2)

  write(0,*)
  write(0,*)"stopping code: --check-panels option detected on command line."
  stop

  return
  end subroutine CheckPanelPair

  end module CheckPanelPairM
