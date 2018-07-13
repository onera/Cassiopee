! ================================================================
  module PrintPolyM
! ================================================================

  implicit none

  private
  public :: PrintPoly

  contains

! ==================================================================
  subroutine PrintPoly(spoly,iu,tag)
! ==================================================================
! print a polygon to a specified unit

  use Types    ,only: polygon
  use UserInput,only: ltri

  implicit none

  type(polygon) :: spoly

  integer,intent(in) :: iu
  integer :: ic,iv

  character(len=*),intent(in) :: tag

  logical :: writec


  continue



  do ic = 1,spoly%num_cont

!   if the "ltri" array is allocated, only write 
!   contours who have encountered errors
    writec = .true.
    if (allocated(ltri)) then
      if (ltri(ic) == 0) writec = .false.
    end if

    if (writec) then
      write(iu,20)trim(tag),ic,spoly%num_cont
      do iv = 1,spoly%cont(ic)%num_vert
        write(iu,30)spoly%cont(ic)%x(iv),    &
                    spoly%cont(ic)%y(iv),    &
                    spoly%cont(ic)%node(iv), &
                    spoly%cont(ic)%edge(iv)
      end do
      write(iu,30)spoly%cont(ic)%x(1),    &
                  spoly%cont(ic)%y(1),    &
                  spoly%cont(ic)%node(1), &
                  spoly%cont(ic)%edge(1)
    end if
  end do

  20  format(1x,'zone t = "',a,' (',i2,' of ',i2,')"')
  30  format(2(1x,es16.9),2(1x,i9))

  return
  end subroutine PrintPoly

  end module PrintPolyM
