! ================================================================
  module CalculatePolygonAreaM
! ================================================================

  implicit none

  private
  public :: CalculatePolygonArea

  contains

! ========================================================
  function CalculatePolygonArea(p,ipan) result (area)
! ========================================================
! Procedure to calculate the area of a polygon
! using the trapezoid rule.  External contours
! (hole=0) add to the area; hole contours (hole=1)
! subtract from the area.  Note that the dummy
! argument p is a pointer, which requires that we
! have an explicit interface for this procedure.

  use IntrType,only: rd
  use Types   ,only: polygon

  implicit none

  type(polygon)      :: p
  real(kind=rd)      :: area
  real(kind=rd)      :: x1,x2,y1,y2,areac
  integer,intent(in) :: ipan
  integer            :: ic,iv1,iv2


  continue


! calculate the area of the polygon
  area = 0.0_rd
  do ic = 1,p%num_cont

    areac = 0.0_rd
    do iv1 = 1,p%cont(ic)%num_vert
      if (iv1 == p%cont(ic)%num_vert) then
        iv2 = 1
      else
        iv2 = iv1 + 1
      end if
      x2 = p%cont(ic)%x(iv2)
      x1 = p%cont(ic)%x(iv1)
      y2 = p%cont(ic)%y(iv2)
      y1 = p%cont(ic)%y(iv1)
      areac = areac + (x1-x2)*(y2+y1)
    end do

    if (p%hole(ic) == 0) then
      area = area + abs(areac)
    else if (p%hole(ic) == 1) then
      area = area - abs(areac)
    else
      write(0,*)"ERROR! invalid value of hole."
      write(0,*)"panel, hole value = ",ipan,p%hole(ic)
      write(0,*)"The area of this polygon is being ignored to allow"
      write(0,*)"the code to continue, but this error MIGHT indicate"
      write(0,*)"memory management issues affecting the results."
!     stop
    end if

  end do
  area = 0.5_rd * area

  return
  end function CalculatePolygonArea

  end module CalculatePolygonAreaM
