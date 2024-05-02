! ================================================================
  module TriQuadM
! ================================================================

  implicit none

  private
  public :: TriQuad

  contains

! ======================================================================
  subroutine TriQuad(spoly)
! ======================================================================
! perform a simple triangulation of a quadrilateral

  use IntrType             ,only: rd
  use NetAlloc             ,only: nalloc
  use ProcessPairProcedures,only: DeallocatePolygon
  use Types                ,only: polygon
  use UserInput            ,only: ttype

  implicit none

  type(polygon),pointer      :: spoly
  real(kind=rd),dimension(4) :: xtmp,ytmp
  real(kind=rd)              :: dx1,dy1,dx2,dy2,cross
  integer,dimension(4)       :: edge_tmp,ie1_tmp,ie2_tmp,ntmp
  integer,dimension(3,2)     :: jv
  integer                    :: ic,iv,jvv


  continue


! check orientation of quad
  do iv = 1,4
    xtmp(iv) = spoly%cont(1)%x(iv)
    ytmp(iv) = spoly%cont(1)%y(iv)
    ntmp(iv) = spoly%cont(1)%node(iv)
    ie1_tmp(iv) = spoly%cont(1)%ie1(iv)
    ie2_tmp(iv) = spoly%cont(1)%ie2(iv)
    edge_tmp(iv) = spoly%cont(1)%edge(iv)
  end do

  dx1   = xtmp(3) - xtmp(1)
  dy1   = ytmp(3) - ytmp(1)
  dx2   = xtmp(4) - xtmp(2)
  dy2   = ytmp(4) - ytmp(2)
  cross = dx1*dy2 - dy1*dx2

  
  if (ttype == "Triangle") then

    if (cross >= 0.0_rd) then
      jv(1,1) = 1
      jv(2,1) = 2
      jv(3,1) = 3
      jv(1,2) = 1
      jv(2,2) = 3
      jv(3,2) = 4
    else
      jv(1,1) = 1
      jv(2,1) = 3
      jv(3,1) = 2
      jv(1,2) = 1
      jv(2,2) = 4
      jv(3,2) = 3
    end if

    call DeallocatePolygon(spoly)

!   create a new subject polygon with two contours
    spoly%num_cont = 2
    allocate(spoly%cont(2)); nalloc = nalloc + 1
    allocate(spoly%hole(2)); nalloc = nalloc + 1

    do ic = 1,2
      spoly%hole(ic) = 0
      spoly%cont(ic)%num_vert = 3

      allocate(spoly%cont(ic)%x(3))   ; nalloc = nalloc + 1
      allocate(spoly%cont(ic)%y(3))   ; nalloc = nalloc + 1
      allocate(spoly%cont(ic)%node(3)); nalloc = nalloc + 1
      allocate(spoly%cont(ic)%ie1(3)) ; nalloc = nalloc + 1
      allocate(spoly%cont(ic)%ie2(3)) ; nalloc = nalloc + 1
      allocate(spoly%cont(ic)%edge(3)); nalloc = nalloc + 1

      do iv = 1,3

        jvv = jv(iv,ic)

        spoly%cont(ic)%x(iv)    = xtmp(jvv)
        spoly%cont(ic)%y(iv)    = ytmp(jvv)
        spoly%cont(ic)%node(iv) = ntmp(jvv)
        spoly%cont(ic)%ie1(iv)  = ie1_tmp(jvv)
        spoly%cont(ic)%ie2(iv)  = ie2_tmp(jvv)
        spoly%cont(ic)%edge(iv) = edge_tmp(jvv)

      end do
    end do

  else if (ttype == "GPC") then

    if (cross > 0.0_rd) then

!     switch nodes 3 & 4 to form GPC tristrips

      xtmp(1) = spoly%cont(1)%x(3)
      spoly%cont(1)%x(3) = spoly%cont(1)%x(4)
      spoly%cont(1)%x(4) = xtmp(1)
  
      xtmp(1) = spoly%cont(1)%y(3)
      spoly%cont(1)%y(3) = spoly%cont(1)%y(4)
      spoly%cont(1)%y(4) = xtmp(1)

      ntmp(1) = spoly%cont(1)%node(3)
      spoly%cont(1)%node(3) = spoly%cont(1)%node(4)
      spoly%cont(1)%node(4) = ntmp(1)

    else

!     switch nodes 1 & 2 to form GPC tristrips

      xtmp(1) = spoly%cont(1)%x(1)
      spoly%cont(1)%x(1) = spoly%cont(1)%x(2)
      spoly%cont(1)%x(2) = xtmp(1)
  
      xtmp(1) = spoly%cont(1)%y(1)
      spoly%cont(1)%y(1) = spoly%cont(1)%y(2)
      spoly%cont(1)%y(2) = xtmp(1)

      ntmp(1) = spoly%cont(1)%node(1)
      spoly%cont(1)%node(1) = spoly%cont(1)%node(2)
      spoly%cont(1)%node(2) = ntmp(1)
    end if

  end if

  return
  end subroutine TriQuad

  end module TriQuadM
