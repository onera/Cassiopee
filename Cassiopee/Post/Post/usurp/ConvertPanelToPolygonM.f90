! ================================================================
  module ConvertPanelToPolygonM
! ================================================================

  implicit none

  private
  public :: ConvertPanelToPolygon

  contains

! ====================================================================
  subroutine ConvertPanelToPolygon(g,x,xo,p,num_vert,node,edge)
! ====================================================================
! This procedure rotates a panel (x) into a facet plane (g,xo) and 
! stores it as a polygon (p).

  use IntrType,only: rd
  use NetAlloc,only: nalloc
  use Types   ,only: polygon

  implicit none

  real(kind=rd),dimension(3,3),intent(in) :: g   !rotation matrix
  real(kind=rd),dimension(3,4),intent(in) :: x   !panel corner points
  real(kind=rd),dimension(3)  ,intent(in) :: xo  !origin

  type(polygon),intent(inout)             :: p

  integer,intent(in),dimension(4)         :: node,edge
  integer,intent(in)                      :: num_vert
  integer                                 :: i,ierr


  continue


  p%num_cont = 1

  allocate(p%hole(p%num_cont),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on hole in ConvertPanelToPolygon"
    stop
  end if
  p%hole(1) = 0

  allocate(p%cont(p%num_cont),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on cont in ConvertPanelToPolygon"
    stop
  end if
  p%cont(1)%num_vert = num_vert

  allocate(p%cont(1)%x(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on x in ConvertPanelToPolygon"
    stop
  end if

  allocate(p%cont(1)%y(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on y in ConvertPanelToPolygon"
    stop
  end if

  allocate(p%cont(1)%node(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on node in ConvertPanelToPolygon"
    stop
  end if

  allocate(p%cont(1)%ie1(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on ie1 in ConvertPanelToPolygon"
    stop
  end if

  allocate(p%cont(1)%ie2(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on ie2 in ConvertPanelToPolygon"
    stop
  end if

  allocate(p%cont(1)%edge(num_vert),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed on edge in ConvertPanelToPolygon"
    stop
  end if

  do i = 1,num_vert
    p%cont(1)%x(i) = g(1,1)*(x(1,i)-xo(1)) &
                   + g(1,2)*(x(2,i)-xo(2)) &
                   + g(1,3)*(x(3,i)-xo(3))
    p%cont(1)%y(i) = g(2,1)*(x(1,i)-xo(1)) &
                   + g(2,2)*(x(2,i)-xo(2)) &
                   + g(2,3)*(x(3,i)-xo(3))
    p%cont(1)%node(i) = node(i)
    p%cont(1)%ie1(i)  = 0
    p%cont(1)%ie2(i)  = 0
    p%cont(1)%edge(i) = edge(i)
  end do

  return
  end subroutine ConvertPanelToPolygon

  end module ConvertPanelToPolygonM
