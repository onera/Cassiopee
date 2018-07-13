! ================================================================
  module InsertNodeM
! ================================================================

  implicit none

  private
  public :: InsertNode

  contains

! ===========================================================
  subroutine InsertNode(rpoly,rc,rv,rpan,xnew,ynew,nnew,enew)
! ===========================================================
! Insert a specified vertex (xnew,ynew,nnew,enew) into an
! existing polygon contour (rpoly, rc) just after a
! specified vertex (rv).  Note that rpoly is a pointer,
! so an explicit interface is required.

  use IntrType ,only: rd
  use NetAlloc ,only: nalloc
  use Types    ,only: polygon,contour
  use UserInput,only: plotpanel

  implicit none

  type(polygon),pointer :: rpoly
  type(contour)         :: temp

  real(kind=rd),intent(in) :: xnew,ynew
  integer,intent(in)       :: rc,rv,rpan,nnew,enew
  integer                  :: jv,num_vert


  continue


  if (rpan == plotpanel) then
    write(*,*)"In InsertNode: ",nnew,enew
  end if

! create a copy of the recipient contour
  temp%num_vert = rpoly%cont(rc)%num_vert
  allocate(temp%x(temp%num_vert))   ; nalloc = nalloc + 1
  allocate(temp%y(temp%num_vert))   ; nalloc = nalloc + 1
  allocate(temp%node(temp%num_vert)); nalloc = nalloc + 1
  allocate(temp%ie1(temp%num_vert)) ; nalloc = nalloc + 1
  allocate(temp%ie2(temp%num_vert)) ; nalloc = nalloc + 1
  allocate(temp%edge(temp%num_vert)); nalloc = nalloc + 1
  do jv = 1,temp%num_vert
    temp%x(jv)    = rpoly%cont(rc)%x(jv)
    temp%y(jv)    = rpoly%cont(rc)%y(jv)
    temp%node(jv) = rpoly%cont(rc)%node(jv)
    temp%ie1(jv)  = rpoly%cont(rc)%ie1(jv)
    temp%ie2(jv)  = rpoly%cont(rc)%ie2(jv)
    temp%edge(jv) = rpoly%cont(rc)%edge(jv)
  end do


! delete the old contour
  deallocate(rpoly%cont(rc)%x)   ; nalloc = nalloc - 1
  deallocate(rpoly%cont(rc)%y)   ; nalloc = nalloc - 1
  deallocate(rpoly%cont(rc)%node); nalloc = nalloc - 1
  deallocate(rpoly%cont(rc)%ie1) ; nalloc = nalloc - 1
  deallocate(rpoly%cont(rc)%ie2) ; nalloc = nalloc - 1
  deallocate(rpoly%cont(rc)%edge); nalloc = nalloc - 1

! create space for the new contour
  num_vert = temp%num_vert + 1
  rpoly%cont(rc)%num_vert = num_vert
  allocate(rpoly%cont(rc)%x(num_vert))   ; nalloc = nalloc + 1
  allocate(rpoly%cont(rc)%y(num_vert))   ; nalloc = nalloc + 1
  allocate(rpoly%cont(rc)%node(num_vert)); nalloc = nalloc + 1
  allocate(rpoly%cont(rc)%ie1(num_vert)) ; nalloc = nalloc + 1
  allocate(rpoly%cont(rc)%ie2(num_vert)) ; nalloc = nalloc + 1
  allocate(rpoly%cont(rc)%edge(num_vert)); nalloc = nalloc + 1

! copy the original contour into the new space
  do jv = 1,rv
    rpoly%cont(rc)%x(jv)    = temp%x(jv)
    rpoly%cont(rc)%y(jv)    = temp%y(jv)
    rpoly%cont(rc)%node(jv) = temp%node(jv)
    rpoly%cont(rc)%ie1(jv)  = temp%ie1(jv)
    rpoly%cont(rc)%ie2(jv)  = temp%ie2(jv)
    rpoly%cont(rc)%edge(jv) = temp%edge(jv)
  end do
  do jv = rv + 1,temp%num_vert
    rpoly%cont(rc)%x(jv+1)    = temp%x(jv)
    rpoly%cont(rc)%y(jv+1)    = temp%y(jv)
    rpoly%cont(rc)%node(jv+1) = temp%node(jv)
    rpoly%cont(rc)%ie1(jv+1)  = temp%ie1(jv)
    rpoly%cont(rc)%ie2(jv+1)  = temp%ie2(jv)
    rpoly%cont(rc)%edge(jv+1) = temp%edge(jv)
  end do

! insert the new node (and new node number)
  rpoly%cont(rc)%x(rv+1)    = xnew
  rpoly%cont(rc)%y(rv+1)    = ynew
  rpoly%cont(rc)%node(rv+1) = nnew
  rpoly%cont(rc)%ie1(rv+1)  = 0
  rpoly%cont(rc)%ie2(rv+1)  = 0
  rpoly%cont(rc)%edge(rv+1) = enew

! delete the temporary copy of the original contour
  deallocate(temp%x)   ; nalloc = nalloc - 1
  deallocate(temp%y)   ; nalloc = nalloc - 1
  deallocate(temp%node); nalloc = nalloc - 1
  deallocate(temp%ie1) ; nalloc = nalloc - 1
  deallocate(temp%ie2) ; nalloc = nalloc - 1
  deallocate(temp%edge); nalloc = nalloc - 1

  return
  end subroutine InsertNode

  end module InsertNodeM
