! ======================================================================
  module UnusedVerticesM
! ======================================================================

  implicit none

  private
  public :: UnusedVertices

  contains

! ======================================================================
  subroutine UnusedVertices(p)
! ======================================================================
! This procedure finds any unused vertices and creates a new numbering
! system that includes only the vertices that are actually used.
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,showcpu
  use VertexData  ,only: num_nodes,new_nodes,num_used,used,idup

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rs) :: time1,time2

  integer :: ic,ierr,ip1,ip2,ipan,iv
  integer :: num_vertices


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)

  num_vertices = num_nodes + new_nodes

  allocate(idup(num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed allocating idup in WriteTriQ"
    write(0,*)"num_vertices = ",num_vertices
    stop
  end if
  idup = 0

  allocate(used(num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed allocating used in WriteTriQ"
    write(0,*)"num_vertices = ",num_vertices
    stop
  end if
  used = .false.


! mark vertices that are actually used
  do ipan = ip1,ip2
    do ic = 1,p(ipan)%poly%num_cont
      do iv = 1,p(ipan)%poly%cont(ic)%num_vert
        used(p(ipan)%poly%cont(ic)%node(iv)) = .true.
      end do
    end do
  end do


! count the number of used vertices and create a new numbering
! system including only vertices that are actually used
  num_used = 0
  do iv = 1,num_vertices
    if (used(iv)) then
      num_used = num_used + 1
      idup(iv) = num_used
    end if
  end do
  

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in UnusedVertices: ',time2-time1
  end if

  return
  end subroutine UnusedVertices

  end module UnusedVerticesM
