! ================================================================
  module LoadPanelsNPHASEM
! ================================================================

  implicit none

  private
  public :: LoadPanelsNPHASE

  contains

! ============================================================
  subroutine LoadPanelsNPHASE(p)
! ============================================================
! Load elements from unstructured NPHASE data file into panels.
! Note that the dummy argument p is a point, so an explicit
! interface is required for this procedure.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: num_panels,num_scalars,icomp
  use Types       ,only: panel,scal
  use UserInput   ,only: cpufmt,showcpu,workdir
  use VertexData  ,only: num_nodes,xv

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rs)         :: time1,time2
  integer,dimension (4) :: node
  integer               :: i,idir,ierr,ipan,pid
  character(len=6)      :: tag


  continue


  call my_cpu_time(time1)


! read the original element data

  open(unit     = 3,                                    &
       file     = trim(workdir)//"nphase_elements.dat", &
       form     = "formatted",                          &
       status   = "old",                                &
       position = "rewind",                             &
       action   = "read")


! read node data
  read(unit=3,fmt=*)num_nodes

  allocate(xv(3,num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(unit=*,fmt=*)"allocation of xv failed in LoadPanelsNPHASE."
    write(unit=*,fmt=*)"num_nodes = ",num_nodes
    stop
  end if

  do i = 1,num_nodes
    read(unit=3,fmt=*)(xv(idir,i),idir=1,3)
  end do


! read element data and store directly into panels
  read(unit=3,fmt=*)num_panels

  do ipan = 1,num_panels
    read(unit=3,fmt=*)p(ipan)%num_vertices,pid, &
             (node(i),i=1,p(ipan)%num_vertices)

    p(ipan)%panel_number = ipan
    p(ipan)%iblank       = 1
    p(ipan)%isurf        = pid + 1
    p(ipan)%icomp        = icomp(p(ipan)%isurf)

    do i = 1,p(ipan)%num_vertices
      p(ipan)%node(i) = node(i)
    end do
  end do

! read scalar data
  allocate(scal(num_scalars,num_panels)); nalloc = nalloc + 1
  do i = 1,num_scalars
  do ipan=1,num_panels
    read(unit=3,fmt=*)tag,scal(i,ipan)
  end do
  end do

  close(unit=3)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in LoadPanelsNPHASE:',time2-time1
  end if

  return
  end subroutine LoadPanelsNPHASE

  end module LoadPanelsNPHASEM
