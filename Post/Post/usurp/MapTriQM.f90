! ================================================================
  module MapTriQM
! ================================================================

  implicit none

  private
  public :: MapTriQ

  contains

! ==================================================
  subroutine MapTriQ(p)
! ==================================================
! Generate a grid.i.triq file from a previously
! generated usurp.map
! 
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,pathname,showcpu,workdir
  use VertexData  ,only: iassoc,new_nodes,num_nodes
  use VertexData  ,only: used,xv2,nqv,num_used,xv,qv

  implicit none

  type(panel),pointer,dimension(:)       :: p
  real(kind=rd),allocatable,dimension(:) :: qt
  real(kind=rd),dimension(3)             :: xt
  real(kind=rd)                          :: wp,wv
  real(kind=rs)                          :: time1,time2
  integer,allocatable,dimension(:)       :: nt
  integer                                :: ierr,iv,num_nodes_prev
  integer                                :: num_elements,num_vertices
  integer                                :: num_used_prev,iqv,ipan,nv
  logical                                :: ex,have_qv


  continue


! initialization
  call my_cpu_time(time1)
  have_qv = allocated(qv)

  inquire(file=trim(workdir)//"usurp.map",exist=ex)
  if (.not.ex) then
    write(0,*)
    write(0,*)"ERROR! MapTriQ could not find an existing usurp.map file."
    stop
  end if

  inquire(file=trim(workdir)//"grid.i.tri",exist=ex)
  if (.not.ex) then
    write(0,*)
    write(0,*)"ERROR! MapTriQ could not find an existing grid.i.tri file."
    stop
  end if


! read the usurp.map info
  open(unit=49,file=trim(workdir)//"usurp.map",form="unformatted", &
       status="old",action="read")
  read(49)num_nodes_prev,new_nodes
  if (num_nodes_prev /= num_nodes) then
    write(0,*)
    write(0,*)"ERROR! usurp.map has num_nodes = ",num_nodes_prev
    write(0,*)"but current assembly has num_nodes = ",num_nodes
    stop
  end if
  num_vertices = num_nodes + new_nodes

  allocate(used(num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate used in MapTriQ."
    write(0,*)"num_vertices = ",num_vertices
    stop
  end if

  allocate(xv2(3,num_nodes+1:num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate xv2 in MapTriQ."
    write(0,*)"num_nodes,num_vertices = ",num_nodes,num_vertices
    stop
  end if

  allocate(iassoc(num_nodes+1:num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate iassoc in MapTriQ."
    write(0,*)"num_nodes,num_vertices = ",num_nodes,num_vertices
    stop
  end if

  read(49)used
  do iv = num_nodes + 1,num_vertices
    if (used(iv)) read(49)xv2(:,iv),iassoc(iv)
  end do
  close(49)


  num_used = 0
  do iv = 1,num_vertices
    if (used(iv)) num_used = num_used + 1
  end do


! read the grid.i.tri file
  open(unit=48,file=trim(workdir)//"grid.i.tri",action="read", &
       status="old")
  read(48,*)num_used_prev,num_elements

  if (num_used_prev /= num_used) then
    write(0,*)
    write(0,*)"ERROR! usurp.map has num_used = ",num_used
    write(0,*)"but grid.i.tri has num_used = ",num_used_prev
    stop
  end if

  open(unit=49,file=trim(pathname)//"grid.i.triq",action="write")
  write(49,*)num_used,num_elements,nqv

  do iv = 1,num_used
    read(48,*)xt(:)
    write(49,*)real(xt(:))
  end do   

  allocate(nt(3)); nalloc = nalloc + 1

  do iv = 1,num_elements
    read(48,*)nt(:)
    write(49,*)nt(:)
  end do

  deallocate(nt) ; nalloc = nalloc - 1
  allocate(nt(1)); nalloc = nalloc + 1

  do iv = 1,num_elements
    read(48,*)nt(1)
    write(49,*)nt(1)
  end do

  deallocate(nt); nalloc = nalloc - 1

  if (have_qv) then

!   write dependent variable data for existing points
    do iv = 1,num_nodes
       if (used(iv)) write(49,*)(real(qv(iqv,iv)),iqv=1,nqv)
    end do

    allocate(qt(nqv)); nalloc = nalloc + 1

!   generate and write dependent variable data for created points
    do iv = num_nodes + 1,num_vertices
      if (used(iv)) then

!       fill in qt using a simple inverse distance interpolation
!       and the corners of the original panel
!       Note: we should probably "i-blank" the corners
        qt(:) = 0.0_rd
        wp = 0.0_rd
        ipan = iassoc(iv)
        do nv = 1,p(ipan)%num_vertices
          wv = (xv2(1,iv) - xv(1,p(ipan)%node(nv)))**2 &
             + (xv2(2,iv) - xv(2,p(ipan)%node(nv)))**2 &
             + (xv2(3,iv) - xv(3,p(ipan)%node(nv)))**2
          if (wv == 0.0_rd) then
            wp = 1.0_rd
            qt(:) = qv(:,p(ipan)%node(nv))
            exit
          else
            wp = wp + 1.0_rd/wv
            qt(:) = qt(:) + qv(:,p(ipan)%node(nv))/wv
          end if
        end do
        qt(:) = qt(:)/wp
        write(49,*)(real(qt(iqv)),iqv=1,nqv)
      end if
    end do

    deallocate(qt); nalloc = nalloc - 1

  else
    write(0,*)
    write(0,*)"ERROR! qv data is not present for grid.i.triq in MapTriQ"
    stop
  end if
  
  close(unit=48)
  close(unit=49)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in MapTriQ: ',time2-time1
  end if

  return
  end subroutine MapTriQ

  end module MapTriQM
