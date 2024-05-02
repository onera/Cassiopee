! ================================================================
  module CreateNewNodeListM
! ================================================================

  implicit none

  private
  public :: CreateNewNodeList

  contains

! ======================================================================
  subroutine CreateNewNodeList(p)
! ======================================================================
! This procedure creates a list of nodes that were created during
! the polygon clipping or triangulation by GPC.
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use RotateBackM ,only: RotateBack
  use Types       ,only: panel,polygon
  use UserInput   ,only: cpufmt,showcpu
  use UserInput   ,only: verbosity
  use VertexData  ,only: nqv,num_nodes,qv,xv,xv2,qv2,new_nodes
  use VertexData  ,only: iassoc,xv2o

  implicit none

  type(panel),pointer,dimension (:) :: p
  real(kind=rd),dimension(3)        :: xt

  real(kind=rd) :: wp,wv
  real(kind=rs) :: time1,time2
  integer       :: ip1,ip2,ipan,ic,iv,jv,nv
  logical       :: have_qv


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  have_qv = allocated(qv)


! count new nodes (node <= 0)
  do ipan = ip1,ip2
  do ic = 1,p(ipan)%poly%num_cont
  do iv = 1,p(ipan)%poly%cont(ic)%num_vert

    if (p(ipan)%poly%cont(ic)%node(iv) <= 0) then

      new_nodes = new_nodes + 1

!     Copy the new node number (num_nodes + new_nodes) into
!     the polygon vertex node number
      p(ipan)%poly%cont(ic)%node(iv) = num_nodes + new_nodes

    end if
  end do !iv
  end do !ic
  end do !ipan

  if (verbosity > 1) then
    write(0,"(1x,a,i6,a)")"CreateNewNodeList found ",new_nodes," new nodes"
  end if


! Create the 3D coordinate data for any new nodes

  allocate(xv2(3,num_nodes+1:num_nodes+new_nodes)) ; nalloc = nalloc + 1
  allocate(xv2o(3,num_nodes+1:num_nodes+new_nodes)); nalloc = nalloc + 1
  allocate(iassoc(num_nodes+1:num_nodes+new_nodes)); nalloc = nalloc + 1

  if (have_qv) then
    allocate(qv2(nqv,num_nodes+1:num_nodes+new_nodes)); nalloc = nalloc + 1
  end if


  do ipan = ip1,ip2
  do ic = 1,p(ipan)%poly%num_cont
  do iv = 1,p(ipan)%poly%cont(ic)%num_vert
    jv = p(ipan)%poly%cont(ic)%node(iv) 
    if (jv > num_nodes) then
!     Rotate the node back into 3D
      call RotateBack(p(ipan)%poly%cont(ic)%x(iv), &
                      p(ipan)%poly%cont(ic)%y(iv), &
                      p(ipan)%g,p(ipan)%xmid,xt)
      xv2(1:3,jv) = xt(1:3)
      xv2o(1:3,jv) = xt(1:3)
      iassoc(jv) = ipan

      if (have_qv) then
!       Fill in qv2 using a simple inverse distance interpolation
!       and the corners of the original panel
!       Note: we should probably "i-blank" the corners
        qv2(:,jv) = 0.0_rd
        wp = 0.0_rd
        do nv = 1,p(ipan)%num_vertices
          wv = (xv2(1,jv) - xv(1,p(ipan)%node(nv)))**2 &
             + (xv2(2,jv) - xv(2,p(ipan)%node(nv)))**2 &
             + (xv2(3,jv) - xv(3,p(ipan)%node(nv)))**2
          if (wv == 0.0_rd) then
            wp = 1.0_rd
            qv2(1:nqv,jv) = qv(1:nqv,p(ipan)%node(nv))
            exit
          else
            wp = wp + 1.0_rd/wv
            qv2(:,jv) = qv2(:,jv) + qv(:,p(ipan)%node(nv))/wv
          end if
        end do
        qv2(:,jv) = qv2(:,jv)/wp
      end if !have_qv

    end if
  end do !iv
  end do !ic
  end do !ipan


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CreateNewNodeList: ',time2-time1
  end if

  return
  end subroutine CreateNewNodeList

  end module CreateNewNodeListM
