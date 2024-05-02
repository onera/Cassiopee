! ================================================================
  module ResolveEdgeIntersectionsM
! ================================================================

  implicit none

  private
  public  :: ResolveEdgeIntersections
  private :: SegmentDistance
  private :: CheckInsertion

  contains

! ======================================================================
  subroutine ResolveEdgeIntersections(p)
! ======================================================================
! This procedure attempts to create a fully watertight surface by
! resolving any hanging nodes that were created during the polygon 
! clipping.
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use EdgeData     ,only: edge
  use IntrType     ,only: rd,rs
  use my_cpu_timeM ,only: my_cpu_time
  use NetAlloc     ,only: nalloc
  use RotateBackM  ,only: RotateBack
  use sort2M       ,only: sort2
  use Types        ,only: panel,polygon
  use UserInput    ,only: cpufmt,showcpu,pathname
  use UserInput    ,only: verbosity
  use VertexData   ,only: nqv,num_nodes,qv,xv,xv2,qv2,xv2o
  use VertexData   ,only: new_nodes,iassoc

  implicit none

  type(panel),pointer,dimension (:)    :: p
  type(panel),pointer                  :: q

  real(kind=rd),dimension(3) :: xt,x1,x2,x3,x4
  real(kind=rd)              :: wp,wv,sc,tc

  real(kind=rs) :: time1,time2,time3,time4,time5,time6
  real(kind=rs) :: bltime,rtime,itime,stime

  integer,allocatable,dimension(:,:) :: epair
  integer :: ip1,ip2,ipan,ic,iv,jv,nv
  integer :: e1,e2,etmp,f1L,f1R,f2L,f2R
  integer :: ierr,i,j,i1,i2
  logical :: have_qv


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  have_qv = allocated(qv)


! To resolve hanging nodes, begin by building a list of the edge pairs
! that intersect to form the hanging nodes.

  call my_cpu_time(time5)


! sweep through all panels to count the number of hanging nodes (node < 0)
  new_nodes = 0
  do ipan = ip1,ip2
  do ic = 1,p(ipan)%poly%num_cont
  do iv = 1,p(ipan)%poly%cont(ic)%num_vert 
    if (p(ipan)%poly%cont(ic)%node(iv) < 0) then
      new_nodes = new_nodes + 1
    end if
  end do
  end do
  end do
  if (new_nodes == 0) return


! create space for a list of the edge pairs
  allocate(epair(2,new_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate epair"
    write(0,*)"new_nodes = ",new_nodes
    stop
  end if


! store the edge pairs associated with any new nodes
  new_nodes = 0
  do ipan = ip1,ip2
  do ic = 1,p(ipan)%poly%num_cont
  do iv = 1,p(ipan)%poly%cont(ic)%num_vert 
    if (p(ipan)%poly%cont(ic)%node(iv) < 0) then
      new_nodes = new_nodes + 1
      e1 = p(ipan)%poly%cont(ic)%ie1(iv)
      if (e1 == 0) then
        write(0,*)"ERROR! intersection edge info e1 = 0"
        write(0,*)"ipan,ic,iv = ",ipan,ic,iv
        stop
      end if
      e2 = p(ipan)%poly%cont(ic)%ie2(iv)
      if (e2 == 0) then
        write(0,*)"ERROR! intersection edge info e2 = 0"
        write(0,*)"ipan,ic,iv = ",ipan,ic,iv
        stop
      end if
      if (e2 < e1) then
        etmp = e2
        e2 = e1
        e1 = etmp
      end if 
      epair(1,new_nodes) = e1
      epair(2,new_nodes) = e2
    end if
  end do
  end do
  end do


! sort the edge pairs
  call my_cpu_time(time3)

  call sort2(new_nodes,epair)

  goto 10
  do i = 1,new_nodes-1
  do j = i+1,new_nodes
    if (epair(1,j) < epair(1,i)) then
      e1 = epair(1,j)
      e2 = epair(2,j)
      epair(1,j) = epair(1,i)
      epair(2,j) = epair(2,i)
      epair(1,i) = e1
      epair(2,i) = e2
    else if (epair(1,j) == epair(1,i)) then
      if (epair(2,j) < epair(2,i)) then
        e2 = epair(2,j)
        epair(2,j) = epair(2,i)
        epair(2,i) = e2
      end if
    end if
  end do
  end do
  10 continue

  call my_cpu_time(time4)
  stime = time4-time3


! remove any duplicate entries
  j = 1
  do i = 2,new_nodes
    if (epair(1,i) /= epair(1,j)) then
      j = j + 1
      epair(1,j) = epair(1,i)
      epair(2,j) = epair(2,i)
    else if (epair(2,i) /= epair(2,j)) then
      j = j + 1
      epair(1,j) = epair(1,i)
      epair(2,j) = epair(2,i)
    end if
  end do
  new_nodes = j
       

! sweep through the panels again: this time, associate with each hanging
! node a unique number representing the edge intersection that caused it
  do ipan = ip1,ip2
  do ic = 1,p(ipan)%poly%num_cont
  do iv = 1,p(ipan)%poly%cont(ic)%num_vert 

    if (p(ipan)%poly%cont(ic)%node(iv) < 0) then

      e1 = p(ipan)%poly%cont(ic)%ie1(iv)
      e2 = p(ipan)%poly%cont(ic)%ie2(iv)
      if (e2 < e1) then
        etmp = e2
        e2 = e1
        e1 = etmp
      end if 

!     search the list for this pair
      i1 = 1
      i2 = new_nodes

      do
        if ((e1 == epair(1,i1)).and.(e2 == epair(2,i1))) then
          j = i1
          exit
        else if ((e1 == epair(1,i2)).and.(e2 == epair(2,i2))) then
          j = i2
          exit
        else
          j = (i1 + i2)/2
        end if
  
        if (e1 < epair(1,j)) then
          i2 = j
        else if (e1 > epair(1,j)) then
          i1 = j
        else
          if (e2 < epair(2,j)) then
            i2 = j
          else if (e2 > epair(2,j)) then
            i1 = j
          else
            exit
          end if
        end if
      end do

      p(ipan)%poly%cont(ic)%node(iv) = num_nodes + j
    end if

  end do
  end do
  end do

  call my_cpu_time(time6)
  bltime = time6-time5-stime



! Allocate array space for the new 3D data
  call my_cpu_time(time3)

  allocate(xv2(3,num_nodes+1:num_nodes+new_nodes)) ; nalloc = nalloc + 1
  allocate(iassoc(num_nodes+1:num_nodes+new_nodes)); nalloc = nalloc + 1
  allocate(xv2o(3,num_nodes+1:num_nodes+new_nodes)); nalloc = nalloc + 1

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
  call my_cpu_time(time4)
  rtime = time4-time3


! Go through the edge intersection list and figure out which 
! faces need to be visited for possible node insertions
  call my_cpu_time(time3)

  if (verbosity >= 2) then
    open(unit=38,file=trim(pathname)//"edge_list")
    write(38,*)"new_nodes = ",new_nodes
  end if

  do i = 1,new_nodes

    e1  = epair(1,i)
    e2  = epair(2,i)

    x1(:) = xv(:,edge(e1)%vertex(1))
    x2(:) = xv(:,edge(e1)%vertex(2))
    x3(:) = xv(:,edge(e2)%vertex(1))
    x4(:) = xv(:,edge(e2)%vertex(2))

    call SegmentDistance(x1,x2,x3,x4,xt,sc,tc) 

    xv2(:,num_nodes + i) = xt(:)

    f1L = edge(e1)%left_face
    if (f1L > 0) then
      q => p(f1L)
      call CheckInsertion(q,e1,sc,num_nodes+i)
    end if

    f1R = edge(e1)%right_face
    if (f1R > 0) then
      q => p(f1R)
      call CheckInsertion(q,e1,sc,num_nodes+i)
    end if

    f2L = edge(e2)%left_face
    if (f2L > 0) then
      q => p(f2L)
      call CheckInsertion(q,e2,tc,num_nodes+i)
    end if

    f2R = edge(e2)%right_face
    if (f2R > 0) then
      q => p(f2R)
      call CheckInsertion(q,e2,tc,num_nodes+i)
    end if

    if (verbosity >= 2) write(38,"(7(1x,i7))")num_nodes+i,e1,e2,f1L,f1R,f2L,f2R

  end do
  if (verbosity >= 2) close(38)
  call my_cpu_time(time4)
  itime = time4-time3


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in ResolveEdgeIntersections: ',time2-time1
    write(0,cpufmt)'  > time to build intersection list: ',bltime
    write(0,cpufmt)'  > time to sort list: ',stime
    write(0,cpufmt)'  > time to build 3D data: ',rtime
    write(0,cpufmt)'  > time to insert connecting nodes: ',itime
  end if


  deallocate(epair); nalloc = nalloc - 1

  return
  end subroutine ResolveEdgeIntersections

! =============================================================================
  subroutine SegmentDistance(xa,xb,xc,xd,xi,sc,tc)
! =============================================================================
! dist3D_segment_to_segment():
! Input: two 3D line segments S1 and S2
! Return: the shortest distance between S1 and S2
! source: http://softsurfer.com/Archive/algorithm_0106/...
!          .../algorithm_0106.htm#dist3D_Segment_to_Segment()
! converted from C to Fortran

  use IntrType,only: rd

  implicit none

  real(kind=rd),intent(in),dimension(3)  :: xa,xb,xc,xd
  real(kind=rd),intent(out),dimension(3) :: xi
  real(kind=rd),intent(out)              :: sc,tc
  real(kind=rd),parameter   :: SMALL_NUM = 1.0e-15

  real(kind=rd) :: ux,uy,uz
  real(kind=rd) :: vx,vy,vz
  real(kind=rd) :: wx,wy,wz
  real(kind=rd) :: a,b,c,d,e,DD
  real(kind=rd) :: sN,sD
  real(kind=rd) :: tN,tD


  continue


  ux = xb(1) - xa(1)
  uy = xb(2) - xa(2)
  uz = xb(3) - xa(3)
  vx = xd(1) - xc(1)
  vy = xd(2) - xc(2)
  vz = xd(3) - xc(3)
  wx = xa(1) - xc(1)
  wy = xa(2) - xc(2)
  wz = xa(3) - xc(3)

  a = ux*ux + uy*uy + uz*uz
  b = ux*vx + uy*vy + uz*vz 
  c = vx*vx + vy*vy + vz*vz
  d = ux*wx + uy*wy + uz*wz
  e = vx*wx + vy*wy + vz*wz

  DD = a*c - b*b
  sD = DD
  tD = DD

! compute the line parameters of the two closest points

  if (DD < SMALL_NUM) then
!   the lines are almost parallel
    sN = 0.0_rd ! force using point P0 on segment S1
    sD = 1.0_rd ! to prevent possible division by 0.0 later
    tN = e
    tD = c
  else
!   get the closest points on the infinite lines
    sN = (b*e - c*d)
    tN = (a*e - b*d)
  end if


  if (sN < 0.0_rd) then 
!   sc < 0 => the s=0 edge is visible
    sN = 0.0_rd
    tN = e
    tD = c
  else if (sN > sD) then 
!   sc > 1 => the s=1 edge is visible
    sN = sD
    tN = e + b
    tD = c
  end if

  if (tN < 0.0_rd) then 
  ! tc < 0 => the t=0 edge is visible
    tN = 0.0_rd
    ! recompute sc for this edge
    if (-d < 0.0_rd) then
      sN = 0.0_rd
    else if (-d > a) then
      sN = sD
    else 
      sN = -d
      sD = a
    end if
  else if (tN > tD) then
!   tc > 1 => the t=1 edge is visible
    tN = tD
    ! recompute sc for this edge
    if ((-d + b) < 0.0_rd) then
      sN = 0.0_rd
    else if ((-d + b) > a) then
      sN = sD
    else
      sN = (-d + b)
      sD = a
    end if
  end if

! finally do the division to get sc and tc
  if (abs(SN) < SMALL_NUM) then
    sc = 0.0_rd
  else
    sc = sN / sD
  end if

  if (abs(tN) < SMALL_NUM) then
    tc = 0.0_rd
  else
    tc = tN / tD
  end if

  if (a < c) then
    xi(1) = xa(1) + sc*ux
    xi(2) = xa(2) + sc*uy
    xi(3) = xa(3) + sc*uz
  else if (c < a) then
    xi(1) = xc(1) + tc*vx
    xi(2) = xc(2) + tc*vy
    xi(3) = xc(3) + tc*vz
  else
    xi(1) = 0.5_rd*(xa(1) + sc*ux + xc(1) + tc*vx)
    xi(2) = 0.5_rd*(xa(2) + sc*uy + xc(2) + tc*vy)
    xi(3) = 0.5_rd*(xa(3) + sc*uz + xc(3) + tc*vz)
  end if


  return
  end subroutine SegmentDistance

! ==============================================================
  subroutine CheckInsertion(q,e2,tc,m)
! ==============================================================
! check to see whether a given point (associated with node m
! at parametric position tc) along edge e2 needs to be inserted
! into panel q.

  use EdgeData   ,only: edge
  use InsertNodeM,only: InsertNode
  use IntrType   ,only: rd
  use RotateBackM,only: RotateBack
  use Types      ,only: panel,polygon
  use UserInput  ,only: plotpanel
  use VertexData ,only: xv

  implicit none

  type(panel),pointer :: q 
  type(polygon),pointer :: spoly

  real(kind=rd),intent(in) :: tc
  real(kind=rd),dimension(3) :: x1,x2,xt1,xt2
  real(kind=rd) :: tc1,tc2,ux,uy,uz,uv,vx,vy,vz,vv
  real(kind=rd) :: ss,xh,yh

  integer,intent(in) :: e2,m
  integer :: ic,iv,nv,span


  continue


  span = q%panel_number


! loop through the vertices of the polygon to see
! whether the node "m" is already present
  do ic = 1,q%poly%num_cont
  do iv = 1,q%poly%cont(ic)%num_vert
    if (q%poly%cont(ic)%node(iv) == m) then
      if (plotpanel == span) then
        write(*,*)"CheckInsertion: m already exists.",m
      end if
      return
    end if
  end do
  end do


! loop through the edges of the polygon to see 
! whether the edge "e2" is present
  do ic = 1,q%poly%num_cont
  do iv = 1,q%poly%cont(ic)%num_vert
    if (q%poly%cont(ic)%edge(iv) == e2) then

      if (plotpanel == span) then
        write(*,*)"CheckInsertion: e2 is present in m ",e2,m
        write(*,*)"check for tc = ",tc
      end if

      if (iv < q%poly%cont(ic)%num_vert) then
        nv = iv + 1
      else 
        nv = 1
      end if

!     we need to figure out whether tc is a member of this edge
      x1(:) = xv(:,edge(e2)%vertex(1))
      x2(:) = xv(:,edge(e2)%vertex(2))

      vx = x2(1)-x1(1)
      vy = x2(2)-x1(2)
      vz = x2(3)-x1(3)
      vv = vx*vx + vy*vy + vz*vz
      if (vv == 0.0_rd) then
        write(0,*)"ERROR! vv = 0.0 in CheckInsertion; unexpected."
        return
      end if

!     Rotate the end points back into 3D
      call RotateBack(q%poly%cont(ic)%x(iv), &
                      q%poly%cont(ic)%y(iv), &
                      q%g,q%xmid,xt1)

      ux  = xt1(1)-x1(1)
      uy  = xt1(2)-x1(2)
      uz  = xt1(3)-x1(3)
      uv  = ux*vx + uy*vy + uz*vz
      tc1 = uv / vv

      call RotateBack(q%poly%cont(ic)%x(nv), &
                      q%poly%cont(ic)%y(nv), &
                      q%g,q%xmid,xt2)

      ux  = xt2(1)-x1(1)
      uy  = xt2(2)-x1(2)
      uz  = xt2(3)-x1(3)
      uv  = ux*vx + uy*vy + uz*vz
      tc2 = uv / vv

      if (((tc1 < tc).and.(tc < tc2)).or. &
          ((tc2 < tc).and.(tc < tc1))) then
!       insert node
        ss = (tc - tc1) / (tc2 - tc1)
        xh = q%poly%cont(ic)%x(iv) + ss*(q%poly%cont(ic)%x(nv) &
                                       - q%poly%cont(ic)%x(iv))
        yh = q%poly%cont(ic)%y(iv) + ss*(q%poly%cont(ic)%y(nv) &
                                       - q%poly%cont(ic)%y(iv))
        spoly => q%poly
        call InsertNode(spoly,ic,iv,span,xh,yh,m,e2)
      end if

    end if
  end do
  end do

  return
  end subroutine CheckInsertion

  end module ResolveEdgeIntersectionsM
