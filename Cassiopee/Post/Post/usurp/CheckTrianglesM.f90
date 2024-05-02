! ================================================================
  module CheckTrianglesM
! ================================================================

  implicit none

  private
  public  :: CheckTriangles
  private :: indexx

  contains

! ======================================================================
  subroutine CheckTriangles(p)
! ======================================================================
! This procedure checks the triangulated surface data for 
! edges belonging to more than one triangle
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use RotateBackM ,only: RotateBack
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,pathname,showcpu

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: stime,etime

  integer,allocatable,dimension(:,:) :: myv
  integer,allocatable,dimension(:) :: edge_element,indx,edge_refs
  integer :: num_elements,num_edges
  integer :: ip1,ip2,ipan,ic,iv,iel,ierr,iedge
  integer :: i,j,k,m

  character(len=17),allocatable,dimension(:) :: edge_name

  logical :: first


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  first = .true.

  num_elements = 0
  do ipan = ip1,ip2
    do ic = 1,p(ipan)%poly%num_cont
      num_elements = num_elements + 1
    end do
  end do

  allocate(myv(4,num_elements),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate myv in CheckTriangles"
    write(0,*)"num_elements = ",num_elements
    stop
  end if


! generate element list
  iel = 1
  do ipan = ip1,ip2
    do ic = 1,p(ipan)%poly%num_cont
      iv = 1
      myv(1,iel) = p(ipan)%poly%cont(ic)%node(iv)
      myv(2,iel) = p(ipan)%poly%cont(ic)%node(iv+1)
      myv(3,iel) = p(ipan)%poly%cont(ic)%node(iv+2)
      myv(4,iel) = ipan
      iel = iel + 1
    end do
  end do


! generate an edge list
  num_edges = 3 * num_elements

  allocate(edge_element(num_edges),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate edge_element in CheckTriangles."
    write(0,*)"num_edges = ",num_edges
    stop
  end if

  allocate(edge_refs(num_edges),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate edge_refs in CheckTriangles."
    write(0,*)"num_edges = ",num_edges
    stop
  end if
  edge_refs = 0

  allocate(edge_name(num_edges),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate edge_name in CheckTriangles."
    write(0,*)"num_edges = ",num_edges
    stop
  end if

  call my_cpu_time(time3)
  iedge = 0
  do iel = 1,num_elements
    do j = 1,3
      k = mod(j,3) + 1
      iedge = iedge + 1

      if (myv(j,iel) < myv(k,iel)) then
        write(edge_name(iedge),"(i8,a1,i8)")myv(j,iel),"-",myv(k,iel)
      else
        write(edge_name(iedge),"(i8,a1,i8)")myv(k,iel),"-",myv(j,iel)
      end if

      edge_element(iedge) = iel

    end do
  end do
  call my_cpu_time(time4)
  etime = time4-time3


! sort the edge list
  allocate(indx(num_edges),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate indx in CheckTriangles"
    write(0,*)"num_edges = ",num_edges
    stop
  end if

  call my_cpu_time(time3)
  call indexx(num_edges,edge_name,indx)
  call my_cpu_time(time4)
  stime = time4 - time3


! check for edges belonging to more than two elements
  do m = 1,2
  j = 1
  do while (j <= num_edges)
    k = j
    EXTENT_LOOP: do
      if (k < num_edges) then
        if (edge_name(indx(k+1)) == edge_name(indx(j))) then
          k = k + 1
          cycle EXTENT_LOOP
        else
          exit EXTENT_LOOP
        end if
      else
        exit EXTENT_LOOP
      end if
    end do EXTENT_LOOP
    if (m == 1) then
      do i = j,k
        edge_refs(indx(i)) = k - j + 1
      end do
    else
      if ((k - j + 1) > 2) then
        if (first) then
          open(unit=38,file=trim(pathname)//"check_triangles.dat")
          first = .false.
        end if
        write(38,"(a,i2,a)")" Edge belongs to ",k-j+1," elements: "// &
                           trim(adjustl(edge_name(indx(j))))
        do i = j,k
          iel = edge_element(indx(i))
          ipan = myv(4,iel)
          write(38,*)"  element,panel: ",iel,ipan,p(ipan)%ratio
          write(38,*)edge_refs((iel-1)*3+1), &
                    edge_refs((iel-1)*3+2), &
                    edge_refs((iel-1)*3+3)
        end do
      end if
    end if
    j = k + 1
  end do
  end do

  
  deallocate(myv,edge_element,edge_name,indx,edge_refs); nalloc = nalloc - 5

  if (.not.first) close(38)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CheckTriangles: ',time2-time1
    write(0,cpufmt)'  > time in sorting',stime
    write(0,cpufmt)'  > time to generate element list ',etime
  end if

  return
  end subroutine CheckTriangles

!     ======================================================================
      subroutine indexx(n,arr,indx)
!     ======================================================================
!     indexes an array arr(1:n),i.e. outputs the array indx(1:n) such that 
!     arr(indx(j)) is in ascending order for j=1,2,..N.  The input n and arr
!     are not changed.

      implicit none

      integer,intent(in)   :: n
      integer,dimension(n) :: indx
      character(len=17),dimension(n),intent(in) :: arr

      integer,parameter :: m=7,nstack=50
      integer,dimension(nstack) :: istack
      integer :: i,indxt,ir,itemp,j,jstack,k,l

      character(len=17) :: a


      continue


      do j=1,n
         indx(j)=j
      end do

      jstack=0
      l=1
      ir=n
  1   if(ir-l.lt.M) then
         do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do i=j-1,1,-1
               if (arr(indx(i)).le.a) goto 2
               indx(i+1)=indx(i)
            enddo
            i=0
  2         indx(i+1)=indxt
         end do
         if(jstack.eq.0) return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)).gt.arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
  3      continue
            i=i+1
         if(arr(indx(i)).lt.a) goto 3
  4      continue
            j=j-1
         if(arr(indx(j)).gt.a)goto 4
         if(j.lt.i)goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
  5      indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
!        push pointers to larger subarray on stack, process smaller 
!        subarray immediately.
         if(jstack > NSTACK) stop 'NSTACK too small in indexx'
         if (ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         end if
      end if
      goto 1

      end subroutine indexx

  end module CheckTrianglesM
