! ================================================================
  module ShrinkVertexListM
! ================================================================

  implicit none

  private
  public  :: ShrinkVertexList
  private :: NodeID

  contains

! ================================================================
  subroutine ShrinkVertexList(p)
! ================================================================
! For structured patch grids, reset the vertex numbering to the
! earliest occurence of a grid point in the global vertex 
! numbering system.  Note that this is done by grid index (i,j,k,n)
! not coordinate (x,y,z), so this routine does NOT merge points
! that are merely coincident in space, except in some special cases
! (axis, periodic, and wake-cut boundaries within a given subset).  
! Also note that the dummy argument is a pointer, so this procedure 
! requires an explicit interface.

  use IntrType    ,only: rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: num_patches,block,i1,i2,j1,j2,k1,k2
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,showcpu
  use VertexData  ,only: xv,num_nodes,alias

  implicit none

  type(panel),pointer,dimension(:) :: p
  type(panel),pointer              :: q

  real(kind=rs)     :: time1,time2
  integer           :: imin,imax,jmin,jmax,kmin,kmax,iv1,iv2,i,j,k
  integer           :: ipan,iv,ip1,ip2,imid,ii,ierr


  continue


  call my_cpu_time(time1)


! allocate memory for the vertex aliases
  allocate(alias(num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate alias in ShrinkVertexList"
    write(0,*)"num_nodes = ",num_nodes
    stop
  end if

! initialize the list
  do i = 1,num_nodes
    alias(i) = i
  end do


! Resolve shared vertices between subsets sharing a common edge
! (split domains)
  do ip1 = 2,num_patches
    do ip2 = ip1-1,1,-1
      if (block(ip1) == block(ip2)) then
        imin = max(i1(ip1),i1(ip2))
        imax = min(i2(ip1),i2(ip2))
        jmin = max(j1(ip1),j1(ip2))
        jmax = min(j2(ip1),j2(ip2))
        kmin = max(k1(ip1),k1(ip2))
        kmax = min(k2(ip1),k2(ip2))
        
        if ((imax >= imin).and.(jmax >= jmin).and.(kmax >= kmin)) then
!         These patches overlap. Set the nodes in patch ip1 back to the
!         node numbers for patch ip2
          if (k1(ip1) == k2(ip1)) then
            k = k1(ip1)
            do j = jmin,jmax
            do i = imin,imax

              iv1 = NodeID(i,j,k,ip1)
              iv2 = NodeID(i,j,k,ip2)
              alias(iv1) = min(alias(iv1),iv2)

            end do
            end do
          else if (i1(ip1) == i2(ip1)) then
            i = i1(ip1)
            do k = kmin,kmax
            do j = jmin,jmax

              iv1 = NodeID(i,j,k,ip1)
              iv2 = NodeID(i,j,k,ip2)
              alias(iv1) = min(alias(iv1),iv2)

            end do
            end do
          else if (j1(ip1) == j2(ip1)) then
            j = j1(ip1)
            do k = kmin,kmax
            do i = imin,imax

              iv1 = NodeID(i,j,k,ip1)
              iv2 = NodeID(i,j,k,ip2)
              alias(iv1) = min(alias(iv1),iv2)

            end do
            end do
          end if
        end if
      end if
    end do
  end do


! SPECIAL BOUNDARY TYPES
  do ip1 = 1,num_patches

!   I-SURFACES
    if (i1(ip1) == i2(ip1)) then
      i = i1(ip1)

!     REDUCE ANY WAKE CUT BOUNDARIES
      if (mod(j2(ip1)-j1(ip1)+1,2) == 1) then
        do k = k1(ip1),k2(ip1),k2(ip1)-k1(ip1)
          imid = (j1(ip1)+j2(ip1))/2
          do ii = j1(ip1),imid-1

            j = ii
            iv1 = NodeID(i,j,k,ip1)

            j = j2(ip1) - ii + j1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

      if (mod(k2(ip1)-k1(ip1)+1,2) == 1) then
        do j = j1(ip1),j2(ip1),j2(ip1)-j1(ip1)
          imid = (k1(ip1)+k2(ip1))/2
          do ii = k1(ip1),imid-1

            k = ii
            iv1 = NodeID(i,j,k,ip1)

            k = k2(ip1) - ii + k1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

!     REDUCE ANY AXIS BOUNDARIES (SINGULARITIES)
      do k = k1(ip1),k2(ip1),k2(ip1)-k1(ip1)
        j = j1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do j = j1(ip1)+1,j2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do
      end do  !k-loop

      do j = j1(ip1),j2(ip1),j2(ip1)-j1(ip1)
        k = k1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do k = k1(ip1)+1,k2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do

      end do  !j-loop

!     REDUCE ANY PERIODIC BOUNDARIES (J-DIRECTION)
      do k = k1(ip1),k2(ip1)

        iv1 = NodeID(i,j1(ip1),k,ip1)
        iv2 = NodeID(i,j2(ip1),k,ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do

!     REDUCE ANY PERIODIC BOUNDARIES (K-DIRECTION)
      do j = j1(ip1),j2(ip1)

        iv1 = NodeID(i,j,k1(ip1),ip1)
        iv2 = NodeID(i,j,k2(ip1),ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do

    end if


!   J-SURFACES
    if (j1(ip1) == j2(ip1)) then
      j = j1(ip1)

!     REDUCE ANY WAKE CUT BOUNDARIES
      if (mod(k2(ip1)-k1(ip1)+1,2) == 1) then
        do i = i1(ip1),i2(ip1),i2(ip1)-i1(ip1)
          imid = (k1(ip1)+k2(ip1))/2
          do ii = k1(ip1),imid-1

            k = ii
            iv1 = NodeID(i,j,k,ip1)

            k = k2(ip1) - ii + k1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

      if (mod(i2(ip1)-i1(ip1)+1,2) == 1) then
        do k = k1(ip1),k2(ip1),k2(ip1)-k1(ip1)
          imid = (i1(ip1)+i2(ip1))/2
          do ii = i1(ip1),imid-1

            i = ii
            iv1 = NodeID(i,j,k,ip1)

            i = i2(ip1) - ii + i1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

!     REDUCE ANY AXIS BOUNDARIES (SINGULARITIES)
      do i = i1(ip1),i2(ip1),i2(ip1)-i1(ip1)
        k = k1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do k = k1(ip1)+1,k2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do
      end do  !i-loop

      do k = k1(ip1),k2(ip1),k2(ip1)-k1(ip1)
        i = i1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do i = i1(ip1)+1,i2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do

      end do  !k-loop

!     REDUCE ANY PERIODIC BOUNDARIES (K-DIRECTION)
      do i = i1(ip1),i2(ip1)

        iv1 = NodeID(i,j,k1(ip1),ip1)
        iv2 = NodeID(i,j,k2(ip1),ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do

!     REDUCE ANY PERIODIC BOUNDARIES (I-DIRECTION)
      do k = k1(ip1),k2(ip1)

        iv1 = NodeID(i1(ip1),j,k,ip1)
        iv2 = NodeID(i2(ip1),j,k,ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do

    end if

  
!   K-SURFACES
    if (k1(ip1) == k2(ip1)) then
      k = k1(ip1)

!     REDUCE ANY WAKE CUT BOUNDARIES
      if (mod(i2(ip1)-i1(ip1)+1,2) == 1) then
        do j = j1(ip1),j2(ip1),j2(ip1)-j1(ip1)
          imid = (i1(ip1)+i2(ip1))/2
          do ii = i1(ip1),imid-1

            i = ii
            iv1 = NodeID(i,j,k,ip1)

            i = i2(ip1) - ii + i1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

      if (mod(j2(ip1)-j1(ip1)+1,2) == 1) then
        do i = i1(ip1),i2(ip1),i2(ip1)-i1(ip1)
          imid = (j1(ip1)+j2(ip1))/2
          do ii = j1(ip1),imid-1

            j = ii
            iv1 = NodeID(i,j,k,ip1)

            j = j2(ip1) - ii + j1(ip1)
            iv2 = NodeID(i,j,k,ip1)

            if ((xv(1,iv2) == xv(1,iv1)).and. &
                (xv(2,iv2) == xv(2,iv1)).and. &
                (xv(3,iv2) == xv(3,iv1))) then

              alias(iv2) = min(alias(iv2),iv1)

            end if
          end do
        end do
      end if

!     REDUCE ANY AXIS BOUNDARIES (SINGULARITIES)
      do j = j1(ip1),j2(ip1),j2(ip1)-j1(ip1)
        i = i1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do i = i1(ip1)+1,i2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do
      end do  !j-loop

      do i = i1(ip1),i2(ip1),i2(ip1)-i1(ip1)
        j = j1(ip1)
        iv1 = NodeID(i,j,k,ip1)
        do j = j1(ip1)+1,j2(ip1)
          iv2 = NodeID(i,j,k,ip1)
          if ((xv(1,iv2) == xv(1,iv1)).and. &
              (xv(2,iv2) == xv(2,iv1)).and. &
              (xv(3,iv2) == xv(3,iv1))) then

            alias(iv2) = min(alias(iv2),iv1)

          end if
        end do

      end do  !i-loop

!     REDUCE ANY PERIODIC BOUNDARIES (I-DIRECTION)
      do j = j1(ip1),j2(ip1)

        iv1 = NodeID(i1(ip1),j,k,ip1)
        iv2 = NodeID(i2(ip1),j,k,ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do

!     REDUCE ANY PERIODIC BOUNDARIES (J-DIRECTION)
      do i = i1(ip1),i2(ip1)

        iv1 = NodeID(i,j1(ip1),k,ip1)
        iv2 = NodeID(i,j2(ip1),k,ip1)

        if ((xv(1,iv2) == xv(1,iv1)).and. &
            (xv(2,iv2) == xv(2,iv1)).and. &
            (xv(3,iv2) == xv(3,iv1))) then

          alias(iv2) = min(alias(iv2),iv1)

        end if
        
      end do


    end if  !k-surface
  end do !patch loop


! condense the list
  do i = 1,num_nodes
    j = alias(i)
    do while (j /= alias(j))
      j = alias(j)
    end do
    alias(i) = j
  end do


! reassign the node names of all panels
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  do ipan = ip1,ip2
    q => p(ipan)
    do iv = 1,q%num_vertices
      q%node(iv) = alias(q%node(iv))
    end do
  end do

      
  deallocate(alias); nalloc = nalloc - 1


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in ShrinkVertexList:',time2-time1
  end if


  return
  end subroutine ShrinkVertexList

! ================================================================
  pure function NodeID(i,j,k,ip1) result (iv1)
! ================================================================
! convert the ijk indices on a structured patch to the 1d vertex
! node number

  use PatchInfo,only: i1,i2,j1,j2,k1,k2
  use Types    ,only: surf

  implicit none

  integer,intent(in) :: i,j,k,ip1
  integer            :: iv1


  continue


  iv1 = surf(ip1)%first_node + (i - i1(ip1)) &
                             + (j - j1(ip1))*(i2(ip1)-i1(ip1)+1) &
                             + (k - k1(ip1))*(j2(ip1)-j1(ip1)+1) &
                                            *(i2(ip1)-i1(ip1)+1)

  return
  end function NodeID

  end module ShrinkVertexListM
