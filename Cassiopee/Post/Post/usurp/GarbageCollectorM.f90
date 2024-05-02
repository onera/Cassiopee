! ================================================================
  module GarbageCollectorM
! ================================================================

  implicit none

  private
  public :: GarbageCollector

  contains

! =================================================
  subroutine GarbageCollector(spoly,pann,mode)
! =================================================
! This subroutine eliminates degenerate contours and
! duplicate vertices from polygons.  Note that the
! tolerance used here is driven by the ability of
! the polygon clipper (GPC) and Shewchuk's Triangle
! to handle degeneracies.

  use IntrType ,only: rd
  use NetAlloc ,only: nalloc
  use Types    ,only: polygon
  use UserInput,only: degen_tol,degen_tol2,watertight

  implicit none

  type(polygon),pointer :: spoly
  real(kind=rd)         :: dist
  integer,intent(in)    :: mode,pann
  integer               :: ic,jc,iv,jv,kv
  logical               :: erase_contour,check_nodes,one_edge


  continue


  check_nodes = (watertight.and.(mode == 2))

! check each contour in the polygon
! (note that contours may be deleted as we go, so
!  we need the range of the do loop to be variable)
  ic = 1
  do while (ic <= spoly%num_cont)

!   see if any points need to be merged
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert)
      jv = iv + 1
      SEARCH1: do while (jv <= spoly%cont(ic)%num_vert)
        dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
             + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2 
        if ((dist <= degen_tol2).or.(check_nodes.and.(spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)))) then
!         do kv = jv,spoly%cont(ic)%num_vert - 1
!         This mod is to properly handle the edge data during a merge
          do kv = iv,spoly%cont(ic)%num_vert - 1
            spoly%cont(ic)%x(kv)    = spoly%cont(ic)%x(kv+1)
            spoly%cont(ic)%y(kv)    = spoly%cont(ic)%y(kv+1)
            spoly%cont(ic)%node(kv) = spoly%cont(ic)%node(kv+1)
            spoly%cont(ic)%ie1(kv)  = spoly%cont(ic)%ie1(kv+1)
            spoly%cont(ic)%ie2(kv)  = spoly%cont(ic)%ie2(kv+1)
            spoly%cont(ic)%edge(kv) = spoly%cont(ic)%edge(kv+1)
          end do
          spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 1
        else
          exit SEARCH1
        end if
      end do SEARCH1
      iv = iv + 1
    end do

!   also test to be sure that the first point is not the same as the last
    if (spoly%cont(ic)%num_vert > 1) then
      iv = 1
      jv = spoly%cont(ic)%num_vert
      dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
           + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2 
      if ((dist <= degen_tol2).or.(check_nodes.and.(spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)))) then
        spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 1
      end if
    end if

!   test for the special case of point pattern ABA
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert - 1)
      jv = iv + 2
      SEARCH2: do while (jv <= spoly%cont(ic)%num_vert)
        dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
             + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2 
        if ((dist <= degen_tol2).or.(check_nodes.and.(spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)))) then
!         do kv = jv-1,spoly%cont(ic)%num_vert - 2
!         This modification is required to properly handle edge data
          do kv = jv-2,spoly%cont(ic)%num_vert - 2
            spoly%cont(ic)%x(kv)    = spoly%cont(ic)%x(kv+2)
            spoly%cont(ic)%y(kv)    = spoly%cont(ic)%y(kv+2)
            spoly%cont(ic)%node(kv) = spoly%cont(ic)%node(kv+2)
            spoly%cont(ic)%ie1(kv)  = spoly%cont(ic)%ie1(kv+2)
            spoly%cont(ic)%ie2(kv)  = spoly%cont(ic)%ie2(kv+2)
            spoly%cont(ic)%edge(kv) = spoly%cont(ic)%edge(kv+2)
          end do
          spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 2
        else
          exit SEARCH2
        end if
      end do SEARCH2
      iv = iv + 1
    end do

!   make sure the ABA pattern doesn't occur at the end of a line
    if (spoly%cont(ic)%num_vert >= 3) then
      iv = 1
      jv = spoly%cont(ic)%num_vert - 1
      dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
           + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2 
      if ((dist <= degen_tol2).or.(check_nodes.and.(spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)))) then
        spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 2
      end if
    end if
   
!   this section needed to be modified to properly handle edge data
!   make sure the ABA pattern doesn't occur at the beginning of a line
    if (spoly%cont(ic)%num_vert >= 4) then
      iv = 2
      jv = spoly%cont(ic)%num_vert
      dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
           + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2 
      if ((dist <= degen_tol2).or.(check_nodes.and.(spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)))) then
        do kv = 1,spoly%cont(ic)%num_vert - 2
          spoly%cont(ic)%x(kv)    = spoly%cont(ic)%x(kv+1)
          spoly%cont(ic)%y(kv)    = spoly%cont(ic)%y(kv+1)
          spoly%cont(ic)%node(kv) = spoly%cont(ic)%node(kv+1)
          spoly%cont(ic)%ie1(kv)  = spoly%cont(ic)%ie1(kv+1)
          spoly%cont(ic)%ie2(kv)  = spoly%cont(ic)%ie2(kv+1)
          spoly%cont(ic)%edge(kv) = spoly%cont(ic)%edge(kv+1)
        end do
        spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 2
      end if
    end if


!   see if the entire contour needs to be erased
    erase_contour = .false.

    if (spoly%cont(ic)%num_vert <= 2) then
       erase_contour = .true.

    else &
    if ((maxval(spoly%cont(ic)%x(1:spoly%cont(ic)%num_vert)) &
        -minval(spoly%cont(ic)%x(1:spoly%cont(ic)%num_vert)) &
       <= degen_tol).or.                                             &
        (maxval(spoly%cont(ic)%y(1:spoly%cont(ic)%num_vert)) &
        -minval(spoly%cont(ic)%y(1:spoly%cont(ic)%num_vert)) &
       <= degen_tol)) then
       erase_contour = .true.
    end if

    if (check_nodes) then
      one_edge = .true.
      do iv = 2,spoly%cont(ic)%num_vert
        if (spoly%cont(ic)%edge(iv) /= spoly%cont(ic)%edge(1)) one_edge = .false.
      end do
      if (one_edge) erase_contour = .true.
    end if

    if (erase_contour) then
      do jc = ic,spoly%num_cont-1
        spoly%hole(jc)          = spoly%hole(jc+1)
        if (spoly%cont(jc)%num_vert &
           /= spoly%cont(jc+1)%num_vert) then
          spoly%cont(jc)%num_vert = spoly%cont(jc+1)%num_vert
          deallocate(spoly%cont(jc)%x); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%x(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
          deallocate(spoly%cont(jc)%y); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%y(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
          deallocate(spoly%cont(jc)%node); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%node(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
          deallocate(spoly%cont(jc)%ie1); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%ie1(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
          deallocate(spoly%cont(jc)%ie2); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%ie2(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
          deallocate(spoly%cont(jc)%edge); nalloc = nalloc - 1
          allocate(spoly%cont(jc)%edge(spoly%cont(jc)%num_vert)); nalloc = nalloc + 1
        end if
        do iv = 1,spoly%cont(jc)%num_vert
          spoly%cont(jc)%x(iv)    = spoly%cont(jc+1)%x(iv)
          spoly%cont(jc)%y(iv)    = spoly%cont(jc+1)%y(iv)
          spoly%cont(jc)%node(iv) = spoly%cont(jc+1)%node(iv)
          spoly%cont(jc)%ie1(iv)  = spoly%cont(jc+1)%ie1(iv)
          spoly%cont(jc)%ie2(iv)  = spoly%cont(jc+1)%ie2(iv)
          spoly%cont(jc)%edge(iv) = spoly%cont(jc+1)%edge(iv)
        end do
      end do
      jc = spoly%num_cont
      deallocate(spoly%cont(jc)%x); nalloc = nalloc - 1
      deallocate(spoly%cont(jc)%y); nalloc = nalloc - 1
      deallocate(spoly%cont(jc)%node); nalloc = nalloc - 1
      deallocate(spoly%cont(jc)%ie1); nalloc = nalloc - 1
      deallocate(spoly%cont(jc)%ie2); nalloc = nalloc - 1
      deallocate(spoly%cont(jc)%edge); nalloc = nalloc - 1
      spoly%num_cont = spoly%num_cont - 1
    else
      ic = ic + 1
    end if

  end do

  return
  end subroutine GarbageCollector

  end module GarbageCollectorM
