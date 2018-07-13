! ================================================================
  module BuildGraphM
! ================================================================

  implicit none

  private
  public :: BuildGraph

  contains

! ================================================================
  subroutine BuildGraph(p)
! ================================================================
! This procedure establishes a graph from the patch overlap table
! and then creates a linear extension of the graph (which forms
! the patch ranking system) and a vertex coloring scheme (for
! visualization in Tecplot)

  use GraphProcedures,only: LinearExtension,ColorGraph,DepthFirst
  use IntrType       ,only: rd,rs
  use my_cpu_timeM   ,only: my_cpu_time
  use NetAlloc       ,only: nalloc
  use PatchInfo      ,only: num_patches,table,color,rankp
  use PatchInfo      ,only: num_panels,overlaps
  use Types          ,only: panel,surf
  use UserInput      ,only: colormap,checkpatches,icpatch1
  use UserInput      ,only: icpatch2,my_bit_size,verbosity
  use UserInput      ,only: cpufmt,showcpu

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rd),allocatable,dimension (:) :: weight
  real(kind=rs)                           :: time1,time2

  integer,allocatable,dimension (:,:) :: edge

  integer :: i,j,k,ip,ipan,ip1,ip2,ierr,ksav
  integer :: kint,kpos,kmax,maxp,misranked,kworst
  integer :: num_edges,dropped


  continue


  call my_cpu_time(time1)


! build overlaps table
  allocate(table(num_patches,num_patches)); nalloc = nalloc + 1
  table = 0

  ip1 = lbound(p,1)
  ip2 = ubound(p,1)

  maxp = ip2-ip1+1
  kmax = (maxp*num_patches)/my_bit_size

! plow straight through the last word, relying on the fact
! that we initialized the extraneous parts to zero earlier
  if (mod(maxp*num_patches,my_bit_size) /= 0) kmax = kmax + 1

  k = 0
  do kint = 1,kmax
    if (overlaps(kint) /= 0) then
      do kpos = 0,my_bit_size-1

        if (btest(overlaps(kint),kpos)) then
          ipan = mod(k,num_panels) + 1
          ip   = k/num_panels + 1
          i    = p(ipan)%isurf
          table(i,ip) = table(i,ip) + 1
        end if
        k = k + 1
      end do
    else
      k = k + my_bit_size
    end if
  end do


  deallocate(overlaps); nalloc = nalloc - 1


! count the number of edges in the digraph
  num_edges = 0
  do i = 1,num_patches
  do j = 1,i-1
    if (table(i,j) + table(j,i) /= 0)  &
       num_edges = num_edges + 1
  end do
  end do


! create array to store the edges of the digraph 
  allocate(edge(2,num_edges)); nalloc = nalloc + 1
  allocate(weight(num_edges)); nalloc = nalloc + 1


! store the edges of the digraph
  k = 0
  do i = 1,num_patches
  do j = 1,i-1
    if (table(i,j) + table(j,i) /= 0) then
      k = k + 1
      if (table(i,j) > table(j,i)) then
        edge(1,k) = i
        edge(2,k) = j
      else
        edge(1,k) = j
        edge(2,k) = i
      end if
      weight(k) = float(max(table(i,j),table(j,i))) &
                / float(max(min(table(i,j),table(j,i)),1))
    end if
  end do
  end do

  if (k /= num_edges) then
    write(0,*)"ERROR! k /= num_edges in BuildGraph"
    stop
  end if


  if (checkpatches) then
    if ((min(icpatch1,icpatch2) < 1).or. &
        (max(icpatch1,icpatch2) > num_patches)) then
      write(0,*)"ERROR!  --check-patches specification out of range."
      write(0,*)"patches requested: ",icpatch1,icpatch2
      write(0,*)"valid range: ",1,num_patches
      stop
    end if
    if (table(icpatch1,icpatch2)+table(icpatch2,icpatch1) /= 0) then
      write(0,"(1x,a,i4,a,i4,a,i7)")"CheckPatches: T(", &
            icpatch1,",",icpatch2,") = ",table(icpatch1,icpatch2)
      write(0,"(1x,a,i4,a,i4,a,i7)")"CheckPatches: T(", &
            icpatch2,",",icpatch1,") = ",table(icpatch2,icpatch1)
      do k = 1,num_edges
        if ((min(edge(1,k),edge(2,k)) == min(icpatch1,icpatch2)).and. &
            (max(edge(1,k),edge(2,k)) == max(icpatch1,icpatch2))) then
          write(0,"(1x,a,i4,a,i4,a,i4)")"CheckPatches: Edge ",k, &
               " relates the patches as ",edge(1,k)," => ",edge(2,k)
          write(0,*)"              The edge weight is ",weight(k)
          ksav = k
        end if
      end do
    else
      write(0,*)"CheckPatches: There is no overlap between the patches."
    end if
  end if


! create a 6-color map of the graph
  if (colormap) then
    allocate(color(num_patches),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! allocate color failed in BuildGraph"
      stop
    end if
    call ColorGraph(color,num_patches,edge,num_edges,6)
  end if


! create a linear extension of the directed graph
  allocate(rankp(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate rankp in BuildGraph failed"
    stop
  end if


  call DepthFirst(rankp,num_patches,edge,num_edges,weight,dropped)
! call LinearExtension(rankp,num_patches,edge,num_edges,weight,dropped)


  if (verbosity > 0) then
    if (dropped > 0) then
      write(0,"(1x,a,i4,a,i4,a)") "LinearExtension dropped ", &
        dropped," of the ",num_edges," original constraints"
    else
      write(0,"(1x,a,i4,a)") "LinearExtension met all ", &
        num_edges," of the original constraints"
    end if
  end if


! print warnings for mis-ranked nodes
  misranked = 0
  do k = 1,num_edges
    if ((rankp(edge(1,k)) < rankp(edge(2,k))).and. &
        (weight(k) /= 1.0)) then
      misranked = misranked + 1
      if (misranked == 1) then
        kworst = k
      else
        if (weight(k) > weight(kworst)) kworst = k
      end if
    end if
  end do
  if ((verbosity > 0).and.(misranked > 0)) then
    write(0,*)"number of misranked patch pairs = ",misranked
    write(0,"(1x,a,f6.2,a,i4,a,i4)")                            &
            "highest weight between misranked patches is ",     &
            real(weight(kworst)),                               &
            " between ",edge(1,kworst)," and ",edge(2,kworst)
  end if


! debugging
  if (checkpatches) then
    if (rankp(edge(1,ksav)) < rankp(edge(2,ksav))) then
      write(0,*)"WARNING! CheckPatches: The patch ranking was reversed!"
    else
      write(0,*)"CheckPatches: The patch ranking was preserved."
    end if
  end if


! free edge memory
  deallocate(weight); nalloc = nalloc - 1
  deallocate(edge)  ; nalloc = nalloc - 1


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in BuildGraph:',time2-time1
  end if


  return
  end subroutine BuildGraph

  end module BuildGraphM
