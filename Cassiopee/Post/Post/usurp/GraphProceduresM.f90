! ========================================================
  module GraphProcedures
! ========================================================
! This module contains procedures for working with a graph


  implicit none

  private

  public  :: LinearExtension,ColorGraph,DepthFirst
  private :: DeleteEdges,DeleteVertex,FindVertex

  contains

! ====================================================================
  subroutine ColorGraph(color,num_vertices,perm_edge, &
                        num_edges,num_colors)
! ====================================================================
! This procedure produces a vertex color scheme for a graph
! using a specified number of colors.

  implicit none

  integer,intent(in) :: num_vertices
  integer,intent(in) :: num_edges
  integer,intent(in) :: num_colors

  integer,intent(out),dimension(num_vertices) :: color
  integer,intent(in),dimension(2,num_edges)   :: perm_edge

  integer,dimension(2,num_edges)  :: edge
  integer,dimension(num_vertices) :: degree,stack
  integer,dimension(0:num_colors) :: used

  integer :: i,j,k,istack,neighbor,nc

  logical,dimension(num_colors) :: available
  logical :: stuck


  continue


! copy the arrays so as not to damage the originals
  edge = perm_edge


! calculate the degree of each vertex
  degree = 0
  do i = 1,num_edges
    do j = 1,2
      k = edge(j,i)
      degree(k) = degree(k) + 1
    end do
  end do

! remove any vertices with degree=num_colors and place it on
! the stack
  i = 1
  istack = num_vertices

  do while (istack > 0)
    stuck = .true.
    do while (i <= num_vertices)
      if ((degree(i) /= -1).and.(degree(i) < num_colors)) then
!       delete this vertex from the graph
        degree(i) = -1
        do j = 1,num_edges
          if ((edge(1,j) == i).or.(edge(2,j) == i)) then
            if (edge(1,j) /= i) then
              neighbor = edge(1,j)
            else
              neighbor = edge(2,j)
            end if
            degree(neighbor) = degree(neighbor) - 1
            edge(:,j) = -1
          end if
        end do
        stuck = .false.
!       add this vertex to the stack
        stack(istack) = i
        istack = istack - 1
      end if
      i = i + 1
    end do
    if (stuck) then
!     Push the vertex with the highest degree to the stack in
!     the hopes that deleting it will free up the rest of the graph.
!     Change the degree to 0 and let the loop take care of the rest.
      degree(maxval(maxloc(degree))) = 0
    end if
    i = 1
  end do

! initialize all nodes as uncolored
  color = 0
  used = 0

! pop vertices from the stack and color them
  color(stack(1)) = 1
  used(1) = 1

  do i = 2,num_vertices

!   pop the next node from the stack
    j = stack(i)

!   assume all colors are available
    available = .true.

!   rule out colors that are not available
    do k = 1,num_edges
      if ((perm_edge(1,k) == j).or.(perm_edge(2,k) == j)) then
        if (perm_edge(1,k) /= j) then
          neighbor = perm_edge(1,k)
        else
          neighbor = perm_edge(2,k)
        end if
        nc = color(neighbor)
        if (nc /= 0) available(nc) = .false.
      end if
    end do

!   assign the least-used available color
    nc       = maxval(minloc(used,mask=available))
    color(j) = nc
    used(nc) = used(nc) + 1

  end do
              
  return
  end subroutine ColorGraph

!   ====================================================================
    subroutine LinearExtension(rankp,max_vertices,perm_edge, &
                               max_edges,perm_weight,dropped)
!   ====================================================================
!   This procedure builds a linear extension of the given weighted, 
!   directed graph and eliminates feedback arc sets if the 
!   graph is not acyclic.  The procedure is a based on
!   a heuristic described by Steven S. Skiena in the CD-ROM
!   version of the book "The Algorithm Design Manual" published
!   by Springer-Verlag, New York, 1997.  The variation here breaks
!   cycles by dropping the arc with the lowest weight. 

    use IntrType,only: rd

    implicit none

    integer,intent(in)                            :: max_vertices
    integer,intent(in)                            :: max_edges

    real(kind=rd),intent(in),dimension(max_edges) :: perm_weight
    real(kind=rd),           dimension(max_edges) :: weight
    real(kind=rd)                                 :: wmin

    integer,intent(in), dimension(2,max_edges)    :: perm_edge
    integer,            dimension(2,max_edges)    :: edge
    integer,intent(out),dimension(max_vertices)   :: rankp
    integer,            dimension(max_vertices)   :: ein,eout,vertex

    integer,intent(out)   :: dropped        ! number of dropped constraints
    integer               :: num_vertices   ! current number of vertices that need ranks
    integer               :: num_edges      ! remaining number of constraints

    integer               :: i,iedge,i1,i2,j
    integer               :: min_deg
    integer               :: low_rank,my_rank,top_rank


    continue


!   Initialization
    dropped      = 0
    ein          = 0
    eout         = 0
    rankp        = 0
    edge         = perm_edge
    weight       = perm_weight
    num_vertices = max_vertices
    num_edges    = max_edges
    do i = 1,num_vertices
       vertex(i) = i
    end do


!   Delete any edges with weight equal to one
    i = 1
    do while (i <= num_edges)
      if (weight(i) == 1.0) then
        do j = i,num_edges-1
          weight(j) = weight(j+1)
          edge(:,j) = edge(:,j+1)
        end do
        num_edges = num_edges - 1
      else
        i = i + 1
      end if
    end do


!   Build the initial in/out statistics
    do i = 1,num_edges
      eout(edge(1,i)) = eout(edge(1,i)) + 1
      ein (edge(2,i)) = ein (edge(2,i)) + 1
    end do


    top_rank = num_vertices
    low_rank = 1


    MAIN_LOOP: do while (num_vertices > 0)

      i = num_vertices
      do while (i >= 1)

        if ((ein(i) == 0).or.(eout(i) == 0)) then
          if (ein(i) == 0) then

!           If there are no edges pointing into this vertex, then this 
!           vertex is a "source" in the remaining arc set and can be assigned 
!           the highest remaining rank

            my_rank  = top_rank
            top_rank = top_rank - 1

          else ! (eout(i) == 0)

!           If there are no edges pointing out of this vertex, then this 
!           vertex is a "sink" in the remaining arc set and can be assigned
!           the lowest remaining rank

            my_rank  = low_rank
            low_rank = low_rank + 1

          end if

          rankp(vertex(i)) = my_rank
          call DeleteEdges(vertex(i),vertex,edge,ein,eout,num_edges, &
                           max_edges,max_vertices,num_vertices,weight)
          call DeleteVertex(i,ein,eout,vertex,max_vertices,num_vertices)
          cycle MAIN_LOOP

        end if

        i = i - 1
      end do


!     If we made it this far, it is because the reminaing arc
!     set is cyclic.  Discard the edge with the lowest weight.

      iedge = 1
      i1 = FindVertex(edge(1,iedge),vertex,max_vertices,num_vertices)
      i2 = FindVertex(edge(2,iedge),vertex,max_vertices,num_vertices)
      wmin    = weight(iedge)
      min_deg = min(eout(i1),ein(i2))

      do i = 2,num_edges
        if (weight(i) < wmin) then
          iedge = i
          i1 = FindVertex(edge(1,i),vertex,max_vertices,num_vertices)
          i2 = FindVertex(edge(2,i),vertex,max_vertices,num_vertices)
          wmin    = weight(i)
          min_deg = min(eout(i1),ein(i2))
        else if (weight(i) == wmin) then
          i1 = FindVertex(edge(1,i),vertex,max_vertices,num_vertices)
          i2 = FindVertex(edge(2,i),vertex,max_vertices,num_vertices)
          if (min(eout(i1),ein(i2)) < min_deg) then
            iedge = i
            min_deg = min(eout(i1),ein(i2))
          end if
        end if
      end do

!     adjust the in/out-degrees of the vertices in the deleted edge
      i1 = FindVertex(edge(1,iedge),vertex,max_vertices,num_vertices)
      i2 = FindVertex(edge(2,iedge),vertex,max_vertices,num_vertices)
      eout(i1) = eout(i1) - 1
      ein(i2)  = ein(i2) - 1

!     delete the edge
      do j = iedge,num_edges-1
         weight(j) = weight(j+1)
         edge(:,j) = edge(:,j+1)
      end do
      num_edges = num_edges - 1

      dropped = dropped + 1

    end do MAIN_LOOP

    if ((num_edges /= 0).or.(num_vertices /= 0)) then
      write(0,*)
      write(0,*)"ERROR! failure in LinearExtension algorithm."
      if (num_edges /= 0) then
        write(0,*)"final value of num_edges = ",num_edges
      end if
      if (num_vertices /= 0) then
        write(0,*)"final value of num_vertices = ",num_vertices
      end if
      stop
    end if

!   if (dropped > 0) then
!     write(0,"(1x,a,i4,a,i4,a)") "LinearExtension dropped ", &
!       dropped," of the ",max_edges," original constraints"
!   else
!     write(0,"(1x,a,i4,a)") "LinearExtension met all ", &
!       max_edges," of the original constraints"
!   end if

    return
    end subroutine LinearExtension

!   ===============================================================
    subroutine DeleteEdges(i,vertex,edge,ein,eout,num_edges, &
                           max_edges,max_vertices,num_vertices, &
                           weight)
!   ===============================================================
!   delete weighted edges containing the vertex "i"

    use IntrType,only: rd

    implicit none

    integer,intent(in) :: i,max_edges,max_vertices,num_vertices

    real(kind=rd),dimension (max_edges),intent(inout) :: weight

    integer,dimension (2,max_edges) ,intent(inout) :: edge
    integer,dimension (max_vertices),intent(inout) :: ein
    integer,dimension (max_vertices),intent(inout) :: eout
    integer,dimension (max_vertices),intent(in   ) :: vertex
    integer                         ,intent(inout) :: num_edges

    integer :: j,k


    continue


    j = 1
    do while (j <= num_edges)
      if ((edge(1,j) == i).or.(edge(2,j) == i)) then
        do k = 1,num_vertices
          if (edge(1,j) == vertex(k)) eout(k) = eout(k) - 1
          if (edge(2,j) == vertex(k)) ein(k) = ein(k) - 1
        end do
        do k = j,num_edges - 1
          edge(1:2,k) = edge(1:2,k+1)
          weight(k)   = weight(k+1)
        end do
        num_edges = num_edges - 1
      else
        j = j + 1
      end if
    end do

    return
    end subroutine DeleteEdges

!   =========================================================================
    subroutine DeleteVertex(i,ein,eout,vertex,max_vertices,num_vertices)
!   =========================================================================
!   delete a vertex from the graph

    implicit none

    integer,intent(in)                            :: max_vertices
    integer,intent(inout),dimension(max_vertices) :: ein,eout,vertex
    integer,intent(inout)                         :: num_vertices
    integer,intent(in)                            :: i
    integer                                       :: j


    continue


    do j = i,num_vertices - 1
      ein(j) = ein(j+1)
      eout(j) = eout(j+1)
      vertex(j) = vertex(j+1)
    end do

    num_vertices = num_vertices - 1

    return
    end subroutine DeleteVertex

!   ========================================================================
    pure integer function FindVertex(ivert,vertex,max_vertices,num_vertices)
!   ========================================================================
!   find the position of vertex "ivert" in the list

    implicit none

    integer,intent(in)                         :: ivert
    integer,intent(in)                         :: num_vertices
    integer,intent(in)                         :: max_vertices
    integer,intent(in),dimension(max_vertices) :: vertex
    integer                                    :: i


    continue


    do i = 1,num_vertices
       if (vertex(i) == ivert) then
         FindVertex = i
         return
       end if
    end do

!   return error flag
    FindVertex = -1

    return
    end function FindVertex

! =====================================================================
  subroutine DepthFirst(my_rank,num_nodes,edge,num_edges,weight, &
                         total_dropped)
! =====================================================================
! create a linear extension of a graph by first grouping the nodes
! into strongly connected components (Tarjan 1972 as found in
! Sedgewick 1983)

  use IntrType,only: rd

  implicit none

  integer,intent(in) :: num_edges,num_nodes

  real(kind=rd),dimension (num_edges)  ,intent(in)  :: weight
  real(kind=rd),dimension (num_edges)               :: sg_weight

  integer,dimension (2,num_edges),intent(in)  :: edge
  integer,dimension (2,num_edges)             :: sg_edge

  integer,dimension (num_nodes)  ,intent(out) :: my_rank
  integer,dimension (num_nodes)               :: list,spot
  integer,dimension (num_nodes)               :: stack,val
  integer,dimension (num_nodes)               :: sg_rank

  integer,intent(out) :: total_dropped
  integer             :: i,m,p
  integer             :: visited,dropped
  integer             :: first,last,members
  integer             :: current_rank,low_rank,top_rank
  integer             :: num_sg_edges


  continue


  current_rank  = 1
  members       = 0
  p             = 1
  val           = 0

  total_dropped = 0
  my_rank       = 0
  visited       = 0

  do i = 1,num_nodes
    if (val(i) == 0) call visit(i,m)
  end do


! process the subgraphs
  first = 1
  do while (first <= num_nodes)

!   determine the subgraph
    current_rank = my_rank(list(first))
    last = first
    SUBGRAPH: do while (last < num_nodes)
       if (my_rank(list(last+1)) == current_rank) then
         last = last + 1
         if (last == num_nodes) exit SUBGRAPH
       else
         exit SUBGRAPH
       end if
    end do SUBGRAPH

!   process the subgraph
    if (last /= first) then
      
      members  = last - first + 1
      low_rank = current_rank
      top_rank = current_rank + members - 1

      num_sg_edges = 0
!     pick out the edges of this subgraph
      do i = 1,num_edges
        if ((my_rank(edge(1,i)) == current_rank).and. &
            (my_rank(edge(2,i)) == current_rank)) then
          num_sg_edges = num_sg_edges + 1
          sg_edge(1,num_sg_edges) = spot(edge(1,i))-first+1
          sg_edge(2,num_sg_edges) = spot(edge(2,i))-first+1
          sg_weight(num_sg_edges) = weight(i)
        end if
      end do

      call LinearExtension(sg_rank,members,sg_edge, &
                                 num_sg_edges,sg_weight,dropped)
      total_dropped = total_dropped + dropped

      do i = 1,members
        my_rank(list(i+first-1)) = sg_rank(i)+current_rank-1
      end do
    end if


!   check the next subgraph
    first = last + 1

  end do

  contains

! ===============================================================
  recursive subroutine visit(k,ret_val)
! ===============================================================
! recursive function to group strongly connected nodes of a graph
! Tarjan 1972 as described by Sedgewick 1983

  implicit none

  integer,intent(in)     :: k
  integer,intent(out)    :: ret_val
  integer                :: j,min_val,m


  continue


  visited  = visited + 1
  val(k)   = visited
  min_val  = visited
  stack(p) = k

  p = p + 1

  do j = 1,num_edges
    if (edge(1,j) == k) then

      if (val(edge(2,j)) == 0) then
        call visit(edge(2,j),m)
        if (m < min_val) min_val = m
      else if (val(edge(2,j)) < min_val) then
        min_val = val(edge(2,j))
      end if

    end if
  end do


  if (min_val == val(k)) then
    do
      p = p - 1
      val(stack(p)) = num_nodes + 1
      my_rank(stack(p)) = current_rank
      list(current_rank + members) = stack(p)
      spot(stack(p)) = current_rank + members
      members = members + 1
      if (stack(p) == k) then
        current_rank = current_rank + members
        members = 0
        exit
      end if
    end do
  end if

  ret_val = min_val

  return
  end subroutine visit

  end subroutine DepthFirst

  end module GraphProcedures
