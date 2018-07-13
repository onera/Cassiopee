! ====================================================
  module EdgeData
! ====================================================
!
! variables:
! num_edges is the total number of edges (including duplicates)

  implicit none
  save

  type,public :: edge_type
    integer :: left_face
    integer :: right_face
    integer,dimension(2) :: vertex
  end type edge_type

  type(edge_type),allocatable,dimension(:),public :: edge

  integer,public :: num_edges=0

  end module EdgeData
