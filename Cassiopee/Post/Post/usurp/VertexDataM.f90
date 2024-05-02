! ====================================================
  module VertexData
! ====================================================
! The array "xv(3,:)" holds all of the grid points
! from the original surface grid. The array "qv(nq,:)" 
! holds all of the flow variables at each original 
! surface grid point. This is intended for use with the
! surface triangulation option in order to produce a
! grid.i.triq file.
!
! variables:
! xv(3,:) stores the grid coordinates at each vertex
! qv(nqv,:) stores the conservative variables at each vertex
! nqv is the number of conservative variables (min(5,nq))
! num_nodes is the total number of vertices (including duplicates)
! num_used is the number of vertices that are used in the triangulation

  use IntrType,only: rd

  implicit none
  save

  real(kind=rd),allocatable,dimension (:,:),public :: xv,xv2,xv2o
  real(kind=rd),allocatable,dimension (:,:),public :: qv,qv2

  integer,allocatable,dimension(:),public :: idup,iassoc
  integer,allocatable,dimension(:),public :: alias

  integer,public :: num_nodes=0
  integer,public :: new_nodes=0
  integer,public :: num_used=0
  integer,public :: nqv
  
  logical,allocatable,dimension(:),public :: used

  end module VertexData
