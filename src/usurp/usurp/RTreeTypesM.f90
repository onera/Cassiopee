! =========================================================
  module RTreeTypes
! =========================================================
! This module contains Fortran derived type definitions
! and parameters associated with the R-tree data structure.

  use IntrType,only: rd
  use Types   ,only: panel

  implicit none

  integer,parameter,public :: ndir = 3
  integer,parameter,public :: MaxEntries =  5
  integer,parameter,public :: MinEntries =  2
  integer,save     ,public :: total_nodes = 0
  logical,save     ,public :: debug = .false.

  type,public :: entry
    type(node) ,pointer :: pn
    type(panel),pointer :: pp
    real(kind=rd),dimension(ndir,2) :: I
  end type entry

  type,public :: entry_ptr
    type(entry),pointer :: pe
  end type entry_ptr

  type,public :: node
    type (node) , pointer                  :: parent
    type (entry_ptr), dimension(:),pointer :: E
    integer                                :: entries
    integer                                :: node_number
    logical                                :: leaf
  end type node

  end module RTreeTypes
