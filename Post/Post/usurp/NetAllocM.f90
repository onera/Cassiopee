! ===================================================
  module NetAlloc
! ===================================================
! variables to keep track of the number of allocated
! arrays, entries, nodes, and entry pointers in order
! to help avoid memory leaks

  implicit none

  integer,save,public :: nalloc = 0
  integer,save,public :: nentry = 0
  integer,save,public :: nnodes = 0
  integer,save,public :: npnten = 0

  end module NetAlloc
