! =============================================================
  module GroupInfo
! =============================================================
! Module to carry arrays that store the integration group data
! for FOMOCO / OVERFLOW operation modes.

  implicit none
  save

  integer,allocatable,dimension(:,:),public :: group_data
  integer,allocatable,dimension(:),public   :: group_size,group_iref

  integer,public :: num_groups

  character(len=40),allocatable,dimension(:),public :: group_name

  end module GroupInfo
