! =====================================================================
  module PatchInfo
! =====================================================================
! Module to carry arrays and parameters that define the surface data
! associated with the structured solvers.  Most of the data structures
! in USURP are defined on a panel-by-panel basis.  Saving the 
! structured patch data allows us to reconstruct the structured grids
! for visualization output and can also be used to organize the
! integrated force and moment data.
!
! num_scalars    is the number of dependent variables in the dataset
! num_patches    is the number of structured grid patches
! num_components is the number of named geometric components

  implicit none
  save

  integer,parameter                  ,public :: max_components = 100
  integer,allocatable,dimension (:)  ,public :: i1,j1,k1,i2,j2,k2
  integer,allocatable,dimension (:)  ,public :: block,icomp,ibdir
  integer,allocatable,dimension (:)  ,public :: rankp
  integer,allocatable,dimension (:)  ,public :: color
  integer,allocatable,dimension (:)  ,public :: overlaps
  integer,allocatable,dimension (:,:),public :: table
  integer,allocatable,dimension (:)  ,public :: ncomp_surfs,refc
  integer                            ,public :: num_components=0
  integer                            ,public :: num_patches=0
  integer                            ,public :: num_scalars=0
  integer                            ,public :: num_panels=0

  logical                            ,public :: firstpass

  character(len=20),allocatable,dimension(:),public :: tag_list

  end module PatchInfo
