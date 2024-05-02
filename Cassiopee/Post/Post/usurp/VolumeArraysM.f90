! =============================================================
  module VolumeArrays
! =============================================================
! Module to carry arrays that store the grid and volumes for
! full blocks, previously passed as arguments to the 
! CalculateVolumes, ReadGrid, and ReadFlow subroutines or else
! created as temporary arrays in those routines.  They have been
! moved here in an attempt to alleviate memory demands on the stack.

  use IntrType,only: rd,rs

  implicit none
  save

  real(kind=rd),allocatable,dimension (:,:,:,:),public :: qd
  real(kind=rd),allocatable,dimension (:,:,:)  ,public :: xb,yb,zb,volv

  real(kind=rs),allocatable,dimension (:,:,:,:),public :: qs
  real(kind=rs),allocatable,dimension (:,:,:)  ,public :: xs,ys,zs

  integer      ,allocatable,dimension (:,:,:)  ,public :: ibv

  end module VolumeArrays
