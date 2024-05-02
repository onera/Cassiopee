! =============================================================
  module FomoArrays
! =============================================================
! Module to carry arrays that store the integrated forces and
! moments separated into pressure and viscous components as
! well as geometric components.

  use IntrType,only: rd

  implicit none
  save

  real(kind=rd),allocatable,dimension (:),public :: cfpx,cfpy,cfpz
  real(kind=rd),allocatable,dimension (:),public :: cfvx,cfvy,cfvz
  real(kind=rd),allocatable,dimension (:),public :: cfmx,cfmy,cfmz
  real(kind=rd),allocatable,dimension (:),public :: cmpx,cmpy,cmpz
  real(kind=rd),allocatable,dimension (:),public :: cmvx,cmvy,cmvz
  real(kind=rd),allocatable,dimension (:),public :: cmmx,cmmy,cmmz
  real(kind=rd),allocatable,dimension (:),public :: ctax,ctay,ctaz
  real(kind=rd),allocatable,dimension (:),public :: weta,cfr

  end module FomoArrays
