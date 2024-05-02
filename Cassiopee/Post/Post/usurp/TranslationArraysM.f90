! =====================================================================
  module TranslationArrays
! =====================================================================
! Module to carry arrays used to pass real and integer data from
! Fortran to the gpc library, which is in written in c.

  use IntrType,only: rd

  implicit none
  save

  integer,parameter             ,public :: asize = 1000
  real(kind=rd),dimension(asize),public :: fltc_array
  integer      ,dimension(asize),public :: intc_array

  end module TranslationArrays
