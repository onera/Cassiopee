! =====================================================================
  module IntrType
! =====================================================================
! module to define single and double precision real types

  implicit none
  save

  integer,parameter,public :: rd = selected_real_kind(p=12,r=37)
  integer,parameter,public :: rs = selected_real_kind(p=6,r=37)

  end module IntrType
