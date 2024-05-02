! =====================================================================
  module TimeKeeper
! =====================================================================
! A module used only to carry the time spent in the gpc library.
! the module is used by the main code and by the ProcessPair
! module.

  use IntrType,only: rs

  implicit none
  save

  real(kind=rs),public :: gpc_time

  end module TimeKeeper
