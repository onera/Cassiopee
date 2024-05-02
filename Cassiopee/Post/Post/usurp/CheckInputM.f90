! ================================================================
  module CheckInputM
! ================================================================

  implicit none

  private
  public :: CheckInput

  contains

! ================================================
  subroutine CheckInput()
! ================================================
! This subroutine checks conflicts in the command
! line input

  use UserInput,only: basis,checkpatches,colormap,plotvel
  use UserInput,only: solution_exists

  implicit none


  continue


  if (basis == "error") then 
    write(0,*)"ERROR! The basis variable was never established."
    write(0,*)"this is a bug in the code but can be avoided by"
    write(0,*)"specifying the basis on the command line."
    stop
  end if

  if (colormap.and.(basis /= "patch")) then
    write(0,*)"WARNING! The --color option is only available"
    write(0,*)"when --basis=patch, so it has been reset to .false."
    colormap = .false.
  end if

  if (checkpatches.and.(basis /= "patch")) then
    write(0,*)"WARNING! The --check-patches option is only available"
    write(0,*)"when --basis=patch, so it will be ignored in this case."
    checkpatches = .false.
  end if

  if (plotvel.and.(.not.solution_exists)) then
    write(0,*)"WARNING! The velocity data is not available, so"
    write(0,*)"the --plot-velocity option will be ignored."
    plotvel = .false.
  end if

  return
  end subroutine CheckInput

  end module CheckInputM
