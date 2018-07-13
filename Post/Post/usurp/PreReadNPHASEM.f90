! ================================================================
  module PreReadNPHASEM
! ================================================================

  implicit none

  private
  public :: PreReadNPHASE

  contains

! ========================================================
  subroutine PreReadNPHASE()
! ========================================================
! Read the NPHASE data file (nphase_elements.dat) in order
! to acquire the number of panels and vertices.

  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: tag_list,max_components,icomp
  use PatchInfo                ,only: num_scalars,num_panels,num_patches
  use StructuredPatchProcedures,only: EchoPatchData
  use UserInput                ,only: workdir
  use VertexData               ,only: num_nodes

  implicit none

  integer             :: i,pid,ierr
  integer             :: j,node(4),nvert


  continue


! initialization
  allocate(tag_list(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate tag_list in PreReadNPHASE"
    write(0,*)"max_components = ",max_components
    stop
  end if
  tag_list(0)    = "untagged"


! read the original vertex data

  open(unit     = 3,                                    &
       file     = trim(workdir)//"nphase_elements.dat", &
       form     = "formatted",                          &
       status   = "old",                                &
       position = "rewind",                             &
       action   = "read")


! read node data
  read(unit=3,fmt=*)num_nodes
  do i = 1,num_nodes
    read(unit=3,fmt=*)
  end do

! read panel data
  read(unit=3,fmt=*)num_panels
  do i = 1,num_panels
    read(unit=3,fmt=*)nvert,pid,(node(j),j=1,nvert)
    num_patches = max(num_patches,pid+1)
  end do

! attempt to read scalar data
  num_scalars = 0
  do
    do i = 1,num_panels
      read(unit=3,fmt=*,iostat=ierr)
      if (ierr /= 0) exit
    end do
    if (ierr == 0) then
      num_scalars = num_scalars + 1
    else
      exit
    end if
  end do
  close(unit=3)

  allocate(icomp(num_patches)); nalloc = nalloc + 1
  icomp = 0


! screen output to confirm results
  call EchoPatchData()


  return
  end subroutine PreReadNPHASE

  end module PreReadNPHASEM
