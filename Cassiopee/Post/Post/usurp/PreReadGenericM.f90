! ================================================================
  module PreReadGenericM
! ================================================================

  implicit none

  private
  public :: PreReadGeneric

  contains

! ======================================================================
  subroutine PreReadGeneric()
! ======================================================================
! Read the boundary condition file (generic.bc) in order
! to acquire the number of panels, patches, and vertices, and to store
! the patch data.
  use PatchInfo                ,only: tag_list
  use PatchInfo                ,only: num_patches,icomp,num_components
  use PatchInfo                ,only: i1,i2,j1,j2,k1,k2,ncomp_surfs,block
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch,EchoPatchData
  use Types                    ,only: surf
  use UserInput                ,only: workdir

  implicit none

  integer             :: ip,n,ierr


  continue

! first read to ascertain number of patches and components
  write(*,*)
  write(unit=*,fmt=*)"Reading generic.bc ..."
  open(unit=3,file=trim(workdir)//"generic.bc",status="old", &
              action="read",position="rewind")
  do
    read(unit=3,fmt=*,iostat=ierr)n
    if (ierr /= 0) exit
    num_patches = num_patches + 1
  end do
  
  num_patches = num_patches/2
  
  call AllocateStructuredPatch()  
  allocate(ncomp_surfs(0:0),stat=ierr)
  allocate(tag_list(0:0),stat=ierr)
! second read to fill in data
  ip = 0
  rewind(unit=3)
  do while (ip < num_patches)
    ip = ip + 1
    read(unit=3,fmt=*,iostat=ierr)block(ip),i1(ip),i2(ip),j1(ip), &
                                            j2(ip),k1(ip),k2(ip)
    if (ierr /= 0) exit

    icomp(ip) = 0
    ncomp_surfs(0) = ncomp_surfs(0) + 1
    call CountStructuredPatch(ip)
  end do
  close(unit=3)

! screen output to confirm results
  call EchoPatchData()

  return
  end subroutine PreReadGeneric

  end module PreReadGenericM
