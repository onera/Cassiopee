! ================================================================
  module PreReadBCinM
! ================================================================

  implicit none

  private
  public :: PreReadBCin

  contains

! ======================================================================
  subroutine PreReadBCin()
! ======================================================================
! Read the boundary condition file (bc.in) for UNCLE / OVER-REL in order
! to acquire the number of panels, patches, and vertices, and to store
! the patch data.

  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: ncomp_surfs,refc,tag_list
  use PatchInfo                ,only: i1,j1,k1,i2,j2,k2,icomp,block
  use PatchInfo                ,only: num_components,max_components
  use PatchInfo                ,only: num_patches
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch,EchoPatchData
  use Types                    ,only: surf
  use UserInput                ,only: workdir

  implicit none

  integer             :: ibc,ic,icindex,ierr,ip,len_str,n
  character(len=132)  :: whole_line
  character(len=20)   :: tag
  logical             :: new_component


  continue


! initialization
  allocate(tag_list(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate tag_list in PreReadBCin"
    write(0,*)"max_components = ",max_components
    stop
  end if
  tag_list(0)    = "untagged"

  allocate(ncomp_surfs(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate ncomp_surfs in PreReadBCin"
    write(0,*)"max_components = ",max_components
    stop
  end if
  ncomp_surfs = 0

  allocate(refc(1:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate refc in PreReadBCin"
    write(0,*)"num_components = ",num_components
    stop
  end if



! first read to ascertain number of patches and components
  write(*,*)
  write(unit=*,fmt=*)"Reading bc.in ..."
  open(unit=3,file=trim(workdir)//"bc.in",status="old", &
              action="read",position="rewind")
  do
    read(unit=3,fmt=*)n,ibc
    if (n == 999) exit
    if ((ibc /= 1).and.(ibc /= 7)) cycle

    num_patches = num_patches + 1

!   check for components
    backspace(unit=3)
    read(unit=3,fmt="(a132)")whole_line
    icindex = index(whole_line,"component=")

    if (icindex /= 0) then 
!     read(unit=whole_line(icindex+10:132),fmt=*)tag
      icindex = icindex + 10
      len_str = len(trim(whole_line))
      tag = adjustl(whole_line(icindex:len_str-1))

!     see if we should add this component to the list
      new_component = .true.
      do ic = 1,num_components
        if (tag == tag_list(ic)) new_component = .false.
      end do
      if (new_component) then
        num_components = num_components + 1
        tag_list(num_components) = tag
        if (num_components > max_components) then
          write(0,*)"ERROR! num_components exceeds max_components."
          write(0,*)"Increase max_components in PatchInfoM.f90"
          stop
        end if
      end if
    end if
  end do


  call AllocateStructuredPatch()


! second read to fill in data
  ip = 0
  rewind(unit=3)
  do
    read(unit=3,fmt=*)n,ibc
    if (n == 999) exit
    if ((ibc /= 1).and.(ibc /= 7)) cycle

    ip = ip + 1
    backspace(unit=3)
    read(unit=3,fmt=*)block(ip),ibc,i1(ip),j1(ip),k1(ip), &
                                    i2(ip),j2(ip),k2(ip)
     

!   check for components
    backspace(unit=3)
    read(unit=3,fmt="(a132)")whole_line
    icindex = index(whole_line,"component=")

    if (icindex /= 0) then 
!     read(unit=whole_line(icindex+10:132),fmt=*)tag
      icindex = icindex + 10
      len_str = len(trim(whole_line))
      tag = adjustl(whole_line(icindex:len_str-1))

      do ic = 1,num_components
        if (tag == tag_list(ic)) then
          icomp(ip) = ic
          ncomp_surfs(ic) = ncomp_surfs(ic) + 1
        end if
      end do
    else
!     untagged
      icomp(ip) = 0
      ncomp_surfs(0) = ncomp_surfs(0) + 1
    end if

    call CountStructuredPatch(ip)

  end do
  close(unit=3)


! screen output to confirm results
  call EchoPatchData()


  return
  end subroutine PreReadBCin

  end module PreReadBCinM
