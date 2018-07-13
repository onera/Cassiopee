! ================================================================
  module PreReadUNCLEMM
! ================================================================

  implicit none

  private
  public :: PreReadUNCLEM

  contains

! ======================================================================
  subroutine PreReadUNCLEM()
! ======================================================================
! Extract the boundary condition info for UNCLE-M from uncle.inp in 
! order to acquire the number of panels, patches, and vertices, and to 
! store the patch data.

  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: i1,i2,j1,j2,k1,k2,block
  use PatchInfo                ,only: ncomp_surfs,icomp,num_patches
  use PatchInfo                ,only: max_components,tag_list
  use PatchInfo                ,only: num_components
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch
  use StructuredPatchProcedures,only: EchoPatchData
  use Types                    ,only: surf
  use UserInput                ,only: workdir

  implicit none

  integer,allocatable,dimension(:) :: nig,njg,nkg

  integer             :: ip,n,ierr,nb,ii,ic
  character(len=132)  :: whole_line
  character(len=132)  :: gridname
  character(len=020)  :: tag
  character(len=010)  :: fileid
  logical             :: ex,new_component


  continue


! initialization
  allocate(tag_list(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate tag_list in PreReadUNCLEM"
    write(0,*)"max_components = ",max_components
    stop
  end if
  tag_list = "untagged"

  allocate(ncomp_surfs(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate ncomp_surfs in PreReadBCin"
    write(0,*)"max_components = ",max_components
    stop
  end if
  ncomp_surfs = 0



! count the number of grid blocks
  nb = 0
  do
    gridname = "GRIDS/uncle.grd"
    write(fileid,"(i10)")nb
    if (nb < 10) gridname = trim(gridname)//"0"
    gridname = trim(gridname)//adjustl(trim(fileid))
    inquire(file=trim(workdir)//trim(gridname),exist=ex)
    if (ex) then
      nb = nb + 1
    else
      exit
    end if
  end do


! store the grid block sizes
  allocate(nig(nb),njg(nb),nkg(nb)); nalloc = nalloc + 3

  do n = 1,nb
    gridname = "GRIDS/uncle.grd"
    if (n - 1 < 10) gridname = trim(gridname)//"0"
    write(fileid,"(i10)")n-1
    gridname = trim(gridname)//adjustl(trim(fileid))
    open(unit=3,file=trim(workdir)//trim(gridname),form="unformatted", &
         status="old",action="read")
    read(3)
    read(3)nig(n),njg(n),nkg(n)
    close(3)
  end do


! first read to ascertain number of patches and components
  write(*,*)
  write(unit=*,fmt=*)"Reading uncle.inp ..."
  open(unit=3,file=trim(workdir)//"uncle.inp",status="old", &
              action="read",position="rewind")
  do
    read(unit=3,fmt="(a132)")whole_line
    if (index(whole_line,"end") == 1) exit
    if (index(whole_line,"wall_patch") == 1) num_patches = num_patches + 1
  end do


  call AllocateStructuredPatch()


! second read to fill in data
  ip = 0
  rewind(unit=3)
  do
    read(unit=3,fmt="(a132)")whole_line

    if (index(whole_line,"end") == 1) exit
    if (index(whole_line,"wall_patch") /= 1) cycle

    ip = ip + 1
    read(unit=3,fmt=*)block(ip),i1(ip),i2(ip),j1(ip), &
                                j2(ip),k1(ip),k2(ip)

    if (i1(ip) == 0) i1(ip) = nig(block(ip))
    if (i2(ip) == 0) i2(ip) = nig(block(ip))
    if (j1(ip) == 0) j1(ip) = njg(block(ip))
    if (j2(ip) == 0) j2(ip) = njg(block(ip))
    if (k1(ip) == 0) k1(ip) = nkg(block(ip))
    if (k2(ip) == 0) k2(ip) = nkg(block(ip))

    ii = index(whole_line,":") + 1
    if (len(trim(whole_line(ii:len(whole_line)))) /= 0) then
      read(whole_line(ii:len(whole_line)),*)ic
      write(tag,"(i20)")ic
      tag = "component "//trim(adjustl(tag))

      icomp(ip) = ic
      ncomp_surfs(ic) = ncomp_surfs(ic) + 1
      tag_list(ic) = tag
      num_components = max(num_components,ic)
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


  deallocate(nkg,njg,nig); nalloc = nalloc - 3

  return
  end subroutine PreReadUNCLEM

  end module PreReadUNCLEMM
