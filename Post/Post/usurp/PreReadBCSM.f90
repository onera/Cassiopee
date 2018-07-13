! ================================================================
  module PreReadBCSM
! ================================================================

  implicit none

  private
  public :: PreReadBCS

  contains

! ===================================================================
  subroutine PreReadBCS()
! ===================================================================
! Read the boundary condition file (.bcs) for CFD-SHIP in 
! order to acquire the number of panels, patches, and vertices, and
! to store the patch data.

  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: i1,i2,j1,j2,k1,k2,block,icomp
  use PatchInfo                ,only: ncomp_surfs,num_patches,tag_list
  use PatchInfo                ,only: max_components
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch
  use StructuredPatchProcedures,only: EchoPatchData
  use Types                    ,only: surf
  use UserInput                ,only: workdir,fnameo

  implicit none

  integer :: i,ip,n,ierr
  integer :: NumBoundarySurfaces,bctype
  integer :: FirstDerivativeFlag
  integer :: i1s,i2s,j1s,j2s,k1s,k2s
  integer :: DonorMesh,NormalDirection


  continue


! initialization
  allocate(tag_list(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate tag_list in PreReadBCS"
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



! first read to ascertain number of patches and components
  write(*,*)
  write(unit=*,fmt=*)"Reading CFD-SHIP BCS file ..."
  open(unit=3, &
       file=trim(workdir)//trim(fnameo)//".bcs", &
       status="old", &
       action="read", &
       position="rewind")
  do
    read(3,*,end=20,err=20)NumBoundarySurfaces
    do i = 1,NumBoundarySurfaces

      read(3,*)bctype

      if ((bctype == 20).or. &
          (bctype == 22).or. &
          (bctype == 24)) num_patches = num_patches + 1
      read(3,*)NormalDirection
      read(3,*)i1s,i2s
      read(3,*)j1s,j2s
      read(3,*)k1s,k2s
      read(3,*)FirstDerivativeFlag

      if (bctype == 30) then
        read(3,*)
      end if

      if ((bctype == 92).or.(bctype == 52)) then
        read(3,*)DonorMesh
        read(3,*)NormalDirection
        read(3,*)i1s,i2s
        read(3,*)j1s,j2s
        read(3,*)k1s,k2s
      end if
    end do
  end do

 20 continue


  call AllocateStructuredPatch()


! second read to fill in data
  ip = 0
  n = 0
  rewind(unit=3)
  do
    read(3,*,end=30,err=30)NumBoundarySurfaces
    n = n + 1
    do i = 1,NumBoundarySurfaces
      read(3,*)bctype
      read(3,*)NormalDirection
      read(3,*)i1s,i2s
      read(3,*)j1s,j2s
      read(3,*)k1s,k2s
      if ((bctype == 20).or.(bctype == 22).or.(bctype == 24)) then
        ip = ip + 1
        i1(ip) = i1s
        i2(ip) = i2s
        j1(ip) = j1s
        j2(ip) = j2s
        k1(ip) = k1s
        k2(ip) = k2s
        block(ip) = n

!       untagged
        icomp(ip) = 0
        ncomp_surfs(0) = ncomp_surfs(0) + 1

!       count panels
        call CountStructuredPatch(ip)

      end if
      read(3,*)FirstDerivativeFlag
      if (bctype == 30) then
        read(3,*)
      end if
      if ((bctype == 92).or.(bctype == 52)) then
        read(3,*)DonorMesh
        read(3,*)NormalDirection
        read(3,*)i1s,i2s
        read(3,*)j1s,j2s
        read(3,*)k1s,k2s
      end if
    end do
  end do
  30 continue
  close(unit=3)


! screen output to confirm results
  call EchoPatchData()


  return
  end subroutine PreReadBCS

  end module PreReadBCSM
