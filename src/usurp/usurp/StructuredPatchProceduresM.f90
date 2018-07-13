! ==============================================================
  module StructuredPatchProcedures
! ==============================================================
! This module groups together the routines having to do with
! setting up the input for structured patch grids.

  implicit none

  public :: AllocateStructuredPatch
  public :: CountStructuredPatch
  public :: EchoPatchData

  contains

! ===================================
  subroutine CountStructuredPatch(ip)
! ===================================
! count the nodes, panels, and edges for the structured patch grids

  use EdgeData   ,only: num_edges
  use PatchInfo  ,only: num_panels,i1,i2,j1,j2,k1,k2
  use Types      ,only: surf
  use VertexData ,only: num_nodes

  implicit none

  integer,intent(in) :: ip


  continue


  surf(ip)%first_panel = num_panels + 1
  surf(ip)%first_node  = num_nodes  + 1
  surf(ip)%first_edge  = num_edges  + 1

  if (k1(ip) == k2(ip)) then
    num_panels = num_panels + (i2(ip)-i1(ip)  )*(j2(ip)-j1(ip)  )
    num_nodes  = num_nodes  + (i2(ip)-i1(ip)+1)*(j2(ip)-j1(ip)+1)
    num_edges  = num_edges  + (i2(ip)-i1(ip)+1)*(j2(ip)-j1(ip)  ) &
                            + (i2(ip)-i1(ip)  )*(j2(ip)-j1(ip)+1)
  else if (i1(ip) == i2(ip)) then
    num_panels = num_panels + (j2(ip)-j1(ip)  )*(k2(ip)-k1(ip)  )
    num_nodes  = num_nodes  + (j2(ip)-j1(ip)+1)*(k2(ip)-k1(ip)+1)
    num_edges  = num_edges  + (j2(ip)-j1(ip)+1)*(k2(ip)-k1(ip)  ) &
                            + (j2(ip)-j1(ip)  )*(k2(ip)-k1(ip)+1)
  else if (j1(ip) == j2(ip)) then
    num_panels = num_panels + (k2(ip)-k1(ip)  )*(i2(ip)-i1(ip)  )
    num_nodes  = num_nodes  + (k2(ip)-k1(ip)+1)*(i2(ip)-i1(ip)+1)
    num_edges  = num_edges  + (k2(ip)-k1(ip)+1)*(i2(ip)-i1(ip)  ) &
                            + (k2(ip)-k1(ip)  )*(i2(ip)-i1(ip)+1)
  end if

  surf(ip)%last_panel = num_panels
  surf(ip)%last_node  = num_nodes
  surf(ip)%last_edge  = num_edges

  return
  end subroutine CountStructuredPatch


! ===================================================================
  subroutine AllocateStructuredPatch()
! ===================================================================
! allocate the basic arrays for the structured patch grids

  use NetAlloc ,only: nalloc
  use PatchInfo,only: i1,i2,j1,j2,k1,k2,block,icomp,num_patches
  use Types    ,only: surf

  implicit none

  integer :: ierr


  continue


  allocate(   i1(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for i1"
    stop
  end if

  allocate(   i2(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for i2"
    stop
  end if

  allocate(   j1(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for j1"
    stop
  end if

  allocate(   j2(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for j2"
    stop
  end if

  allocate(   k1(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for k1"
    stop
  end if

  allocate(   k2(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for k2"
    stop
  end if

  allocate(block(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for block"
    stop
  end if

  allocate(icomp(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for icomp"
    stop
  end if

  allocate(surf(num_patches),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for surf"
    stop
  end if

  return
  end subroutine AllocateStructuredPatch

! ==================================================================
  subroutine EchoPatchData()
! ==================================================================
! write summary patch info to standard out

  use EdgeData  , only: num_edges
  use PatchInfo , only: num_panels,num_patches,num_scalars
  use PatchInfo , only: num_components,tag_list
  use VertexData, only: num_nodes

  implicit none

  integer :: ic


  continue


  write(unit=*,fmt="(1x,a,i8)")"number of surface panels = ",num_panels
  write(unit=*,fmt="(1x,a,i8)")"number of vertex nodes = ",num_nodes
  write(unit=*,fmt="(1x,a,i8)")"number of edges = ",num_edges
  write(unit=*,fmt="(1x,a,i2)")"number of scalars = ",num_scalars

  write(unit=*,fmt="(1x,a,i2)") &
       "number of named geometric components = ",num_components

  do ic = 1,num_components
    write(unit=*,fmt="(1x,i2,1x,a)")ic," "//trim(tag_list(ic))
  end do

  write(unit=*,fmt=*)

  end subroutine EchoPatchData

  end module StructuredPatchProcedures
