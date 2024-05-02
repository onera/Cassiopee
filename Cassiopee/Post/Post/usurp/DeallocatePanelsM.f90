! ================================================================
  module DeallocatePanelsM
! ================================================================

  implicit none

  private
  public :: DeallocatePanels

  contains

! ========================================================
  subroutine DeallocatePanels(p)
! ========================================================
! note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use EdgeData             ,only: edge
  use FomoArrays           ,only: cfpx,cfpy,cfpz,cmpx,cmpy,cmpz
  use FomoArrays           ,only: cfvx,cfvy,cfvz,cmvx,cmvy,cmvz
  use FomoArrays           ,only: cfmx,cfmy,cfmz,cmmx,cmmy,cmmz
  use FomoArrays           ,only: ctax,ctay,ctaz,weta,cfr
  use GroupInfo            ,only: group_data,group_iref
  use GroupInfo            ,only: group_name,group_size
  use IntrType             ,only: rs
  use my_cpu_timeM         ,only: my_cpu_time
  use NetAlloc             ,only: nalloc
  use PatchInfo            ,only: block,color,ibdir,icomp,rankp
  use PatchInfo            ,only: refc,ncomp_surfs,tag_list
  use PatchInfo            ,only: i1,i2,j1,j2,k1,k2
  use PriPairs             ,only: ipri
  use ProcessPairProcedures,only: DeallocatePolygon
  use Types                ,only: panel,polygon,scal,surf
  use UserInput            ,only: cpufmt,refq,showcpu
  use VertexData           ,only: xv,xv2,qv,qv2,iassoc,xv2o

  implicit none

  type(panel),pointer,dimension (:) :: p
  type(polygon),pointer             :: q
  real(kind=rs)                     :: time1,time2
  integer                           :: ip1,ip2,ipan


  continue


  call my_cpu_time(time1)

! panel data
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  do ipan = ip1,ip2
    q => p(ipan)%poly
    call DeallocatePolygon(q)
  end do
  deallocate(p); nalloc = nalloc - 1


! vertex data
  if (allocated(xv)) then
    deallocate(xv); nalloc = nalloc - 1
  end if
  if (allocated(qv)) then
    deallocate(qv); nalloc = nalloc - 1
  end if
  if (allocated(xv2)) then
    deallocate(xv2); nalloc = nalloc - 1
  end if
  if (allocated(xv2o)) then
    deallocate(xv2o); nalloc = nalloc - 1
  end if
  if (allocated(qv2)) then
    deallocate(qv2); nalloc = nalloc - 1
  end if
  if (allocated(iassoc)) then
    deallocate(iassoc); nalloc = nalloc - 1
  end if
  if (allocated(scal)) then
    deallocate(scal); nalloc = nalloc - 1
  end if


! patch info
  if (allocated(i1)) then
    deallocate(i1); nalloc = nalloc - 1
  end if
  if (allocated(i2)) then
    deallocate(i2); nalloc = nalloc - 1
  end if
  if (allocated(j1)) then
    deallocate(j1); nalloc = nalloc - 1
  end if
  if (allocated(j2)) then
    deallocate(j2); nalloc = nalloc - 1
  end if
  if (allocated(k1)) then
    deallocate(k1); nalloc = nalloc - 1
  end if
  if (allocated(k2)) then
    deallocate(k2); nalloc = nalloc - 1
  end if
  if (allocated(ibdir)) then
    deallocate(ibdir); nalloc = nalloc - 1
  end if
  if (allocated(icomp)) then
    deallocate(icomp); nalloc = nalloc - 1
  end if
  if (allocated(block)) then
    deallocate(block); nalloc = nalloc - 1
  end if
  if (allocated(color)) then
    deallocate(color); nalloc = nalloc - 1
  end if
  if (allocated(rankp)) then
    deallocate(rankp); nalloc = nalloc - 1
  end if
  if (allocated(refc)) then
    deallocate(refc); nalloc = nalloc - 1
  end if
  if (allocated(ncomp_surfs)) then
    deallocate(ncomp_surfs); nalloc = nalloc - 1
  end if
  if (allocated(tag_list)) then
    deallocate(tag_list); nalloc = nalloc - 1
  end if


! patch info (in module Types)
  if (allocated(surf)) then
    deallocate(surf); nalloc = nalloc - 1
  end if


! edge info (in module EdgeData)
  if (allocated(edge)) then
    deallocate(edge); nalloc = nalloc - 1
  end if


! fomo arrays 
  if (allocated(cfpx))  then
    deallocate(cfpx); nalloc = nalloc - 1
  end if
  if (allocated(cfpy))  then
    deallocate(cfpy); nalloc = nalloc - 1
  end if
  if (allocated(cfpz))  then
    deallocate(cfpz); nalloc = nalloc - 1
  end if
  if (allocated(cfvx))  then
    deallocate(cfvx); nalloc = nalloc - 1
  end if
  if (allocated(cfvy))  then
    deallocate(cfvy); nalloc = nalloc - 1
  end if
  if (allocated(cfvz))  then
    deallocate(cfvz); nalloc = nalloc - 1
  end if
  if (allocated(cmpx))  then
    deallocate(cmpx); nalloc = nalloc - 1
  end if
  if (allocated(cmpy))  then
    deallocate(cmpy); nalloc = nalloc - 1
  end if
  if (allocated(cmpz))  then
    deallocate(cmpz); nalloc = nalloc - 1
  end if
  if (allocated(cmvx))  then
    deallocate(cmvx); nalloc = nalloc - 1
  end if
  if (allocated(cmvy))  then
    deallocate(cmvy); nalloc = nalloc - 1
  end if
  if (allocated(cmvz))  then
    deallocate(cmvz); nalloc = nalloc - 1
  end if
  if (allocated(weta))  then
    deallocate(weta); nalloc = nalloc - 1
  end if
  if (allocated(ctax))  then
    deallocate(ctax); nalloc = nalloc - 1
  end if
  if (allocated(ctay))  then
    deallocate(ctay); nalloc = nalloc - 1
  end if
  if (allocated(ctaz))  then
    deallocate(ctaz); nalloc = nalloc - 1
  end if
  if (allocated(cfmx))  then
    deallocate(cfmx); nalloc = nalloc - 1
  end if
  if (allocated(cfmy))  then
    deallocate(cfmy); nalloc = nalloc - 1
  end if
  if (allocated(cfmz))  then
    deallocate(cfmz); nalloc = nalloc - 1
  end if
  if (allocated(cmmx))  then
    deallocate(cmmx); nalloc = nalloc - 1
  end if
  if (allocated(cmmy))  then
    deallocate(cmmy); nalloc = nalloc - 1
  end if
  if (allocated(cmmz))  then
    deallocate(cmmz); nalloc = nalloc - 1
  end if
  if (allocated(cfr))  then
    deallocate(cfr); nalloc = nalloc - 1
  end if

! GroupInfo
  if (allocated(group_name)) then
    deallocate(group_name); nalloc = nalloc - 1
  end if
  if (allocated(group_size)) then
    deallocate(group_size); nalloc = nalloc - 1
  end if
  if (allocated(group_iref)) then
    deallocate(group_iref); nalloc = nalloc - 1
  end if
  if (allocated(group_data)) then
    deallocate(group_data); nalloc = nalloc - 1
  end if

! Priority Pairs
  if (allocated(ipri)) then
    deallocate(ipri); nalloc = nalloc - 1
  end if

! UserInfo
  if (allocated(refq)) then
    deallocate(refq); nalloc = nalloc - 1
  end if


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in DeallocatePanels:',time2-time1
  end if


  return
  end subroutine DeallocatePanels

  end module DeallocatePanelsM
