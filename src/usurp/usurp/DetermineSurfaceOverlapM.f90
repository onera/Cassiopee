! ================================================================
  module DetermineSurfaceOverlapM
! ================================================================

  implicit none

  private
  public :: DetermineSurfaceOverlap

  contains

! =========================================================
  subroutine DetermineSurfaceOverlap(p,run_gpc)
! =========================================================
! Check for surface overlap using a two-layer sytem of boxes
! where the first layer checks whether two "patches" overlap
! and (if so) the second layer checks all panel pair 
! combinations within those two patches.  This method serves
! as a "safe" alternative to the R-tree-based method.  Both
! methods should return the same number of intersecting boxes
! but may return them in a different order.

  use Intrtype             ,only: rs
  use Types                ,only: panel,surf
  use OverlappedM          ,only: Overlapped
  use PatchInfo            ,only: num_patches
  use ProcessPairProcedures,only: times_called,ProcessPair
  use TimeKeeper           ,only: gpc_time
  use UserInput            ,only: disjoin
  use VertexData           ,only: xv

  implicit none

  type(panel),pointer,dimension (:) :: p
  type(panel),pointer               :: p1,p2

  real(kind=rs)      :: time43
  integer            :: ip,ip1,ip2,ipan1,ipan2,idir
  logical            :: ran_gpc
  logical,intent(in) :: run_gpc


  continue


  nullify(p1,p2)

  gpc_time = 0.0_rs
  time43   = 0.0_rs


! compute a bounding box for every patch
  do ip = 1,num_patches
    do idir = 1,3
      surf(ip)%box(idir,1) = &
          minval(xv(idir,surf(ip)%first_node:surf(ip)%last_node))
      surf(ip)%box(idir,2) = &
          maxval(xv(idir,surf(ip)%first_node:surf(ip)%last_node))
    end do
  end do


! check every pair of patches for overlap
  do ip1 = 1,num_patches-1
    do ip2 = ip1+1,num_patches

      if (overlapped(surf(ip1)%box,surf(ip2)%box)) then
!        write(*,*)"patches ",ip1,ip2

!        the patches might overlap, so test the facets
         do ipan1 = surf(ip1)%first_panel,surf(ip1)%last_panel
           do ipan2 = surf(ip2)%first_panel,surf(ip2)%last_panel

             if ((overlapped(p(ipan1)%bound,p(ipan2)%bound)).and. &
                ((p(ipan1)%icomp == p(ipan2)%icomp).or.(.not.disjoin))) then

               p1 => p(ipan1)
               p2 => p(ipan2)

               times_called = times_called + 1
               call ProcessPair(p1,p2,time43,run_gpc,ran_gpc)

               if (ran_gpc) gpc_time = gpc_time + time43
             end if

           end do
         end do

       end if

     end do
  end do

  nullify(p1,p2)

  return
  end subroutine DetermineSurfaceOverlap

  end module DetermineSurfaceOverlapM
