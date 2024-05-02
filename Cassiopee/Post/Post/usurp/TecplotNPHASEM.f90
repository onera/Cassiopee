! ================================================================
  module TecplotNPHASEM
! ================================================================

  implicit none

  private
  public :: TecplotNPHASE

  contains

! ======================================================================
  subroutine TecplotNPHASE(p)
! ======================================================================
! Procedure to convert polygons to triangles for any panels that have
! been partially eliminated, and to write the triangles and associated
! scalar data to the Tecplot output file.  This procedure differs from
! "TecplotTriangles" in that it attempts to eliminate duplicate nodes
! between triangles, thus generating a "more watertight" surface.  The
! surface is not strictly watertight due to "hanging nodes" that are
! created in some panels.
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: num_scalars
  use RotateBackM ,only: RotateBack
  use Types       ,only: panel,scal
  use UserInput   ,only: cpufmt,merge_tol2,pathname,showcpu,ttype
  use VertexData  ,only: xv,num_nodes

  implicit none

  type(panel),pointer,dimension (:) :: p

  real(kind=rd),allocatable,dimension (:,:) :: xv2
  real(kind=rd),dimension(3)                :: xt
  real(kind=rd)                             :: dist2
  real(kind=rs)                             :: time1,time2

  integer,dimension(4) :: node
  integer              :: ipan,ip1,ip2,ic,iv,idir,i,j
  integer              :: merged_nodes,ntri,mergept
  integer              :: num_elements,num_vertices,ierr
  integer              :: max_vertices
  character(len=80)    :: zonestuff


  continue


! initialization
  call my_cpu_time(time1)
  num_elements = 0
  merged_nodes = 0
  max_vertices = num_nodes
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)


  do ipan = ip1,ip2

    if (p(ipan)%itri == 1) then
!     this panel was already triangulated in CreateTriangles

!     count the vertices and elements from the triangulation
      do ic = 1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          num_elements = num_elements + p(ipan)%poly%cont(ic)%num_vert-2
        else if (ttype == "Triangle") then
          num_elements = num_elements + 1
        end if

!       count only new vertices
        do iv = 1,p(ipan)%poly%cont(ic)%num_vert
          if (p(ipan)%poly%cont(ic)%node(iv) <= 0) then

!           try one last time to merge this with a pre-existing point
            call RotateBack(p(ipan)%poly%cont(ic)%x(iv), &
                            p(ipan)%poly%cont(ic)%y(iv), &
                            p(ipan)%g,p(ipan)%xmid,xt)
            mergept = 0
            MERGE1: do i = 1,num_nodes
              dist2 = (xv(1,i)-xt(1))**2 &
                    + (xv(2,i)-xt(2))**2 &
                    + (xv(3,i)-xt(3))**2
              if (dist2 <= merge_tol2) then
                mergept = i
                exit MERGE1
              end if
            end do MERGE1

            if (mergept /= 0) then
              p(ipan)%poly%cont(ic)%node(iv) = mergept
              merged_nodes = merged_nodes + 1
            else
              max_vertices = max_vertices + 1
            end if
          end if
        end do
      end do

    else if (p(ipan)%ratio >= 0.999) then

!     this is an unbroken element, which adds no new vertices
      num_elements = num_elements + 1
    end if
  end do


! allocate space for new vertices
  allocate(xv2(3,num_nodes+1:max_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)'ERROR! allocate failed for xv2 in TecplotNPHASE.'
    write(0,*)'size was ',max_vertices-num_nodes
  end if


! generate triangle vertices
  num_vertices = num_nodes
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic = 1,p(ipan)%poly%num_cont
      do iv = 1,p(ipan)%poly%cont(ic)%num_vert

        if (p(ipan)%poly%cont(ic)%node(iv) <= 0) then

          call RotateBack(p(ipan)%poly%cont(ic)%x(iv), &
                          p(ipan)%poly%cont(ic)%y(iv), &
                          p(ipan)%g,p(ipan)%xmid,xt)

!         try to merge this point with one of the other new points
          mergept = 0
          MERGE2: do j = num_nodes + 1,num_vertices
            dist2 = (xv2(1,j)-xt(1))**2 &
                  + (xv2(2,j)-xt(2))**2 &
                  + (xv2(3,j)-xt(3))**2
            if (dist2 <= merge_tol2) then
              mergept = j
              exit MERGE2
            end if
          end do MERGE2

          if (mergept /= 0) then
            merged_nodes = merged_nodes + 1
            p(ipan)%poly%cont(ic)%node(iv) = mergept
          else
            num_vertices = num_vertices + 1
            xv2(1:3,num_vertices) = xt(1:3)
            p(ipan)%poly%cont(ic)%node(iv) = num_vertices
          end if

        end if

      end do
      end do
    end if
  end do


! open Tecplot file and write header information
  open(unit=42,file=trim(pathname)//"usurp-nphase.dat",action="write")
  write(42,*)'variables = x,y,z,iblank,ratio,ipan,component'
  do i = 1,num_scalars
    write(42,'(a2,i2.2)')',q',i
  end do
  write(42,10)num_vertices,num_elements
  10 format(1x,'zone n=',i6,',e=',i6, &
            ',zonetype=fequadrilateral,datapacking=block')
  zonestuff = "  varlocation=([4-  ]=cellcentered)"
  write(zonestuff(19:20),'(i2.2)')7+num_scalars
  write(42,*)zonestuff 


! write vertices in block format
  do idir = 1,3

!   write vertices from original list
    do iv = 1,num_nodes
      write(42,*)xv(idir,iv)
    end do

!   write vertices from newly generated triangles
    do iv = num_nodes+1,num_vertices
      write(42,*)xv2(idir,iv)
    end do

  end do


! generate iblank data for each element
! (all triangles have iblank=1)
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      if (ttype == "GPC") then
        write(42,*)((1,iv=1,p(ipan)%poly%cont(ic)%num_vert-2), &
                     ic=1,p(ipan)%poly%num_cont)
      else if (ttype == "Triangle") then
        write(42,*)(1,ic=1,p(ipan)%poly%num_cont)
      end if
    else if (p(ipan)%ratio >= 0.999) then
      write(42,*)p(ipan)%iblank
    end if
  end do


! generate ratio data for each element
! (all triangles have ratio=1)
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      if (ttype == "GPC") then
        write(42,*)((1.0,iv=1,p(ipan)%poly%cont(ic)%num_vert-2), &
                       ic=1,p(ipan)%poly%num_cont)
      else if (ttype == "Triangle") then
        write(42,*)(1.0,ic=1,p(ipan)%poly%num_cont)
      end if
    else if (p(ipan)%ratio >= 0.999) then
      write(42,*)p(ipan)%ratio
    end if
  end do


! generate panel number data for each element
! (all triangles have the panel number of their parent)
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      if (ttype == "GPC") then
        write(42,*)((ipan,iv=1,p(ipan)%poly%cont(ic)%num_vert-2), &
                        ic=1,p(ipan)%poly%num_cont)
      else if (ttype == "Triangle") then
        write(42,*)(ipan,ic=1,p(ipan)%poly%num_cont)
      end if
    else if (p(ipan)%ratio >= 0.999) then
      write(42,*)ipan
    end if
  end do


! generate component number data for each element
! (all triangles have the component number of their parent)
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      if (ttype == "GPC") then
        write(42,*)((p(ipan)%icomp,                 &
                   iv=1,p(ipan)%poly%cont(ic)%num_vert-2), &
                   ic=1,p(ipan)%poly%num_cont)
      else if (ttype == "Triangle") then
        write(42,*)(p(ipan)%icomp,                 &
                   ic=1,p(ipan)%poly%num_cont)
      end if
    else if (p(ipan)%ratio >= 0.999) then
      write(42,*)p(ipan)%icomp
    end if
  end do


! generate scalar data for each element
! (all triangles have the scalar data of their parent)
  do i = 1,num_scalars
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      if (ttype == "GPC") then
        write(42,*)((scal(i,ipan), &
                   iv=1,p(ipan)%poly%cont(ic)%num_vert-2), &
                   ic=1,p(ipan)%poly%num_cont)
      else if (ttype == "Triangle") then
        write(42,*)(scal(i,ipan), &
                   ic=1,p(ipan)%poly%num_cont)
      end if
    else if (p(ipan)%ratio >= 0.999) then
      write(42,*)scal(i,ipan)
    end if
  end do
  end do


! generate element list
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic = 1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          ntri = p(ipan)%poly%cont(ic)%num_vert-2
        else
          ntri = 1
        end if
        do iv = 1,ntri
          node(1) = abs(p(ipan)%poly%cont(ic)%node(iv))
          if (mod(iv,2) == 1) then
            node(2) = abs(p(ipan)%poly%cont(ic)%node(iv+1))
            node(3) = abs(p(ipan)%poly%cont(ic)%node(iv+2))
            node(4) = abs(p(ipan)%poly%cont(ic)%node(iv+2))
          else
            node(3) = abs(p(ipan)%poly%cont(ic)%node(iv+1))
            node(4) = abs(p(ipan)%poly%cont(ic)%node(iv+1))
            node(2) = abs(p(ipan)%poly%cont(ic)%node(iv+2))
          end if

          write(42,*)node
        end do
      end do
    else if (p(ipan)%ratio >= 0.999) then
      if (p(ipan)%num_vertices == 3) then
        write(42,*)(p(ipan)%node(iv),iv=1,p(ipan)%num_vertices), &
                    p(ipan)%node(3)
      else
        write(42,*)(p(ipan)%node(iv),iv=1,p(ipan)%num_vertices)
      end if
    end if
  end do
  write(0,*)'number of merged nodes = ',merged_nodes

  close(unit=42)

  deallocate(xv2); nalloc = nalloc - 1

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in TecplotNPHASE:',time2-time1
  end if


  return
  end subroutine TecplotNPHASE

  end module TecplotNPHASEM
