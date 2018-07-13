! ================================================================
  module TecplotTrianglesM
! ================================================================

  implicit none

  private
  public :: TecplotTriangles

  contains

! ======================================================================
  subroutine TecplotTriangles(p)
! ======================================================================
! For any panels that have been partially eliminated, write the 
! triangles and associated scalar data to the Tecplot output file.
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: num_patches
  use Types       ,only: panel,scal
  use UserInput   ,only: colormap,cpufmt,plotvel,showcpu
  use UserInput   ,only: solution_exists,tecbinary,ttype
  use VertexData  ,only: xv,xv2,num_nodes,new_nodes

  implicit none

  type(panel),pointer,dimension (:)       :: p

  real(kind=rs),allocatable,dimension (:) :: v2
  real(kind=rs)                           :: time1,time2

  integer,allocatable,dimension (:,:)     :: n2
  integer,allocatable,dimension (:)       :: varloc
  integer,allocatable,dimension (:)       :: index

  integer :: i,ierr,ipan,ip1,ip2,ic,iv,idir
  integer :: num_elements,num_vertices,ntri
#ifdef USE_TECPLOT
  integer,parameter :: visdouble=0
  integer,parameter :: zero=0
  integer :: TecDat100,TecEnd100,TecFil100,TecNod100,TecZne100
#endif


  continue


  call my_cpu_time(time1)

  if (ttype /= "none") then

! initialization
  num_elements = 0
  num_vertices = 0
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)


! create an array to keep track of which vertices are used
  allocate(index(1:num_nodes+new_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)'ERROR! allocate failed for index in TecplotTriangles.'
    write(0,*)'num_nodes = ',num_nodes
    write(0,*)'new_nodes = ',new_nodes
    stop
  end if
  index = 0


! count any triangles and keep track of which vertices are needed
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic = 1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          num_elements = num_elements + p(ipan)%poly%cont(ic)%num_vert-2
        else if (ttype == "Triangle") then
          num_elements = num_elements + 1
        end if
        do iv = 1,p(ipan)%poly%cont(ic)%num_vert
          i = p(ipan)%poly%cont(ic)%node(iv)
          if ((i <= 0).or.(i > num_nodes + new_nodes)) then
            write(0,*)"ERROR! invalid node value in TecplotTriangles"
            write(0,*)"ipan,ic,iv,node = ",ipan,ic,iv,i
            stop
          else
            index(i) = 1
          end if
        end do
      end do
    end if
  end do


! add up the total number of vertices used and assign new values to them
  num_vertices = 0
  do iv = 1,num_nodes + new_nodes
    if (index(iv) == 1) then
      num_vertices = num_vertices + 1
      index(iv) = num_vertices
    end if
  end do


! if there are no triangles, close the Tecplot file and return
  if (num_elements == 0) then
    if (tecbinary) then
#ifdef USE_TECPLOT
      ierr = TecFil100(1)
      ierr = Tecend100()
      if (colormap) close(unit=43)
#endif
    else
      close(unit=42)
    end if
    return
  end if


! generate Tecplot zone header information
  if (solution_exists) then
    if (plotvel) then
      allocate(varloc(12)); nalloc = nalloc + 1
      varloc = (/ 1,1,1,0,0,0,0,0,0,0,0,0 /)
    else
      allocate(varloc(8)); nalloc = nalloc + 1
      varloc = (/ 1,1,1,0,0,0,0,0 /)
    end if
  else
    allocate(varloc(7)); nalloc = nalloc + 1
    varloc = (/ 1,1,1,0,0,0,0 /)
  end if

  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecFil100(1)
    ierr = Teczne100 ( " "//char(0), &
                        2, &              ! <--- zonetype=fetriangle
                        num_vertices,num_elements,0, &
                        0,0,0, &
                        1,0,0, &
                        varloc, &
                        %val(zero), &
                        0)
    if (colormap) then
      write(43,"(a,i4,a)") &
        "$!Field [",num_patches+1,"] Mesh{Color = Black}"
      write(43,"(a,i4,a)") &
        "$!Field [",num_patches+1,"] Shade{Color = Black}"
    end if
#endif
  else

    write(42,10)num_vertices,num_elements
 10 format(1x,'zone n=',i6,',e=',i6,',zonetype=fetriangle,datapacking=block')
    if (solution_exists) then
      if (plotvel) then
        write(42,*) "  varlocation=([4-12]=cellcentered)"
      else
        write(42,*) "  varlocation=([4-8]=cellcentered)"
      end if
    else
      write(42,*) "  varlocation=([4-7]=cellcentered)"
    end if
    if (colormap) write(42,*) "  c=black"
  end if


! set v2 to store the vertex locations
  allocate(v2(num_vertices),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for v2 in TecplotTriangles"
    write(0,*)"num_vertices = ",num_vertices
    stop
  end if


! copy the needed vertex coordinates
  do idir = 1,3
    do iv = 1,num_nodes
      i = index(iv)
      if (i /= 0) v2(i) = xv(idir,iv)
    end do
    do iv = num_nodes + 1,num_nodes+new_nodes
      i = index(iv)
      if (i /= 0) v2(i) = xv2(idir,iv)
    end do
    
!   write the triangle vertices in block format
    if (tecbinary) then
#ifdef USE_TECPLOT
      ierr = TecDat100(num_vertices,v2,visdouble)
#endif
    else
      write(42,*)v2
    end if
  end do

! reset v2 to store the element data
  deallocate(v2); nalloc = nalloc - 1
  allocate(v2(num_elements),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate failed for v2 in TecplotTriangles"
    write(0,*)"num_elements = ",num_elements
    stop
  end if

  if (solution_exists) then
! generate pressure data for each element
! (each triangle has the pressure of its parent)
    i = 0
    do ipan = ip1,ip2
      if (p(ipan)%itri == 1) then
        do ic=1,p(ipan)%poly%num_cont
          if (ttype == "GPC") then
            ntri = p(ipan)%poly%cont(ic)%num_vert-2
          else
            ntri = 1
          end if
          do iv=1,ntri
            i = i + 1
            v2(i) = 2.0*scal(2,ipan)
          end do
        end do
      end if
    end do

    if (tecbinary) then
#ifdef USE_TECPLOT
      ierr = TecDat100(num_elements,v2,visdouble)
#endif
    else
      write(42,*)v2
    end if
  end if


! generate iblank data for each element
! (all triangles have iblank=1)
  v2 = 1.0
  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecDat100(num_elements,v2,visdouble)
#endif
  else
    write(42,*)v2
  end if


! generate ratio data for each element
! (all triangles have ratio=1)
! v2 = 1.0    !===> still equal to 1.0 because of previous iblank data
  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecDat100(num_elements,v2,visdouble)
#endif
  else
    write(42,*)v2
  end if


! generate panel number data for each element
! (each triangle has the panel number of its parent)
  i = 0
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic=1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          ntri = p(ipan)%poly%cont(ic)%num_vert-2
        else
          ntri = 1
        end if
        do iv=1,ntri
          i = i + 1
          v2(i) = ipan
        end do
      end do
    end if
  end do
  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecDat100(num_elements,v2,visdouble)
#endif
  else
    write(42,*)v2
  end if


! generate component number data for each element
! (each triangle has the component number of its parent)
  i = 0
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic=1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          ntri = p(ipan)%poly%cont(ic)%num_vert-2
        else
          ntri = 1
        end if
        do iv=1,ntri
          i = i + 1
          v2(i) = p(ipan)%icomp
        end do
      end do
    end if
  end do
  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecDat100(num_elements,v2,visdouble)
#endif
  else
    write(42,*)v2
  end if


  if (plotvel) then
    do idir = 1,4
! generate velocity and skin friction data for each element
! (each triangle shares the value of its parent)
    i = 0
    do ipan = ip1,ip2
      if (p(ipan)%itri == 1) then
        do ic=1,p(ipan)%poly%num_cont
          if (ttype == "GPC") then
            ntri = p(ipan)%poly%cont(ic)%num_vert-2
          else
            ntri = 1
          end if
          do iv=1,ntri
            i = i + 1
            v2(i) = scal(8+idir,ipan)
          end do
        end do
      end if
    end do

    if (tecbinary) then
#ifdef USE_TECPLOT
      ierr = TecDat100(num_elements,v2,visdouble)
#endif
    else
      write(42,*)v2
    end if

    end do !idir
  end if


! we are finished with v2
  deallocate(v2); nalloc = nalloc - 1


! create space for the element list
  allocate(n2(3,num_elements),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate n2 in TecplotTriangles"
    write(0,*)"num_elements = ",num_elements
    stop
  end if


! fill in the element list
  i = 0
  do ipan = ip1,ip2
    if (p(ipan)%itri == 1) then
      do ic = 1,p(ipan)%poly%num_cont
        if (ttype == "GPC") then
          ntri = p(ipan)%poly%cont(ic)%num_vert-2
        else
          ntri = 1
        end if
        do iv = 1,ntri
          i = i + 1
          n2(1,i) = index(p(ipan)%poly%cont(ic)%node(iv))
          if (mod(iv,2) == 1) then
            n2(2,i) = index(p(ipan)%poly%cont(ic)%node(iv+1))
            n2(3,i) = index(p(ipan)%poly%cont(ic)%node(iv+2))
          else
            n2(3,i) = index(p(ipan)%poly%cont(ic)%node(iv+1))
            n2(2,i) = index(p(ipan)%poly%cont(ic)%node(iv+2))
          end if
        end do
      end do
    end if
  end do

  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecNod100(n2)
#endif
  else
    do iv = 1,num_elements
      write(42,*)n2(1,iv),n2(2,iv),n2(3,iv)
    end do
  end if

  deallocate(n2)    ; nalloc = nalloc - 1
  deallocate(varloc); nalloc = nalloc - 1
  deallocate(index) ; nalloc = nalloc - 1

  end if !ttype /= none


  if (tecbinary) then
#ifdef USE_TECPLOT
    ierr = TecEnd100()
    if (colormap) close(unit=43)
#endif
  else
    close(unit=42)
  end if


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in TecplotTriangles:',time2-time1
  end if


  return
  end subroutine TecplotTriangles

  end module TecplotTrianglesM
