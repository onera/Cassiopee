! ================================================================
  module WriteTriQM
! ================================================================

  implicit none

  private
  public :: WriteTriQ

  contains

! ======================================================================
  subroutine WriteTriQ(p)
! ======================================================================
! This procedure writes the triangulated surface data to a
! grid.i.tri or grid.i.triq file, with a debugging copy written
! in a Tecplot format.
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: panel
  use UserInput   ,only: cpufmt,full_surface,pathname,plotpanel,showcpu
  use UserInput   ,only: tecbinary,ttype,verbosity
  use VertexData  ,only: nqv,num_nodes,qv,xv,xv2,qv2,new_nodes,xv2o
  use VertexData  ,only: num_used,used,idup,iassoc

  implicit none

  type(panel),pointer,dimension (:) :: p

! variables used to debug flipped triangles
! real(kind=rd),dimension(3,3) :: xx
! real(kind=rd) :: v1x,v2x,v3x,v3
! real(kind=rd) :: v1y,v2y,v3y
! real(kind=rd) :: v1z,v2z,v3z,dot

  real(kind=rs),allocatable,dimension(:) :: v2
  real(kind=rs)                          :: time1,time2

  integer,allocatable,dimension(:,:) :: n2
  integer,allocatable,dimension(:)   :: varloc
  integer                            :: ip1,ip2,ipan,ic,iv,idir,ierr,i
  integer                            :: num_vertices,num_elements,ntri
! integer                            :: flipped=0

! integer :: ii  !debug flipped triangles

  character(len=132) :: tecvariables
  character(len=80)  :: zonestuff

  logical :: have_qv


  continue


! initialization
  call my_cpu_time(time1)
  have_qv = allocated(qv)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)


! assume that EVERY panel has been triangulated
! (this may be a problem if an error occurred in trying to 
! triangulate the polygon)

  num_elements = 0
  do ipan = ip1,ip2
    if (ttype == "GPC") then
      ntri = 0
      do ic = 1,p(ipan)%poly%num_cont
        ntri = ntri + p(ipan)%poly%cont(ic)%num_vert - 2
      end do
    else
      ntri = p(ipan)%poly%num_cont
    end if
    num_elements = num_elements + ntri
  end do

  num_vertices = num_nodes + new_nodes

! open wetted surface triangulation format
  if (have_qv) then
    open(unit=45,file=trim(pathname)//"grid.i.triq",action="write")
    write(45,*)num_used,num_elements,nqv
  else
    open(unit=45,file=trim(pathname)//"grid.i.tri",action="write")
    write(45,*)num_used,num_elements
  end if
  
! write the vertices to the grid.i.tri file in point format
  do iv = 1,num_nodes
     if (used(iv)) write(45,*)(real(xv(idir,iv)),idir=1,3)
  end do
  do iv = num_nodes+1,num_vertices
     if (used(iv)) write(45,*)(real(xv2(idir,iv)),idir=1,3)
  end do


! create space for the element list
  allocate(n2(3,num_elements),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate n2 in WriteTriQ"
    write(0,*)"num_elements = ",num_elements
    stop
  end if


! fill in the element list
  i = 0
  do ipan = ip1,ip2
    do ic = 1,p(ipan)%poly%num_cont
      if (ttype == "GPC") then
        ntri = p(ipan)%poly%cont(ic)%num_vert - 2
      else
        ntri = 1
      end if
      do iv = 1,ntri
        i = i + 1
        n2(1,i) = idup(p(ipan)%poly%cont(ic)%node(iv))
        if (mod(iv,2) == 1) then
          n2(2,i) = idup(p(ipan)%poly%cont(ic)%node(iv+1))
          n2(3,i) = idup(p(ipan)%poly%cont(ic)%node(iv+2))
        else
          n2(3,i) = idup(p(ipan)%poly%cont(ic)%node(iv+1))
          n2(2,i) = idup(p(ipan)%poly%cont(ic)%node(iv+2))
        end if
      end do
    end do
  end do
  if (i /= num_elements) then
    write(0,*)"ERROR! in element count in WriteTriQ"
    write(0,*)"i,num_elements = ",i,num_elements
    stop
  end if


  do iv = 1,num_elements
    write(45,*)n2(:,iv)
  end do


! generate component number data for each element
! (all triangles have the component number of their parent)
  do ipan = ip1,ip2
    if (ttype == "GPC") then
      ntri = 0
      do ic = 1,p(ipan)%poly%num_cont
        ntri = ntri + p(ipan)%poly%cont(ic)%num_vert-2
      end do
    else
      ntri = p(ipan)%poly%num_cont
    end if
    do iv = 1,ntri
    write(45,*)p(ipan)%icomp
    end do
  end do


  if (have_qv) then
!   write dependent variable data
    do iv = 1,num_nodes
       if (used(iv)) write(45,*)(real(qv(idir,iv)),idir=1,nqv)
    end do
    do iv = num_nodes+1,num_vertices
       if (used(iv)) write(45,*)(real(qv2(idir,iv)),idir=1,nqv)
    end do
  end if

  close(unit=45)


! ===================================
! WRITE A TECPLOT VERSION OF THE FILE
! ===================================

! if (verbosity >= 1) then
!   open Tecplot file and write header information

    if (have_qv) then
      tecvariables = "x,y,z,q5,panel"
    else
      tecvariables = "x,y,z,panel"
    end if

    if (tecbinary) then


    else
      open(unit=44,file=trim(pathname)//"usurp-triangles.dat", &
           action="write")
      write(44,*)"variables = "//trim(tecvariables)
      write(44,10)num_used,num_elements
 10   format(1x,'zone n=',i8,',e=',i8, &
                ',zonetype=fetriangle,datapacking=block')
      if (have_qv) then
        zonestuff = "  varlocation=([5]=cellcentered)"
      else
        zonestuff = "  varlocation=([4]=cellcentered)"
      end if
      write(44,*)zonestuff 
    end if !tec_binary

!   set v2 to store the vertex locations
    allocate(v2(num_used),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! allocate failed for v2 in WriteTriQ"
      write(0,*)"num_used = ",num_used
      stop
    end if


!   copy the needed vertex coordinates
    do idir = 1,3
      do iv = 1,num_nodes
        if (used(iv)) v2(idup(iv)) = xv(idir,iv)
      end do
      do iv = num_nodes + 1,num_vertices
        if (used(iv)) v2(idup(iv)) = xv2(idir,iv)
      end do

      if (tecbinary) then

      else 
        write(44,*)v2
      end if
    end do

    if (have_qv) then
!     write dependent variable to the Tecplot file in block format
!     write variables from original list
      idir = min(5,nqv)
      do iv = 1,num_nodes
        if (used(iv)) v2(idup(iv)) = qv(idir,iv)
      end do
      do iv = num_nodes+1,num_vertices
        if (used(iv)) v2(idup(iv)) = qv2(idir,iv)
      end do

      if (tecbinary) then

      else 
        write(44,*)v2
      end if

    end if !have_qv


!   reset v2 to store the element data
    deallocate(v2); nalloc = nalloc - 1
    allocate(v2(num_elements),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! allocated failed for v2 in WriteTriQ"
      write(0,*)"num_elements = ",num_elements
      stop
    end if

!   generate component number data for each element
!   (all triangles have the component number of their parent)
    i = 0
    do ipan = ip1,ip2

      if (ttype == "GPC") then
        ntri = 0
        do ic = 1,p(ipan)%poly%num_cont
          ntri = ntri + p(ipan)%poly%cont(ic)%num_vert-2
        end do
      else
        ntri = p(ipan)%poly%num_cont
      end if

      do iv = 1,ntri
        i = i + 1

!       check the dot product with the underlying panel
!       do idir = 1,3
!       do ii=1,num_vertices
!         if (n2(idir,i) == idup(ii)) then
!           if (ii <= num_nodes) then
!             xx(1:3,idir) = xv(1:3,ii)
!           else
!             xx(1:3,idir) = xv2(1:3,ii)
!           end if
!           exit
!         end if
!       end do
!       end do

!       v1x = xx(1,2) - xx(1,1)
!       v1y = xx(2,2) - xx(2,1)
!       v1z = xx(3,2) - xx(3,1)
!       v2x = xx(1,3) - xx(1,1)
!       v2y = xx(2,3) - xx(2,1)
!       v2z = xx(3,3) - xx(3,1)
!       v3x = v1y*v2z - v2y*v1z
!       v3y = v1z*v2x - v2z*v1x
!       v3z = v1x*v2y - v2x*v1y

!       dot = p(ipan)%xn(1)*v3x+p(ipan)%xn(2)*v3y+p(ipan)%xn(3)*v3z
!       v3 = sqrt(v3x**2 + v3y**2 + v3z**2)
!       if (ipan == plotpanel) write(0,*)"tri ",i,v3x/v3,v3y/v3,v3z/v3

!       if (dot < 0) then
!         if (ipan == plotpanel) then
!           write(0,*)"ERROR! triangle is flipped."
!           write(0,*)"ipan,iv = ",ipan,iv
!           do ii=1,num_vertices
!             if (n2(1,i) == idup(ii)) write(0,*)ii
!           end do
!           do ii=1,num_vertices
!             if (n2(2,i) == idup(ii)) write(0,*)ii
!           end do
!           do ii=1,num_vertices
!             if (n2(3,i) == idup(ii)) write(0,*)ii
!           end do
!         end if
!         flipped = flipped + 1
!         stop
!         if (flipped == 10) stop
!       end if

        v2(i) = p(ipan)%panel_number
!       v2(i) = p(ipan)%component

!       if (dot > 0.0) then
!         v2(i) = 1.0
!       else if (dot < 0.0) then
!         v2(i) = -1.0
!       else
!         v2(i) = 0.0
!       end if

      end do
    end do
    if (i /= num_elements) then
      write(0,*)"ERROR! in element count in WriteTriQ"
      write(0,*)"i,num_elements = ",i,num_elements
      stop
    end if
      
    if (tecbinary) then

    else
      write(44,*)v2
    end if

!   we are finished with v2
    deallocate(v2); nalloc = nalloc - 1

    if (tecbinary) then

    else
      do iv = 1,num_elements
        write(44,*)n2(:,iv)
      end do
    end if !tecbinary

    if (tecbinary) then

    else
      close(unit=44)
    end if

! end if !verbosity


! output a map of the triangulation
  if (full_surface) then
    open(unit=49,file=trim(pathname)//"usurp.map", &
         form="unformatted",action="write")
    write(49)num_nodes,new_nodes
    write(49)used
    do iv = num_nodes + 1,num_vertices
      if (used(iv)) write(49)xv2o(:,iv),iassoc(iv)
    end do
    close(49)
  end if


  deallocate(used); nalloc = nalloc - 1
  deallocate(idup); nalloc = nalloc - 1
  deallocate(n2)  ; nalloc = nalloc - 1


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in WriteTriQ: ',time2-time1
  end if

! if (verbosity >= 1) write(0,*)"number of reversed normals = ",flipped

  return
  end subroutine WriteTriQ

  end module WriteTriQM
