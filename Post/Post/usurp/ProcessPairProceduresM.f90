! ==============================================================
  module ProcessPairProcedures
! ==============================================================
! This module groups together the routines that directly process
! two overlapping panels by passing the data to the general
! polygon clipping (gpc) library


  implicit none

  integer,public :: times_called

  public  :: ProcessPair
  public  :: gpc_f_translator
  public  :: triangle_translator_f
  public  :: DeallocatePolygon
  private :: SetRank

  contains

! ==============================================================
  subroutine ProcessPair(p1,p2,time43,run_gpc,ran_gpc)
! ==============================================================
! This subroutine makes the final determination as to whether
! to go foward with the polygon boolean difference operation
! depending on the dot product of the two panel surface normals.
! It calls the gpc library if final clearance is obtained.
! Note that the dummy arguments p1 and p2 are both pointers,
! which means that an explicit interface is required.  Containing
! the procedure in a module satisfies this requirement.

  use ConvertPanelToPolygonM,only: ConvertPanelToPolygon
  use GarbageCollectorM     ,only: GarbageCollector
  use IntrType              ,only: rd,rs
  use my_cpu_timeM          ,only: my_cpu_time
  use NetAlloc              ,only: nalloc
  use PatchInfo             ,only: table,overlaps,firstpass
  use PatchInfo             ,only: num_panels,rankp
  use PrintPolyM            ,only: PrintPoly
  use PriPairs              ,only: ipri
  use Types                 ,only: panel,polygon
  use UserInput             ,only: basis,my_bit_size,pathname,plotpanel
  use UserInput             ,only: dotlimit,trap_clip
  use UserInput             ,only: overflow,use_priority,verbosity
  use VertexData            ,only: xv

  implicit none

  type(panel),pointer          :: p1,p2
  type(panel),pointer          :: ipans,ipanc
  type(polygon),pointer        :: subjpoly,q
  type(polygon),target         :: clippoly

  real(kind=rd),dimension(3,4) :: x
! real(kind=rd)                :: CalculatePolygonArea
  real(kind=rd)                :: dot
  real(kind=rs)                :: time3,time4
  real(kind=rs),intent(out)    :: time43

  integer                      :: rank1,rank2
  integer                      :: ic,is1,is2
  integer                      :: iv,idir
  integer                      :: k,kint,kpos
  logical,intent(in)           :: run_gpc
  logical,intent(out)          :: ran_gpc


  continue


  nullify(ipans,ipanc,subjpoly)
  ran_gpc = .false.


! ensure that panels are on the same side of the surface
  dot = p1%xn(1)*p2%xn(1) + p1%xn(2)*p2%xn(2) + p1%xn(3)*p2%xn(3)
  if (dot <= dotlimit) return

! if (dot <= 0.25) return !included angle must be less than 75 degrees
! if (dot <= 0.00) return !included angle must be less than 90 degrees


! skip if a non-communicating priority pair has been set
  if (overflow.and.use_priority) then
    is1 = p1%isurf
    is2 = p2%isurf
    if (ipri(is1,is2) < 0) return
  end if


! first pass: just build table
  if (firstpass) then

!   panel 1 overlaps surface 2
    k = (p2%isurf - 1)*num_panels + p1%panel_number
    kint = (k-1)/my_bit_size + 1
    kpos = mod(k-1,my_bit_size)
    overlaps(kint) = ibset(overlaps(kint),kpos)

!   panel 2 overlaps surface 1
    k = (p1%isurf - 1)*num_panels + p2%panel_number
    kint = (k-1)/my_bit_size + 1
    kpos = mod(k-1,my_bit_size)
    overlaps(kint) = ibset(overlaps(kint),kpos)

    return
  end if


! determine panel precedence based on iblank
  rank1 = SetRank(p1%iblank)
  rank2 = SetRank(p2%iblank)


! if iblank has not settled precedence, base precedence on 
! either panel area (basis = "panel") or the overlap table
! (basis = "patch")
  if (rank1 == rank2) then
    if (basis == "panel") then
      if (p1%area_ref < p2%area_ref) then
        rank1 = 1
        rank2 = 0
      else if (p2%area_ref <= p1%area_ref) then
        rank1 = 0
        rank2 = 1
      end if
    else ! if (basis == "patch") then
      is1 = p1%isurf
      is2 = p2%isurf
      rank1 = rankp(is1)
      rank2 = rankp(is2)
    end if

!   over-rule all of this if a priority pair has been set
    if (overflow.and.use_priority) then
      is1 = p1%isurf
      is2 = p2%isurf
      if (ipri(is1,is2) > 0) then
        rank1 = 1
        rank2 = 0
      else if (ipri(is2,is1) > 0) then
        rank1 = 0
        rank2 = 1
      end if
    end if
  end if


! set which polygon is the "subject" and which is the "clipper"
  if (rank2 > rank1) then
    ipanc => p2
    ipans => p1
  else if (rank1 > rank2) then
    ipanc => p1
    ipans => p2
  else
    write(0,*)'ERROR! equal panel ranks in ProcessPair'
    write(0,*)'1:',1,p1%iblank,p1%area_ref,rank1
    write(0,*)'2:',2,p2%iblank,p2%area_ref,rank2
    stop
  end if


! if the subject polygon already has zero weight, do not bother
  if ((ipans%ratio == 0.0_rd).or.(ipanc%ratio == 0.0_rd)) then
    return
  end if


! generate clip polygon
! if (CalculatePolygonArea(ipanc%poly,ipanc) < ipanc%ratio) then

!   project ipanc%poly into the ipans polygon plane
!   clippoly%num_cont = ipanc%poly%num_cont

!   allocate(clippoly%cont(clippoly%num_cont))
!   nalloc = nalloc + 1
!   allocate(clippoly%hole(clippoly%num_cont))
!   nalloc = nalloc + 1

!   do ic = 1,clippoly%num_cont
!     clippoly%hole(ic)          = ipanc%poly%hole(ic)
!     clippoly%cont(ic)%num_vert = ipanc%poly%cont(ic)%num_vert

!     allocate(clippoly%cont(ic)%x(clippoly%cont(ic)%num_vert))
!     nalloc = nalloc + 1
!     allocate(clippoly%cont(ic)%y(clippoly%cont(ic)%num_vert))
!     nalloc = nalloc + 1
!     allocate(clippoly%cont(ic)%node(clippoly%cont(ic)%num_vert))
!     nalloc = nalloc + 1

!     do iv = 1,clippoly%cont(ic)%num_vert
!       do idir = 1,3
!         x(idir,1) = ipanc%g(1,idir)*ipanc%poly%cont(ic)%x(iv) &
!                   + ipanc%g(2,idir)*ipanc%poly%cont(ic)%y(iv) &
!                   + ipanc%xmid(idir)
!       end do
!       clippoly%cont(ic)%x(iv) = ipans%g(1,1)*(x(1,1)-ipans%xmid(1)) &
!                               + ipans%g(1,2)*(x(2,1)-ipans%xmid(2)) &
!                               + ipans%g(1,3)*(x(3,1)-ipans%xmid(3))
!       clippoly%cont(ic)%y(iv) = ipans%g(2,1)*(x(1,1)-ipans%xmid(1)) &
!                               + ipans%g(2,2)*(x(2,1)-ipans%xmid(2)) &
!                               + ipans%g(2,3)*(x(3,1)-ipans%xmid(3))
!       clippoly%cont(ic)%node(iv) = ipanc%poly%cont(ic)%node(iv)
!     end do
!   end do
! else
    do iv = 1,ipanc%num_vertices
    do idir = 1,3
      x(idir,iv) = xv(idir,ipanc%node(iv))
    end do
    end do
    call ConvertPanelToPolygon(ipans%g,x,ipans%xmid,clippoly, &
                               ipanc%num_vertices,ipanc%node, &
                               ipanc%edge)
! end if


! calculate the polygon boolean difference
  if (ipanc%poly%num_cont > 0) then
    subjpoly => ipans%poly
  
!   if (verbosity >= 3) then
    if (trap_clip) then
      open(unit=90,file=trim(pathname)//"scratch_clipping.dat", &
           action="write")
      write(90,*)"variables = x,y,node,edge"
      do ic = 1,subjpoly%num_cont
        write(90,*)'zone T = "Subject ',ipans%panel_number,ic,'"'
        do iv = 1,subjpoly%cont(ic)%num_vert
          write(90,*)subjpoly%cont(ic)%x(iv), &
                     subjpoly%cont(ic)%y(iv), &
                     subjpoly%cont(ic)%node(iv), &
                     subjpoly%cont(ic)%edge(iv)
        end do
      end do
      do ic = 1,clippoly%num_cont
        write(90,*)'zone T = "Clipper ',ipanc%panel_number,ic,'"'
        do iv = 1,clippoly%cont(ic)%num_vert
          write(90,*)clippoly%cont(ic)%x(iv), &
                     clippoly%cont(ic)%y(iv), &
                     clippoly%cont(ic)%node(iv), &
                     clippoly%cont(ic)%edge(iv)
        end do
      end do
      close(90)
    end if

    call my_cpu_time(time3)
    if (run_gpc) call gpc_f_translator(subjpoly,clippoly,1, &
                                       ipans%panel_number,  &
                                       ipanc%panel_number)
    call my_cpu_time(time4)

!   if (ipans%panel_number == plotpanel) &
!     call PrintPoly(subjpoly,6,"before Garbage")

!   The subroutine call was added to attempt garbage removal
!   on the resulting polygon out of GPC, i.e. to remove "degenerate"
!   contours and merge points in close proximity. 
!   It has proved necessary in curing a segmentation faults in
!   several cases, and examination of the polygons during operation
!   shows many degenerate cases and repeated vertices.  The other question,
!   of course, is to what extent it impacts execution time and whether it
!   can be done more efficiently.
    call GarbageCollector(subjpoly,ipans%panel_number,1)

    ran_gpc = .true.
    time43  = time4 - time3

    if (ipans%panel_number == plotpanel) &
      call PrintPoly(subjpoly,6,"after Garbage")

  end if

! clean-up clip polygon
! note: this is important to avoid a memory leak
  q => clippoly
  call DeallocatePolygon(q)

  nullify(ipans,ipanc,subjpoly)
  
  return
  end subroutine ProcessPair

! =================================================================
  subroutine gpc_f_translator (spoly, cpoly, operation, span, cpan)
! =================================================================
! Procedure to translate data between Fortran and c.
! Note that the dummy argument spoly is a pointer, so
! an explicit interface is required for this procedure.
! containing the procedure within a module satisfies that
! requirement.

  use InsertNodeM       ,only: InsertNode
  use IntrType          ,only: rd,rs
  use my_cpu_timeM      ,only: my_cpu_time
  use NetAlloc          ,only: nalloc
  use PrintPolyM        ,only: PrintPoly
  use TranslationArrays ,only: intc_array,fltc_array,asize
  use Types             ,only: polygon
  use UserInput         ,only: gpcf_time,never_skip_clip,pathname
  use UserInput         ,only: plotpanel,watertight

  implicit none

  type(polygon),pointer    :: spoly,q
  type(polygon),intent(in) :: cpoly
  type(polygon),target     :: opoly
! type(polygon),pointer    :: cpoly

! real(kind=rd)            :: dist2,match_tol,match_tol2
  real(kind=rd)            :: x1,x2,y1,y2,xmid,ymid,dist,distmin
  real(kind=rd)            :: t0_num,t0_den,t0,t0sav
  real(kind=rd)            :: rn,theta,delta,twopi
  real(kind=rs)            :: time1,time2

  integer,intent(in)       :: operation,span,cpan
  integer                  :: i,ic,ierr,iv,j,jc,jv
  integer                  :: nattempts,nextv
  integer                  :: number_floats,number_integers
  integer                  :: n1,n2,nsav
  integer,save             :: slide_number = 0

  character(len=80)        :: slide_file,fileid

  EXTERNAL gpc_c_translator_clip
  EXTERNAL gpc_c_translator_strip
  EXTERNAL toto
  
  continue


! initialization
  delta      = 1.0d-08
! match_tol  = 0.0_rd
! match_tol2 = match_tol**2
  twopi      = 2.0_rd*acos(-1.0_rd)


! check the amount of integer data required to be passed to gpc
  number_integers = 3 + 2 * spoly%num_cont + 2 * cpoly%num_cont
  if (number_integers > asize) then
    write(0,*)"ERROR! number_integers > asize: increase asize to ",asize
    write(0,*)"while working on panel ",span
    stop
  end if

! copy the subject polygon for future reference
  opoly%num_cont = spoly%num_cont
  allocate(opoly%cont(opoly%num_cont))
  nalloc = nalloc + 1
  allocate(opoly%hole(opoly%num_cont))
  nalloc = nalloc + 1

  do ic = 1,opoly%num_cont
    opoly%hole(ic)          = spoly%hole(ic)
    opoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert

    allocate(opoly%cont(ic)%x(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1
    allocate(opoly%cont(ic)%y(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1
    allocate(opoly%cont(ic)%node(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1
    allocate(opoly%cont(ic)%ie1(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1
    allocate(opoly%cont(ic)%ie2(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1
    allocate(opoly%cont(ic)%edge(opoly%cont(ic)%num_vert))
    nalloc = nalloc + 1

    do iv = 1,opoly%cont(ic)%num_vert
      opoly%cont(ic)%x(iv)    = spoly%cont(ic)%x(iv)
      opoly%cont(ic)%y(iv)    = spoly%cont(ic)%y(iv)
      opoly%cont(ic)%node(iv) = spoly%cont(ic)%node(iv)
      opoly%cont(ic)%ie1(iv)  = spoly%cont(ic)%ie1(iv)
      opoly%cont(ic)%ie2(iv)  = spoly%cont(ic)%ie2(iv)
      opoly%cont(ic)%edge(iv) = spoly%cont(ic)%edge(iv)
    end do
  end do


  CLIP_LOOP: do nattempts = 1,5

    if ((span == plotpanel).and.(operation==1)) then
      slide_number = slide_number + 1
      write(fileid,*)slide_number
      if (slide_number < 10) then
        slide_file = "slide0"//trim(adjustl(fileid))//".dat"
      else
        slide_file = "slide"//trim(adjustl(fileid))//".dat"
      end if
      open(unit=98,file=trim(pathname)//trim(slide_file),action="write")
      write(98,*)"variables = x,y,node,edge"
      write(fileid,*)span
      call PrintPoly(spoly,98,"subject "//trim(adjustl(fileid)))
      write(fileid,*)cpan
      call PrintPoly(cpoly,98,"clip "//trim(adjustl(fileid)))
!      call flush(98)
    end if


!   fill an array of integers for the two polygons
    i = 1
    intc_array(i) = asize
    i = i + 1
    intc_array(i) = spoly%num_cont
    i = i + 1
    do ic = 1,spoly%num_cont
      intc_array(i)   = spoly%cont(ic)%num_vert
      intc_array(i+1) = spoly%hole(ic)
      i = i + 2
    end do
    intc_array(i) = cpoly%num_cont
    i = i + 1
    do ic = 1,cpoly%num_cont
      intc_array(i)   = cpoly%cont(ic)%num_vert
      intc_array(i+1) = cpoly%hole(ic)
      i = i + 2
    end do


!   check the amount of real data required to be passed to gpc
    number_floats = 0
    do ic = 1,spoly%num_cont
      number_floats = number_floats + 2*spoly%cont(ic)%num_vert
    end do
    do ic = 1,cpoly%num_cont
      number_floats = number_floats + 2*cpoly%cont(ic)%num_vert
    end do
    if (number_floats > asize) then
      write(0,*)"ERROR! number_floats > asize: increase asize to ",asize
      write(0,*)"while working on panel ",span
      stop
    end if


!   fill in the real array
    i = 1
    do ic = 1,spoly%num_cont
    do iv = 1,spoly%cont(ic)%num_vert
      fltc_array(i+0) = spoly%cont(ic)%x(iv)
      fltc_array(i+1) = spoly%cont(ic)%y(iv)
      i = i + 2
    end do
    end do
    do ic = 1,cpoly%num_cont
    do iv = 1,cpoly%cont(ic)%num_vert
      if (nattempts == 1) then
        fltc_array(i+0) = cpoly%cont(ic)%x(iv)
        fltc_array(i+1) = cpoly%cont(ic)%y(iv)
      else
        call random_number(rn)
        theta = rn*twopi
        fltc_array(i+0) = cpoly%cont(ic)%x(iv) + delta*cos(theta)
        fltc_array(i+1) = cpoly%cont(ic)%y(iv) + delta*sin(theta)
      end if
      i = i + 2
    end do
    end do


!   call the c translator
    call my_cpu_time(time1)

    if (operation == 1) then
       call gpc_c_translator_clip(intc_array, fltc_array)
    else
       call gpc_c_translator_strip(intc_array, fltc_array)
    end if

    call my_cpu_time(time2)
    gpcf_time(1) = gpcf_time(1) + (time2-time1)

    ierr = intc_array(1)

    if (ierr < -2) then
      write(0,*)"ERROR! gpc_c_translator needs asize increased to ",-ierr
      stop
    end if

    if ((ierr == -1).or.(ierr == -2)) then
      write(0,"(1x,a,i7)")"WARNING! problem in gpc_polygon_clip "// &
                          "for panel ",span
      if (ierr == -1) then
        write(0,*)"clipping is being skipped to avoid a fatal error in gpc."
      else if (ierr == -2) then
        write(0,*)"(invalid value of hole)"
        write(0,*)"clipping is being skipped to avoid later problems."
      end if
      q => opoly
      call DeallocatePolygon(q)
    end if


    if (.not.never_skip_clip) then
      if (ierr < 0) then
        write(0,*)
        return
      else
        exit CLIP_LOOP
      end if
    else ! never_skip_clip
      if (ierr < 0) then
        if (nattempts < 5) then
          write(0,"(1x,a,i1)")"try again, attempt ",nattempts+1
        else
          write(0,*)"WARNING! giving up after 5 attempts"
          write(0,*)
          return
        end if
      else
        if (nattempts > 1) then
          write(0,*)"succeeded"
          write(0,*)
        end if
        exit CLIP_LOOP
      end if
    end if

  end do CLIP_LOOP


! erase the previous subject polygon memory
  call my_cpu_time(time1)
  call DeallocatePolygon(spoly)
  call my_cpu_time(time2)
  gpcf_time(2) = gpcf_time(2) + (time2-time1)


! reassign the data associated with the subject polygon
  call my_cpu_time(time1)
  i = 2
  j = 1
  spoly%num_cont = intc_array(i)
  i = i + 1

  allocate(spoly%cont(spoly%num_cont),stat=ierr)
  nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate spoly%cont"
    stop
  end if

  allocate(spoly%hole(spoly%num_cont),stat=ierr)
  nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate spoly%hole"
    stop
  end if

  do ic = 1,spoly%num_cont
    spoly%cont(ic)%num_vert = intc_array(i)
    spoly%hole(ic)          = intc_array(i+1)
    i = i + 2

    allocate(spoly%cont(ic)%x(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%x"
      stop
    end if

    allocate(spoly%cont(ic)%y(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%y"
      stop
    end if

    allocate(spoly%cont(ic)%node(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%node"
      stop
    end if

    allocate(spoly%cont(ic)%ie1(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%ie1"
      stop
    end if

    allocate(spoly%cont(ic)%ie2(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%ie2"
      stop
    end if

    allocate(spoly%cont(ic)%edge(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%edge"
      stop
    end if

    do iv = 1,spoly%cont(ic)%num_vert
      spoly%cont(ic)%x(iv) = fltc_array(j)
      spoly%cont(ic)%y(iv) = fltc_array(j+1)
      j = j + 2
!     initialize the node and edge numbers to zero for now
      spoly%cont(ic)%node(iv) = 0
      spoly%cont(ic)%ie1(iv)  = 0
      spoly%cont(ic)%ie2(iv)  = 0
      spoly%cont(ic)%edge(iv) = 0
    end do
  end do
  call my_cpu_time(time2)
  gpcf_time(3) = gpcf_time(3) + (time2-time1)


! double-check the amount of integer data being passed back from gpc
  number_integers = 2 + 2 * spoly%num_cont
  if (number_integers > asize) then
    write(0,*)"ERROR! number_integers > asize: increase asize to ",asize
    write(0,*)"while working on panel ",span
    stop
  end if


! double-check the amount of float data being passed back from gpc
  number_floats = 0
  do ic = 1,spoly%num_cont
    number_floats = number_floats + 2*spoly%cont(ic)%num_vert
  end do
  if (number_floats > asize) then
    write(0,*)"ERROR! number_floats > asize: increase asize to ",asize
    write(0,*)"while working on panel ",span
    stop
  end if


! The clipping is done and we have a new subject polygon.  While we
! still have the original panels, check and see whether any
! of the subject vertices were "inherited" from the original polygons.
! Keep match_tol very small (or zero) so that we do not mistakenly
! inherit any nodes that are merely nearby.

  do ic = 1,spoly%num_cont
  VERTEX_LOOP: do iv = 1,spoly%cont(ic)%num_vert

    if (operation == 1) then

!     compare to the clipping polygon
      do jc = 1,cpoly%num_cont
      do jv = 1,cpoly%cont(jc)%num_vert

!       allow some margin of error:
!       dist2 = (spoly%cont(ic)%x(iv) - cpoly%cont(jc)%x(jv))**2 &
!             + (spoly%cont(ic)%y(iv) - cpoly%cont(jc)%y(jv))**2
!       if (dist2 <= match_tol2) then
!       allow no margin of error:
        if ((spoly%cont(ic)%x(iv) == cpoly%cont(jc)%x(jv)).and. &
            (spoly%cont(ic)%y(iv) == cpoly%cont(jc)%y(jv))) then
          spoly%cont(ic)%node(iv) = cpoly%cont(jc)%node(jv)
          spoly%cont(ic)%ie1(iv)  = cpoly%cont(jc)%ie1(jv)
          spoly%cont(ic)%ie2(iv)  = cpoly%cont(jc)%ie2(jv)
          cycle VERTEX_LOOP
        end if

      end do
      end do
    end if


    do jc = 1,opoly%num_cont
      do jv = 1,opoly%cont(jc)%num_vert

!       compare each vertex in the new subject polygon to 
!       each vertex in the original subject polygon to see 
!       whether there are any matches
!       allow some margin of error:
!       dist2 = (spoly%cont(ic)%x(iv) - opoly%cont(jc)%x(jv))**2 &
!             + (spoly%cont(ic)%y(iv) - opoly%cont(jc)%y(jv))**2
!       if (dist2 <= match_tol2) then
!       allow no margin of error:
        if ((spoly%cont(ic)%x(iv) == opoly%cont(jc)%x(jv)).and. &
            (spoly%cont(ic)%y(iv) == opoly%cont(jc)%y(jv))) then
          spoly%cont(ic)%node(iv) = opoly%cont(jc)%node(jv)
          spoly%cont(ic)%ie1(iv)  = opoly%cont(jc)%ie1(jv)
          spoly%cont(ic)%ie2(iv)  = opoly%cont(jc)%ie2(jv)
          cycle VERTEX_LOOP
        end if

      end do
    end do

  end do VERTEX_LOOP


  if (watertight) then
!   The GPC library removes "superfluous" vertices from horizontal edges,
!   but USURP needs these in order to create watertight surfaces.  Test
!   whether either the clip polygon or original subject polygon has any 
!   vertices embedded in any horizontal edges of the new subject polygon, 
!   and insert them if it does.
    iv = 1
    do while (iv <= spoly%cont(ic)%num_vert)
      if (iv == spoly%cont(ic)%num_vert) then
        nextv = 1
      else
        nextv = iv + 1
      end if
      if (spoly%cont(ic)%y(iv) == spoly%cont(ic)%y(nextv)) then
!       horizontal edge
!       test the clip polygon
        do jc = 1,cpoly%num_cont
        do jv = 1,cpoly%cont(jc)%num_vert
          if (cpoly%cont(jc)%y(jv) == spoly%cont(ic)%y(iv)) then
            if ((cpoly%cont(jc)%x(jv) > min(spoly%cont(ic)%x(iv), &
                                            spoly%cont(ic)%x(nextv))).and. &
                (cpoly%cont(jc)%x(jv) < max(spoly%cont(ic)%x(iv), &
                                            spoly%cont(ic)%x(nextv)))) then
!              insert a new (known) node into the subject polygon
               call InsertNode(spoly,ic,iv,span, &
                               cpoly%cont(jc)%x(jv), &
                               cpoly%cont(jc)%y(jv), &
                               cpoly%cont(jc)%node(jv), &
                                0)
!                              cpoly%cont(jc)%edge(jv))
             end if
          end if
        end do
        end do

!       test the original polygon
        do jc = 1,opoly%num_cont
        do jv = 1,opoly%cont(jc)%num_vert
          if (opoly%cont(jc)%y(jv) == spoly%cont(ic)%y(iv)) then
            if ((opoly%cont(jc)%x(jv) > min(spoly%cont(ic)%x(iv), &
                                            spoly%cont(ic)%x(nextv))).and. &
                (opoly%cont(jc)%x(jv) < max(spoly%cont(ic)%x(iv), &
                                            spoly%cont(ic)%x(nextv)))) then
!              insert a new (known) node into the subject polygon
               call InsertNode(spoly,ic,iv,span, &
                               opoly%cont(jc)%x(jv), &
                               opoly%cont(jc)%y(jv), &
                               opoly%cont(jc)%node(jv), &
!                              opoly%cont(jc)%edge(jv))
                                0)
             end if
          end if
        end do
        end do
  
      end if
      iv = iv + 1
    end do
  end if !watertight


  VERTEX_LOOP2: do iv = 1,spoly%cont(ic)%num_vert

    if (operation == 1) then

      if (spoly%cont(ic)%node(iv) == 0) then

!       Assume any remaining free (hanging) nodes were created as a result
!       of the clipping panel breaking an edge.

        spoly%cont(ic)%node(iv) = -cpan
        spoly%cont(ic)%ie1(iv)  = 0
        spoly%cont(ic)%ie2(iv)  = 0

!       If a watertight surface is desired, do not allow new nodes to be
!       created too near existing nodes of the clipping polygon, because
!       these are often the result of an edge being projected onto two
!       adjacent clipping panels (non-unique projection).

        if (watertight) then

          xmid = spoly%cont(ic)%x(iv)
          ymid = spoly%cont(ic)%y(iv)
          distmin = 1.0e+20

          do jc = 1,cpoly%num_cont
            do jv = 1,cpoly%cont(jc)%num_vert
              x1 = cpoly%cont(jc)%x(jv)
              y1 = cpoly%cont(jc)%y(jv)
              n1 = cpoly%cont(jc)%node(jv)
              if (jv < cpoly%cont(jc)%num_vert) then
                x2 = cpoly%cont(jc)%x(jv+1)
                y2 = cpoly%cont(jc)%y(jv+1)
                n2 = cpoly%cont(jc)%node(jv+1)
              else
                x2 = cpoly%cont(jc)%x(1)
                y2 = cpoly%cont(jc)%y(1)
                n2 = cpoly%cont(jc)%node(1)
              end if
              t0_num = (xmid-x1)*(x2-x1) + (ymid-y1)*(y2-y1)
              t0_den = (  x2-x1)*(x2-x1) + (  y2-y1)*(y2-y1)
              if (t0_den /= 0.0_rd) t0 = t0_num / t0_den
              t0 = min(max(0.0_rd,t0),1.0_rd)
              x2 = x1 + t0*(x2-x1)
              y2 = y1 + t0*(y2-y1)
              dist = (xmid-x2)**2 + (ymid-y2)**2
              if (dist < distmin) then
                distmin = dist
                if (t0 < 0.5_rd) then
                  t0sav = t0
                  nsav = n1
                else
                  t0sav = 1.0_rd - t0
                  nsav = n2
                end if
              end if
            end do
          end do

!         Define what we mean by "too close" as a ratio
!         of the edge length (e.g. 0.001 = 0.1%)

          if (t0sav < 0.001_rd) then
            if (iv > 1) then
              n1 = spoly%cont(ic)%node(iv-1)
            else
              n1 = spoly%cont(ic)%node(spoly%cont(ic)%num_vert)
            end if
            if (iv < spoly%cont(ic)%num_vert) then
              n2 = spoly%cont(ic)%node(iv+1)
            else
              n2 = spoly%cont(ic)%node(1)
            end if

!           if (nsav == 2401) then
!             write(0,*)"iv,n1,n2,nsav = ",iv,n1,n2,nsav
!             write(0,*)"ic = ",ic
!             write(0,*)"num_vert = ",spoly%cont(ic)%num_vert
!             do jv = 1,spoly%cont(ic)%num_vert
!               write(0,*)jv,spoly%cont(ic)%node(jv)
!             end do
!           end if

            if ((n1 == nsav).or.(n2 == nsav)) then
!             leave the node named zero and deal with it later
              spoly%cont(ic)%node(iv) = 0
            else 
!             name the node according to the node of the clipping polygon
              spoly%cont(ic)%node(iv) = nsav
            end if
          end if

        end if !watertight
      end if !node(iv) = 0
    end if !clipping operation=1

  end do VERTEX_LOOP2


! For watertight surfaces, we may elect to delete nodes that were
! introduced too near to existing vertices on the clipping polygon 
  if ((operation == 1).and.(watertight)) then
    iv = 1
    do while (iv <= spoly%cont(ic)%num_vert)
      if (spoly%cont(ic)%node(iv) == 0) then
        do jv = iv,spoly%cont(ic)%num_vert - 1
          spoly%cont(ic)%node(jv) = spoly%cont(ic)%node(jv+1)
          spoly%cont(ic)%x(jv) = spoly%cont(ic)%x(jv+1)
          spoly%cont(ic)%y(jv) = spoly%cont(ic)%y(jv+1)
          spoly%cont(ic)%ie1(jv) = spoly%cont(ic)%ie1(jv+1)
          spoly%cont(ic)%ie2(jv) = spoly%cont(ic)%ie2(jv+1)
          spoly%cont(ic)%edge(jv) = spoly%cont(ic)%edge(jv+1)
        end do
        spoly%cont(ic)%num_vert = spoly%cont(ic)%num_vert - 1
      else
        iv = iv + 1
      end if
    end do
  end if

  end do



! While we still have the original panels, reconstruct the
! edge data from the original polygons.

  do ic = 1,spoly%num_cont
  EDGE_LOOP: do iv = 1,spoly%cont(ic)%num_vert

!   calculate the midpoint of each edge
    x1 = spoly%cont(ic)%x(iv)
    y1 = spoly%cont(ic)%y(iv)
    if (iv < spoly%cont(ic)%num_vert) then
      x2 = spoly%cont(ic)%x(iv+1)
      y2 = spoly%cont(ic)%y(iv+1)
    else
      x2 = spoly%cont(ic)%x(1)
      y2 = spoly%cont(ic)%y(1)
    end if
    xmid = 0.5_rd * (x1 + x2)
    ymid = 0.5_rd * (y1 + y2)
    distmin = 1.0e+20

    do jc = 1,opoly%num_cont
      do jv = 1,opoly%cont(jc)%num_vert
        x1 = opoly%cont(jc)%x(jv)
        y1 = opoly%cont(jc)%y(jv)
        if (jv < opoly%cont(jc)%num_vert) then
          x2 = opoly%cont(jc)%x(jv+1)
          y2 = opoly%cont(jc)%y(jv+1)
        else
          x2 = opoly%cont(jc)%x(1)
          y2 = opoly%cont(jc)%y(1)
        end if
        t0_num = (xmid-x1)*(x2-x1) + (ymid-y1)*(y2-y1)
        t0_den = (  x2-x1)*(x2-x1) + (  y2-y1)*(y2-y1)
        if (t0_den /= 0.0_rd) t0 = t0_num / t0_den
        t0 = min(max(0.0_rd,t0),1.0_rd)
        x2 = x1 + t0*(x2-x1)
        y2 = y1 + t0*(y2-y1)
        dist = (xmid-x2)**2 + (ymid-y2)**2
        if (dist < distmin) then
          distmin = dist
          spoly%cont(ic)%edge(iv) = opoly%cont(jc)%edge(jv)
        end if
      end do
    end do

    if (operation == 1) then

!     compare to the clipping polygon
      do jc = 1,cpoly%num_cont
      do jv = 1,cpoly%cont(jc)%num_vert
        x1 = cpoly%cont(jc)%x(jv)
        y1 = cpoly%cont(jc)%y(jv)
        if (jv < cpoly%cont(jc)%num_vert) then
          x2 = cpoly%cont(jc)%x(jv+1)
          y2 = cpoly%cont(jc)%y(jv+1)
        else
          x2 = cpoly%cont(jc)%x(1)
          y2 = cpoly%cont(jc)%y(1)
        end if
        t0_num = (xmid-x1)*(x2-x1) + (ymid-y1)*(y2-y1)
        t0_den = (  x2-x1)*(x2-x1) + (  y2-y1)*(y2-y1)
        if (t0_den /= 0.0_rd) t0 = t0_num / t0_den
        t0 = min(max(0.0_rd,t0),1.0_rd)
        x2 = x1 + t0*(x2-x1)
        y2 = y1 + t0*(y2-y1)
        dist = (xmid-x2)**2 + (ymid-y2)**2
        if (dist < distmin) then
          distmin = dist
          spoly%cont(ic)%edge(iv) = cpoly%cont(jc)%edge(jv)
        end if
      end do
      end do

    end if

  end do EDGE_LOOP
  end do


! put in the edge info for any newly formed intersection points
  if (watertight) then
    do ic = 1,spoly%num_cont
    do iv = 1,spoly%cont(ic)%num_vert
      if (spoly%cont(ic)%node(iv) < 0) then
        if (spoly%cont(ic)%ie2(iv) == 0) then
          spoly%cont(ic)%ie2(iv) = spoly%cont(ic)%edge(iv)
        end if
        if (spoly%cont(ic)%ie1(iv) == 0) then
          if (iv > 1) then
            jv = iv - 1
          else
            jv = spoly%cont(ic)%num_vert
          end if
          spoly%cont(ic)%ie1(iv) = spoly%cont(ic)%edge(jv)
        end if
      end if
    end do
    end do
  end if


  if ((span == plotpanel).and.(operation == 1)) then
    call PrintPoly(spoly,98,"result")
    close(98)
  end if


! clean up copy of subject polygon
  q => opoly
  call DeallocatePolygon(q)

  return
  end subroutine gpc_f_translator

! ===========================================================
  subroutine triangle_translator_f (spoly, pann)
! ==========================================================
! Procedure to translate data between Fortran and c.
! Note that the dummy argument spoly is a pointer, so
! an explicit interface is required for this procedure.

  use IntrType         ,only: rs
  use my_cpu_timeM     ,only: my_cpu_time
  use NetAlloc         ,only: nalloc
  use TranslationArrays,only: intc_array, fltc_array, asize
  use Types            ,only: polygon
  use UserInput        ,only: gpcf_time,ltri,pathname,verbosity
  use UserInput        ,only: trap_tri

  implicit none

  type(polygon),pointer :: spoly

  real(kind=rs)         :: time1,time2

  integer,allocatable,dimension(:) :: saved_node
  integer,intent(in)    :: pann
  integer               :: current_point,first_point
  integer               :: i,ic,ierr,iv
  integer               :: number_floats,number_integers
  integer               :: number_nodes,number_segments,number_holes
  integer               :: num_corners


  continue


  if (trap_tri) then
    open(unit=90,file=trim(pathname)//"scratch_triangle.dat",action="write")
    write(90,*)"variables = x,y,node,edge"
  end if


! Count the number of vertices in the polygon
  number_nodes = 0
  do ic = 1,spoly%num_cont
!   if (pann == 123846) then
!     write(*,*)ic,ltri(ic),spoly%cont(ic)%num_vert
!   end if
    if (ltri(ic) == 0) then
      number_nodes = number_nodes + spoly%cont(ic)%num_vert
    end if
  end do
! if (pann == 123846) then
!   write(*,*)"number_nodes = ",number_nodes
! end if

  allocate(saved_node(number_nodes))
  nalloc = nalloc + 1

! CAUTION: hole logic not yet coded
  number_holes = 0


! Each contour has number of segments equal to number of vertices
! So the total number of segments is equal to the total number of
! vertices.
  number_segments = number_nodes


! Check the amount of integer data required to be passed to Triangle.
! For Triangle, integer data consists of the following:
! (1) dimension of integer array
! (2) panel number
! (3) total number of nodes
! (4) total number of segments
! (5) total number of holes
! (6) 2 times the number of segments (segment endpoint data)
  number_integers = 5 + 2*number_segments


  if (number_integers > asize) then
    write(0,*)"ERROR! number_integers > asize: increase asize to ",asize
    write(0,*)"while working on panel ",pann
    write(0,*)"subroutine triangle_translator_f"
    stop
  end if


! fill an array of integers for the polygon
  i = 1
  intc_array(i) = asize
  i = i + 1
  intc_array(i) = pann
  i = i + 1
  intc_array(i) = number_nodes
  i = i + 1
  intc_array(i) = number_segments
  i = i + 1
  intc_array(i) = number_holes
  i = i + 1
  first_point = 0
  do ic = 1,spoly%num_cont
    if (ltri(ic) == 0) then
      current_point = first_point
      do iv = 1,spoly%cont(ic)%num_vert
        intc_array(i)   = current_point
        if (iv < spoly%cont(ic)%num_vert) then
          intc_array(i+1) = current_point + 1
        else
          intc_array(i+1) = first_point
        end if
        i = i + 2
        current_point = current_point + 1
      end do
      first_point = first_point + spoly%cont(ic)%num_vert
    end if
  end do
! if (pann == 123846) then
!   write(*,*)intc_array(1:number_integers)
! end if


! Check the amount of real data required to be passed to Triangle.
! For Triangle, integer data consists of the following:
! (1) x,y coordinates of each node (2*number_nodes)
! (2) x,y point in each hole (2*number_holes)
! CAUTION: hole logic not yet coded
  number_floats = 2*(number_nodes + number_holes)


  if (number_floats > asize) then
    write(0,*)"ERROR! number_floats > asize: increase asize to ",asize
    write(0,*)"while working on panel ",pann
    stop
  end if


! fill in the real array
  i = 1
  do ic = 1,spoly%num_cont
    if (ltri(ic) == 0) then
!     if (verbosity >= 3) write(90,*)'zone T = "',pann,ic,'"'
      if (trap_tri) write(90,*)'zone T = "',pann,ic,'"'
      do iv = 1,spoly%cont(ic)%num_vert
        fltc_array(i  ) = spoly%cont(ic)%x(iv)
        fltc_array(i+1) = spoly%cont(ic)%y(iv)
!       if (verbosity >= 3) then
        if (trap_tri) then
          write(90,*)spoly%cont(ic)%x(iv), &
                     spoly%cont(ic)%y(iv), &
                     spoly%cont(ic)%node(iv), &
                     spoly%cont(ic)%edge(iv)
        end if
        saved_node((i+1)/2) = spoly%cont(ic)%node(iv)
        i = i + 2
      end do
!     if (verbosity >= 3) then
      if (trap_tri) then
        write(90,*)spoly%cont(ic)%x(1), &
                   spoly%cont(ic)%y(1), &
                   spoly%cont(ic)%node(1), &
                   spoly%cont(ic)%edge(1)
      end if
    end if
  end do
  if (trap_tri) close(90)


! call the c translator
  call my_cpu_time(time1)

  call triangle_translator_c(intc_array, fltc_array)

  call my_cpu_time(time2)
  gpcf_time(1) = gpcf_time(1) + (time2-time1)

  ierr = intc_array(1)
  if (ierr < 0) then
    write(0,*)"ERROR! triangle_translator_c has returned with an error"
    write(0,*)"and may need asize to be increased to ",-ierr
    stop
  end if

  num_corners = intc_array(3)
  if (num_corners /= 3) then
    write(0,*)"ERROR! triangle_translator_c has returned with "
    write(0,*)"num_corners not equal to 3: num_corners = ",num_corners
    write(0,*)"Check options passed to Triangle."
    stop
  end if


! erase the previous subject polygon memory
  call my_cpu_time(time1)
  call DeallocatePolygon(spoly)
  call my_cpu_time(time2)
  gpcf_time(2) = gpcf_time(2) + (time2-time1)


! reassign the data associated with the subject polygon
  call my_cpu_time(time1)
  spoly%num_cont = intc_array(2)

  allocate(spoly%cont(spoly%num_cont),stat=ierr)
  nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate spoly%cont"
    stop
  end if

  allocate(spoly%hole(spoly%num_cont),stat=ierr)
  nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! failed to allocate spoly%hole"
    stop
  end if

  i = 4

  if (maxval(intc_array(i:i+intc_array(2)*intc_array(3)-1)) &
       > number_nodes - 1) then
    write(0,*)"WARNING!  Triangle has failed on panel ",pann
  end if

  do ic = 1,spoly%num_cont
    spoly%cont(ic)%num_vert = num_corners
    spoly%hole(ic)          = 0

    allocate(spoly%cont(ic)%x(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%x"
      stop
    end if

    allocate(spoly%cont(ic)%y(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%y"
      stop
    end if

    allocate(spoly%cont(ic)%node(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%node"
      stop
    end if

    allocate(spoly%cont(ic)%ie1(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%ie1"
      stop
    end if

    allocate(spoly%cont(ic)%ie2(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%ie2"
      stop
    end if

    allocate(spoly%cont(ic)%edge(spoly%cont(ic)%num_vert),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! failed to allocate spoly%cont%edge"
      stop
    end if

    do iv = 1,spoly%cont(ic)%num_vert
      if (intc_array(i)+1 > number_nodes) then
!       write(0,*)"WARNING!  Triangle has failed on panel ",pann
!       Allow the code to proceed by setting the illegitimate node
!       number back to the first node of the polygon.........?????
        intc_array(i) = 0
      end if
      spoly%cont(ic)%x(iv)    = fltc_array(2*intc_array(i)+1)
      spoly%cont(ic)%y(iv)    = fltc_array(2*intc_array(i)+2)
      spoly%cont(ic)%node(iv) = saved_node(intc_array(i)+1)
      spoly%cont(ic)%ie1(iv) = 0
      spoly%cont(ic)%ie2(iv) = 0
      spoly%cont(ic)%edge(iv) = 0
      i = i + 1
    end do
  end do
  call my_cpu_time(time2)
  gpcf_time(3) = gpcf_time(3) + (time2-time1)


! double-check the amount of integer data being passed back from gpc
  number_integers = 3 + spoly%num_cont * num_corners
  if (number_integers > asize) then
    write(0,*)"ERROR! number_integers > asize: increase asize to ",asize
    write(0,*)"while working on panel ",pann
    stop
  end if


! double-check the amount of float data being passed back from gpc
! (not needed: output is identical to input)

  deallocate(saved_node)
  nalloc = nalloc - 1

  return
  end subroutine triangle_translator_f

! ===========================================
  subroutine DeallocatePolygon(spoly)
! ===========================================
! erase the previous subject polygon memory

  use NetAlloc,only: nalloc
  use Types   ,only: polygon

  implicit none

  type(polygon),pointer :: spoly

  integer :: ic


  continue


  do ic = 1,spoly%num_cont
    if (associated(spoly%cont(ic)%x)) then
      deallocate(spoly%cont(ic)%x)
      nalloc = nalloc - 1
    end if
    if (associated(spoly%cont(ic)%y)) then
      deallocate(spoly%cont(ic)%y)
      nalloc = nalloc - 1
    end if
    if (associated(spoly%cont(ic)%node)) then
      deallocate(spoly%cont(ic)%node)
      nalloc = nalloc - 1
    end if
    if (associated(spoly%cont(ic)%ie1)) then
      deallocate(spoly%cont(ic)%ie1)
      nalloc = nalloc - 1
    end if
    if (associated(spoly%cont(ic)%ie2)) then
      deallocate(spoly%cont(ic)%ie2)
      nalloc = nalloc - 1
    end if
    if (associated(spoly%cont(ic)%edge)) then
      deallocate(spoly%cont(ic)%edge)
      nalloc = nalloc - 1
    end if
  end do
  if (associated(spoly%cont)) then
    deallocate(spoly%cont)
    nalloc = nalloc - 1
  end if
  if (associated(spoly%hole)) then
    deallocate(spoly%hole)
    nalloc = nalloc - 1
  end if
  spoly%num_cont = 0


  return
  end subroutine DeallocatePolygon

! ======================================================================
  function SetRank(iblank) result (rank)
! ======================================================================
! Set the precedence of panels based on their i-blank value

  implicit none

  integer,intent(in) :: iblank
  integer            :: rank


  continue


  if (iblank == 1) then
    rank = 4
  else if (iblank < 0) then
    rank = 3
  else if (iblank == 101) then
    rank = 2
  else if (iblank == 0) then
    rank = 1
  else
    write(0,*)'ERROR! unrecognized iblank value in SetRank, iblank = ',iblank
    stop
  end if

  return
  end function SetRank

  end module ProcessPairProcedures
