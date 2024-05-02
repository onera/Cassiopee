! ================================================================
  module CreateTrianglesM
! ================================================================

  implicit none

  private
  public  :: CreateTriangles
  private :: SplitContours,SelfIntersection

  contains

! ======================================================================
  subroutine CreateTriangles(p)
! ======================================================================
! Convert designated polygons into triangles.
!
! Note that this procedure contains a pointer as the
! dummy argument p, so this procedure requires an explicit interface.

  use GarbageCollectorM    ,only: GarbageCollector
  use IntrType             ,only: rd,rs
  use my_cpu_timeM         ,only: my_cpu_time
  use NetAlloc             ,only: nalloc
  use PrintPolyM           ,only: PrintPoly
  use ProcessPairProcedures,only: gpc_f_translator,triangle_translator_f
  use ProcessPairProcedures,only: DeallocatePolygon
  use TriQuadM             ,only: TriQuad
  use Types                ,only: panel,polygon
  use UserInput            ,only: verbosity,pathname,plotpanel,showcpu
  use UserInput            ,only: full_surface,ttype,ltri,never_skip_tri
  use UserInput            ,only: cpufmt,trap_tri

  implicit none

  type(panel),pointer,dimension (:) :: p
  type(polygon),pointer             :: subjpoly
  type(polygon)                     :: clippoly

  real(kind=rd) :: reset_max
  real(kind=rd) :: triangulated_min
  real(kind=rd) :: max_weight_limit
  real(kind=rd) :: min_weight_limit
  real(kind=rd) :: x1,x2,y1,y2,areac

  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: gtime,stime,ttime,utime

  integer :: ic,ip1,ip2,ipan,iv1,iv2
  integer :: num_holes
  integer :: gone,smallest
  integer :: too_big,too_small,triangulated
  integer :: whole,npskipped
  integer :: num_degen_area
  integer :: ncskipped

  character(len=32) :: tag

  logical,parameter :: check_areas = .true.
  logical           :: unclipped,lop


  continue


! initialization
  call my_cpu_time(time1)
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)
  clippoly%num_cont = 0
  nullify(clippoly%hole)
  nullify(clippoly%cont)
  gtime = 0.0_rs
  stime = 0.0_rs
  ttime = 0.0_rs
  utime = 0.0_rs
  reset_max = 0.0_rd
  triangulated_min = 2.0_rd
  lop = .false.
  num_holes = 0
  gone = 0
  too_big = 0
  too_small = 0
  triangulated = 0
  whole = 0
  npskipped = 0
  num_degen_area = 0
  ncskipped = 0


! set the bounds of what we want to include in the
! partial triangulation for Tecplot
  max_weight_limit = 0.999_rd
  min_weight_limit = 0.001_rd

  do ipan = ip1,ip2

    p(ipan)%itri = 0
    subjpoly => p(ipan)%poly

    if (p(ipan)%ratio == 0.0_rd) then
!     WARNING!!!!!!!!!!!! The next line causes a memory leak
!     subjpoly%num_cont = 0
!     write(*,*)"ratio = 0 for ipan ",ipan,p(ipan)%iblank
!     call PrintPoly(subjpoly,0,"Delete this:")
      call DeallocatePolygon(subjpoly)
    end if

    if (p(ipan)%panel_number == plotpanel) &
      call PrintPoly(subjpoly,6,"Last garbage")
     
    call my_cpu_time(time3)
    call GarbageCollector(subjpoly,p(ipan)%panel_number,1)
!   if (never_skip_tri) then
!     call SplitContours(subjpoly,ipan)
      call SplitContours(subjpoly)
      call GarbageCollector(subjpoly,ipan,2)
!   end if
    call my_cpu_time(time4)
    gtime = gtime + (time4-time3)

    if (subjpoly%num_cont == 0) then
      reset_max = max(reset_max,p(ipan)%ratio)
      p(ipan)%ratio = 0.0_rd
    end if

    if (((full_surface).or. &
       ((p(ipan)%ratio > min_weight_limit).and.  &
        (p(ipan)%ratio < max_weight_limit))).and. &
      (subjpoly%num_cont > 0)) then

      if ((p(ipan)%ratio > min_weight_limit).and. &
          (p(ipan)%ratio < max_weight_limit)) then
!       include this in the partial Tecplot triangulation
        p(ipan)%itri = 1
      end if

      triangulated = triangulated + 1

      if (p(ipan)%ratio < triangulated_min) then
        triangulated_min = p(ipan)%ratio
        smallest = p(ipan)%panel_number
      end if

      do ic = 1,subjpoly%num_cont
        num_holes = num_holes + subjpoly%hole(ic)
      end do

      if (ipan == plotpanel) then
        open(unit=96,file=trim(pathname)//"triangulate.dat")
        write(96,*)"variables = x,y,node,edge"
        call PrintPoly(subjpoly,96,"triangulate")
        close(96)
      end if

      if ((p(ipan)%ratio == 1.0_rd).and. &
          (p(ipan)%poly%num_cont == 1).and. &
          (p(ipan)%poly%hole(1) == 0).and. &
          (p(ipan)%poly%cont(1)%num_vert == 4)) then
        unclipped = .true.
      else
        unclipped = .false.
      end if

!     perform triangulation
      if (ttype == "GPC") then

        if (unclipped) then
          call my_cpu_time(time3)
          call TriQuad(subjpoly)
          call my_cpu_time(time4)
          utime = utime + (time4-time3)
        else
          call my_cpu_time(time3)
          call gpc_f_translator(subjpoly,clippoly,2,ipan,0)
          call my_cpu_time(time4)
          ttime = ttime + (time4-time3)
        end if

      else if (ttype == "Triangle") then

        if (unclipped) then
          call my_cpu_time(time3)
          call TriQuad(subjpoly)
          call my_cpu_time(time4)
          utime = utime + (time4-time3)
        else

          allocate(ltri(subjpoly%num_cont)); nalloc = nalloc + 1
          ltri=0

          call my_cpu_time(time3)
          if (.not.never_skip_tri) then

            if (check_areas) then
!           throw out any contours that are "too small"
            do ic = 1,subjpoly%num_cont
              if (ltri(ic) > 0) cycle
              areac = 0.0_rd
              do iv1 = 1,subjpoly%cont(ic)%num_vert
                if (iv1 == subjpoly%cont(ic)%num_vert) then
                  iv2 = 1
                else
                  iv2 = iv1 + 1
                end if
                x2 = subjpoly%cont(ic)%x(iv2)
                x1 = subjpoly%cont(ic)%x(iv1)
                y2 = subjpoly%cont(ic)%y(iv2)
                y1 = subjpoly%cont(ic)%y(iv1)
                areac = areac + 0.5_rd*(x1-x2)*(y2+y1)
              end do
!             if (ipan == 123846) then
!               write(*,*)"area ",ic,areac
!             end if
              if (abs(areac) < 1.0d-10) then
                num_degen_area = num_degen_area + 1
                ltri(ic) = 1
              end if
            end do
            end if !check_areas

            call SelfIntersection(subjpoly,ipan)

          end if !.not.never_skip_tri
          call my_cpu_time(time4)
          stime = stime + (time4-time3)


          if (maxval(ltri) > 0) then
!           at least one contour in this polygon will be skipped
            do ic = 1,subjpoly%num_cont
              if (ltri(ic) > 0) ncskipped = ncskipped + 1
            end do
            npskipped = npskipped + 1
!           if (verbosity >= 1) then
!             write(0,*)"Skip triangulation, self-intersecting "// &
!                     "polygon, panel = ",ipan
!           end if
            if (verbosity >= 2) then
              write(tag,"(a6,i7)")"panel ",ipan
              if (.not.lop) then
                open(unit=91, &
                file=trim(pathname)//"skipped.dat",action="write")
                write(91,*)"variables = x,y,node,edge"
                lop = .true.
              end if
              call PrintPoly(subjpoly,91,trim(tag))
            end if
          end if

          if (minval(ltri) == 0) then
!           at least one contour in this polygon will be triangulated
            call my_cpu_time(time3)
            call triangle_translator_f(subjpoly,ipan)
            call my_cpu_time(time4)
            ttime = ttime + (time4-time3)
          end if

          deallocate(ltri); nalloc = nalloc - 1

        end if !clipped

      else !GPC/Triangle

        write(0,*)"ERROR! invalid value of ttype ",trim(ttype)
        stop
      end if

      if (ipan == plotpanel) then
        open(unit=96,file=trim(pathname)//"triangulated.dat")
        write(96,*)"variables = x,y,node,edge"
        call PrintPoly(subjpoly,96,"triangulated")
        close(96)
      end if

    else !do not triangulate this polygon

      if (p(ipan)%ratio == 1.0_rd) then
        whole = whole + 1
      else if ((p(ipan)%ratio == 0.0_rd).or. &
               (subjpoly%num_cont == 0)) then
        gone = gone + 1
      else if ((p(ipan)%ratio >=  max_weight_limit).and. &
               (p(ipan)%ratio < 1.0_rd)) then
        too_big = too_big + 1
      else if ((p(ipan)%ratio <= min_weight_limit).and. &
               (p(ipan)%ratio > 0.0_rd)) then
        too_small = too_small + 1
      end if

    end if
  
  end do

  if (verbosity >= 2) then
    write(0,*)
    write(0,*)"From CreateTriangles: "
    write(0,*)"total panels                             = ",ip2-ip1+1
    if (.not.full_surface) then
      write(0,*)"panels that remained unclipped           = ",whole
      write(0,*)"panels that were completely removed      = ",gone
      write(0,*)"panels too big to bother triangulating   = ",too_big
      write(0,*)"panels too small to bother triangulating = ",too_small
    end if
    write(0,*)"triangulated panels                      = ",triangulated
    write(0,*)"number of holes in triangulated panels   = ",num_holes
    write(0,*)"smallest panel that was triangulated     = ",smallest
    write(0,*)"smallest weight that was triangulated    = ",triangulated_min
    write(0,*)"largest weight that was reset to zero    = ",reset_max
  end if

  if ((num_holes > 0).and.(ttype == "Triangle")) then
    write(0,*)
    write(0,"(a,i5,a)")" WARNING! The triangulated polygons included ", &
                         num_holes," holes."
    write(0,*)"The Triangle library is not yet set-up to handle holes."
    write(0,*)"Consider using --ttype=GPC instead."
  end if

  if ((npskipped > 0).and.(verbosity < 2)) then
    write(0,*)
    if (ncskipped == 1) then
      write(0,"(a,i5,a)")" WARNING! USURP skipped the triangulation of ", &
      ncskipped," contour"
    else
      write(0,"(a,i5,a)")" WARNING! USURP skipped the triangulation of ", &
      ncskipped," contours"
    end if
    if (npskipped == 1) then
      write(0,"(a,i5,a)")" on ",npskipped, &
         " polygon in hopes of avoiding a crash by Triangle."
    else
      write(0,"(a,i5,a)")" on ",npskipped, &
         " polygons in hopes of avoiding a crash by Triangle."
    end if

    if (num_degen_area == 1) then
      write(0,"(a,i5,a)")" (",num_degen_area, &
          " of these contours was deemed too small.)"
    else if (num_degen_area > 1) then
      write(0,"(a,i5,a)")" (",num_degen_area, &
          " of these contours were deemed too small.)"
    end if

    write(0,*)"Use --never-skip=tri to force triangulation, or"
    write(0,"(a,i1,a)")" Increase to --verbose=",verbosity+1, &
                       " for more details."
    write(0,*)
  end if

  if (trap_tri) then
    open(unit=90,file=trim(pathname)//"scratch_triangle.dat")
    close(90,status="delete")
  end if

  inquire(file=trim(pathname)//"skipped.dat",opened=lop)
  if (lop) close(91)

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CreateTriangles: '  ,time2-time1
    write(0,cpufmt)'  > time in GarbageCollector: ' ,gtime
    write(0,cpufmt)'  > time in SelfIntersection: ' ,stime
    write(0,cpufmt)'  > time in '//trim(ttype)//': ',ttime
    write(0,cpufmt)'  > time in TriQuad: '          ,utime
  end if

  return
  end subroutine CreateTriangles

! ====================================================
  subroutine SplitContours(spoly)
! subroutine SplitContours(spoly,pann)
! ====================================================
! This procedure checks for self-intersecting
! polygon contours, which cause problems in Triangle,
! and attempts to split the contour at the 
! point of self-intersection.
! Note that the dummy argument spoly is a pointer
! so this procedure requires an explicit interface.

  use IntrType             ,only: rd
  use NetAlloc             ,only: nalloc
  use ProcessPairProcedures,only: DeallocatePolygon
  use Types                ,only: polygon

  implicit none

  type(polygon),pointer :: spoly,q
  type(polygon),target  :: opoly
  real(kind=rd)         :: xtmp,ytmp
! integer,intent(in)    :: pann
  integer               :: ic,jc,iv,jv,kv,kv2,nco
  integer               :: ntmp,ie1_tmp,ie2_tmp,edge_tmp


  continue


! check each contour in the polygon to see whether any self-
! intersections occur (note that contours may later need to be 
! split as we go, so we need the range of the do loop to be 
! variable).  In this routine, check only for two points with
! the same x,y location.

  ic = 1
  ICLOOP: do while (ic <= spoly%num_cont)

!   see if any points are revisited
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert)
      jv = iv + 1
      do while (jv <= spoly%cont(ic)%num_vert)
        if ((spoly%cont(ic)%x(iv) == spoly%cont(ic)%x(jv)).and. &
            (spoly%cont(ic)%y(iv) == spoly%cont(ic)%y(jv))) then

!         write(0,*)
!         write(0,*)"WARNING! repeated vertex in panel ",pann
!         write(0,*)"SplitContour indices: ",iv,jv

          do while (iv /= 1)
            xtmp = spoly%cont(ic)%x(1)
            ytmp = spoly%cont(ic)%y(1)
            ntmp = spoly%cont(ic)%node(1)
            ie1_tmp = spoly%cont(ic)%ie1(1)
            ie2_tmp = spoly%cont(ic)%ie2(1)
            edge_tmp = spoly%cont(ic)%edge(1)
            do kv = 1,spoly%cont(ic)%num_vert-1
              spoly%cont(ic)%x(kv)    = spoly%cont(ic)%x(kv+1)
              spoly%cont(ic)%y(kv)    = spoly%cont(ic)%y(kv+1)
              spoly%cont(ic)%node(kv) = spoly%cont(ic)%node(kv+1)
              spoly%cont(ic)%ie1(kv)  = spoly%cont(ic)%ie1(kv+1)
              spoly%cont(ic)%ie2(kv)  = spoly%cont(ic)%ie2(kv+1)
              spoly%cont(ic)%edge(kv) = spoly%cont(ic)%edge(kv+1)
            end do
            kv = spoly%cont(ic)%num_vert
            spoly%cont(ic)%x(kv) = xtmp
            spoly%cont(ic)%y(kv) = ytmp
            spoly%cont(ic)%node(kv) = ntmp
            spoly%cont(ic)%ie1(kv) = ie1_tmp
            spoly%cont(ic)%ie2(kv) = ie2_tmp
            spoly%cont(ic)%edge(kv) = edge_tmp
            iv = iv - 1
            jv = jv - 1
          end do

!         write(0,*)"Split contour; ipan,ic,iv = ",pann,ic,jv

!         split the contour while copying the entire polygon into opoly

          nco = spoly%num_cont + 1
          opoly%num_cont = nco

          allocate(opoly%cont(nco)); nalloc = nalloc + 1
          allocate(opoly%hole(nco)); nalloc = nalloc + 1

          do jc = 1,spoly%num_cont

            opoly%hole(jc) = spoly%hole(jc) 
            if (jc /= ic) then
              ntmp = spoly%cont(jc)%num_vert
            else
              ntmp = jv - 1
            end if
            opoly%cont(jc)%num_vert = ntmp

            allocate(opoly%cont(jc)%x(ntmp))   ; nalloc = nalloc + 1
            allocate(opoly%cont(jc)%y(ntmp))   ; nalloc = nalloc + 1
            allocate(opoly%cont(jc)%node(ntmp)); nalloc = nalloc + 1
            allocate(opoly%cont(jc)%ie1(ntmp)) ; nalloc = nalloc + 1
            allocate(opoly%cont(jc)%ie2(ntmp)) ; nalloc = nalloc + 1
            allocate(opoly%cont(jc)%edge(ntmp)); nalloc = nalloc + 1

            do kv = 1,ntmp
              opoly%cont(jc)%x(kv)    = spoly%cont(jc)%x(kv)
              opoly%cont(jc)%y(kv)    = spoly%cont(jc)%y(kv)
              opoly%cont(jc)%node(kv) = spoly%cont(jc)%node(kv)
              opoly%cont(jc)%ie1(kv)  = spoly%cont(jc)%ie1(kv)
              opoly%cont(jc)%ie2(kv)  = spoly%cont(jc)%ie2(kv)
              opoly%cont(jc)%edge(kv) = spoly%cont(jc)%edge(kv)
            end do

          end do

          ntmp = spoly%cont(ic)%num_vert - jv + 1
          opoly%cont(nco)%num_vert = ntmp
          opoly%hole(nco) = spoly%hole(ic)
          allocate(opoly%cont(nco)%x(ntmp))   ; nalloc = nalloc + 1
          allocate(opoly%cont(nco)%y(ntmp))   ; nalloc = nalloc + 1
          allocate(opoly%cont(nco)%node(ntmp)); nalloc = nalloc + 1
          allocate(opoly%cont(nco)%ie1(ntmp)) ; nalloc = nalloc + 1
          allocate(opoly%cont(nco)%ie2(ntmp)) ; nalloc = nalloc + 1
          allocate(opoly%cont(nco)%edge(ntmp)); nalloc = nalloc + 1
          kv2 = jv
          do kv = 1,ntmp
            opoly%cont(nco)%x(kv)    = spoly%cont(ic)%x(kv2)
            opoly%cont(nco)%y(kv)    = spoly%cont(ic)%y(kv2)
            opoly%cont(nco)%node(kv) = spoly%cont(ic)%node(kv2)
            opoly%cont(nco)%ie1(kv)  = spoly%cont(ic)%ie1(kv2)
            opoly%cont(nco)%ie2(kv)  = spoly%cont(ic)%ie2(kv2)
            opoly%cont(nco)%edge(kv) = spoly%cont(ic)%edge(kv2)
            kv2 = kv2 + 1
          end do

          call DeallocatePolygon(spoly)

!         copy opoly back into spoly
          spoly%num_cont = opoly%num_cont
          allocate(spoly%hole(spoly%num_cont)); nalloc = nalloc + 1
          allocate(spoly%cont(spoly%num_cont)); nalloc = nalloc + 1
          do jc = 1,opoly%num_cont
            spoly%hole(jc)          = opoly%hole(jc)
            spoly%cont(jc)%num_vert = opoly%cont(jc)%num_vert
            ntmp = spoly%cont(jc)%num_vert
            allocate(spoly%cont(jc)%x(ntmp))   ; nalloc = nalloc + 1
            allocate(spoly%cont(jc)%y(ntmp))   ; nalloc = nalloc + 1
            allocate(spoly%cont(jc)%node(ntmp)); nalloc = nalloc + 1
            allocate(spoly%cont(jc)%ie1(ntmp)) ; nalloc = nalloc + 1
            allocate(spoly%cont(jc)%ie2(ntmp)) ; nalloc = nalloc + 1
            allocate(spoly%cont(jc)%edge(ntmp)); nalloc = nalloc + 1
            do kv = 1,spoly%cont(jc)%num_vert
              spoly%cont(jc)%x(kv)    = opoly%cont(jc)%x(kv)
              spoly%cont(jc)%y(kv)    = opoly%cont(jc)%y(kv)
              spoly%cont(jc)%node(kv) = opoly%cont(jc)%node(kv)
              spoly%cont(jc)%ie1(kv)  = opoly%cont(jc)%ie1(kv)
              spoly%cont(jc)%ie2(kv)  = opoly%cont(jc)%ie2(kv)
              spoly%cont(jc)%edge(kv) = opoly%cont(jc)%edge(kv)
            end do
          end do

          q => opoly
          call DeallocatePolygon(q)

          cycle ICLOOP
        end if
        jv = jv + 1
      end do
      iv = iv + 1
    end do

    ic = ic + 1
  end do ICLOOP

  return
  end subroutine SplitContours

! ====================================================
  subroutine SelfIntersection(spoly,pann)
! ====================================================
! This procedure checks for self-intersecting
! polygon contours, which cause problems in Triangle.
! Note that the dummy argument spoly is a pointer
! so this procedure requires an explicit interface.

  use IntrType ,only: rd
  use Types    ,only: polygon
  use UserInput,only: verbosity,full_surface,ltri

  implicit none

  type(polygon),pointer :: spoly
  real(kind=rd)         :: ax,bx,cx,dx,ay,by,cy,dy,rr,ss
  real(kind=rd)         :: denom,dist
  integer,intent(in)    :: pann
  integer               :: ic,iv,jv
  logical,parameter     :: check_colinear = .false.


  continue


! check each contour in the polygon to see whether any self-
! intersections occur (note that contours may later need to be 
! split as we go, so we need the range of the do loop to be 
! variable)

  ic = 1
  ICLOOP: do while (ic <= spoly%num_cont)

    if (ltri(ic) > 0) then
      ic = ic + 1
      cycle ICLOOP
    end if

!   check for repeated node numbers only when using --full-surface or
!   --watertight (note that full_surface=.true. for watertight as well)
    if (full_surface) then


!   SEE IF ANY NODES ARE REVISITED
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert)
      jv = iv + 1
      do while (jv <= spoly%cont(ic)%num_vert)
        if (spoly%cont(ic)%node(iv) == spoly%cont(ic)%node(jv)) then
          if (verbosity >= 2) then
            write(0,*)
            write(0,*)"WARNING! repeated node in panel ",pann
            write(0,*)iv,jv,spoly%cont(ic)%node(iv)
            if ((spoly%cont(ic)%x(iv) /= spoly%cont(ic)%x(jv)).or. &
                (spoly%cont(ic)%y(iv) /= spoly%cont(ic)%y(jv))) then
              write(0,*)"x,y values are not equal!!!"
              write(0,*)spoly%cont(ic)%x(iv),spoly%cont(ic)%y(iv)
              write(0,*)spoly%cont(ic)%x(jv),spoly%cont(ic)%y(jv)
              dist = (spoly%cont(ic)%x(iv) - spoly%cont(ic)%x(jv))**2 &
                   + (spoly%cont(ic)%y(iv) - spoly%cont(ic)%y(jv))**2
              write(0,*)"dist (squared) = ",dist
            end if
          end if
          ltri(ic) = 1
          ic = ic + 1
          cycle ICLOOP
        end if
        jv = jv + 1
      end do
      iv = iv + 1
    end do
    end if  !full_surface


!   SEE IF ANY POINTS ARE REVISITED
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert)
      jv = iv + 1
      do while (jv <= spoly%cont(ic)%num_vert)
        if ((spoly%cont(ic)%x(iv) == spoly%cont(ic)%x(jv)).and. &
            (spoly%cont(ic)%y(iv) == spoly%cont(ic)%y(jv))) then
          if (verbosity >= 2) then
            write(0,*)
            write(0,*)"WARNING! repeated vertex in panel ",pann
            write(0,*)iv,jv
          end if
          ltri(ic) = 1
          ic = ic + 1
          cycle ICLOOP
        end if
        jv = jv + 1
      end do
      iv = iv + 1
    end do


!   SEE IF ANY NON-ADJACENT SEGMENTS INTERSECT
    iv = 1
    do while (iv < spoly%cont(ic)%num_vert)
      ax = spoly%cont(ic)%x(iv)
      ay = spoly%cont(ic)%y(iv)
      bx = spoly%cont(ic)%x(iv+1)
      by = spoly%cont(ic)%y(iv+1)

      jv = iv + 2
      do while (jv <= spoly%cont(ic)%num_vert)
        cx = spoly%cont(ic)%x(jv)
        cy = spoly%cont(ic)%y(jv)
        if (jv < spoly%cont(ic)%num_vert) then
          dx = spoly%cont(ic)%x(jv+1)
          dy = spoly%cont(ic)%y(jv+1)
        else if (iv == 1) then
          jv = jv + 1
          cycle
        else
          dx = spoly%cont(ic)%x(1)
          dy = spoly%cont(ic)%y(1)
        end if

        denom = (bx-ax)*(dy-cy) - (dx-cx)*(by-ay)
        rr    = (dx-cx)*(ay-cy) - (ax-cx)*(dy-cy)
        ss    = (cx-ax)*(by-ay) - (bx-ax)*(cy-ay)


!       CHECK FOR OVERLAPPING COLINEAR SEGMENTS
        if ((check_colinear).and.(abs(denom) < 1.0d-10)) then
          if ((abs(rr) < 1.0d-10).or.(abs(ss) < 1.0d-10)) then
            if ((min(ax,bx) <= max(cx,dx)).and. &
                (max(ax,bx) >= min(cx,dx)).and. &
                (min(ay,by) <= max(cy,dy)).and. &
                (max(ay,by) >= min(cy,dy))) then
              if (verbosity >= 2) then
                write(0,*)
                write(0,*)"WARNING! overlapping colinear "// &
                          "segments in panel ",pann
                write(0,*)"contour = ",ic
                write(0,*)"iv,jv = ",iv,jv
                write(0,*)rr,ss,denom
              end if
              ltri(ic) = 1
              ic = ic + 1
              cycle ICLOOP
            end if
          end if
        end if


!       CHECK FOR INTERSECTING SEGMENTS
        if (abs(denom) >= 1.0d-10) then
          rr = rr / denom
          ss = ss / denom
          if ((0.0_rd <= ss).and.(ss <= 1.0_rd).and. &
              (0.0_rd <= rr).and.(rr <= 1.0_rd)) then
            if (verbosity >= 2) then
              write(0,*)
              write(0,*)"WARNING! segments intersect in panel ",pann
              write(0,*)iv,jv,rr,ss,denom
            end if
            ltri(ic) = 1
            ic = ic + 1
            cycle ICLOOP
          end if
        end if

        jv = jv + 1
      end do
      iv = iv + 1
    end do

    ic = ic + 1
  end do ICLOOP

  return
  end subroutine SelfIntersection

  end module CreateTrianglesM
