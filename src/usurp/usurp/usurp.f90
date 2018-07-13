! ----------------------------------------------------------------
!
! USURP: Unique Surfaces Using Ranked Polygons
! Copyright 2005-2007 The Pennsylvania State University
! Author: David Boger
!
! ----------------------------------------------------------------
!
! David Boger
! Applied Research Laboratory
! Penn State University
! phone:  (215) 682-4044 (Warminster, PA)
!         (814) 863-3055 (State College, PA)
! e-mail: dab143@only.arl.psu.edu (work)
!         dab143@psu.edu (personal)
!
! Integration of hydrodynamic forces and moments on overset surfaces.
! USURP was inspired by Larry Wigton's presentation of "Polymixsur" 
! at the 7th Overset Symposium in Huntingdon Beach, CA in 2004.
! USURP also relies on the General Polygon Clipping library (Version 2.32)
! by Alan Murta for polygon Boolean operations and triangulation and
! "Triangle" (Version 1.6) by Jonathan Shewchuk for triangulation.
! Note that Triangle is freely available but copyrighted by its author,
! Jonathan Shewchuk, and may not be sold or included in commercial products
! without a license.  Similarly, the GPC library is free for non-commerical
! use.  Anyone wishing to use the GPC library in support of a commercial
! product should email gpc@cs.man.ac.uk.
!
! USURP calls the gpc and Triangle libraries via "translator" routines, one in
! Fortran that wraps the subject and clipping polygon structures into
! an integer array and float array of fixed length, and a second in
! c that takes these arrays, converts them back into polygon structures,
! and calls the c libraries.  These routines also handle the information
! moving in the other direction.  The calls could be made from 
! Fortran to the c libraries directly, but this would involve Fortran 
! descriptors, and the resulting code would not be portable.
!
! Reference:
! Boger, D.A., and Dreyer, J.J., "Prediction of Hydrodynamic Forces and 
! Moments for Underwater Vehicles Using Overset Grids," AIAA-2006-1148, 
! 44th AIAA Aerospace Sciences Meeting and Exhibit, Reno, Nevada, 
! Jan. 9-12, 2006.
!
! compile with "make"  
! usage: usurp --help (see user's guide for details)
! the input path should be the CFD main working directory; 
! the input path is the current directory by default.
!
! ====================================================
  program usurp
! ====================================================
  use BuildGraphM               , only: BuildGraph
  use CalculateForcesAndMomentsM, only: CalculateForcesAndMoments
  use CalculatePolygonAreaM     , only: CalculatePolygonArea 
  use CalculateProjectedAreasM  , only: CalculateProjectedAreas
  use CheckInputM               , only: CheckInput
  use CheckPanelPairM           , only: CheckPanelPair
  use CheckRatiosM              , only: CheckRatios
  use CheckTrianglesM           , only: CheckTriangles
  use CommandLineM              , only: CommandLine
  use CreateNewNodeListM        , only: CreateNewNodeList
  use CreateTrianglesM          , only: CreateTriangles
  use CutOutStringsM            , only: CutOutStrings
  use DeallocatePanelsM         , only: DeallocatePanels
  use DefinePanelsM             , only: DefinePanels
  use DetermineSurfaceOverlapM  , only: DetermineSurfaceOverlap
  use DPLRM                     , only: PreReadDPLR,StorePatchesDPLR
  use ExaminePrioritiesM        , only: ExaminePriorities
  use GetUserInputM             , only: GetUserInput
  use IntrType                  , only: rd,rs
  use LoadPanelsM               , only: LoadPanels
  use LoadPanelsNPHASEM         , only: LoadPanelsNPHASE
  use MapTriQM                  , only: MapTriQ
  use my_cpu_timeM              , only: my_cpu_time
  use NetAlloc                  , only: nalloc,nentry,nnodes,npnten
  use OutputForcesAndMomentsM   , only: OutputForcesAndMoments
  use OutputOVERFLOWM           , only: OutputOVERFLOW
  use PatchInfo                 , only: num_panels,num_patches
  use PatchInfo                 , only: table,firstpass,overlaps
  use PreReadBCinM              , only: PreReadBCin
  use PreReadBCSM               , only: PreReadBCS
  use PreReadGenericM           , only: PreReadGeneric
  use PreReadNPHASEM            , only: PreReadNPHASE
  use PreReadOVERFLOWM          , only: PreReadOVERFLOW
  use PreReadUNCLEMM            , only: PreReadUNCLEM
  use ProcessPairProcedures     , only: times_called
  use ReadPanelWeightsM         , only: ReadPanelWeights
  use ResolveEdgeIntersectionsM , only: ResolveEdgeIntersections
  use RTreeTypes                , only: debug,node
  use RTreeProcedures           , only: BuildRTree,CheckForOverlapWithin
  use RTreeProcedures           , only: DeleteTree,PrintTree
  use ShrinkVertexListM         , only: ShrinkVertexList
  use StorePatchesCFDSHIPM      , only: StorePatchesCFDSHIP
  use StorePatchesGenericM      , only: StorePatchesGeneric
  use StorePatchesOVERFLOWM     , only: StorePatchesOVERFLOW
  use StorePatchesUNCLEM        , only: StorePatchesUNCLE
  use TecplotNPHASEM            , only: TecplotNPHASE
  use TecplotTrianglesM         , only: TecplotTriangles
  use TimeKeeper                , only: gpc_time
  use Types                     , only: panel
  use UnusedVerticesM           , only: UnusedVertices
  use UserInput                 , only: basis,cfd_ship,checkpanels,cpufmt
  use UserInput                 , only: debugcl,degen_tol,degen_tol2
  use UserInput                 , only: dplr,dumptree,allow_fringe
  use UserInput                 , only: freememory,full_surface,generic
  use UserInput                 , only: gpcf_time,keep_weights
  use UserInput                 , only: merge_tol,merge_tol2,my_bit_size
  use UserInput                 , only: notree,nphase,over_rel,overflow
  use UserInput                 , only: pathname,showcpu,solution_exists
  use UserInput                 , only: ttype,unclem,watertight,verbosity
  use UserInput                 , only: trap_clip,use_map,use_priority
  use WriteGridIBIM             , only: WriteGridIBI
  use WritePanelWeightsM        , only: WritePanelWeights
  use WritePatchesM             , only: WritePatches
  use WriteTriQM                , only: WriteTriQ
  use WriteMixsurFMPM           , only: WriteMixsurFMP

  implicit none

  type(node) ,pointer                 :: root
  type(panel),pointer,dimension (:)   :: p
  real(kind=rd)                       :: pvalue
  real(kind=rs)                       :: time1,time2
  integer,dimension(2)                :: seed
  integer                             :: ierr,ipan,osize
  logical                             :: run_gpc


  continue


  write(*,"(/,1x,a,/)")"This is USURP, version 2.38"


! initialization
  nullify(root)
  my_bit_size = bit_size(my_bit_size)
  pvalue  = 0.0_rd
  seed(1) = 1239
  seed(2) = 9863
  call random_seed(put=seed)


! set defined tolerances

! note: larger values of degen_tol help avoid GPC/Triangle crashes by cleaning
! up degenerate polygons, but smaller values are needed to create watertight 
! surfaces. Many GPC output degeneracies can be avoided by increasing gpc_eps
! (which encourages the merging of collinear, overlapping edges), but 
! increasing gpc_eps causes GPC to be less stable.
  degen_tol  = 1.0e-10
  degen_tol2 = degen_tol**2

! note: merge_tol is only used in one essentially obsolete routine
  merge_tol  = 1.0e-5   
  merge_tol2 = merge_tol**2


! get user input
  call CommandLine()
  call GetUserInput()
  debug = debugcl

! read the CFD solver boundary conditions
  call my_cpu_time(time1)
  if (over_rel) then
    call PreReadBCin()     !UNCLE / OVER-REL
  else if (cfd_ship) then
    call PreReadBCS()      !CFD-SHIP
  else if (overflow) then
    call PreReadOVERFLOW() !MIXSUR / FOMOCO / OVERFLOW
  else if (nphase) then
    call PreReadNPHASE()   !NPHASE
  else if (unclem) then
    call PreReadUNCLEM()   !UNCLE-M
  else if (generic) then
    call PreReadGeneric()  !Generic
  else if (dplr) then
    call PreReadDPLR()     !DPLR
  else
    write(0,*)"ERROR! flow solver could not be determined."
    stop
  end if
  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in ReadBC: ',time2-time1
  end if

! echo the MIXSUR input stream to mixsur.fmp
  if (overflow) call WriteMixsurFMP()
  
! get the grid data for the bcpatches
  if (over_rel .or. unclem) then
    call StorePatchesUNCLE()     !UNCLE / OVER-REL
  else if (cfd_ship) then
    call StorePatchesCFDSHIP()   !CFD-SHIP
  else if (overflow) then
    call StorePatchesOVERFLOW()  !MIXSUR / FOMOCO / OVERFLOW
  else if (generic) then
    call StorePatchesGeneric()   !Generic
  else if (dplr) then
    call StorePatchesDPLR()      !DPLR
  end if

  call CheckInput()

! extra output to support OVERFLOW mode for OVERFLOW users
  if (overflow) call WriteGridIBI()


! allocate memory for the panels
  call my_cpu_time(time1)
  allocate(p(num_panels),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(unit=0,fmt=*)"ERROR! allocation of p failed."
    write(unit=0,fmt=*)"num_panels = ",num_panels
    stop
  end if
  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time to allocate panel memory:',time2-time1
  end if
     

! define panel corners
  if (nphase) then
    call LoadPanelsNPHASE(p)
  else
    call LoadPanels(p)
  end if

  if ((full_surface).and.(.not.nphase)) then
    call ShrinkVertexList(p)
  end if

! calculate the geometric data for each panel
  call DefinePanels(p)

! check for cut-out strings
  call CutOutStrings(p)


! if user wants to use existing weights, read them here and then skip
! ahead to the force and moment integration
  if (keep_weights) then
    call ReadPanelWeights(p)

  else !if .not.keep_weights

!   calculate the initial polygon areas
!   note: we already calculated areas in a manner consistent with UNCLE
!   which is useful for integrating the forces, but now we need to 
!   calculate the area in a manner consistent with USURP, which will
!   be useful for calculating ratios of the initial and final projected
!   polygons.
    call my_cpu_time(time1)
    do ipan = 1,num_panels
      if (p(ipan)%iblank == 0) then
        p(ipan)%ratio = 0.0_rd
      else if ((p(ipan)%iblank /= 1).and. &
            (over_rel.or.unclem).and.(.not.allow_fringe)) then
        p(ipan)%ratio = 0.0_rd
      else
        p(ipan)%ratio = CalculatePolygonArea(p(ipan)%poly,ipan)
      end if
    end do
    if (showcpu) then
      call my_cpu_time(time2)
      write(0,cpufmt)'cpu time in initial areas:',time2-time1
    end if


!   debug option: check whether two specified panels overlap
    if (checkpanels) then
      call CheckPanelPair(p)
    end if

    if (basis == "patch") then
      osize = (num_panels*num_patches)/my_bit_size
      if (mod(num_panels*num_patches,my_bit_size) /= 0) osize = osize + 1
      allocate(overlaps(osize),stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! allocate failed for overlaps."
        write(0,*)"num_panels, num_patches = ",num_panels,num_patches
        write(0,*)"my_bit_size = ",my_bit_size
        write(0,*)"osize = ",osize
        stop
      end if
      overlaps = 0
    end if

!   build (and print) an R-tree
    if (.not.notree) then

      call my_cpu_time(time1)
      call BuildRTree( root, p)
      if (showcpu) then
        call my_cpu_time(time2)
        write(0,cpufmt)'cpu time in BuildRTree:',time2-time1
      end if

      if (dumptree) then
!       output the r-tree
        call my_cpu_time(time1)
        open(unit=3,file="tree.dat")
        write(3,*)"variables = x,y,z"
        call PrintTree( root ,1 ,1, pvalue)
        close(unit=3)
        if (showcpu) then
          call my_cpu_time(time2)
          write(0,cpufmt)'cpu time in PrintTree:',time2-time1
        end if
      end if
    end if


!   Clip overlapping panels. Find overlapping panels using either 
!   structured surface bounding boxes or an R-Tree search.
!   First pass:  Build overlap info if basis = "patch"
!   Second pass: Perform polygon clipping

    firstpass = .true.

    TWO_PASS_LOOP: do

      if (basis == "panel") firstpass = .false.

      times_called = 0
      gpc_time     = 0
      gpcf_time    = 0
      run_gpc      = .not.firstpass

      if (run_gpc.and.trap_clip) then
        write(0,*)"WARNING! Enabling --trap=clip will slow ", &
                  "this portion of the code."
      end if

      call my_cpu_time(time1)
      if (notree) then
        call DetermineSurfaceOverlap(p,run_gpc)
      else
        call CheckForOverlapWithin(root,run_gpc)
      end if
      if (showcpu) then
        call my_cpu_time(time2)
        if (notree) then
          write(0,cpufmt)'cpu time in DetermineSurfaceOverlap:', &
                        time2-time1-gpc_time
        else
          write(0,cpufmt)'cpu time in CheckForOverlapWithin:', &
                        time2-time1-gpc_time
        end if
        if (run_gpc) then
          write(0,cpufmt)'cpu time in gpc:',gpc_time
          write(0,cpufmt)'  > time in gpc_c: ',gpcf_time(1)
          write(0,cpufmt)'  > time in dealloc: ',gpcf_time(2)
          write(0,cpufmt)'  > time in alloc: ',gpcf_time(3)
        end if
      end if

      if (verbosity >= 1) then
        if (notree) then
          write(0,'(1x,a,i7,a)')'DetermineSurfaceOverlap found ', &
             times_called,' valid overlapping bounding boxes'
        else
          write(0,'(1x,a,i7,a)')'CheckForOverlapWithin found ', &
             times_called,' valid overlapping bounding boxes'
        end if
      end if

      if (firstpass) then
        call BuildGraph(p)
        if (use_priority.and.(verbosity >= 1)) call ExaminePriorities()
      else
        exit TWO_PASS_LOOP
      end if

      firstpass = .false.
    end do TWO_PASS_LOOP


!   memory clean-up
    if (freememory.and.(.not.notree)) then
      call my_cpu_time(time1)
      call DeleteTree(root)
      deallocate(root); nalloc = nalloc - 1
      nnodes = nnodes - 1
      if (showcpu) then
        call my_cpu_time(time2)
        write(0,cpufmt)'cpu time in DeleteTree:',time2-time1
      end if
    end if
    if (basis == "patch") then
      deallocate(table); nalloc = nalloc - 1
    end if

    if (trap_clip) then
      open(unit=90,file=trim(pathname)//"scratch_clipping.dat")
      close(unit=90,status="delete")
    end if


!   calculate the final polygon areas and ratios
    call my_cpu_time(time1)
    do ipan = 1,num_panels
      if (p(ipan)%iblank == 0) then
        p(ipan)%ratio = 0.0_rd
      else
        if (p(ipan)%ratio == 0.0) then
!         write(*,*)"warning! zero area for panel ",ipan
!         stop
        else
          p(ipan)%ratio = CalculatePolygonArea(p(ipan)%poly,ipan) &
                        / p(ipan)%ratio
        end if
      end if
    end do
    if (showcpu) then
      call my_cpu_time(time2)
      write(0,cpufmt)'cpu time in final areas:',time2-time1
    end if

  end if !keep_weights


! check the area ratios
  call CheckRatios(p)


  if (solution_exists) then

!   calculate the forces and moments acting on all panels
    call CalculateForcesAndMoments(p)


!   output the forces and moments
    if (overflow) then
      call OutputOVERFLOW()
    else 
      call OutputForcesAndMoments()
    end if
 
  end if


! calculate wetted and projected areas for all panels
  call CalculateProjectedAreas(p)


  if (keep_weights) then

!   write a TriQ file using usurp.map
    if (solution_exists.and.use_map) call MapTriQ(p)
  else !if .not.keep_weights

!   write the panel weights to a data file
    call WritePanelWeights(p)

!   write the patches to a Tecplot file
    if (.not.nphase) call WritePatches(p)

!   create 3D coordinates for any new vertices created during clipping
!   For ttype == "Triangle" (watertight or not) this must be done prior
!   to CreateTriangle so that new nodes don't appear to be repeated nodes
!   in SelfIntersect
    write(*,*) ' 9'

    if (ttype == "Triangle") then
      if (watertight) then
        call ResolveEdgeIntersections(p)
      else
        call CreateNewNodeList(p)
      end if
    end if

!   triangulate whatever needs to be triangulated
    if (ttype /= "none") then
      if (verbosity >= 1) write(0,*)"Calling CreateTriangles"
      call CreateTriangles(p)
      if (verbosity >= 1) write(0,*)"Finished CreateTriangles"
    end if


!   for GPC, the new node numbering must occur after the triangulation
!   step because GPC introduces new nodes
    if (ttype == "GPC") then
      call CreateNewNodeList(p)
    end if


    if (full_surface) call UnusedVertices(p)
    if (watertight.and.(verbosity >= 2)) call CheckTriangles(p)
    if (full_surface) call WriteTriQ(p)

    write(*,*) ' 10'

!   write triangles to a Tecplot file
    if (nphase) then
      call TecplotNPHASE(p)
    else
      call TecplotTriangles(p)
    end if

  end if

  write(*,*) ' 11'
      

! free memory
  if (freememory) then

    call DeallocatePanels(p)

    if (abs(nalloc)+abs(nentry)+abs(nnodes)+abs(npnten) /= 0) then
      write(0,*)
      write(0,*)"WARNING! possible memory leak detected."
      write(0,*)"Be sure to deallocate all allocated memory."
      write(0,*)
      write(0,*)"net allocations   : ",nalloc
      write(0,*)"net entries       : ",nentry
      write(0,*)"net nodes         : ",nnodes
      write(0,*)"net entry pointers: ",npnten
    end if
  end if

  end program usurp
