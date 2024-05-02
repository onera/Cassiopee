SUBROUTINE usurp( nzone, nit, njt, nkt, iblank, ncellmax, coord, nptsmax,&
     ratio)

  use BuildGraphM ,only: BuildGraph
  use CalculateForcesAndMomentsM, only: CalculateForcesAndMoments
  use CalculatePolygonAreaM,only: CalculatePolygonArea
  use CalculateProjectedAreasM, only: CalculateProjectedAreas 
  use CheckInputM   ,only: CheckInput
  use CheckPanelPairM ,only: CheckPanelPair
  use CheckRatiosM ,only: CheckRatios
  use CheckTrianglesM           , only: CheckTriangles
  use CommandLineM              , only: CommandLine
  use CreateNewNodeListM        , only: CreateNewNodeList
  use CreateTrianglesM          , only: CreateTriangles
  use CutOutStringsM,only: CutOutStrings
  use DeallocatePanelsM         , only: DeallocatePanels
  use DefinePanelsM ,only: DefinePanels
  use DetermineSurfaceOverlapM , only: DetermineSurfaceOverlap
  use EdgeData  , only: num_edges
  use ExaminePrioritiesM , only: ExaminePriorities
  use IntrType      ,only: rd,rs
  use LoadPanelsM ,only: LoadPanels
  use LoadPanelsNPHASEM ,only: LoadPanelsNPHASE
  use MapTriQM,only: MapTriQ
  use NetAlloc                  , only: nalloc,nentry,nnodes,npnten
  use OutputForcesAndMomentsM   , only: OutputForcesAndMoments
  use PatchInfo     ,only: num_panels,num_patches,icomp,num_components,num_scalars
  use PatchInfo     ,only: i1,i2,j1,j2,k1,k2,ncomp_surfs,block
  use PatchInfo     ,only: tag_list, overlaps, firstpass, table
  use ProcessPairProcedures ,only: times_called
  use ReadPanelWeightsM,only: ReadPanelWeights
  use ResolveEdgeIntersectionsM , only: ResolveEdgeIntersections
  use RTreeProcedures ,only: DeleteTree, BuildRTree,PrintTree,CheckForOverlapWithin
  use RTreeTypes    ,only: debug,node
  use Types         ,only: panel, scal,surf
  use ShrinkVertexListM ,only: ShrinkVertexList
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch,EchoPatchData
  use TecplotNPHASEM ,only: TecplotNPHASE
  use TecplotTrianglesM ,only: TecplotTriangles
  use TimeKeeper ,only: gpc_time
  use UnusedVerticesM           , only: UnusedVertices
  use UserInput ,only: gpcf_time, verbosity
  use UserInput ,only: merge_tol,merge_tol2,my_bit_size
  use UserInput ,only: debugcl,degen_tol,degen_tol2
  use UserInput ,only: disjoin, default_disjoin, default_basis
  use UserInput ,only: allow_fringe, colormap,debugcl
  use UserInput ,only: basis, checkpanels, checkpatches
  use UserInput ,only: dumptree,freememory, over_rel, unclem
  use UserInput ,only: full_surface,icpanel1,icpanel2,icpatch1,icpatch2
  use UserInput ,only: ignore_solution,nphase,keep_weights,use_map
  use UserInput ,only: notree,trap_clip,use_priority,pathname,solution_exists
  use UserInput ,only: ttype,watertight
  use VertexData ,only: num_nodes,xv
  use WritePanelWeightsM        , only: WritePanelWeights
  use WritePatchesM             , only: WritePatches
  use WriteTriQM                , only: WriteTriQ
!

  IMPLICIT NONE

! IN
  integer*4   :: ncellmax, nptsmax
  integer*4, dimension(nzone) :: nit, njt, nkt    
  integer*4, dimension(ncellmax) :: iblank
  real(kind=rd),dimension(nptsmax,3) :: coord
! LOCAL
  type(node) ,pointer               :: root
  type(panel),pointer,dimension (:) :: p
  real(kind=rd)                     :: pvalue
  integer,dimension(2)              :: seed
  integer                           :: ierr,ipan,osize
  logical                           :: run_gpc
  integer                           :: nzone
  integer,allocatable,dimension (:) :: nig,njg,nkg
  integer                           :: ip,n
  real(kind=rd),allocatable,dimension (:,:,:) :: xb,yb,zb
  integer,allocatable,dimension (:,:,:) :: ibv
  integer ::i,j,k, iv, imax, jmax, kmax, ncell, nnode, nnode2
  integer :: dir1,dir2,dir3, idir
  integer,dimension (3) :: i1s,i1e
! OUT
  real(kind=rd),dimension(ncellmax) :: ratio

  ! initialization
  nullify(root)
  my_bit_size = bit_size(my_bit_size)
  pvalue  = 0.0_rd
  seed(1) = 1239
  seed(2) = 9863
  !call random_seed(put=seed)
  call random_seed()

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

!------------------------------------------------------------------------------
! Command line arguments
!------------------------------------------------------------------------------
  call CommandLine()

!------------------------------------------------------------------------------
! valeurs par defaut
!------------------------------------------------------------------------------

  if (default_disjoin) disjoin = .false.
  if (default_basis)   basis = "panel"
  
!------------------------------------------------------------------------------
! convert input data (connectivity, iblank)
!------------------------------------------------------------------------------
  num_patches = nzone 

  call AllocateStructuredPatch()
  allocate(ncomp_surfs(0:0),stat=ierr)
  allocate(tag_list(0:0),stat=ierr)

! read imin, imax etc for each block 
  ip = 0
  allocate(nig(num_patches),njg(num_patches),nkg(num_patches)); nalloc = nalloc + 3

  do while (ip < num_patches)
    ip = ip + 1
    nig(ip) = nit(ip)
    njg(ip) = njt(ip)
    nkg(ip) = nkt(ip)

    block(ip)  = ip
    i1(ip) = 1 
    i2(ip) = nig(ip)
    j1(ip) = 1
    j2(ip) = njg(ip)
    k1(ip) = 1
    k2(ip) = nkg(ip)
    icomp(ip) = 0
    ncomp_surfs(0) = ncomp_surfs(0) + 1
    call CountStructuredPatch(ip)

  end do
  ! screen output to confirm results
  call EchoPatchData()

  !---------------------------
  ! Read the grid connectivity
  !---------------------------
  ! assume there is no flow data
  allocate(xv(3,num_nodes)); nalloc = nalloc + 1
  allocate(scal(1,num_panels)); nalloc = nalloc + 1
  ncell = 1
  nnode = 1
  
  do n = 1,num_patches
     allocate( xb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
     allocate( yb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
     allocate( zb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
     imax = nig(n)
     jmax = njg(n)
     kmax = nkg(n)
     
     if ( nig(n) == 1 ) then 
        imax = 2         
     else if ( njg(n) == 1 ) then
        jmax = 2
     else
        kmax = 2    
     endif
     allocate(ibv(2:imax,2:jmax,2:kmax)); nalloc = nalloc + 1
         
     !read iblank array
     do k = 2, kmax
     do j = 2, jmax
     do i = 2, imax
         ibv(i,j,k) = iblank(ncell)
         ncell = ncell+1
     enddo
     enddo
     enddo

     ! grid connectivity
     
     do k = 1, nkg(n)
     do j = 1, njg(n)
     do i = 1, nig(n)
        xb(i,j,k) = coord(nnode,1) 
        yb(i,j,k) = coord(nnode,2)
        zb(i,j,k) = coord(nnode,3)                      
        nnode = nnode + 1    
     enddo
     enddo
     enddo

     ! fill in the grid for each patch
     do ip = 1,num_patches  
        if (block(ip) == n) then
           ! store the grid points for the patch
           iv = surf(ip)%first_node
           do k = k1(ip),k2(ip)
           do j = j1(ip),j2(ip)
           do i = i1(ip),i2(ip)
              xv(1,iv) = xb(i,j,k)
              xv(2,iv) = yb(i,j,k)
              xv(3,iv) = zb(i,j,k)
      
              iv = iv + 1
           end do
           end do
           end do

           ! figure out which direction is constant
           if (i1(ip) == i2(ip)) then
              dir1 = 1
              dir2 = 2
              dir3 = 3
           else if (j1(ip) == j2(ip)) then
              dir1 = 2
              dir2 = 3
              dir3 = 1
           else if (k1(ip) == k2(ip)) then
              dir1 = 3
              dir2 = 1
              dir3 = 2
           else
              write(0,*)"ERROR!  no constant direction on patch ",ip
              stop
           end if

           ! store the i-blank variable
           i1s(1) = i1(ip)
           i1e(1) = i2(ip)
           i1s(2) = j1(ip)
           i1e(2) = j2(ip)
           i1s(3) = k1(ip)
           i1e(3) = k2(ip)
           if (i1s(dir2) == i1e(dir2)) then
              write(0,*)"error in patch ",ip
              stop
           end if
           if (i1s(dir3) == i1e(dir3)) then
              write(0,*)"error in patch ",ip
              stop
           end if
           
           i1s(dir2) = i1s(dir2) + 1
           i1s(dir3) = i1s(dir3) + 1
           
           ipan = surf(ip)%first_panel
           
           do k = i1s(3),i1e(3)
           do j = i1s(2),i1e(2)
           do i = i1s(1),i1e(1)
              
              if (dir1 == 1) then
                 if (i == 1) then
                    scal(1,ipan) = ibv(2,j,k)
                 else
                    scal(1,ipan) = ibv(nig(n),j,k)
                 end if
              else if (dir1 == 2) then
                 if (j == 1) then
                    scal(1,ipan) = ibv(i,2,k)
                 else
                    scal(1,ipan) = ibv(i,njg(n),k)               
                 end if
              else if (dir1 == 3) then
                 if (k == 1) then
                    scal(1,ipan) = ibv(i,j,2)
                 else
                    scal(1,ipan) = ibv(i,j,nkg(n))
                 end if
              else
                 write(0,*)"ERROR! invalid dir1 in StorePatchesGeneric"
                 stop
              end if

              
              ipan = ipan + 1
           end do
           end do
           end do
        
        endif
     enddo
     
     deallocate(ibv,zb,yb,xb); nalloc = nalloc - 4
  
   enddo
   deallocate(nkg,njg,nig); nalloc = nalloc - 3

!----------------------------------------------------------------------------------
!  continue usurp 
!----------------------------------------------------------------------------------
  call CheckInput()

! allocate memory for the panels
  allocate(p(num_panels),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(unit=0,fmt=*)"ERROR! allocation of p failed."
    write(unit=0,fmt=*)"num_panels = ",num_panels
    stop
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
     do ipan = 1,num_panels
        if (p(ipan)%iblank == 0) then
           p(ipan)%ratio = 0.0_rd
        else if ((p(ipan)%iblank /= 1).and. &
             (over_rel.or.unclem).and.(.not.allow_fringe)) then
           p(ipan)%ratio = 0.0_rd
        else
           p(ipan)%ratio = CalculatePolygonArea(p(ipan)%poly,ipan)
        endif
     enddo
     
     !debug option: check whether two specified panels overlap
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

     ! build (and print) an R-tree
     if (.not.notree) then
        
        call BuildRTree( root, p)
      
        if (dumptree) then
           ! output the r-tree
           open(unit=3,file="tree.dat")
           write(3,*)"variables = x,y,z"
           call PrintTree( root ,1 ,1, pvalue)
           close(unit=3)
        end if
     endif
             
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
        
        if (notree) then
           call DetermineSurfaceOverlap(p,run_gpc)
        else             
           call CheckForOverlapWithin(root,run_gpc)
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

     end do TWO_PASS_LOOP

     ! memory clean-up
     freememory = .true.
     if (freememory.and.(.not.notree)) then
        call DeleteTree(root)
        deallocate(root); nalloc = nalloc - 1
        nnodes = nnodes - 1
     end if
     
     if (basis == "patch") then
        deallocate(table); nalloc = nalloc - 1
     end if
     
     if (trap_clip) then
        open(unit=90,file=trim(pathname)//"scratch_clipping.dat")
        close(unit=90,status="delete")
     end if

     ! calculate the final polygon areas and ratios
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
     
  end if !keep_weights
  
  ! check the area ratios
  call CheckRatios(p)
  
  if (solution_exists) then
!   calculate the forces and moments acting on all panels
    call CalculateForcesAndMoments(p)
    !   output the forces and moments
    call OutputForcesAndMoments()
  end if
  
! calculate wetted and projected areas for all panels
  call CalculateProjectedAreas(p)

  

!!$  if (keep_weights) then
!!$     !   write a TriQ file using usurp.map
!!$     if (solution_exists.and.use_map) call MapTriQ(p)
!!$  else !if .not.keep_weights
!!$     
!!$     !   write the panel weights to a data file
!!$     call WritePanelWeights(p)
!!$     
!!$!   write the patches to a Tecplot file
!!$    if (.not.nphase) call WritePatches(p)
!!$
!!$!   create 3D coordinates for any new vertices created during clipping
!!$!   For ttype == "Triangle" (watertight or not) this must be done prior
!!$!   to CreateTriangle so that new nodes don't appear to be repeated nodes
!!$!   in SelfIntersect
!!$
!!$
!!$    if (ttype == "Triangle") then
!!$      if (watertight) then
!!$        call ResolveEdgeIntersections(p)
!!$      else
!!$        call CreateNewNodeList(p)
!!$      end if
!!$    end if
!!$!   triangulate whatever needs to be triangulated
!!$    if (ttype /= "none") then
!!$      if (verbosity >= 1) write(0,*)"Calling CreateTriangles"
!!$      call CreateTriangles(p)
!!$      if (verbosity >= 1) write(0,*)"Finished CreateTriangles"
!!$    end if
!!$
!!$
!!$!   for GPC, the new node numbering must occur after the triangulation
!!$!   step because GPC introduces new nodes
!!$    if (ttype == "GPC") then
!!$      call CreateNewNodeList(p)
!!$    end if
!!$
!!$
!!$    if (full_surface) call UnusedVertices(p)
!!$    if (watertight.and.(verbosity >= 2)) call CheckTriangles(p)
!!$    if (full_surface) call WriteTriQ(p)
!!$
!!$
!!$
!!$!   write triangles to a Tecplot file
!!$    if (nphase) then
!!$      call TecplotNPHASE(p)
!!$    else
!!$      call TecplotTriangles(p)
!!$    end if
!!$
!!$ end if
 
! convert for output : all data are defined in nodes
 nnode = 1
 nnode2 = 1

 do ip = 1, num_patches
    if (i1(ip) == i2(ip)) then
       nit(ip) = 1
       njt(ip) = j2(ip)-j1(ip)+1
       nkt(ip) = k2(ip)-k1(ip)+1
    else if (j1(ip) == j2(ip)) then
       nit(ip) = i2(ip)-i1(ip)+1
       njt(ip) = 1
       nkt(ip) = k2(ip)-k1(ip)+1
    else if (k1(ip) == k2(ip)) then
       nit(ip) = i2(ip)-i1(ip)+1
       njt(ip) = j2(ip)-j1(ip)+1
       nkt(ip) = 1
    end if


    do iv=surf(ip)%first_node,surf(ip)%last_node
       do idir = 1,3
          coord(nnode,idir)  = xv(idir,iv)
       enddo
       nnode = nnode + 1
    enddo
    do ipan = surf(ip)%first_panel,surf(ip)%last_panel 
       ratio(nnode2) = p(ipan)%ratio
       nnode2 = nnode2 + 1    
    enddo
 enddo
 
 ! free memory
 if (freememory) then

    call DeallocatePanels(p)
    
!!$    if (abs(nalloc)+abs(nentry)+abs(nnodes)+abs(npnten) /= 0) then
!!$      write(0,*)
!!$      write(0,*)"WARNING! possible memory leak detected."
!!$      write(0,*)"Be sure to deallocate all allocated memory."
!!$      write(0,*)
!!$      write(0,*)"net allocations   : ",nalloc
!!$      write(0,*)"net entries       : ",nentry
!!$      write(0,*)"net nodes         : ",nnodes
!!$      write(0,*)"net entry pointers: ",npnten
!!$    end if
  end if

  num_nodes = 0
  num_components=0
  num_patches=0
  num_scalars=0
  num_panels=0
  num_edges = 0
  nalloc = 0
  return 
END SUBROUTINE usurp
