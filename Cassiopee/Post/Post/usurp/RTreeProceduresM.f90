! =====================================================================
  module RTreeProcedures
! =====================================================================
! This module contains all of the procedures associated with building
! and searching an R-tree (Guttman, 1984).  They are defined in this
! module for two reasons.  First, most pass pointers through their arguments
! and so must have explicit interfaces, and containing them within this
! module satisfies that requirement.  Second, containing them within this
! module limits their scope so that R-tree procedure names can't clash
! with procedure names in the main program.

  implicit none

  public  :: BuildRTree,PrintTree,CheckForOverlapWithin,DeleteTree
  public  :: SearchTree

  private :: CheckForOverlapBetween,CreateNode,FindOverlap
  private :: Insert,AdjustTree,LinearSplitNode,RemoveEntry
  private :: ChooseLeaf,LinearPickSeeds,QuadraticPickSeeds,PickNext

  contains

!   =========================================================
    subroutine BuildRTree(root,p)
!   =========================================================
!   subroutine BuildRTree (and related procedures)
!   David Boger
!   November 2004
! 
!   Fortran 90 implementation of Antonin Guttman's paper
!   "R-trees: A dynamic index structure for spatial searching"
!   1984 ACM
!
    use NetAlloc  ,only: nalloc,nentry
    use Types     ,only: panel
    use RTreeTypes,only: node,entry,debug,ndir

    implicit none

    type (node) , pointer :: root
    type (entry), pointer :: temp_entry
    type (panel), pointer, dimension (:) :: p

    integer :: idir,ip,ip1,ip2
 

    continue


    nullify(root,temp_entry)


!   set up the r-tree root
    call CreateNode(root)
    root%leaf = .true.

    ip1 = lbound(p,1)
    ip2 = ubound(p,1)


!   insert each polygon into the tree
    do ip = ip1,ip2


!     stow each panel into an entry
      allocate(temp_entry)
      nalloc = nalloc + 1
      nentry = nentry + 1
      temp_entry%pp => p(ip)
      nullify(temp_entry%pn)
      do idir = 1,ndir
        temp_entry%I(idir,1) = p(ip)%bound(idir,1)
        temp_entry%I(idir,2) = p(ip)%bound(idir,2)
      end do


!     insert each entry into the tree
      if (debug) then
        write(unit=*,fmt=*)
        write(unit=*,fmt=*)"**********************************"
        write(unit=*,fmt=*)"calling insert to insert panel ",ip
      end if
      call Insert(root,temp_entry)

    end do

    return
    end subroutine BuildRTree

!   ======================================================
    recursive subroutine CheckForOverlapWithin (T,run_gpc)
!   ======================================================

    use IntrType             ,only: rs
    use OverlappedM          ,only: Overlapped
    use PatchInfo            ,only: block
    use ProcessPairProcedures,only: times_called,ProcessPair
    use RTreeTypes           ,only: node,debug
    use TimeKeeper           ,only: gpc_time
    use UserInput            ,only: disjoin

    implicit none

    type (node), pointer :: T
    real(kind=rs)        :: time43
    integer              :: i,j
!   logical              :: overlapped
    logical              :: ran_gpc
    logical,intent(in)   :: run_gpc


    continue


    if (debug) write(unit=*,fmt=*)"check for overlap within node ", &
        T%node_number

    if (.not.T%leaf) then

!     test for overlap within the sub-nodes of a tree
      do i = 1,T%entries
        call CheckForOverlapWithin(T%E(i)%pe%pn,run_gpc)
      end do

!     test for overlap among the sub-nodes of a tree
      do i = 1,T%entries-1
      do j = i+1,T%entries
        if (overlapped(T%E(i)%pe%I,T%E(j)%pe%I)) then
          call CheckForOverlapBetween(T%E(i)%pe%pn,T%E(j)%pe%pn,run_gpc)
        end if
      end do
      end do

    else

!     test for overlap within a leaf
      do i = 1,T%entries-1
      do j = i+1,T%entries

        if (overlapped(T%E(i)%pe%I,T%E(j)%pe%I).and. &
            (T%E(i)%pe%pp%isurf /= T%E(j)%pe%pp%isurf).and. &
            (block(T%E(i)%pe%pp%isurf) /= block(T%E(j)%pe%pp%isurf))) then

          if ((T%E(i)%pe%pp%icomp == T%E(j)%pe%pp%icomp) &
             .or.(.not.disjoin)) then

            times_called = times_called + 1

            call ProcessPair(T%E(i)%pe%pp,T%E(j)%pe%pp, &
                             time43,run_gpc,ran_gpc)

            if (ran_gpc) gpc_time = gpc_time + time43

          end if
        end if
      end do
      end do
    end if

    return
    end subroutine CheckForOverlapWithin

!   =========================================================
    recursive subroutine CheckForOverlapBetween(n,nn, run_gpc)
!   =========================================================
    use IntrType             ,only: rs
    use RTreeTypes           ,only: node,debug
    use ProcessPairProcedures,only: times_called,ProcessPair
    use TimeKeeper           ,only: gpc_time
    use OverlappedM          ,only: Overlapped
    use PatchInfo            ,only: block
    use UserInput            ,only: disjoin

    implicit none

    type(node),pointer    :: n,nn
    real(kind=rs)         :: time43
    integer               :: i,j
!   logical               :: overlapped
    logical               :: ran_gpc
    logical,intent(in)    :: run_gpc


    continue


    if (N%node_number == NN%node_number) stop

    if (debug) write(unit=*,fmt=*)"check for overlap between ", &
                   N%node_number,NN%node_number

    do i = 1,N%entries
    do j = 1,NN%entries

      if (overlapped(N%E(i)%pe%i,NN%E(j)%pe%i)) then

        if (debug) write(unit=*,fmt=*)"node", &
             N%node_number,"e",i,"overlaps ",NN%node_number,j

        if (N%leaf) then
          if (N%E(i)%pe%pp%isurf /= NN%E(j)%pe%pp%isurf) then
          if (block(N%E(i)%pe%pp%isurf) /= block(NN%E(j)%pe%pp%isurf)) then
            if ((N%E(i)%pe%pp%icomp == NN%E(j)%pe%pp%icomp) &
                  .or.(.not.disjoin)) then

              times_called = times_called + 1

              call ProcessPair(N%E(i)%pe%pp,NN%E(j)%pe%pp,time43, &
                               run_gpc,ran_gpc)

              if (ran_gpc) gpc_time = gpc_time + time43

            end if
          end if
          end if
        else
          call CheckForOverlapBetween(N%E(i)%pe%pn,NN%E(j)%pe%pn,run_gpc)
        end if
      end if
    end do
    end do

    return
    end subroutine CheckForOverlapBetween

!   ==========================================================
    recursive subroutine PrintTree(T,level,level_entry,pvalue)
!   ==========================================================
!   purpose: print the tree out to a Tecplot file

    use IntrType  ,only: rd
    use RTreeTypes,only: node,debug,ndir,MaxEntries

    implicit none

    type (node), pointer     :: T

    real(kind=rd),intent(in) :: pvalue
    real(kind=rd)            :: cvalue

    integer,intent(in)       :: level,level_entry
    integer                  :: i,unum

    character(len=80)        :: tag1
    character(len=20)        :: myfmt10,c1,c2,c3,c4


    continue


    c1 = "zone t = "
    c2 = "object, l="
    c3 = "node, l="
    c4 = "e="

    unum = 30 + level

    myfmt10="(a,i3,a,2i6,a1)"

    if (T%leaf) then
      if (debug) write(unit=*,fmt=*)"node ",T%node_number, &
                      " is a leaf with ",T%entries," entries"
      do i = 1,T%entries

        if (level == 1) then 
          cvalue = real(i) 
        else    
          cvalue = pvalue-0.5/float((MaxEntries+1)**(level-2)) &
                 + float(i)/float((MaxEntries+1)**(level-1))
        end if  

!       if (debug) write(unit=*,fmt=*)"printing panel ", &
!                  T%E(i)%pe%pp%panel_number
        tag1 = trim(c1)//char(34)//trim(c2)
        write(unit=unum,fmt=myfmt10)trim(tag1),level,trim(c4), &
                                 T%E(i)%pe%pp%panel_number,i,char(34)
        write(unit=unum,fmt=*)"i=2"
        write(unit=unum,fmt=*)"j=2"
        if (ndir == 3) write(unit=unum,fmt=*)"k=2"
        write(unit=unum,fmt=*)"f=block"
        write(unit=unum,fmt=*)"varlocation=([4]=cellcentered)"

!       write block corners (block format)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(1,1), &
                              T%E(i)%pe%pp%bound(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(1,1), &
                              T%E(i)%pe%pp%bound(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(1,1), &
                              T%E(i)%pe%pp%bound(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(1,1), &
                              T%E(i)%pe%pp%bound(1,2)

        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(2,1), &
                              T%E(i)%pe%pp%bound(2,1)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(2,2), &
                              T%E(i)%pe%pp%bound(2,2)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(2,1), &
                              T%E(i)%pe%pp%bound(2,1)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(2,2), &
                              T%E(i)%pe%pp%bound(2,2)

        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(3,1), &
                              T%E(i)%pe%pp%bound(3,1)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(3,1), &
                              T%E(i)%pe%pp%bound(3,1)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(3,2), &
                              T%E(i)%pe%pp%bound(3,2)
        write(unit=unum,fmt=*)T%E(i)%pe%pp%bound(3,2), &
                              T%E(i)%pe%pp%bound(3,2)
        write(unit=unum,fmt=*)cvalue

      end do
    else
      if (debug) write(unit=*,fmt=*)"node ",T%node_number, &
                 " is a node with ",T%entries," entries"
      do i = 1,T%entries

        if (level == 1) then 
          cvalue = real(i) 
        else    
          cvalue = pvalue-0.5/float((MaxEntries+1)**(level-2)) &
                 + float(i)/float((MaxEntries+1)**(level-1))
        end if  

        tag1 = trim(c1)//char(34)//trim(c3)
        write(unit=unum,fmt=myfmt10)trim(tag1),level,trim(c4), &
                        level_entry,i,char(34)
        if (level == 1) write(unit=unum,fmt=*)"c=black"
        write(unit=unum,fmt=*)"i=2"
        write(unit=unum,fmt=*)"j=2"
        if (ndir == 3) write(unit=unum,fmt=*)"k=2"

        write(unit=unum,fmt=*)"f=block"
        write(unit=unum,fmt=*)"varlocation=([4]=cellcentered)"

!       write block corners
        write(unit=unum,fmt=*)T%E(i)%pe%I(1,1),T%E(i)%pe%I(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%I(1,1),T%E(i)%pe%I(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%I(1,1),T%E(i)%pe%I(1,2)
        write(unit=unum,fmt=*)T%E(i)%pe%I(1,1),T%E(i)%pe%I(1,2)

        write(unit=unum,fmt=*)T%E(i)%pe%I(2,1),T%E(i)%pe%I(2,1)
        write(unit=unum,fmt=*)T%E(i)%pe%I(2,2),T%E(i)%pe%I(2,2)
        write(unit=unum,fmt=*)T%E(i)%pe%I(2,1),T%E(i)%pe%I(2,1)
        write(unit=unum,fmt=*)T%E(i)%pe%I(2,2),T%E(i)%pe%I(2,2)

        write(unit=unum,fmt=*)T%E(i)%pe%I(3,1),T%E(i)%pe%I(3,1)
        write(unit=unum,fmt=*)T%E(i)%pe%I(3,1),T%E(i)%pe%I(3,1)
        write(unit=unum,fmt=*)T%E(i)%pe%I(3,2),T%E(i)%pe%I(3,2)
        write(unit=unum,fmt=*)T%E(i)%pe%I(3,2),T%E(i)%pe%I(3,2)
        write(unit=unum,fmt=*)cvalue

        call PrintTree(T%E(i)%pe%pn,level+1,i,cvalue)
      end do
    end if

    return
    end subroutine PrintTree

!   =============================================
    subroutine Insert(T,E)
!   =============================================
!   purpose: insert a new entry E into an R-Tree

    use NetAlloc  ,only: nalloc,nentry
    use RTreeTypes,only: node,entry,debug,ndir,MaxEntries

    implicit none

    type (node),  pointer :: T,L,LL
    type (entry), pointer :: E
    integer               :: i,idir,istat


    continue


    nullify(L)
    nullify(LL)

!   (I1)[find position for new record]:
!   invoke ChooseLeaf to select a leaf node L in which to place E
    if (debug) write(unit=*,fmt=*)"Insert: calling ChooseLeaf"
    if (debug) write(unit=*,fmt=*)"T is node ",T%node_number
    if (debug) write(unit=*,fmt=*)"E points to ",E%pp%panel_number
    call ChooseLeaf(T,E,L)

!   (I2)[add record to leaf node]: if L has room for another
!   entity, install E.
    if (L%entries < MaxEntries) then
      if (debug) write(unit=*,fmt=*)"Insert: node ",L%node_number, &
                 " has room to install entry"
      L%entries = L%entries + 1
      if (debug) write(unit=*,fmt=*)"Insert: node ",L%node_number, &
                 " now has ",L%entries, " entries"
      L%E(L%entries)%pe => E

!   otherwise invoke splitnode to obtain L and LL containing E
!   and all the old entries of L
    else
      if (debug) then
        write(unit=*,fmt=*)"Insert: node ",L%node_number, &
                 " does not have room to install entry"
        write(unit=*,fmt=*)"Insert: calling LinearSplitNode"
      end if
      call LinearSplitNode(L,LL,E)
      if (.not.associated(LL)) then
        write(unit=0,fmt=*)"ERROR! LL not associated after LinearSplitNode"
        stop
      end if
    end if

!   (I3)[propagate changes upward]: invoke AdjustTree on L,
!   also passing LL if a split was performed
    if (debug) write(unit=*,fmt=*)"Insert: calling AdjustTree"
    call AdjustTree(L,LL)

    if (associated(L%parent)) then
      write(unit=0,fmt=*)"ERROR! L is not the root after AdjustTree"
      stop
    end if

!   (I4)[grow tree taller]: if node split propagation caused
!   the root to split, create a new root whose children are
!   the two resulting nodes
    if ((.not.associated(L%parent)).and.(associated(LL))) then

      call CreateNode(T)

      T%entries = 2

!     create an entry (node pointer and interval) for L
      allocate(T%E(1)%pe,stat=istat)
      nalloc = nalloc + 1
      nentry = nentry + 1
      if (istat /= 0) stop
      T%E(1)%pe%pn => L
      nullify(T%E(1)%pe%pp)
      do idir = 1,ndir
        T%E(1)%pe%I(idir,1) = L%E(1)%pe%I(idir,1)
        T%E(1)%pe%I(idir,2) = L%E(1)%pe%I(idir,2)
        do i = 2,L%entries
          T%E(1)%pe%I(idir,1) = min(T%E(1)%pe%I(idir,1), &
                                    L%E(i)%pe%I(idir,1))
          T%E(1)%pe%I(idir,2) = max(T%E(1)%pe%I(idir,2), &
                                    L%E(i)%pe%I(idir,2))
        end do
      end do


!     create an entry (node pointer and interval) for LL
      allocate(T%E(2)%pe,stat=istat)
      nalloc = nalloc + 1
      nentry = nentry + 1
      if (istat /= 0) stop
      T%E(2)%pe%pn => LL
      nullify(T%E(2)%pe%pp)
      do idir = 1,ndir
        T%E(2)%pe%I(idir,1) = LL%E(1)%pe%I(idir,1)
        T%E(2)%pe%I(idir,2) = LL%E(1)%pe%I(idir,2)
        do i = 2,LL%entries
          T%E(2)%pe%I(idir,1) = min(T%E(2)%pe%I(idir,1), &
                                   LL%E(i)%pe%I(idir,1))
          T%E(2)%pe%I(idir,2) = max(T%E(2)%pe%I(idir,2), &
                                   LL%E(i)%pe%I(idir,2))
        end do
      end do

      T%leaf = .false.
      nullify(T%parent)
      L%parent => T
      LL%parent => T

      if (debug) then
        write(unit=*,fmt=*)"Insert: the root node (node ", &
                     L%node_number,") has been split."
        write(unit=*,fmt=*) &
         "Insert: a new root T has been allocated (node",T%node_number,")"
        write(unit=*,fmt=*)"Insert: T now has ",T%entries," entries."
        write(unit=*,fmt=*)"Insert: L (node ", &
                     L%node_number,") is in the first entry of T"
        write(unit=*,fmt=*)"Insert: LL (node ", &
                    LL%node_number,") is in the second entry of T"
        write(unit=*,fmt=*)"Insert: T is not a leaf."
        write(unit=*,fmt=*)"Insert: T has no parent."
        write(unit=*,fmt=*)"Insert: T (node ", &
        T%node_number,") is the parent of L (node",L%node_number,")"
        write(unit=*,fmt=*)"Insert: T (node ", &
        T%node_number,") is the parent of LL (node",LL%node_number,")"
      end if

    else
      if (debug) write(unit=*,fmt=*)"Insert: the root has not been split."
    end if

    if (debug) write(unit=*,fmt=*)"Insert: complete. return to main code."

    return
    end subroutine Insert

!   ==========================================================
    recursive subroutine AdjustTree(L,LL)
!   ==========================================================
!   ascend from a leaf node L to the root, adjusting covering
!   rectangles and propagating node splits as necessary

    use NetAlloc  ,only: nalloc,nentry
    use RTreeTypes,only: node,entry,debug,ndir,MaxEntries

    implicit none

    type (node), pointer  :: L,LL,N,NN,P,PP
    type (entry), pointer :: ENN,EN
    integer               :: idir,i,istat,ien


    continue


    nullify(N,NN,P,PP,EN,ENN)


!   (AT1)[initialize]: set N=L.  if L was split previously,
!   set NN to be the resulting second node.
    N => L
    if (associated(LL)) then
      if (debug) write(unit=*,fmt=*) &
         "adjusttree: LL is associated; point NN to LL."
      NN => LL
    else
      if (debug) write(unit=*,fmt=*)"Adjusttree: LL is not associated."
    end if


!   (AT2)[check if done]: if N is the root, stop.
    if (.not.associated(N%parent)) then
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: N is the root; nothing to adjust"
      return
    else
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: N is not the root; adjust the covering rectangle."
    end if


!   (AT3)[adjust covering rectangle in parent entry]:
!   let P be the parent node of N
    P => N%parent
 
    if (debug) write(unit=*,fmt=*) "AdjustTree: I (",N%node_number, &
        ") am one of the ",P%entries," entries of node ",P%node_number

    i = 0
    do ien = 1,P%entries
      if (P%E(ien)%pe%pn%node_number == N%node_number) then
        i = ien
        exit
      end if
    end do
    if (i == 0) stop

    EN => P%E(i)%pe
    if (debug) then
      write(unit=*,fmt=*)"AdjustTree: my entry number is ",i
      write(unit=*,fmt=*)"AdjustTree: EN is now pointing to P%E(i)%pe"
    end if


!   adjust EN%I so that it tightly encloses all entry rectangles in N.
    if (debug) write(unit=*,fmt=*)"AdjustTree: adjust EN%I"
    do idir = 1,ndir
      EN%I(idir,1) = N%E(1)%pe%I(idir,1)
      EN%I(idir,2) = N%E(1)%pe%I(idir,2)
      do i = 2,N%entries
        EN%I(idir,1) = min(EN%I(idir,1),N%E(i)%pe%I(idir,1))
        EN%I(idir,2) = max(EN%I(idir,2),N%E(i)%pe%I(idir,2))
      end do
    end do
    if (debug) then
      write(unit=*,fmt=*)"AdjustTree: finished adjusting EN%I"
      write(unit=*,fmt=*)"AdjustTree: EN%I(1,1:2) = ",EN%I(1,1),EN%I(1,2)
    end if


!   (AT4)[propagate node split upward]: if N has a partner NN
!   resulting from an earlier split,
    if (associated(NN)) then
      if (debug) write(unit=*,fmt=*) &
         "adjusttree: NN is associated; create a new entry ENN"
!     create a new entry ENN 
      allocate(ENN,stat=istat)
      nalloc = nalloc + 1
      nentry = nentry + 1
      if (istat /= 0) stop
!     with ENN%P pointing to NN 
      ENN%PN => NN
      nullify(ENN%PP)
!     and ENN%I enclosing all rectangles in NN.
      do idir = 1,ndir
        ENN%I(idir,1) = NN%E(1)%pe%I(idir,1)
        ENN%I(idir,2) = NN%E(1)%pe%I(idir,2)
        do i = 2,NN%entries
          ENN%I(idir,1) = min(ENN%I(idir,1),NN%E(i)%pe%I(idir,1))
          ENN%I(idir,2) = max(ENN%I(idir,2),NN%E(i)%pe%I(idir,2))
        end do
      end do
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: ENN%I(1,1:2) = ",ENN%I(1,1),ENN%I(1,2)
         
!     add ENN to P if there is room.
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: add ENN to P if there is room."
      if (P%entries < MaxEntries) then
        if (debug) write(unit=*,fmt=*)"AdjustTree: P has room for ENN."
        P%entries = P%entries + 1
        if (debug) write(unit=*,fmt=*)"AdjustTree: P (node ", &
           P%node_number,") now has ",P%entries," entries."
        P%E(P%entries)%pe => ENN
      else
!       otherwise, invoke SplitNode to produce P and PP containing
!       ENN and all P's old entries
        if (debug) write(unit=*,fmt=*)"AdjustTree: P (node ", &
           P%node_number,") does not have room for ENN."
        if (debug) write(unit=*,fmt=*)"AdjustTree: calling LinearSplitNode"
        call LinearSplitNode(P,PP,ENN)
      end if
    else
      if (debug) write(unit=*,fmt=*)"AdjustTree: NN is not associated"
    end if


!   (AT5)[move up to next level.] set N=P 
    if (debug) write(unit=*,fmt=*)"AdjustTree: move up to next level; L => P"
    L => P
!   and set NN=PP if a split occurred.  
    if (associated(PP)) then
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: PP is associated; LL => PP"
      LL => PP
    else
      if (debug) write(unit=*,fmt=*) &
         "AdjustTree: PP is not associated; LL => NULL()"
      nullify(LL)
    end if

!   repeat from (AT2)
    if (debug) write(unit=*,fmt=*)"AdjustTree: calling AdjustTree"
    call AdjustTree(L,LL)


    return
    end subroutine AdjustTree

!   ==========================================================
    recursive subroutine LinearSplitNode(L,LL,E)
!   ==========================================================
!   purpose: divide a set of m+1 index entries into two groups

    use IntrType  ,only: rd
    use NetAlloc  ,only: nalloc,nnodes,npnten
    use RTreeTypes,only: node,entry,debug,MaxEntries,MinEntries

    implicit none

    type (node), pointer  :: L,LL,TEMP
    type (entry), pointer :: E
    real(kind=rd)         :: dL, dLL,areaL,areaLL
    integer               :: i,i1,i2,istat


    continue


    nullify(TEMP)


!   before getting started, define a temporary node that
!   contains all m+1 index entries from the original two nodes

    allocate(TEMP,stat=istat)
    nalloc = nalloc + 1
    nnodes = nnodes + 1

    if (istat /= 0) stop

    allocate(TEMP%E(MaxEntries + 1),stat=istat)
    nalloc = nalloc + 1
    npnten = npnten + 1

    if (istat /= 0) stop

    TEMP%entries = L%entries + 1
    TEMP%parent  => L%parent
    TEMP%leaf  = L%leaf

    if (debug) then
      write(unit=*,fmt=*)"LinearSplitNode: created temp node with ", &
                          TEMP%entries," entries"
      if (associated(TEMP%parent)) then
        write(unit=*,fmt=*)"LinearSplitNode: temp parent node = ", &
             TEMP%parent%node_number
      else
        write(unit=*,fmt=*)"LinearSplitNode: temp is at root level."
      end if
      write(unit=*,fmt=*)"LinearSplitNode: temp leaf = ",TEMP%leaf
    end if

    if (TEMP%entries /= MaxEntries + 1) then
      write(unit=0,fmt=*)"ERROR! miscount in LinearSplitNode"
      write(unit=0,fmt=*)"TEMP%entries = ",TEMP%entries
      write(unit=0,fmt=*)"MaxEntries = ",MaxEntries
      write(unit=0,fmt=*)"L%entries = ",L%entries
      stop
    end if


!   define temp node
    do i = 1,L%entries
      TEMP%E(i)%pe => L%E(i)%pe
    end do
    TEMP%E(TEMP%entries)%pe => E

    do i = 1,MaxEntries
      nullify(L%E(i)%pe)
    end do

      
!   (S1)[pick first entry for each group]:  apply algorithm
!   PickSeeds to choose two entries to be the first
!   elements of the groups.  

!   The initial logic in LinearPickSeeds was problematic and allowed
!   i1 = i2 if one box was inside another.
!   if (debug) write(unit=*,fmt=*) &
!       "LinearSplitNode: calling LinearPickSeeds"
!   call LinearPickSeeds(TEMP,i1,i2)
!   if (debug) write(unit=*,fmt=*) &
!       "LinearSplitNode: LinearPickSeeds has chosen entries ",i1,i2

    if (debug) write(unit=*,fmt=*) &
        "Splitnode: calling QuadraticPickSeeds"
    call QuadraticPickSeeds(TEMP,i1,i2)
    if (debug) write(unit=*,fmt=*) &
       "SplitNode: QuadraticPickSeeds has chosen entries ",i1,i2

    if ((max(i1,i2) > TEMP%entries).or.(min(i1,i2) < 1)) then
      write(unit=0,fmt=*) &
          "ERROR! entries out of range, TEMP%entries = ",TEMP%entries
      stop
    end if
    if (i1 == i2) then
      write(unit=0,fmt=*) &
         "ERROR! PickSeeds has chosen the same seed twice!"
      stop
    end if

!   assign each to a group.
    if (debug) write(unit=*,fmt=*) &
       "SplitNode: assign first entry to L (node ",L%node_number,")"
    L%entries = 1
    L%E(1)%pe => TEMP%E(i1)%pe
    L%parent => TEMP%parent
    L%leaf = TEMP%leaf
    if (.not.L%leaf) L%E(1)%pe%pn%parent => L

!   create a new node LL whose parent and leaf status is the same as L
    call CreateNode(LL)

    if (debug) write(unit=*,fmt=*) &
        "LinearSplitNode: assign second entry to LL (new node ", &
        LL%node_number,")"

    LL%entries = 1
    LL%E(1)%pe => TEMP%E(i2)%pe
    LL%parent => TEMP%parent
    LL%leaf = TEMP%leaf
    if (.not.LL%leaf) LL%E(1)%pe%pn%parent => LL

!   remove the first two entries from the temp node
    if (debug) write(unit=*,fmt=*)"LinearSplitNode: remove entry ", &
       max(i1,i2)," from temp"
    call RemoveEntry(TEMP, max(i1,i2))
    if (debug) write(unit=*,fmt=*)"LinearSplitNode: remove entry ", &
       min(i1,i2)," from temp"
    call RemoveEntry(TEMP, min(i1,i2))

    do

!     (S2)[check if done]: if all entries have been assigned, stop.
      if (TEMP%entries == 0) then
        deallocate(TEMP%E)
        nalloc = nalloc - 1
        npnten = npnten - 1

        deallocate(TEMP)
        nalloc = nalloc - 1
        nnodes = nnodes - 1

        if (debug) then
          write(unit=*,fmt=*)"LinearSplitNode: TEMP is out of entries; return"
          write(unit=*,fmt=*)"LinearSplitNode: L (node ", &
             L%node_number,") now has ",L%entries," entries"
          write(unit=*,fmt=*)"LinearSplitNode: LL (node ", &
             LL%node_number,") now has ",LL%entries," entries"
        end if

        return
      end if

!     if one group has so few entries that all the rest must be 
!     assigned to it in order for it to have the minimum number m, 
!     assign them and stop.
      if (TEMP%entries + L%entries == MinEntries) then
!       assign the rest of the entries to L
        if (debug) write(unit=*,fmt=*) &
           "LinearSplitNode: assign the rest of the entries to L"
        do i = 1,TEMP%entries
          L%E(i + L%entries)%pe => TEMP%E(i)%pe
          if (.not.L%leaf) L%E(i + L%entries)%pe%pn%parent => L
        end do
        L%entries = L%entries + TEMP%entries

        deallocate(TEMP%E)
        nalloc = nalloc - 1
        npnten = npnten - 1

        deallocate(TEMP)
        nalloc = nalloc - 1
        nnodes = nnodes - 1

        return
      end if
 
      if (TEMP%entries + LL%entries == minentries) then
!       assign the rest of the entries to LL
        if (debug) write(unit=*,fmt=*) &
           "LinearSplitNode: assign the rest of the entries to LL"
        do i = 1,TEMP%entries
          LL%E(i + LL%entries)%pe => TEMP%E(i)%pe
          if (.not.LL%leaf) LL%E(i + LL%entries)%pe%pn%parent => LL
        end do
        LL%entries = LL%entries + TEMP%entries

        deallocate(TEMP%E)
        nalloc = nalloc - 1
        npnten = npnten - 1

        deallocate(TEMP)
        nalloc = nalloc - 1
        nnodes = nnodes - 1

        return
      end if

!     (S3)[select entry to assign]: invoke algorithm PickNext
!     to choose the next entry to assign.  
      if (debug) write(unit=*,fmt=*)"LinearSplitNode: calling PickNext"
      call PickNext(TEMP,L,LL,i,dL,dLL,areaL,areaLL)
      if (debug) write(unit=*,fmt=*) &
         "LinearSplitNode: PickNext selected entry",i
      if ((i > TEMP%entries).or.(i < 1)) then
        write(unit=0,fmt=*)"ERROR! PickNext chose invalid entry ",i
        write(unit=0,fmt=*)"when TEMP%entries = ",TEMP%entries
        stop
      end if

!     add it to the group whose covering rectangle will have 
!     to be enlarged least to accommodate it.  
      if (dL < dLL) then
        if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to L"
        L%entries = L%entries + 1
        L%E(L%entries)%pe => TEMP%E(i)%pe
        if (.not.L%leaf) L%E(L%entries)%pe%pn%parent => L
      else if (dLL < dL) then
        if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to LL"
        LL%entries = LL%entries + 1
        LL%E(LL%entries)%pe => TEMP%E(i)%pe
        if (.not.LL%leaf) LL%E(LL%entries)%pe%pn%parent => LL
      else
!       resolve ties by adding the entry to
!       the group with smaller area
        if (areal < areall) then
          if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to l"
          L%entries = L%entries + 1
          L%E(L%entries)%pe => TEMP%E(i)%pe
          if (.not.L%leaf) L%E(L%entries)%pe%pn%parent => L
        else if (areall < areal) then
          if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to LL"
          LL%entries = LL%entries + 1
          LL%E(LL%entries)%pe => TEMP%E(i)%pe
          if (.not.LL%leaf) LL%E(LL%entries)%pe%pn%parent => LL
        else
!         then to the one with fewer entries, 
          if (L%entries < LL%entries) then
            if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to L"
            L%entries = L%entries + 1
            L%E(L%entries)%pe => TEMP%E(i)%pe
            if (.not.L%leaf) L%E(L%entries)%pe%pn%parent => L
          else if (LL%entries < L%entries) then
            if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to LL"
            LL%entries = LL%entries + 1
            LL%E(LL%entries)%pe => TEMP%E(i)%pe
            if (.not.LL%leaf) LL%E(LL%entries)%pe%pn%parent => LL
          else
!           then to either.
            if (debug) write(unit=*,fmt=*)"LinearSplitNode: add to L"
            L%entries = L%entries + 1
            L%E(L%entries)%pe => TEMP%E(i)%pe
            if (.not.L%leaf) L%E(L%entries)%pe%pn%parent => L
          end if
        end if
      end if
      if (debug) write(unit=*,fmt=*)"LinearSplitNode: remove entry ",i, &
                " from temp"
      call RemoveEntry(temp, i)

!     repeat from (s2)
    end do

    return
    end subroutine LinearSplitNode 

!   ================================================
    subroutine RemoveEntry(n, i1)
!   ================================================
!   purpose: remove one entry from a node

    use RTreeTypes,only: node
 
    implicit none

    type(node),pointer :: N
    integer,intent(in) :: i1
    integer            :: i


    continue


    do i = i1,N%entries-1
      N%E(i)%pe => N%E(i+1)%pe
    end do

    nullify(N%E(N%entries)%pe)
    N%entries = N%entries - 1

    return
    end subroutine RemoveEntry

!   ====================================================
    subroutine QuadraticPickSeeds(TEMP,i1,i2)
!   ====================================================
!   purpose: select two entries to be the first elements
!   of the groups

    use IntrType  ,only: rd
    use RTreeTypes,only: node,ndir

    implicit none

    type (node), pointer :: temp
    integer,intent(out) :: i1,i2
    integer :: i,j,idir
    real(kind=rd) :: xmin,xmax,areaj,area1,area2,d,dmax


    continue


    do i = 1,TEMP%entries-1
    do j = i+1,TEMP%entries

!     (PS1)[calculate inefficiency of grouping entries together]
!     for each pair of entries E1 and E2, compose a rectangle j
!     including E1%I and E2%I.
      areaj = 1.0
      area1 = 1.0
      area2 = 1.0
      do idir = 1,ndir
        xmin = min(TEMP%E(i)%pe%I(idir,1),TEMP%E(j)%pe%I(idir,1))
        xmax = max(TEMP%E(i)%pe%I(idir,2),TEMP%E(j)%pe%I(idir,2))
        areaj = areaj * (xmax-xmin)
        area1 = area1 * (TEMP%E(i)%pe%I(idir,2)-TEMP%E(i)%pe%I(idir,1))
        area2 = area2 * (TEMP%E(j)%pe%I(idir,2)-TEMP%E(j)%pe%I(idir,1))
      end do

!     (PS2)[choose the most wasteful pair]: choose the pair with
!     the largest d
      d = areaj - area1 - area2
      if ((i == 1).and.(j == 2)) then
        dmax = d
        i1 = i
        i2 = j
      else
        if (d > dmax) then
          dmax = d
          i1 = i
          i2 = j
        end if
      end if
    end do
    end do

    return
    end subroutine QuadraticPickSeeds

!   ====================================================
    subroutine LinearPickSeeds(TEMP,i1,i2)
!   ====================================================
!   purpose: select two entries to be the first elements
!   of the groups

    use IntrType  ,only: rd
    use RTreeTypes,only: node,ndir,debug

    implicit none

    type (node), pointer :: temp

    real(kind=rd),dimension(ndir) :: width
    real(kind=rd) :: xmin,xmax,maxsep,tmpsep,highest_low,lowest_high

    integer,intent(out) :: i1,i2
    integer :: idir,i,i1_tmp,i2_tmp


    continue


!   (LPS1)[find extreme rectangles along all dimensions]:
!   along each dimension, find the entry whose rectangle
!   has the highest low side, and the one with the lowest
!   high side.  record the separation.
!   DAB -----> defer this to (LPS3)


!   (LPS2)[adjust for shape of the rectangle cluster]:
!   normalize the spearations by dividing the width of the
!   entire set along the corresponding dimension.
    do idir = 1,ndir
!     initialize the extremes using the new entry
      xmin = TEMP%E(1)%pe%I(idir,1)
      xmax = TEMP%E(1)%pe%I(idir,2)
!     update the extremes using all entries in the current node
      do i = 2,TEMP%entries
        xmin = min(xmin,TEMP%E(i)%pe%I(idir,1))
        xmax = max(xmax,TEMP%E(i)%pe%I(idir,2))
      end do
      width(idir) = xmax - xmin
      if (debug) write(unit=*,fmt=*) &
         "LinearPickSeeds: width(",idir,") = ",width(idir)
    end do


!   (LPS3)[select the most extreme pair]:  choose the pair
!   with the greatest normalized separation along any dimension.
    do idir = 1,ndir

!     find the entry whose rectangle has the highest low side
      i1_tmp = 1
      highest_low = TEMP%E(1)%pe%I(idir,1)
      do i = 2,TEMP%entries
        if (TEMP%E(i)%pe%I(idir,1) > highest_low) then
          i1_tmp = i
          highest_low = TEMP%E(i)%pe%I(idir,1)
          if (debug) write(unit=*,fmt=*) &
             "LinearPickSeeds: reset i1_tmp to ",i1_tmp
        end if
      end do

!     find the entry whose rectangle has the lowest high side
      if (i1_tmp /= 1) then
        i2_tmp = 1
      else
        i2_tmp = 2
      end if
      lowest_high = TEMP%E(i2_tmp)%pe%I(idir,2)
      do i = 2,TEMP%entries
        if (TEMP%E(i)%pe%I(idir,2) < lowest_high) then
          if (i /= i1_tmp) then
            i2_tmp = i
            lowest_high = TEMP%E(i)%pe%I(idir,2)
            if (debug) write(unit=*,fmt=*) &
               "LinearPickSeeds: reset i2_tmp to ",i2_tmp
          end if
        end if
      end do
         
!     record the normalized separation    
      tmpsep = (highest_low - lowest_high) / width(idir)

      if (idir == 1) then
        maxsep = tmpsep
        i1 = i1_tmp
        i2 = i2_tmp
      else
!       choose the pair with the greatest normalized separation
        if (tmpsep > maxsep) then
          i1 = i1_tmp
          i2 = i2_tmp
          if (debug) write(unit=*,fmt=*) &
             "LinearPickSeeds: reset i1,i2 to ",i1,i2
        end if
      end if

    end do


    return
    end subroutine LinearPickSeeds

!   ===========================================================
    recursive subroutine ChooseLeaf(T,E,N)
!   ===========================================================
!   purpose: select a leaf node in which to place a new entry E

    use IntrType  ,only: rd
    use RTreeTypes,only: node,entry,debug,ndir

    implicit none

    type (node), pointer          :: T,N,F
    type (entry), pointer         :: E
    real(kind=rd),dimension(ndir) :: xmin1,xmax1,xmin2,xmax2
    real(kind=rd)                 :: d,dmin
    real(kind=rd)                 :: area1,area2
    integer                       :: i,imin,idir


    continue


    nullify(F)

!   (CL1)[initialize]: set N to be the root node
    if (debug) write(unit=*,fmt=*) &
       "ChooseLeaf: point N to T (node ",T%node_number,")"
    N => T

!   (CL2)[leaf check]: if N is a leaf, return N
    if (N%leaf) then
      if (debug) write(unit=*,fmt=*) &
         "ChooseLeaf: N (node ",N%node_number,") is a leaf; return"
      return
    end if

!   (CL3)[choose subtree]: If N is not a leaf, let F be the
!   entry in N whose rectangle F%I needs least enlargement to
!   include E%pe%I.  Resolve ties by choosing the entry with the
!   rectangle of smallest area.
    if (debug) write(unit=*,fmt=*) &
       "ChooseLeaf: N (node ",N%node_number,") is not a leaf."

    dmin = 1.0e+20
    imin = 0
    do i = 1,N%entries
!     compute the original area of each entry in N
      area1 = 1.0
      area2 = 1.0
      do idir = 1,ndir
        xmin1(idir) = N%E(i)%pe%I(idir,1)
        xmax1(idir) = N%E(i)%pe%I(idir,2)
        area1 = area1 * (xmax1(idir)-xmin1(idir))
        xmin2(idir) = min(xmin1(idir),E%I(idir,1))
        xmax2(idir) = max(xmax1(idir),E%I(idir,2))
        area2 = area2 * (xmax2(idir)-xmin2(idir))
      end do
      d = area2 - area1
      if (d < dmin) then
        dmin = d
        imin = i
      end if
    end do
    if (imin == 0) then
      write(unit=0,fmt=*)"ERROR! no min area found." 
      stop
    end if

    F => N%E(imin)%pe%pn

    if (debug) then
      write(unit=*,fmt=*) &
         "ChooseLeaf: rectangle ",imin," needs least adjustment"
      write(unit=*,fmt=*) &
         "ChooseLeaf: go to F (node ",F%node_number,")"
    end if

!   (CL4)[descend util a leaf is reached.] set N to be the child node
!   pointed to by Fp and repeat from CL2
    if (debug) write(unit=*,fmt=*)"ChooseLeaf: calling ChooseLeaf"
    call ChooseLeaf(F,E,N)

    return
    end subroutine ChooseLeaf

!   =================================================================
    subroutine PickNext(TEMP,N1,N2,imax,d1i,d2i,area1a,area2a)
!   =================================================================
!   purpose: select one remaining entry for classification in a group

    use IntrType  ,only: rd
    use RTreeTypes,only: ndir,node

    implicit none

    type(node),pointer            :: TEMP,N1,N2

    real(kind=rd),intent(out)     :: d1i,d2i,area1a,area2a

    real(kind=rd),dimension(ndir) :: xmin1,xmax1,xmin2,xmax2
    real(kind=rd),dimension(ndir) :: xmin,xmax

    real(kind=rd)                 :: area1b,area2b
    real(kind=rd)                 :: d1,d2,dmax

    integer,intent(out)           :: imax
    integer                       :: i,idir


    continue


!   calculate the initial area of the covering rectangle of group 1
    area1a = 1.0
    do idir = 1,ndir
      xmin1(idir) = N1%E(1)%pe%I(idir,1)
      xmax1(idir) = N1%E(1)%pe%I(idir,2)
      do i = 2,N1%entries
        xmin1(idir) = min(xmin1(idir),N1%E(i)%pe%I(idir,1))
        xmax1(idir) = max(xmax1(idir),N1%E(i)%pe%I(idir,2))
      end do
      area1a = area1a * (xmax1(idir)-xmin1(idir))
    end do


!   calculate the initial area of the covering rectangle of group 2
    area2a = 1.0
    do idir = 1,ndir
      xmin2(idir) = N2%E(1)%pe%I(idir,1)
      xmax2(idir) = N2%E(1)%pe%I(idir,2)
      do i = 2,N2%entries
        xmin2(idir) = min(xmin2(idir),N2%E(i)%pe%I(idir,1))
        xmax2(idir) = max(xmax2(idir),N2%E(i)%pe%I(idir,2))
      end do
      area2a = area2a * (xmax2(idir)-xmin2(idir))
    end do


!   (PN1)[determine cost of putting each entry in each group.]

    do i = 1,TEMP%entries

!     for each entry E not yet in a group, calculate d1 = the
!     area increase required in the covering rectangle of group 1
!     to include EI.
      area1b = 1.0
      do idir = 1,ndir
        xmin(idir) = min(xmin1(idir),TEMP%E(i)%pe%I(idir,1))
        xmax(idir) = max(xmax1(idir),TEMP%E(i)%pe%I(idir,2))
        area1b = area1b * (xmax(idir)-xmin(idir))
      end do
      d1 = area1b - area1a

!     calculate d2 similarly for group 2.
      area2b = 1.0
      do idir = 1,ndir
        xmin(idir) = min(xmin2(idir),TEMP%E(i)%pe%I(idir,1))
        xmax(idir) = max(xmax2(idir),TEMP%E(i)%pe%I(idir,2))
        area2b = area2b * (xmax(idir)-xmin(idir))
      end do
      d2 = area2b - area2a
         
!     (PN2)[find entry with greatest preference for one group.]
!     choose any entry with the maximum difference between d1 and d2.
      if ((i == 1).or.(abs(d2-d1) > dmax)) then
        dmax = abs(d2-d1)
        imax = i
        d1i = d1
        d2i = d2
      end if

    end do

    return
    end subroutine PickNext

!   ===============================================
    subroutine CreateNode(NewNode)
!   ===============================================
!   procedure to create a new node of the R-Tree

    use NetAlloc  ,only: nalloc,nnodes,npnten
    use RTreeTypes,only: node,MaxEntries,total_nodes

    implicit none

    type(node),pointer :: NewNode
    integer :: i,istat


    continue


    allocate(NewNode,stat=istat)
    nalloc = nalloc + 1
    nnodes = nnodes + 1
    if (istat /= 0) stop

    allocate(NewNode%E(MaxEntries),stat=istat)
    nalloc = nalloc + 1
    npnten = npnten + 1

    if (istat /= 0) stop

    do i = 1,MaxEntries
      nullify(NewNode%E(i)%pe)
    end do

    total_nodes = total_nodes + 1
    NewNode%node_number = total_nodes

    nullify(NewNode%parent)
    NewNode%entries = 0

    return
    end subroutine CreateNode

!   ===============================================
    recursive subroutine DeleteTree(TopNode)
!   ===============================================
!   Recursive procedure to delete an R-Tree

    use NetAlloc  ,only: nalloc,nentry,npnten,nnodes
    use RTreeTypes,only: node,MaxEntries

    implicit none

    type(node),pointer :: TopNode,P
    integer :: i


    continue


    nullify(P)


    if (TopNode%leaf) then

!     entries contain pointers to panels -- don't delete the panels
      do i = 1,MaxEntries
        if (associated(TopNode%E(i)%pe)) then
          nullify(TopNode%E(i)%pe%pp)
          nullify(TopNode%E(i)%pe%pn)
          deallocate(TopNode%E(i)%pe)
          nalloc = nalloc - 1
          nentry = nentry - 1
        end if
      end do
    else

!     entries contain pointers to nodes that must be deleted
      do i = 1,MaxEntries
        if (associated(TopNode%E(i)%pe)) then
          nullify(TopNode%E(i)%pe%pp)
          if (associated(TopNode%E(i)%pe%pn)) &
             call DeleteTree(TopNode%E(i)%pe%pn)
          deallocate(TopNode%E(i)%pe)
          nalloc = nalloc - 1
          nentry = nentry - 1
        end if
      end do
    end if

    deallocate(TopNode%E)
    nalloc = nalloc - 1
    npnten = npnten - 1

    TopNode%entries = 0

    if (associated(TopNode%parent)) then
      P => TopNode%parent
      deallocate(TopNode)
      nalloc = nalloc - 1
      nnodes = nnodes - 1
      TopNode => P
      nullify(P)
    end if

    return
    end subroutine DeleteTree

!   =================================================================
    subroutine SearchTree(p,root)
!   =================================================================
!   construct a table counting the number of panel bounding boxes
!   in each surface that overlap any panel bounding boxes in another

    use NetAlloc  ,only: nalloc
    use PatchInfo ,only: num_patches,table
    use RTreeTypes,only: node
    use Types     ,only: panel

    implicit none

    type (node) , pointer               :: root
    type (panel), pointer, dimension(:) :: p
    type (panel), pointer               :: q

    integer :: ipan,ip1,ip2,ip,iq,ierr

    logical,dimension(num_patches) :: list


    continue


    ip1 = lbound(p,1)
    ip2 = ubound(p,1)

    allocate(table(num_patches,num_patches),stat=ierr)
    nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR! allocation of table failed."
      write(0,*)"num_patches = ",num_patches
      stop
    end if

    table = 0

    do ipan = ip1,ip2

      list =  .false.
      q    => p(ipan)
      iq   =  q%isurf

      call FindOverlap( q,root,list,num_patches )

      do ip = 1,num_patches
        if (list(ip)) table(iq,ip) = table(iq,ip) + 1
      end do

    end do

    return
    end subroutine SearchTree

!   =========================================================
    recursive subroutine FindOverlap(q,T,list,num_patches)
!   =========================================================
!   count the number of panels in a given patch that
!   overlap panels in a second patch

    use OverlappedM,only: Overlapped
    use RTreeTypes,only: node
    use Types     ,only: panel

    implicit none

    type (node) , pointer :: T
    type (panel), pointer :: q

    integer,intent(in) :: num_patches
    integer            :: i,ipq,ipt

    logical,dimension(num_patches) :: list
!   logical                        :: overlapped


    continue


    if (.not.T%leaf) then

!     if T is not a leaf, then see if q overlaps 
!     any of the entries in T

      do i = 1,T%entries
        call FindOverlap(q,T%E(i)%pe%pn,list,num_patches)
      end do

    else

!     T is a leaf. See if q overlaps any of the panels in T

      do i = 1,T%entries
        if (overlapped(q%bound,T%E(i)%pe%i)) then
          ipq = q%isurf             ! patch number for panel q
          ipt = T%E(i)%pe%pp%isurf  ! patch number for panel T%E(i)%pe%pp
          if (ipq /= ipt) then
            list(ipt) = .true.
          end if
        end if
      end do

    end if

    return 
    end subroutine FindOverlap

  end module RTreeProcedures
