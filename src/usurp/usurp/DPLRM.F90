! ====================================================
  module DPLRM
! ====================================================

  implicit none
  save

  integer,allocatable,dimension(:,:),private :: jbc
  integer,allocatable,dimension(:),private :: mtx,mty,mtz

  integer,private :: nblk,igrid

  character(len=200),private :: gname

  private :: rdagrd,rdwrthead,read_dplr_input
#ifdef FXDR_HOME
  private :: swaporder
#endif
  public  :: PreReadDPLR,StorePatchesDPLR

 
  contains

! ====================================================
  subroutine PreReadDPLR()
! ====================================================
! generate USURP data for DPLR

  use NetAlloc                 ,only: nalloc
  use PatchInfo                ,only: i1,i2,j1,j2,k1,k2,block
  use PatchInfo                ,only: ncomp_surfs,icomp,num_patches
  use PatchInfo                ,only: max_components,tag_list
  use StructuredPatchProcedures,only: AllocateStructuredPatch
  use StructuredPatchProcedures,only: CountStructuredPatch
  use StructuredPatchProcedures,only: EchoPatchData


  implicit none

  integer :: n,iface,nig,njg,nkg,ip,ibc,ierr


  continue


! initialization
  allocate(tag_list(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate tag_list in PreReadUNCLEM"
    write(0,*)"max_components = ",max_components
    stop
  end if
  tag_list(0)    = "untagged"

  allocate(ncomp_surfs(0:max_components),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! Failed to allocate ncomp_surfs in PreReadBCin"
    write(0,*)"max_components = ",max_components
    stop
  end if
  ncomp_surfs = 0


! read the DPLR input file to ascertain wall boundaries
  write(*,*)
  write(*,*)"Waiting for DPLR input data to be fed from standard input ..."
  call read_dplr_input()

! count wall patches
  do n = 1,nblk
    do iface = 1,6
      ibc = jbc(iface,n)     

      select case(ibc)
!       pick out any type of wall boundary condition
        case(8:10,25:27,30,40,49,70:71,110,125,130,140)
          num_patches = num_patches + 1
      end select

    end do
  end do


  call AllocateStructuredPatch()


! store wall patches
  ip = 0
  do n = 1,nblk
    nig = mtx(n) - 1
    njg = mty(n) - 1
    nkg = mtz(n) - 1
    do iface = 1,6
      ibc = jbc(iface,n)     

      select case(ibc)
!       pick out any type of wall boundary condition
        case(8:10,25:27,30,40,49,70:71,110,125,130,140)
          ip = ip + 1
          block(ip) = n
          i1(ip) = 1
          j1(ip) = 1
          k1(ip) = 1
          i2(ip) = nig
          j2(ip) = njg
          k2(ip) = nkg
          select case(iface)
            case(1)
              i2(ip) = 1
            case(2)
              i1(ip) = nig
            case(3)
              j2(ip) = 1
            case(4)
              j1(ip) = njg
            case(5)
              k2(ip) = 1
            case(6)
              k1(ip) = nkg
          end select

          call CountStructuredPatch(ip)

          icomp(ip) = 0
          ncomp_surfs(0) = ncomp_surfs(0) + 1
      end select
    end do
  end do

! screen output to confirm results
  call EchoPatchData()

  deallocate(mtz,mty,mtx); nalloc = nalloc - 3
  deallocate(jbc); nalloc = nalloc - 1

  return
  end subroutine PreReadDPLR

! ====================================================
  subroutine read_dplr_input()
! ====================================================
! Read the CFD input deck for the DPLR code.  Derived from the
! DPLR code, written by  Michael Wright (mjwright@mail.arc.nasa.gov)

  use NetAlloc,only: nalloc

  implicit none

  integer :: igp,inm
  integer :: ntx,nty,ntz,iconr,isim,ifree,initi
  integer :: i,ierr,ios,m,n

  character(len=200) :: errmsg
  character(len=  5) :: psuf


  continue


! initialize
  ios   = 0


! skip down and read the grid filename
  errmsg = "problem reading input filenames"
  do i = 1,3
    read(*,*,err=501,iostat=ios)
  end do
  read(*,*,err=501,iostat=ios) gname ! grid file


! skip the rest of the input filenames
  do i = 1,4
    read(*,*,err=501,iostat=ios)
  end do
  501 continue
  if (ios > 0) goto 312

 
! read the nblk line
  errmsg = "problem reading nblk line"
  read(*,*,err=502,iostat=ios)
  read(*,*,err=502,iostat=ios)
  read(*,*,err=502,iostat=ios) nblk,igrid
  502 continue
  if (ios > 0) goto 312


! assign grid file suffix
  select case(igrid)
    case(1)
!     native unformatted parallel grid file
      psuf = ".pgrd"
    case(11)
#ifdef USE_FXDR
!     XDR parallel grid file
      psuf = ".pgrx"
#else
      write(0,*)
      write(0,*)"ERROR! USURP is not instrumented to use XDR format"
      write(0,*)"because the FXDR_HOME environment variable"
      write(0,*)"has not been set to point to the FXDR library."
      write(0,*)"See the USURP User's Manual for more information."
      stop
#endif
    case(21)
!     ASCII parallel grid file
      psuf = ".pgra"
    case default
      write(0,111) "igrid",igrid
      write(0,*)"ERROR! Invalid grid file type selected"
      stop
  end select

! attach the suffix if it is needed
  igp = len_trim(gname)
  inm = max(1,igp-4)
  if (gname(inm:igp) /= psuf) gname = trim(gname)//psuf

  write(0,*)"DPLR grid file = "//trim(gname)


! check that the grid file can be opened
  open(55, file=gname, status="old", iostat=ios, err=510, action="read")
  close(55)
  510 continue
  if (ios /= 0) then
    write(0,113) trim(gname)
    write(0,*)"ERROR! Invalid input file specified."
    stop
  end if


! skip down to the block-specific flags
  errmsg = "problem skipping down to block-specific flags"
  do i = 1,43
    read(*,*,err=520,iostat=ios)
  end do
  520 continue
  if (ios > 0) goto 312


! allocate memory for block-specific data
  allocate(mtx(nblk),mty(nblk),mtz(nblk), stat=ierr); nalloc = nalloc + 3
  if (ierr /= 0) then
    write(0,*)
    write(0,*)"ERROR! failed to allocate mtx,mty,mtz in read_dplr_input"
    write(0,*)"nblk = ",nblk
    stop
  end if

  allocate(jbc(6,nblk),stat=ierr)
  if (ierr /= 0) then
    write(0,*)
    write(0,*)"ERROR! failed to allocate jbc in read_dplr_input"
    write(0,*)"nblk = ",nblk
    stop
  end if


! read block-specific data
  do n = 1,nblk

!   read block ntx line
    errmsg = "problem reading ntx line"
    do i = 1,6
      read(*,*,err=530,iostat=ios)
    end do
    read(*,*,err=530,iostat=ios) ntx,nty,ntz,iconr,isim,ifree,initi

!   compute grid sizes
    mtx(n) = ntx + 2
    mty(n) = nty + 2
    mtz(n) = ntz + 2

!   skip over the flux lines
    errmsg = "problem skipping flux lines"
    do i = 1,12
      read(*,*,err=530,iostat=ios)
    end do

!   read the boundary conditions
    errmsg = "problem reading jbc line"
    do i = 1,3
      read(*,*,err=530,iostat=ios)
    end do
    read(*,*,err=530,iostat=ios)jbc(:,n)

  end do
  530 continue
  if (ios > 0) goto 312


! come here on an error during input deck read
  312 continue


! report any errors
  if (ios > 0) then
    write(0,126) trim(errmsg)
    write(0,*)"ERROR! read_dplr_input failed."
    stop
  end if



  101 format(a200)
  111 format(/" ERROR! ",a,"= ",i4," is not a valid choice")
  113 format(/" ERROR! Problem opening grid file ",a)
  126 format(/" ERROR! ",a," in CFD input file")

  return
  end subroutine read_dplr_input

! =================================================================
  subroutine rdagrd(iug,n,ni,nj,nk)
! =================================================================
! Read one block of grid coordinates from the DPLR archival grid 
! file. Note that the loop over number of master blocks is 
! performed in the routine that calls this one, and the grid file
! header information has already been read. This subroutine is
! derived from a subroutine written for DPLR by Michael Wright 
! (mjwright@mail.arc.nasa.gov)

#ifdef USE_FXDR
  use IntrType    ,only: rd
  use NetAlloc    ,only: nalloc
#endif
  use VolumeArrays,only: xb,yb,zb

  implicit none

  integer,intent(in) :: iug,n,ni,nj,nk

#ifdef USE_FXDR
  real(kind=rd),allocatable,dimension(:,:,:) :: xtmp
  integer :: ierr
  integer :: isz,ixdrdmat,ixdrdouble
  integer :: mlx,mly,mlz
#endif


  continue


! write(*,115) n,ni-1,nj-1,nk-1
 
  select case(igrid)

    case(1)
!     native unformatted
      read(iug)xb
      read(iug)yb
      read(iug)zb

#ifdef USE_FXDR
    case(11)
!     XDR binary

!     allocate memory for grid coordinates: note that the DPLR
!     convention is to order the array kji instead of ijk
!     mlx, mly, mlz is the size of the grid array (with phantoms)
      mlx = ni + 2
      mly = nj + 2
      mlz = nk + 2

      allocate(xtmp(mlz,mly,mlx), stat=ierr); nalloc = nalloc + 1
      if (ierr /= 0) then
        write(0,*)"ERROR! failed to allocate xtmp"
        write(0,*)"mlz,mly,mlx = ",mlz,mly,mlx
        stop
      end if

      isz = mlx*mly*mlz

      ierr = ixdrdmat(iug, isz, xtmp)
      call swaporder(xtmp,xb,mlx,mly,mlz)
      ierr = ixdrdmat(iug, isz, xtmp)
      call swaporder(xtmp,yb,mlx,mly,mlz)
      ierr = ixdrdmat(iug, isz, xtmp)
      call swaporder(xtmp,zb,mlx,mly,mlz)

      deallocate(xtmp); nalloc = nalloc - 1
#endif

    case(21)
!     ASCII (formatted)
      read(iug,*)xb
      read(iug,*)yb
      read(iug,*)zb

    case default
      write(0,*)
      write(0,*)"ERROR! unknown grid format in rdagrd"
      write(0,*)"igrid = ",igrid
      stop

  end select

  115  format ('# --> Reading block ',i3,': grid cell size ',2(i4,'X'),i4)

  return
  end subroutine rdagrd

! =================================================================
  subroutine rdwrthead(iu,mmt)
! =================================================================
! Read the header information from DPLR grid files.
! This subroutine is derived from one coding written for DPLR by
! Michael Wright (mjwright@mail.arc.nasa.gov).
! The XDR binary formatted file type is read  using the
! FXDR library by David Pierce of UCSD, which can be obtained
! from http://meteora.ucsd.edu:80/~pierce/fxdr_home_page.html.
! For the XDR option to be available in USURP, the FXDR_HOME
! environment variable must be set.

  use IntrType,only: rd

  implicit none

  real(kind=rd) :: qvers

  integer,intent(out) :: iu

  integer,dimension(3,nblk) :: mmt
  integer,dimension(2) :: nitot,nin2
  integer,dimension(8) :: inta2
  integer :: ierr,ni,m,n,nbout,iotyp,nblkg,ier,idum
  integer :: id,nn,idc,npc,nblk2,nbout2
#ifdef USE_FXDR
  integer :: initxdr,ixdrint,ixdrdouble,ixdrs
#endif


  continue


  select case(igrid)

    case(1)

      iu = 21

      open(iu,file=gname,form="unformatted",status="old", &
              iostat=ierr,action="read")
      if (ierr /= 0) then
        write(0,*)
        write(0,*)"ERROR! problem opening grid file in rdwrthead"
        stop
      end if

      read(iu)iotyp,qvers

      if (iotyp >= 0) then
        iotyp = -1
        qvers = 2.27
      end if

      if (qvers >= 3.00) read(iu)

      read(iu) nblkg,nitot(1),nitot(2)

      if (nblkg /= nblk) then
        write(0,*)
        write(0,*)"ERROR! nblkg /= nblk in rdwrthead"
        stop
      end if

      do ni = 1,nitot(1)+nitot(2)
         read(iu)
         read(iu)
      end do

      read(iu) ((mmt(m,n), m=1,3), n=1,nblk)

      if (iotyp == -1) nbout = 1
      if (iotyp == -2) read(iu) nbout

      do id = 0,nbout-1

!       read(iu,*)nn,nblk2,idc,npc,nbout2,nin2(1,id),nin2(2,id), &
!               mmt(1,3),mmt(2,3),mmt(3,3)
!       read(iu,*)mmt(1,1),mmt(2,1),mmt(3,1),mmt(1,2),mmt(2,2),mmt(3,2)
        read(iu)nn,nblk2,idc,npc,nbout2,nin2(1),nin2(2), &
                idum,idum,idum
        read(iu)idum,idum,idum,idum,idum,idum

        if (nbout2 /= nbout) then
          write(0,*)"ERROR! nbout /= nbout2"
          write(0,*)"nbout,nbout2 = ",nbout,nbout2
          stop
        end if

        do n = 1,nin2(1)+nin2(2)
           read(iu)(inta2(m), m=1,8)
        end do

      end do

#ifdef USE_FXDR
    case(11)

!     open the file and get a unit file 
      ixdrs = initxdr(gname, 'r', .FALSE.)
      iu = ixdrs

      ier  = ixdrint(ixdrs, iotyp)
      ier  = ixdrdouble(ixdrs, qvers)

      if (iotyp >= 0) then
        iotyp = -1
        qvers = 2.27
      end if

      if (qvers < 2.31) then
        write(0,*)
        write(0,*)"ERROR! no logic for XDR file versions 2.31 or older"
        stop
      end if

      if (qvers >= 2.31) ier  = ixdrint(ixdrs, idum)
      if (qvers >= 3.00) ier  = ixdrint(ixdrs, idum)

      ier = ixdrint(ixdrs, nblkg)

      if (nblkg /= nblk) then
        write(0,*)
        write(0,*)"ERROR! nblkg /= nblk in xdrheader"
        write(0,*)"nblk,nblkg = ",nblk,nblkg
        stop
      end if

      ier = ixdrint(ixdrs, nitot(1))
      ier = ixdrint(ixdrs, nitot(2))

      do n = 1,nitot(1)+nitot(2)
      do m = 1,8
        ier = ixdrint(ixdrs, idum)
        ier = ixdrint(ixdrs, idum)
      end do
      end do

      do n = 1,nblk
      do m = 1,3
        ier = ixdrint(ixdrs, mmt(m,n))
      end do
      end do

      if (iotyp == -1) nbout = 1
      if (iotyp == -2) ier = ixdrint(ixdrs, nbout)

      do id = 0,nbout-1

        ier = ixdrint(ixdrs, nn)
        ier = ixdrint(ixdrs, nblk2)
        ier = ixdrint(ixdrs, idc)
        ier = ixdrint(ixdrs, npc)
        ier = ixdrint(ixdrs, nbout2)

        if (nbout2 /= nbout) then
          write(0,*)"ERROR! nbout /= nbout2"
          write(0,*)"nbout,nbout2 = ",nbout,nbout2
          stop
        end if

        ier = ixdrint(ixdrs,nin2(1))
        ier = ixdrint(ixdrs,nin2(2))

        do n = 1,3
        do m = 1,3
!         ier = ixdrint(ixdrs, mmt(m,n))
          ier = ixdrint(ixdrs, idum)
        end do
        end do

        do m = 1,8
        do n = 1,nin2(1)+nin2(2)
          ier = ixdrint(ixdrs, inta2(m))
        end do
        end do

      end do

      write(*,*)"rdwrthead ended successfully for igrid=11"
#endif

    case(21)

      iu = 21

      open(iu,file=gname,form="formatted",status="old", &
           iostat=ierr,action="read")
      if (ierr /= 0) then
        write(0,*)
        write(0,*)"ERROR! problem opening grid file in rdwrthead"
        stop
      end if

      read(iu,*)iotyp,qvers

      if (iotyp >= 0) then
        iotyp = -1
        qvers = 2.27
      end if

      if (qvers >= 3.00) read(iu,*)

      read(iu,*) nblkg,nitot(1),nitot(2)

      if (nblkg /= nblk) then
        write(0,*)
        write(0,*)"ERROR! nblkg /= nblk in rdwrthead"
        stop
      end if

      do ni = 1,nitot(1)+nitot(2)
         read(iu,*)
         read(iu,*)
      end do

      read(iu,*) ((mmt(m,n), m=1,3), n=1,nblk)

      if (iotyp == -1) nbout = 1
      if (iotyp == -2) read(iu,*) nbout

      do id = 0,nbout-1

!       read(iu,*)nn,nblk2,idc,npc,nbout2,nin2(1,id),nin2(2,id), &
!               mmt(1,3),mmt(2,3),mmt(3,3)
!       read(iu,*)mmt(1,1),mmt(2,1),mmt(3,1),mmt(1,2),mmt(2,2),mmt(3,2)
        read(iu,*)nn,nblk2,idc,npc,nbout2,nin2(1),nin2(2), &
                idum,idum,idum
        read(iu,*)idum,idum,idum,idum,idum,idum

        if (nbout2 /= nbout) then
          write(0,*)"ERROR! nbout /= nbout2"
          write(0,*)"nbout,nbout2 = ",nbout,nbout2
          stop
        end if

        do n = 1,nin2(1)+nin2(2)
           read(iu,*)(inta2(m), m=1,8)
        end do

      end do
 
    case default

      write(0,*)
      write(0,*)"ERROR! invalid grid format in rdwrthead"
      write(0,*)"igrid = ",igrid
      stop

  end select

  return
  end subroutine rdwrthead

! ========================================
  subroutine swaporder(xtmp,x,mlx,mly,mlz)
! ========================================
! reverse the index order of a 3D array

  use IntrType,only: rd

  implicit none

  integer,intent(in) :: mlx,mly,mlz

  real(kind=rd),intent(in),dimension(mlz,mly,mlx) :: xtmp
  real(kind=rd),intent(out),dimension(mlx,mly,mlz) :: x

  integer :: i,j,k


  continue


  do k = 1,mlz
  do j = 1,mly
  do i = 1,mlx
    x(i,j,k) = xtmp(k,j,i)
  end do
  end do
  end do

  return
  end subroutine swaporder

! ========================================================
  subroutine StorePatchesDPLR()
! ========================================================
! Read in the grid and iblank data for the surface
! patches using DPLR-formatted grid files and the SUGGAR
! DCI file

  use IntrType          ,only: rs
  use my_cpu_timeM      ,only: my_cpu_time
  use NetAlloc          ,only: nalloc
  use Types             ,only: scal,surf
  use PatchInfo         ,only: block,i1,i2,j1,j2,k1,k2
  use PatchInfo         ,only: num_panels,num_patches,num_scalars
  use StorePatchesUNCLEM,only: ReadDCIFile
  use UserInput         ,only: dci_file_name
  use UserInput         ,only: overset,showcpu
  use UserInput         ,only: cpufmt
  use UserInput         ,only: workdir
  use VertexData        ,only: num_nodes,xv
  use VolumeArrays      ,only: ibv,xb,yb,zb

  implicit none

  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: utime,dtime,ltime

  integer,dimension(3,nblk) :: mmt
  integer,dimension (3)     :: i1s,i1e,i2s,i2e
  integer                   :: ierr,ip,ipan,iug,n,nmax
  integer                   :: ni,nj,nk,i,j,k,iv
  integer                   :: dir1,dir2,dir3
#ifdef USE_FXDR
  integer                   :: ixdrclose
#endif

  logical,allocatable,dimension (:) :: used


  continue


! initialization
  call my_cpu_time(time1)
  utime = 0
  dtime = 0
  ltime = 0

  allocate(xv(3,num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)'ERROR! allocate of xv failed in StorePatchesUNCLE'
    write(0,*)'num_nodes = ',num_nodes
    stop
  end if

! for DPLR, gen_dirt.dci must be in the work directory in order
! for the iblank data to be used
  dci_file_name = "gen_dirt.dci"
  inquire(file=trim(workdir)//trim(dci_file_name),exist=overset)
  if (overset) then
    open(unit     = 41, &
         file     = trim(workdir)//trim(dci_file_name), &
         action   = "read", &
         status   = "old",  &
         position = "rewind")
  end if

! count the iblank variable as the only scalar
  num_scalars = 1

  allocate(scal(num_scalars,num_panels)); nalloc = nalloc + 1


! ascertain the largest block number that we need to use
  if (num_patches == 0) then
    write(0,*)
    write(0,*)"ERROR! num_patches = 0 in StorePatchesDPLR"
    stop
  else
    nmax = maxval(block(1:num_patches))
    if ((nmax <= 0).or.(nmax > nblk)) then
      write(0,*)
      write(0,*)"ERROR! invalid value of nmax; nmax = ",nmax
      stop
    end if
  end if


! determine which blocks are actually used
  allocate(used(nmax)); nalloc = nalloc + 1
  used = .false.
  do ip = 1,num_patches
    used(block(ip)) = .true.
  end do


! Open the grid file and read the header information from it.
! Return iug (unit number), nblk (number of master blocks),
! and mmt (array of grid dimensions).  Note that mmt is 
! the number of actual vertices plus one.
  call rdwrthead(iug,mmt)


! loop over the blocks (only up to the highest block 
! that is actually used)
  do n = 1,nmax

!   read the current block from the grid file
    call my_cpu_time(time3)

    ni = mmt(1,n) - 1
    nj = mmt(2,n) - 1
    nk = mmt(3,n) - 1

!   note that the DPLR grid files use arrays sized to hold
!   grid vertex phantom nodes
    allocate(xb(0:ni+1,0:nj+1,0:nk+1)); nalloc = nalloc + 1
    allocate(yb(0:ni+1,0:nj+1,0:nk+1)); nalloc = nalloc + 1
    allocate(zb(0:ni+1,0:nj+1,0:nk+1)); nalloc = nalloc + 1

    call rdagrd(iug,n,ni,nj,nk)

    call my_cpu_time(time4)
    utime = utime + (time4-time3)

    if (used(n)) then

      allocate(ibv(2:ni,2:nj,2:nk)); nalloc = nalloc + 1

      if (overset) then
!       read the i-blanking info from the dci file
        call my_cpu_time(time3)
        call ReadDCIFile(ibv,ni,nj,nk,n)
        call my_cpu_time(time4)
        dtime = dtime + (time4-time3)
      else
        ibv = 1
      end if


      call my_cpu_time(time3)
!     load the variables into the surface patches

      do ip = 1,num_patches

        if (block(ip) == n) then
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

!         figure out which direction is constant
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

!         copy the volume scalars to the surface
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
            i2s(1) = i
            i2s(2) = j
            i2s(3) = k

            i2e(1) = i
            i2e(2) = j
            i2e(3) = k

!           note: i2e is the interior cell, i2s is the ghost cell
            if (i2e(dir1) == 1) then
              i2e(dir1) = i2e(dir1) + 1
            else
              i2s(dir1) = i2s(dir1) + 1
            end if

            scal(1,ipan) = ibv(i2e(1),i2e(2),i2e(3))
            ipan = ipan + 1
          end do
          end do
          end do

        end if ! block(ip) == n
      end do ! ip = 1,num_patches

      call my_cpu_time(time4)
      ltime = ltime + (time4-time3)
      deallocate(ibv); nalloc = nalloc - 1

    end if ! used(n)

    deallocate(zb,yb,xb); nalloc = nalloc - 3

  end do ! n = 1,nmax

  if (overset) close(unit=41)
  deallocate(used); nalloc = nalloc - 1

! close the grid file
#ifdef USE_FXDR
  if (igrid == 11) then
    ierr = ixdrclose(iug)
  else
#endif
    close(iug)
#ifdef USE_FXDR
  end if
#endif

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in StorePatchesDPLR: ',time2-time1
    write(0,cpufmt)'  > time in grid I/O: ',utime
    if (overset) write(0,cpufmt)'  > time in DCI  I/O: ',dtime
    write(0,cpufmt)'  > time in storage : ',ltime
  end if

  return
  end subroutine StorePatchesDPLR

  end module DPLRM
