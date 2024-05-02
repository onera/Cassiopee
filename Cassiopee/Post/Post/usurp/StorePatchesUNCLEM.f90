! ================================================================
  module StorePatchesUNCLEM
! ================================================================

  implicit none

  private
  public  :: StorePatchesUNCLE
  public  :: ReadDCIFile
  private :: ReadEddyViscosityFile,ReadDistanceFile,ReadRestartFile

  contains


! ===========================================================
  subroutine StorePatchesUNCLE()
! ===========================================================
! Read in the grid and dependent variables for the surface
! patches using over-rel-formatted data files

  use CalculateVolumesM,only: CalculateVolumes
  use IntrType         ,only: rd,rs
  use my_cpu_timeM     ,only: my_cpu_time
  use NetAlloc         ,only: nalloc
  use Types            ,only: scal,surf
  use PatchInfo        ,only: block,i1,i2,j1,j2,k1,k2
  use PatchInfo        ,only: num_panels,num_patches,num_scalars
  use UserInput        ,only: cpufmt,dci_file_name,ignore_solution
  use UserInput        ,only: iwfn,overset,plotvel,prefix,showcpu
  use UserInput        ,only: solution_exists,unclem,workdir
  use VertexData       ,only: num_nodes,xv
  use VolumeArrays     ,only: ibv,volv,xb,yb,zb

  implicit none

  real(kind=rd),allocatable,dimension (:,:,:) :: pv,uv,vv,wv,ev
  real(kind=rd),allocatable,dimension (:,:,:) :: dwcv
  real(kind=rd),parameter                     :: one = 1.0_rd
  real(kind=rd),parameter                     :: half = 0.5_rd

  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: utime,dtime,ltime

  integer,dimension (3) :: i1s,i1e,i2s,i2e
  integer               :: ierr,ip,ipan,n,nmax,nn
  integer               :: ni,nj,nk,i,j,k,iv
  integer               :: dir1,dir2,dir3
  integer               :: ncyc

  character(len=80) :: gname,fname,ename,dname
  character(len= 5) :: fileid

  logical,allocatable,dimension (:) :: used
  logical :: ex


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

  if (overset) then
    inquire(file=trim(workdir)//trim(dci_file_name),exist=ex)
    if (.not.ex) then
      write(*,*)
      write(*,*)"file not found: ",trim(workdir)//trim(dci_file_name)
      write(*,*)"enter dci_file_name, relative to ",trim(workdir)
      read(*,'(a80)')dci_file_name
    end if
    open(unit     = 41, &
         file     = trim(workdir)//trim(dci_file_name), &
         action   = "read", &
         status   = "old",  &
         position = "rewind")
  end if


! see if the flow solution is present
  if (.not.unclem) then
    inquire(file=trim(workdir)//"IO_FILES/"//prefix//"-10", &
            exist=solution_exists)
    if (ignore_solution) solution_exists = .false.
  end if

  if (solution_exists) then
    if (plotvel) then
      num_scalars = 12
    else
      if (iwfn == 1) then
        num_scalars = 11
      else
        num_scalars = 7
      end if
    end if
  else
    num_scalars = 1
  end if
  allocate(scal(num_scalars,num_panels)); nalloc = nalloc + 1


! ascertain the largest block number that we need to use
  if (num_patches == 0) then
    write(*,*)"ERROR! num_patches = 0"
    stop
  end if
  nmax = block(1)
  do ip = 2,num_patches
    nmax = max(nmax,block(ip))
  end do


! determine which blocks are actually used
  allocate(used(nmax)); nalloc = nalloc + 1
  used = .false.
  do ip = 1,num_patches
    used(block(ip)) = .true.
  end do

    
  do n = 1,nmax
    if (used(n)) then

      if (unclem) then
        write(unit=fileid,fmt="(i5)")n-1
        if (n-1 < 10) then
          gname = trim(workdir)//"GRIDS/uncle.grd0"//adjustl(trim(fileid))
        else
          gname = trim(workdir)//"GRIDS/uncle.grd"//adjustl(trim(fileid))
        end if
      else
        write(unit=fileid,fmt="(i5)")n*10
        gname = trim(workdir)//"IO_FILES/G"//adjustl(trim(fileid))
      end if


      call my_cpu_time(time3)
!     read the grid file
      open(unit     = 3,             &
           file     = gname,         &
           form     = "unformatted", &
           action   = "read",        &
           status   = "old",         &
           position = "rewind")
      read(unit=3)nn
      read(unit=3)ni,nj,nk
      allocate(xb(ni,nj,nk)); nalloc = nalloc + 1
      allocate(yb(ni,nj,nk)); nalloc = nalloc + 1
      allocate(zb(ni,nj,nk)); nalloc = nalloc + 1
      read(unit=3)xb,yb,zb
      close(unit=3)


      if (solution_exists) then

!       calculate the cell volumes
        allocate(volv(ni+1,nj+1,nk+1)); nalloc = nalloc + 1
        call CalculateVolumes(ni,nj,nk)


        if (iwfn == 1) then
!         read the distance file
          allocate(dwcv(ni,nj,nk)); nalloc = nalloc + 1
          dname = trim(workdir)//"IO_FILES/D"//adjustl(fileid)
          call ReadDistanceFile(dname,dwcv,ni,nj,nk)
        else
          dname = " "
        end if


!       read the restart file
        allocate(  pv(0:ni+2,0:nj+2,0:nk+2)); nalloc = nalloc + 1
        allocate(  uv(0:ni+2,0:nj+2,0:nk+2)); nalloc = nalloc + 1
        allocate(  vv(0:ni+2,0:nj+2,0:nk+2)); nalloc = nalloc + 1
        allocate(  wv(0:ni+2,0:nj+2,0:nk+2)); nalloc = nalloc + 1
        fname = trim(workdir)//"IO_FILES/"//prefix//"-"//adjustl(fileid)
        call ReadRestartFile(fname,pv,uv,vv,wv,ni,nj,nk,ncyc)
        if (n == 1) then
          if (prefix == "RMS") then
            write(*,'(1x,a30,i7)') &
              'Reading '//prefix//' files with navg = ',ncyc
          else
            write(*,'(1x,a30,i7)') &
              'Reading '//prefix//' files with ncyc = ',ncyc
          end if
        end if


!       read the eddy viscosity file
        ename = trim(workdir)//"IO_FILES/MUT-"//adjustl(fileid)
        allocate(ev(0:ni+2,0:nj+2,0:nk+2)); nalloc = nalloc + 1
        inquire(file=ename,exist=ex)
        if (ex) then
          call ReadEddyViscosityFile(ename,ev,ni,nj,nk)
        else
          ename = " "
          ev = 1.0_rd
        end if

      else
        fname = " "
        ename = " "
        dname = " "
      end if !solution_exists

      call my_cpu_time(time4)
      utime = utime + (time4-time3)


      call my_cpu_time(time3)
!     read the i-blanking info from the dci file
      allocate(ibv(2:ni,2:nj,2:nk)); nalloc = nalloc + 1
      if (overset) then
        call ReadDCIFile(ibv,ni,nj,nk,n)
      else
        ibv = 1
      end if
      call my_cpu_time(time4)
      dtime = dtime + (time4-time3)



      call my_cpu_time(time3)
!     load the variables into the surface patches

      do ip = 1,num_patches

        if (block(ip) == n) then

!         store the grid points for the patch
!         xv(1,surf(ip)%first_node:surf(ip)%last_node) =  &
!           reshape( xb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                   (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
!         xv(2,surf(ip)%first_node:surf(ip)%last_node) =  &
!           reshape( yb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                   (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
!         xv(3,surf(ip)%first_node:surf(ip)%last_node) =  &
!           reshape( zb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                   (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
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

!         store ave p,del u,del v,del w,ave e, vol, wall distance
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

            if (solution_exists) then

              scal(2,ipan) = half*(pv(i,j,k)+pv(i2e(1),i2e(2),i2e(3)))
              scal(3,ipan) = uv(i2e(1),i2e(2),i2e(3))-uv(i2s(1),i2s(2),i2s(3))
              scal(4,ipan) = vv(i2e(1),i2e(2),i2e(3))-vv(i2s(1),i2s(2),i2s(3))
              scal(5,ipan) = wv(i2e(1),i2e(2),i2e(3))-wv(i2s(1),i2s(2),i2s(3))
              scal(6,ipan) = max(half*(ev(i,j,k)+ev(i2e(1),i2e(2),i2e(3))),one)
              scal(7,ipan) = volv(i2e(1),i2e(2),i2e(3))

            end if


            if (solution_exists.and.((iwfn == 1).or.plotvel)) then
              if (iwfn == 1) then
                scal( 8,ipan) = dwcv(i2e(1),i2e(2),i2e(3))
              else
                scal( 8,ipan) = 0.0
              end if
              scal( 9,ipan) =   uv(i2e(1),i2e(2),i2e(3))
              scal(10,ipan) =   vv(i2e(1),i2e(2),i2e(3))
              scal(11,ipan) =   wv(i2e(1),i2e(2),i2e(3))
            end if
 
            ipan = ipan + 1
          end do
          end do
          end do

        end if ! block(ip) == n
      end do ! ip = 1,num_patches

      call my_cpu_time(time4)
      ltime = ltime + (time4-time3)
      deallocate(ibv); nalloc = nalloc - 1
      if (solution_exists) then
        deallocate(ev,wv,vv,uv,pv); nalloc = nalloc - 5
      end if
      if (allocated(dwcv)) then
        deallocate(dwcv); nalloc = nalloc - 1
      end if
      if (allocated(volv)) then
        deallocate(volv); nalloc = nalloc - 1
      end if
      deallocate(zb,yb,xb); nalloc = nalloc - 3

    end if ! used(n)
  end do ! n = 1,nmax

  if (overset) close(unit=41)
  deallocate(used); nalloc = nalloc - 1

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in StorePatchesUNCLE: ',time2-time1
    write(0,cpufmt)'  > time in UNCLE I/O: ',utime
    write(0,cpufmt)'  > time in DCI   I/O: ',dtime
    write(0,cpufmt)'  > time in storage  : ',ltime
  end if

  return
  end subroutine StorePatchesUNCLE

! ===========================================================
  subroutine ReadDistanceFile(dname,disw,ni,nj,nk)
! ===========================================================
! This procedure reads the cell-centered distance-to-the-wall
! which is used in the wall function formulation in OVER-REL
! mode.

  use IntrType,only: rd

  implicit none

  integer,intent(in)           :: ni,nj,nk

  real(kind=rd),intent(out),dimension(ni,nj,nk) :: disw

  integer                      :: ni_temp,nj_temp,nk_temp
  character(len=80),intent(in) :: dname


  continue


  open(unit=3,file=dname,form     = "unformatted", &
                         action   = "read",        &
                         status   = "old",         &
                         position = "rewind")
  read(3)ni_temp,nj_temp,nk_temp
  read(3)disw
  close(3)

  return
  end subroutine ReadDistanceFile

! ==================================================================
  subroutine ReadEddyViscosityFile(ename,emut,ni,nj,nk)
! ==================================================================
! This procedure reads the eddy viscosity files for UNCLE / OVER-REL
! which is used in the viscous force integration (especially for
! cases using wall functions or surface roughness).

  use IntrType,only: rd

  implicit none

  integer,intent(in)                            :: ni,nj,nk
  real(kind=rd),intent(out),dimension(0:ni+2,0:nj+2,0:nk+2) :: emut
  character(len=80),intent(in)                  :: ename


  continue


  open(unit=3,file=ename,form     ="unformatted", &
                         action   ="read",        &
                         status   ="old",         &
                         position ="rewind")
  read(3)emut
  emut = emut + 1.0  !change turbulent viscosity to effective viscosity
  close(3)

  end subroutine ReadEddyViscosityFile

! ===========================================================
  subroutine ReadRestartFile(fname,p,u,v,w,ni,nj,nk,ncyc)
! ===========================================================
! Read the UNCLE / OVER-REL restart files containing the
! flowfield data needed for the force and moment integration.

  use IntrType,only: rd

  implicit none

  integer,intent(in) :: ni,nj,nk

  real(kind=rd),dimension(0:ni+2,0:nj+2,0:nk+2),intent(out) :: p
  real(kind=rd),dimension(0:ni+2,0:nj+2,0:nk+2),intent(out) :: u
  real(kind=rd),dimension(0:ni+2,0:nj+2,0:nk+2),intent(out) :: v
  real(kind=rd),dimension(0:ni+2,0:nj+2,0:nk+2),intent(out) :: w
  real(kind=rd),dimension(0:ni+2,0:nj+2,0:nk+2) :: qt1,qt2

  character(len=80),intent(in) :: fname
  integer,intent(out) :: ncyc


  continue


  open(unit=3,file=fname,form="unformatted", &
                         action="read", &
                         status="old", &
                         position="rewind")
  read(3)ncyc
  read(3)p,u,v,w,qt1,qt2
  close(3)

  return
  end subroutine ReadRestartFile

! ===============================================================
  subroutine ReadDCIFile(ibv,ni,nj,nk,n)
! ===============================================================
! This procedure reads a portion of the domain connectivity
! information (DCI) file generated by SUGGAR, which stores
! the i-blanking information for overset grids.

  implicit none

  integer,intent(in)  :: n,ni,nj,nk
  integer,intent(out) :: ibv(2:ni,2:nj,2:nk)
  integer,save        :: nsofar = 0
  integer             :: ntag,i

  character(len=132)  :: tag



  continue


  skip_through: do

    read(unit=41,fmt="(a)")tag

    if (index(tag,":general_cell_centered_drt_pairs") /= 0) then
      read(unit=tag(33:132),fmt=*)ntag
      do i = 1,ntag
        read(unit=41,fmt=*)
      end do
    else if (index(tag,":general_cell_centered_drt_iblank") /= 0) then
      read(unit=tag(34:132),fmt=*)ntag
      nsofar = nsofar + 1
      if (nsofar == n) then
        read(unit=41,fmt=*)ibv
        exit skip_through
      else
        do i = 1,ntag
          read(unit=41,fmt=*)
        end do
      end if
    end if

  end do skip_through


  return
  end subroutine ReadDCIFile

  end module StorePatchesUNCLEM
