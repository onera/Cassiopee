! ================================================================
  module StorePatchesGenericM
! ================================================================

  implicit none

  private
  public :: StorePatchesGeneric

  contains


! ===========================================================
  subroutine StorePatchesGeneric()
! ===========================================================
! Read in the grid for the surface
! patches from a formatted, plot3d file

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: scal,surf
  use PatchInfo   ,only: i1,i2,j1,j2,k1,k2,block
  use PatchInfo   ,only: num_panels,num_patches,num_scalars
  use UserInput   ,only: cpufmt,fgrid,overset,showcpu
  use UserInput   ,only: workdir
  use VertexData  ,only: num_nodes,xv

  implicit none

  real(kind=rd),allocatable,dimension (:,:,:) :: xb,yb,zb
  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: utime,dtime,ltime

  integer,allocatable,dimension (:,:,:) :: ibv
  integer,allocatable,dimension (:) :: nig,njg,nkg
  integer,dimension (3) :: i1s,i1e

  integer :: ierr,ip,ipan,n,nb
  integer :: i,j,k,iv
  integer :: dir1,dir2,dir3
  logical :: ex


  continue

! initialization
  call my_cpu_time(time1)
  utime = 0.0_rs
  dtime = 0.0_rs
  ltime = 0.0_rs
  allocate(xv(3,num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)"ERROR! allocate of xv failed in StorePatchesGeneric"
    write(0,*)"num_nodes = ",num_nodes
    stop
  end if


! look for the iblank file
  inquire(file  = trim(workdir)//"generic.ib",exist = ex)

  if (.not.ex) then
    write(*,*)
    write(*,*)"file not found: ",trim(workdir)//"generic.ib"
    write(*,*)"assuming no overset."
    overset = .false.
  else
    write(*,*)
    write(*,*)"file found: ",trim(workdir)//"generic.ib"
    write(*,*)"assuming overset."
    overset = .true.
  end if

! open the iblank file
  if (overset) then
    open(unit     = 41, &
         file     = trim(workdir)//"generic.ib", &
         action   = "read", &
         status   = "old",  &
         position = "rewind")
  end if


! assume there is no flow data
  allocate(scal(1,num_panels)); nalloc = nalloc + 1


  fgrid = "generic.grd"

! open the grid file
  open ( unit    = 3,                          &
         file    = trim(workdir)//trim(fgrid), &
         form    = "formatted",                &
         status  = "old",                      &
         action  = "read",                     &
         position= "rewind"                    )


  read(unit=3,fmt=*)nb
  allocate(nig(nb),njg(nb),nkg(nb)); nalloc = nalloc + 3
  read(unit=3,fmt=*)(nig(n),njg(n),nkg(n),n=1,nb)
  do n = 1,nb
     
!   allocate space for the grid and i-blank variable
    allocate( xb(nig(n),njg(n),nkg(n)))      ; nalloc = nalloc + 1
    allocate( yb(nig(n),njg(n),nkg(n)))      ; nalloc = nalloc + 1
    allocate( zb(nig(n),njg(n),nkg(n)))      ; nalloc = nalloc + 1
!    allocate(ibv(2:nig(n),2:njg(n),2:2))     ; nalloc = nalloc + 1
    if ( nig(n) == 1 ) then 
       allocate(ibv(2:2,2:njg(n),2:nkg(n)))     ; nalloc = nalloc + 1
    else if (njg(n) == 1 ) then
       allocate(ibv(2:nig(n),2:2,2:nkg(n)))     ; nalloc = nalloc + 1
    else
       allocate(ibv(2:nig(n),2:njg(n),2:2))     ; nalloc = nalloc + 1
    endif

!   store i-blank variable for next block
    call my_cpu_time(time3)
    if (overset) then
!     read i-blank variable
      read(41,*)ibv
    else
      ibv = 1
    end if

    call my_cpu_time(time4)
    dtime = dtime + (time4-time3)


!   read the next grid block
    call my_cpu_time(time3)
    read(unit=3,fmt=*)xb,yb,zb
    call my_cpu_time(time4)
    utime = utime + (time4-time3)


!   fill in the grid for each patch
    call my_cpu_time(time3)
    do ip = 1,num_patches
      if (block(ip) == n) then


!       store the grid points for the patch
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


!       figure out which direction is constant
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


!       store the i-blank variable
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

      end if
    end do
    call my_cpu_time(time4)
    ltime = ltime + (time4-time3)
    deallocate(ibv,zb,yb,xb); nalloc = nalloc - 4

  end do

  close(unit=3)
  if (overset) close(unit=41)
  deallocate(nkg,njg,nig); nalloc = nalloc - 3

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)"cpu time in StorePatchesGeneric: ",time2-time1
    write(0,cpufmt)"  > time in grid I/O: "  ,utime
    write(0,cpufmt)"  > time in iblank I/O: ",dtime
    write(0,cpufmt)"  > time in storage  : " ,ltime
  end if

  return
  end subroutine StorePatchesGeneric

  end module StorePatchesGenericM
