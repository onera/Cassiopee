! ================================================================
  module StorePatchesCFDSHIPM
! ================================================================

  implicit none

  private
  public  :: StorePatchesCFDSHIP
  private :: InterpolateIBlank

  contains


! ===========================================================
  subroutine StorePatchesCFDSHIP()
! ===========================================================
! Read in the grid for the surface
! patches from a formatted, plot3d file

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use Types       ,only: scal,surf
  use PatchInfo   ,only: block,i1,i2,j1,j2,k1,k2,num_panels,num_patches
  use UserInput   ,only: cpufmt,fgrid,fnameo,overset,showcpu
  use UserInput   ,only: workdir
  use VertexData  ,only: num_nodes,xv

  implicit none

  real(kind=rd),allocatable,dimension (:,:,:) :: xb,yb,zb
  real(kind=rd) :: dumm
  real(kind=rs) :: time1,time2,time3,time4
  real(kind=rs) :: utime,dtime,ltime

  integer,allocatable,dimension (:,:,:) :: ibv
  integer,allocatable,dimension (:) :: nig,njg,nkg
  integer,dimension (3) :: i1s,i1e

  integer :: ierr,ip,ipan,n,nb
  integer :: i,j,k,iv
  integer :: ibv1,ibv2,ibv3,ibv4
  integer :: dir1,dir2,dir3
  integer :: ibpnts,iipnts,iieptr,iisptr,imax,jmax,kmax
  integer :: idum
  logical :: ex


  continue


! initialization
  call my_cpu_time(time1)
  utime = 0.0_rs
  dtime = 0.0_rs
  ltime = 0.0_rs

  allocate(xv(3,num_nodes),stat=ierr); nalloc = nalloc + 1
  if (ierr /= 0) then
    write(0,*)'ERROR! allocate of xv failed in StorePatchesCFDSHIP'
    write(0,*)'num_nodes = ',num_nodes
    stop
  end if


! look for the pegasus file
  inquire(file  = trim(workdir)//trim(fnameo)//'.xintout', &
          exist = ex)

  if (.not.ex) then
    write(*,*)
    write(*,*)"file not found: ",trim(workdir)//trim(fnameo)//'.xintout'
    write(*,*)"assuming no overset."
    overset = .false.
  else
    write(*,*)
    write(*,*)"file found: ",trim(workdir)//trim(fnameo)//'.xintout'
    write(*,*)"assuming overset."
    overset = .true.
  end if


! open the pegasus file
  if (overset) then
    open(unit     = 41, &
         file     = trim(workdir)//trim(fnameo)//'.xintout', &
         action   = "read", &
         status   = "old",  &
         position = "rewind")
  end if


! assume there is no flow data
  allocate(scal(1,num_panels)); nalloc = nalloc + 1


! open the grid file
  open ( unit    = 3,                          &
         file    = trim(workdir)//trim(fgrid), &
         form    = 'formatted',                &
         status  = 'old',                      &
         action  = 'read',                     &
         position= 'rewind'                    )


  read(unit=3,fmt=*)nb
  allocate(nig(nb),njg(nb),nkg(nb)); nalloc = nalloc + 3
! do n = 1,nb
!   read(unit=3,fmt=*)nig(n),njg(n),nkg(n)
! end do
  read(unit=3,fmt=*)(nig(n),njg(n),nkg(n),n=1,nb)

  do n = 1,nb


!   allocate space for the grid and i-blank variable
    allocate( xb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate( yb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate( zb(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1
    allocate(ibv(nig(n),njg(n),nkg(n))); nalloc = nalloc + 1


!   store i-blank variable for next block
    call my_cpu_time(time3)

    if (overset) then
      read(unit=41,fmt=*)ibpnts,iipnts,iieptr,iisptr,imax,jmax,kmax

      if (imax /= nig(n)) then
        write(0,*)'ERROR! xintout file indices do not match grid'
        write(0,*)'block ',n,' imax in xintout = ',imax,' ni = ',nig(n)
        stop
      end if
      if (jmax /= njg(n)) then
        write(0,*)'ERROR! xintout file indices do not match grid'
        write(0,*)'block ',n,' jmax in xintout = ',jmax,' nj = ',njg(n)
        stop
      end if
      if (kmax /= nkg(n)) then
        write(0,*)'ERROR! xintout file indices do not match grid'
        write(0,*)'block ',n,' kmax in xintout = ',kmax,' nk = ',nkg(n)
        stop
      end if

!     skip over interpolation stencils
      read(unit=41,fmt=*)(idum,i=iisptr,iieptr), &
                         (idum,i=iisptr,iieptr), &
                         (idum,i=iisptr,iieptr), &
                         (dumm,i=iisptr,iieptr), &
                         (dumm,i=iisptr,iieptr), &
                         (dumm,i=iisptr,iieptr)
      read(unit=41,fmt=*)(idum,i=1,ibpnts), &
                         (idum,i=1,ibpnts), &
                         (idum,i=1,ibpnts), &
                         (idum,i=1,ibpnts)


!     read i-blank variable
      read(unit=41,fmt=*)ibv
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
!       xv(1,surf(ip)%first_node:surf(ip)%last_node) =  &
!         reshape( xb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                 (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
!       xv(2,surf(ip)%first_node:surf(ip)%last_node) =  &
!         reshape( yb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                 (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
!       xv(3,surf(ip)%first_node:surf(ip)%last_node) =  &
!         reshape( zb(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip)),&
!                 (/ surf(ip)%last_node-surf(ip)%first_node+1 /))
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
            ibv1 = ibv(i,j  ,k  )
            ibv2 = ibv(i,j-1,k  )
            ibv3 = ibv(i,j  ,k-1)
            ibv4 = ibv(i,j-1,k-1)
          else if (dir1 == 2) then
            ibv1 = ibv(i  ,j,k  )
            ibv2 = ibv(i-1,j,k  )
            ibv3 = ibv(i  ,j,k-1)
            ibv4 = ibv(i-1,j,k-1)
          else if (dir1 == 3) then
            ibv1 = ibv(i  ,j  ,k)
            ibv2 = ibv(i-1,j  ,k)
            ibv3 = ibv(i  ,j-1,k)
            ibv4 = ibv(i-1,j-1,k)
          else
            write(0,*)'ERROR! invalid dir1 in StorePatchesCFDSHIP'
            stop
          end if

          scal(1,ipan) = InterpolateIBlank(ibv1,ibv2,ibv3,ibv4)

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
    write(0,cpufmt)'cpu time in StorePatchesCFDSHIP: ',time2-time1
    write(0,cpufmt)'  > time in grid I/O: '   ,utime
    write(0,cpufmt)'  > time in XINTOUT I/O: ',dtime
    write(0,cpufmt)'  > time in storage  : '  ,ltime
  end if

  return
  end subroutine StorePatchesCFDSHIP

! =============================================================
  pure function InterpolateIBLANK(in1,in2,in3,in4) result (ic)
! =============================================================
! Function to move the iblank variable from corners to
! the panel center.  For now, assume iblank can only be
! equal to 0 or 1, so we can take the center to be the
! minimum of the corners.

  implicit none

  integer,intent(in) :: in1,in2,in3,in4
  integer            :: iv1,iv2,iv3,iv4
  integer            :: ic


  continue


  ic = min(in1,in2,in3,in4)
  return


! logic from here down is for the more complicated case
! where iblank can include negative values (fringe) and
! flags for oprhans (101)

  iv1 = max(in1,-1)
  iv2 = max(in2,-1)
  iv3 = max(in3,-1)
  iv4 = max(in4,-1)

! if all values are equal, return that value
  if ((iv1 == iv2).and.(iv2 == iv3).and.(iv3 == iv4)) then
    ic = iv1
    return
  end if

! if any value is 101 (orphan), return that value
  if (max(iv1,iv2,iv3,iv4) == 101) then
    ic = 101
    return
  end if

! if any value is zero (hole), return that value
  if (iv1*iv2*iv3*iv4 == 0) then
    ic = 0
    return
  end if

! we must have a mix of fringe and interior cells, return -1
  ic = -1
  return

  end function InterpolateIBLANK

  end module StorePatchesCFDSHIPM
