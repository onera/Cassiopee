! ================================================================
  module WriteGridIBIM
! ================================================================

  implicit none

  private
  public :: WriteGridIBI

  contains

! =================================================================
  subroutine WriteGridIBI()
! =================================================================
! write a grid.ibi file for use in OVERFLOW mode.  

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use PatchInfo   ,only: i1,i2,j1,j2,k1,k2,ibdir,num_patches
  use Types       ,only: surf
  use UserInput   ,only: showcpu,cpufmt,double,pathname
  use VertexData  ,only: xv

  implicit none

  real(kind=rd),allocatable,dimension(:,:) :: xd,yd,zd
  real(kind=rs)                            :: time1,time2

  integer,allocatable,dimension(:,:) :: ibi
  integer,allocatable,dimension(:)   :: njo,nko

  integer :: ip,iv,j,k,ierr
  integer :: kbeg,kend,inc


  continue


  call my_cpu_time(time1)

  allocate(njo(num_patches),nko(num_patches)); nalloc = nalloc + 2

  do ip = 1,num_patches
    if (abs(ibdir(ip)) == 1) then
      njo(ip) = j2(ip) - j1(ip) + 1
      nko(ip) = k2(ip) - k1(ip) + 1
    else if (abs(ibdir(ip)) == 2) then
      njo(ip) = k2(ip) - k1(ip) + 1
      nko(ip) = i2(ip) - i1(ip) + 1
    else if (abs(ibdir(ip)) == 3) then
      njo(ip) = i2(ip) - i1(ip) + 1
      nko(ip) = j2(ip) - j1(ip) + 1
    end if
  end do

  open(unit    = 7,                          &
       file    = trim(pathname)//"grid.ibi", &
       form    = "unformatted",              &
#ifdef CONVERT_BIG_ENDIAN
       convert = "big_endian",               &
#endif
       action  = "write"                     )

  write(7)num_patches
  write(7)(njo(ip),nko(ip),1,ip=1,num_patches)

  do ip = 1,num_patches

    allocate( xd(njo(ip),nko(ip)),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR!  failed to allocate xd."
      stop
    end if

    allocate( yd(njo(ip),nko(ip)),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR!  failed to allocate yd."
      stop
    end if

    allocate( zd(njo(ip),nko(ip)),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR!  failed to allocate zd."
      stop
    end if

    allocate(ibi(njo(ip),nko(ip)),stat=ierr); nalloc = nalloc + 1
    if (ierr /= 0) then
      write(0,*)"ERROR!  failed to allocate ibi."
      stop
    end if


!   Setting the integration i-blanks to one should be okay because
!   of the fact that we are going to use the USURP panel weights in the
!   integration.  Assuming they are equal to one instead of writing
!   the original values avoids having to store the data (which is not
!   otherwise used) or having to re-read the original grid file. 

    ibi = 1

!   The following section is needed to store the grid points in a
!   temporary array in order to reverse the order of the points for
!   the case where ibdir = 2

    iv = surf(ip)%first_node
    if ((abs(ibdir(ip)) == 1).or.(abs(ibdir(ip)) == 3)) then

!     copy the array into temporary storage
      do k = 1,nko(ip)
      do j = 1,njo(ip)
        xd(j,k) = xv(1,iv)
        yd(j,k) = xv(2,iv)
        zd(j,k) = xv(3,iv)
        iv = iv + 1
      end do
      end do

    else if (abs(ibdir(ip)) == 2) then

!     reverse the order of the array storage
      do j = 1,njo(ip)
      do k = 1,nko(ip)
        xd(j,k) = xv(1,iv)
        yd(j,k) = xv(2,iv)
        zd(j,k) = xv(3,iv)
        iv = iv + 1
      end do
      end do

    end if


!   Reverse the column order for "max" surfaces (ibdir < 0)
    if (ibdir(ip) > 0) then
      kbeg = 1
      kend = nko(ip)
      inc  = 1
    else
      kbeg = nko(ip)
      kend = 1
      inc  = -1
    end if

        
    if (double) then
      write(7)((xd(j,k),j=1,njo(ip)),k=kbeg,kend,inc), &
              ((yd(j,k),j=1,njo(ip)),k=kbeg,kend,inc), &
              ((zd(j,k),j=1,njo(ip)),k=kbeg,kend,inc), &
               ibi
    else
      write(7)((real(xd(j,k)),j=1,njo(ip)),k=kbeg,kend,inc), &
              ((real(yd(j,k)),j=1,njo(ip)),k=kbeg,kend,inc), &
              ((real(zd(j,k)),j=1,njo(ip)),k=kbeg,kend,inc), &
               ibi
    end if

    deallocate(ibi,zd,yd,xd); nalloc = nalloc - 4

  end do

  close(7)

  deallocate(nko,njo); nalloc = nalloc - 2


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time to write grid.ibi: ',time2-time1
  end if


  return
  end subroutine WriteGridIBI

  end module WriteGridIBIM
