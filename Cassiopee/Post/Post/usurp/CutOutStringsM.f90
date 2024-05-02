! ================================================================
  module CutOutStringsM
! ================================================================

  implicit none

  private
  public  :: CutOutStrings
  private :: intercept2

  contains

! ======================================================
  subroutine CutOutStrings(p)
! ======================================================
! This procedure redefines the component number of panels
! based on their relationship with provided strings.
! Each string is specified along with a component number
! marking the interior and exterior of the string.
! Note that because the dummy argument is a pointer,
! this procedure must have an explicit interface.
!
! Revisions:
! 01/02/2007: Extended from jmin patches only to all types

  use IntrType    ,only: rd,rs
  use my_cpu_timeM,only: my_cpu_time
  use NetAlloc    ,only: nalloc
  use OverlappedM ,only: Overlapped
  use PatchInfo   ,only: num_patches,i1,i2,j1,j2,k1,k2
  use Types       ,only: panel,surf
  use UserInput   ,only: cpufmt,showcpu,workdir
  use VertexData  ,only: xv

  implicit none
  
  type(panel),pointer,dimension(:) :: p

  real(kind=rd),allocatable,dimension(:,:,:) :: xp,yp,zp
  real(kind=rd),allocatable,dimension(:)     :: xstr,ystr,zstr
  real(kind=rd),dimension(3,2)               :: str_box

  real(kind=rd) :: xp1,yp1,xp2,yp2,xp3,yp3,xp4,yp4,xs1,ys1,xs2,ys2,ss
  real(kind=rd) :: cross,dxv,dyv,dxs,dys
  real(kind=rs) :: time1,time2

  integer,allocatable,dimension(:,:,:)   :: mc,mp
  integer,allocatable,dimension(:,:)     :: mstack

  integer :: nstr,i,ierr,ip,idir,iseg,iv,ipan,iflag,in_stack,j,k
  integer :: iinc,jinc,kinc,nn
  integer :: ip1,jp1,kp1,ip2,jp2,kp2
  integer :: ip3,jp3,kp3,ip4,jp4,kp4
  integer :: new_interior,mpmax,mpmin,new_exterior
  logical :: ex


  continue


  call my_cpu_time(time1)

  inquire(file=trim(workdir)//"strings.dat",exist=ex)
  if (.not.ex) return

! compute a bounding box for every patch
  do ip = 1,num_patches
    do idir = 1,3
      surf(ip)%box(idir,1) = &
          minval(xv(idir,surf(ip)%first_node:surf(ip)%last_node))
      surf(ip)%box(idir,2) = &
          maxval(xv(idir,surf(ip)%first_node:surf(ip)%last_node))
    end do
  end do

  open(unit=3,file=trim(workdir)//"strings.dat", &
       action="read",status="old")

  STRING_LOOP : do

    read(3,*,iostat=ierr)nstr,new_interior,new_exterior
    if (ierr /= 0) exit STRING_LOOP

!   string allocation
    allocate(xstr(nstr),ystr(nstr),zstr(nstr)); nalloc = nalloc + 3

    do i = 1,nstr
      read(3,*)xstr(i),ystr(i),zstr(i)
    end do

!   compare each patch with the string
    PATCH_LOOP: do ip = 1,num_patches

      iinc = 1
      jinc = 1
      kinc = 1
      if (i1(ip) == i2(ip)) then
        iinc = 0
      else if (j1(ip) == j2(ip)) then
        jinc = 0
      else if (k1(ip) == k2(ip)) then
        kinc = 0
      else
        write(0,*)
        write(0,*)"ERROR! no constant direction for ip = ",ip
        write(0,*)"i1,i2 = ",i1(ip),i2(ip)
        write(0,*)"j1,j2 = ",j1(ip),j2(ip)
        write(0,*)"k1,k2 = ",k1(ip),k2(ip)
        write(0,*)"reported by procedure CutOutStrings"
        stop
      end if

!     patch allocation
      allocate(xp(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip))); nalloc = nalloc + 1
      allocate(yp(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip))); nalloc = nalloc + 1
      allocate(zp(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip))); nalloc = nalloc + 1
      allocate(mp(i1(ip):i2(ip),j1(ip):j2(ip),k1(ip):k2(ip))); nalloc = nalloc + 1

      allocate(mc(i1(ip)+iinc:i2(ip),j1(ip)+jinc:j2(ip),k1(ip)+kinc:k2(ip))); nalloc = nalloc + 1
      nn = (k2(ip)-k1(ip)+1-kinc)*(j2(ip)-j1(ip)+1-jinc)*(i2(ip)-i1(ip)+1-iinc)
      allocate(mstack(nn,3)); nalloc = nalloc + 1

!     patch initialization
      iv = surf(ip)%first_node
      do k = k1(ip),k2(ip)
      do j = j1(ip),j2(ip)
      do i = i1(ip),i2(ip)
        xp(i,j,k) = xv(1,iv)
        yp(i,j,k) = xv(2,iv)
        zp(i,j,k) = xv(3,iv)
        iv = iv + 1
      end do
      end do
      end do
      mp = 0
      in_stack = 0

      SEGMENT_LOOP: do iseg = 2,nstr

        str_box(1,1) = min(xstr(iseg),xstr(iseg-1))-1e-2
        str_box(1,2) = max(xstr(iseg),xstr(iseg-1))+1e-2
        str_box(2,1) = min(ystr(iseg),ystr(iseg-1))-1e-2
        str_box(2,2) = max(ystr(iseg),ystr(iseg-1))+1e-2
        str_box(3,1) = min(zstr(iseg),zstr(iseg-1))-1e-2
        str_box(3,2) = max(zstr(iseg),zstr(iseg-1))+1e-2

        if (.not.overlapped(surf(ip)%box,str_box)) cycle SEGMENT_LOOP

!       The segment intersects this patch; switch to panel level
        ipan = surf(ip)%first_panel - 1
        do k = k1(ip)+kinc,k2(ip)
        do j = j1(ip)+jinc,j2(ip)
        INNER_PANEL_LOOP: do i = i1(ip)+iinc,i2(ip)
          ipan = ipan + 1

!         test whether the panel intersects the segment
          if (.not.overlapped(str_box,p(ipan)%bound)) cycle INNER_PANEL_LOOP

          ip1 = i-iinc ; jp1 = j-jinc ; kp1 = k-kinc
          ip3 = i      ; jp3 = j      ; kp3 = k

          if (i1(ip) == i2(ip)) then
            if (i1(ip) == 1) then
              ip2 = i ; jp2 = j   ; kp2 = k-1
              ip4 = i ; jp4 = j-1 ; kp4 = k
            else
              ip2 = i ; jp2 = j-1 ; kp2 = k
              ip4 = i ; jp4 = j   ; kp4 = k-1
            end if
          else if (j1(ip) == j2(ip)) then
            if (j1(ip) == 1) then
              ip2 = i-1 ; jp2 = j ; kp2 = k
              ip4 = i   ; jp4 = j ; kp4 = k-1
            else
              ip2 = i   ; jp2 = j ; kp2 = k-1
              ip4 = i-1 ; jp4 = j ; kp4 = k
            end if
          else if (k1(ip) == k2(ip)) then
            if (k1(ip) == 1) then
              ip2 = i   ; jp2 = j-1 ; kp2 = k
              ip4 = i-1 ; jp4 = j   ; kp4 = k
            else
              ip2 = i-1 ; jp2 = j   ; kp2 = k
              ip4 = i   ; jp4 = j-1 ; kp4 = k
            end if
          end if

!         rotate the corners down
          xp1 = p(ipan)%g(1,1)*(xp(ip1,jp1,kp1)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(yp(ip1,jp1,kp1)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zp(ip1,jp1,kp1)-p(ipan)%xmid(3))
          yp1 = p(ipan)%g(2,1)*(xp(ip1,jp1,kp1)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(yp(ip1,jp1,kp1)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zp(ip1,jp1,kp1)-p(ipan)%xmid(3))
          xp2 = p(ipan)%g(1,1)*(xp(ip2,jp2,kp2)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(yp(ip2,jp2,kp2)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zp(ip2,jp2,kp2)-p(ipan)%xmid(3))
          yp2 = p(ipan)%g(2,1)*(xp(ip2,jp2,kp2)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(yp(ip2,jp2,kp2)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zp(ip2,jp2,kp2)-p(ipan)%xmid(3))
          xp3 = p(ipan)%g(1,1)*(xp(ip3,jp3,kp3)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(yp(ip3,jp3,kp3)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zp(ip3,jp3,kp3)-p(ipan)%xmid(3))
          yp3 = p(ipan)%g(2,1)*(xp(ip3,jp3,kp3)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(yp(ip3,jp3,kp3)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zp(ip3,jp3,kp3)-p(ipan)%xmid(3))
          xp4 = p(ipan)%g(1,1)*(xp(ip4,jp4,kp4)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(yp(ip4,jp4,kp4)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zp(ip4,jp4,kp4)-p(ipan)%xmid(3))
          yp4 = p(ipan)%g(2,1)*(xp(ip4,jp4,kp4)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(yp(ip4,jp4,kp4)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zp(ip4,jp4,kp4)-p(ipan)%xmid(3))

          xs1 = p(ipan)%g(1,1)*(xstr(iseg-1)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(ystr(iseg-1)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zstr(iseg-1)-p(ipan)%xmid(3))
          ys1 = p(ipan)%g(2,1)*(xstr(iseg-1)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(ystr(iseg-1)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zstr(iseg-1)-p(ipan)%xmid(3))
          xs2 = p(ipan)%g(1,1)*(xstr(iseg)-p(ipan)%xmid(1)) &
              + p(ipan)%g(1,2)*(ystr(iseg)-p(ipan)%xmid(2)) &
              + p(ipan)%g(1,3)*(zstr(iseg)-p(ipan)%xmid(3))
          ys2 = p(ipan)%g(2,1)*(xstr(iseg)-p(ipan)%xmid(1)) &
              + p(ipan)%g(2,2)*(ystr(iseg)-p(ipan)%xmid(2)) &
              + p(ipan)%g(2,3)*(zstr(iseg)-p(ipan)%xmid(3))

          dxs = xs2 - xs1
          dys = ys2 - ys1

          call intercept2(xs1,ys1,xs2,ys2,xp1,yp1,xp2,yp2,ss,iflag)
          if (mod(iflag,2) == 1) then
            dxv = xp1 - xs1
            dyv = yp1 - ys1
            cross = dxs*dyv - dxv*dys
            if (cross > 0.0) then
              mp(ip1,jp1,kp1) = 1
            else if (cross < 0.0) then
              mp(ip1,jp1,kp1) = -1
            end if
          end if

          call intercept2(xs1,ys1,xs2,ys2,xp2,yp2,xp3,yp3,ss,iflag)
          if (mod(iflag,2) == 1) then
            dxv = xp2 - xs1
            dyv = yp2 - ys1
            cross = dxs*dyv - dxv*dys
            if (cross > 0.0) then
              mp(ip2,jp2,kp2) = 1
            else if (cross < 0.0) then
              mp(ip2,jp2,kp2) = -1
            end if
          end if

          call intercept2(xs1,ys1,xs2,ys2,xp3,yp3,xp4,yp4,ss,iflag)
          if (mod(iflag,2) == 1) then
            dxv = xp3 - xs1
            dyv = yp3 - ys1
            cross = dxs*dyv - dxv*dys
            if (cross > 0.0) then
              mp(ip3,jp3,kp3) = 1
            else if (cross < 0.0) then
              mp(ip3,jp3,kp3) = -1
            end if
          end if

          call intercept2(xs1,ys1,xs2,ys2,xp4,yp4,xp1,yp1,ss,iflag)
          if (mod(iflag,2) == 1) then
            dxv = xp4 - xs1
            dyv = yp4 - ys1
            cross = dxs*dyv - dxv*dys
            if (cross > 0.0) then
              mp(ip4,jp4,kp4) = 1
            else if (cross < 0.0) then
              mp(ip4,jp4,kp4) = -1
            end if
          end if

        end do INNER_PANEL_LOOP
        end do
        end do

      end do SEGMENT_LOOP

!     sweep cells to identify boundary/interior/exterior
      do k = k1(ip)+kinc,k2(ip)
      do j = j1(ip)+jinc,j2(ip)
      do i = i1(ip)+iinc,i2(ip)

        !set corners
        ip1 = i-iinc ; jp1 = j-jinc ; kp1 = k-kinc
        ip3 = i      ; jp3 = j      ; kp3 = k

        if (i1(ip) == i2(ip)) then
          if (i1(ip) == 1) then
            ip2 = i ; jp2 = j   ; kp2 = k-1
            ip4 = i ; jp4 = j-1 ; kp4 = k
          else
            ip2 = i ; jp2 = j-1 ; kp2 = k
            ip4 = i ; jp4 = j   ; kp4 = k-1
          end if
        else if (j1(ip) == j2(ip)) then
          if (j1(ip) == 1) then
            ip2 = i-1 ; jp2 = j ; kp2 = k
            ip4 = i   ; jp4 = j ; kp4 = k-1
          else
            ip2 = i   ; jp2 = j ; kp2 = k-1
            ip4 = i-1 ; jp4 = j ; kp4 = k
          end if
        else if (k1(ip) == k2(ip)) then
          if (k1(ip) == 1) then
            ip2 = i   ; jp2 = j-1 ; kp2 = k
            ip4 = i-1 ; jp4 = j   ; kp4 = k
          else
            ip2 = i-1 ; jp2 = j   ; kp2 = k
            ip4 = i   ; jp4 = j-1 ; kp4 = k
          end if
        end if

        mpmax = max(mp(ip1,jp1,kp1),mp(ip2,jp2,kp2),mp(ip3,jp3,kp3),mp(ip4,jp4,kp4))
        mpmin = min(mp(ip1,jp1,kp1),mp(ip2,jp2,kp2),mp(ip3,jp3,kp3),mp(ip4,jp4,kp4))

        if ((mpmax == 1).and.(mpmin == -1)) then
!         boundary cell
          mc(i,j,k) = 0
        else if (mpmax == 1) then
!         interior cell
          mc(i,j,k) = 1
          in_stack = in_stack + 1
          mstack(in_stack,1) = i
          mstack(in_stack,2) = j
          mstack(in_stack,3) = k
        else if (mpmin == -1) then
!         exterior cell
          mc(i,j,k) = -1
          in_stack = in_stack + 1
          mstack(in_stack,1) = i
          mstack(in_stack,2) = j
          mstack(in_stack,3) = k
        else
!         unknown
          mc(i,j,k) = -2
        end if
      end do
      end do
      end do

!     use i-blanking to get rid of any "off" cells
      ipan = surf(ip)%first_panel
      do k = k1(ip)+kinc,k2(ip)
      do j = j1(ip)+jinc,j2(ip)
      do i = i1(ip)+iinc,i2(ip)
        if (p(ipan)%iblank /= 1) mc(i,j,k) = -3
        ipan = ipan + 1
      end do
      end do
      end do

!     flood fill
      do while (in_stack > 0)
        i = mstack(in_stack,1)
        j = mstack(in_stack,2)
        k = mstack(in_stack,3)
        in_stack = in_stack - 1

!       plus-k
        if (k < k2(ip)) then
          if (mc(i,j,k+1) == -2) then
            mc(i,j,k+1) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i
            mstack(in_stack,2) = j
            mstack(in_stack,3) = k+1
          end if
        end if

!       minus-k
        if (k > k1(ip) + 1) then
          if (mc(i,j,k-1) == -2) then
            mc(i,j,k-1) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i
            mstack(in_stack,2) = j
            mstack(in_stack,3) = k-1
          end if
        end if

!       plus-j
        if (j < j2(ip)) then
          if (mc(i,j+1,k) == -2) then
            mc(i,j+1,k) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i
            mstack(in_stack,2) = j+1
            mstack(in_stack,3) = k
          end if
        end if

!       minus-j
        if (j > j1(ip) + 1) then
          if (mc(i,j-1,k) == -2) then
            mc(i,j-1,k) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i
            mstack(in_stack,2) = j-1
            mstack(in_stack,3) = k
          end if
        end if

!       plus-i
        if (i < i2(ip)) then
          if (mc(i+1,j,k) == -2) then
            mc(i+1,j,k) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i+1
            mstack(in_stack,2) = j
            mstack(in_stack,3) = k
          end if
        end if

!       minus-i
        if (i > i1(ip) + 1) then
          if (mc(i-1,j,k) == -2) then
            mc(i-1,j,k) = mc(i,j,k)
            in_stack = in_stack + 1
            mstack(in_stack,1) = i-1
            mstack(in_stack,2) = j
            mstack(in_stack,3) = k
          end if
        end if

      end do
            

!     change the component number of the interior/boundary/exterior
      ipan = surf(ip)%first_panel
      do k = k1(ip)+kinc,k2(ip)
      do j = j1(ip)+jinc,j2(ip)
      do i = i1(ip)+iinc,i2(ip)
        if (mc(i,j,k) == 1) then
          p(ipan)%icomp = new_interior
        else if (mc(i,j,k) == -1) then
          p(ipan)%icomp = new_exterior
        else if (mc(i,j,k) == 0) then
          p(ipan)%icomp = new_exterior
        end if
        ipan = ipan + 1
      end do
      end do
      end do

!     patch deallocation
      deallocate(mstack,mc,mp,zp,yp,xp); nalloc = nalloc - 6

    end do PATCH_LOOP

!   string deallocation
    deallocate(xstr,ystr,zstr); nalloc = nalloc - 3

  end do STRING_LOOP

  close(3)


  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in CutOutStrings:',time2-time1
  end if

  return
  end subroutine CutOutStrings

! =========================================================
  subroutine intercept2(ax,ay,bx,by,cx,cy,dx,dy,ss,iflag)
! =========================================================

  use IntrType,only: rd

  implicit none

  real(kind=rd),intent (in)  :: ax,ay,bx,by,cx,cy,dx,dy
  real(kind=rd),intent (out) :: ss
  real(kind=rd)              :: denom,rr

  integer,intent (out) :: iflag


  continue


  denom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx)
  rr    = ((ay-cy)*(dx-cx) - (ax-cx)*(dy-cy))
  ss    = ((ay-cy)*(bx-ax) - (ax-cx)*(by-ay))

  if (denom == 0.0) then
    if ((rr == 0.0).or.(ss == 0.0)) then
      iflag = 3
    else
      iflag = 2
    end if
    return
  end if

  rr = rr / denom
  ss = ss / denom

  if ((0.0 <= ss).and.(ss <= 1.0).and. &
      (0.0 <= rr).and.(rr <= 1.0)) then
    iflag = 1
  else
    iflag = 0
  end if

  return
  end subroutine intercept2

  end module CutOutStringsM
