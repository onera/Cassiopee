! ================================================================
  module DefinePanelsM
! ================================================================

  implicit none

  private
  public :: DefinePanels

  contains

! ==============================================
  subroutine DefinePanels(p)
! ==============================================
! This procedure calculates most of the basic
! attributes of a panel given its vertices.
! Note that because the dummy argument is a pointer,
! this procedure must have an explicit interface.

  use ConvertPanelToPolygonM,only: ConvertPanelToPolygon
  use IntrType              ,only: rd,rs
  use my_cpu_timeM          ,only: my_cpu_time
  use Types                 ,only: panel
  use UserInput             ,only: cpufmt,min_factor,min_size
  use UserInput             ,only: showcpu,inflate
  use VertexData            ,only: xv

  implicit none

  type(panel),pointer,dimension(:) :: p
  type(panel),pointer              :: q

  real(kind=rd),dimension(3,4)     :: x
  real(kind=rd),dimension(4)       :: xx,yy,zz
  real(kind=rd),dimension(4)       :: xc,yc,zc
  real(kind=rd),dimension(3)       :: xdim
  real(kind=rd)                    :: diag1,atmp
  real(kind=rd)                    :: dx1,dy1,dz1,dx2,dy2,dz2
  real(kind=rd)                    :: tol
  real(kind=rs)                    :: time1,time2,time3,time4,time43
  integer                          :: iv,idir,ipan,ip1,ip2


  continue


! initialization
  call my_cpu_time(time1)
  time43 = 0.0_rs
  ip1 = lbound(p,1)
  ip2 = ubound(p,1)


  do ipan = ip1,ip2
    q => p(ipan)

    do iv = 1,q%num_vertices
    do idir = 1,3
      x(idir,iv) = xv(idir,q%node(iv))
      xx(iv) = xv(1,q%node(iv))
      yy(iv) = xv(2,q%node(iv))
      zz(iv) = xv(3,q%node(iv))
    end do
    end do

    if (q%num_vertices == 3) then
      xx(4) = xx(3)
      yy(4) = yy(3)
      zz(4) = zz(3)
    end if


!   calculate the panel bounding box
!   q%bound = reshape((/minval(xx),minval(yy),minval(zz), &
!                       maxval(xx),maxval(yy),maxval(zz)/),(/3,2/))
    q%bound(1,1) = minval(xx)
    q%bound(2,1) = minval(yy)
    q%bound(3,1) = minval(zz)
    q%bound(1,2) = maxval(xx)
    q%bound(2,2) = maxval(yy)
    q%bound(3,2) = maxval(zz)


!   increase the bounding box dimensions to a minimum size
    do idir = 1,3
      xdim(idir) = q%bound(idir,2) - q%bound(idir,1)
    end do

    tol = max(min_size,min_factor*sqrt(xdim(1)**2+xdim(2)**2+xdim(3)**2))

    do idir = 1,3

!     if (xdim(idir) < tol) then
!       atmp = 0.5_rd * (q%bound(idir,1) + q%bound(idir,2))
!       q%bound(idir,1) = atmp - 0.5_rd * tol
!       q%bound(idir,2) = atmp + 0.5_rd * tol
!     end if

      xdim(idir) = max(xdim(idir),tol)*inflate
      atmp = 0.5_rd * (q%bound(idir,1) + q%bound(idir,2))
      q%bound(idir,1) = atmp - 0.5_rd * xdim(idir)
      q%bound(idir,2) = atmp + 0.5_rd * xdim(idir)

    end do
    
!   calculate the unit normal based on the diagonals
    if (q%num_vertices == 4) then
      dx1 = xx(3)-xx(1)
      dy1 = yy(3)-yy(1)
      dz1 = zz(3)-zz(1)
      dx2 = xx(4)-xx(2)
      dy2 = yy(4)-yy(2)
      dz2 = zz(4)-zz(2)
    else
      dx1 = xx(2)-xx(1)
      dy1 = yy(2)-yy(1)
      dz1 = zz(2)-zz(1)
      dx2 = xx(3)-xx(1)
      dy2 = yy(3)-yy(1)
      dz2 = zz(3)-zz(1)
    end if
    q%xn(1) = 0.5*(dy1*dz2 - dy2*dz1)
    q%xn(2) = 0.5*(dz1*dx2 - dz2*dx1)
    q%xn(3) = 0.5*(dx1*dy2 - dx2*dy1)
    q%area_ref = sqrt(q%xn(1)**2 + q%xn(2)**2 + q%xn(3)**2)
    if (q%area_ref /= 0.0) q%xn = q%xn / q%area_ref


!   calculate the midpoint of the panel
    if (q%num_vertices == 4) then
      q%xmid(1) = 0.25*(xx(1)+xx(2)+xx(3)+xx(4))
      q%xmid(2) = 0.25*(yy(1)+yy(2)+yy(3)+yy(4))
      q%xmid(3) = 0.25*(zz(1)+zz(2)+zz(3)+zz(4))
    else
      q%xmid(1) = (xx(1)+xx(2)+xx(3))/3.0_rd
      q%xmid(2) = (yy(1)+yy(2)+yy(3))/3.0_rd
      q%xmid(3) = (zz(1)+zz(2)+zz(3))/3.0_rd
    end if


    if (q%num_vertices == 4) then

!     calculate the midpoint of the side of each panel
      xc(1) = 0.5*(xx(1)+xx(2))
      yc(1) = 0.5*(yy(1)+yy(2))
      zc(1) = 0.5*(zz(1)+zz(2))
      xc(2) = 0.5*(xx(2)+xx(3))
      yc(2) = 0.5*(yy(2)+yy(3))
      zc(2) = 0.5*(zz(2)+zz(3))
      xc(3) = 0.5*(xx(3)+xx(4))
      yc(3) = 0.5*(yy(3)+yy(4))
      zc(3) = 0.5*(zz(3)+zz(4))
      xc(4) = 0.5*(xx(4)+xx(1))
      yc(4) = 0.5*(yy(4)+yy(1))
      zc(4) = 0.5*(zz(4)+zz(1))


!     calculate the unit normal based on the side midpoints
      dx1 = xc(3)-xc(1)
      dy1 = yc(3)-yc(1)
      dz1 = zc(3)-zc(1)
      dx2 = xc(4)-xc(2)
      dy2 = yc(4)-yc(2)
      dz2 = zc(4)-zc(2)
      q%xn(1) = 0.5*(dy1*dz2 - dy2*dz1)
      q%xn(2) = 0.5*(dz1*dx2 - dz2*dx1)
      q%xn(3) = 0.5*(dx1*dy2 - dx2*dy1)
      atmp = sqrt(q%xn(1)**2 + q%xn(2)**2 + q%xn(3)**2)
      if (atmp /= 0.0) q%xn = q%xn / atmp
    end if

    diag1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1)


!   define the rotation matrix into the facet plane
    if (diag1 /= 0.0) then
      q%g(1,1) = dx1/diag1
      q%g(1,2) = dy1/diag1
      q%g(1,3) = dz1/diag1
    else
      q%g(1,1) = dx1
      q%g(1,2) = dy1
      q%g(1,3) = dz1
    end if
    q%g(3,1) = q%xn(1)
    q%g(3,2) = q%xn(2)
    q%g(3,3) = q%xn(3)
    q%g(2,1) = q%g(3,2)*q%g(1,3)-q%g(1,2)*q%g(3,3)
    q%g(2,2) = q%g(3,3)*q%g(1,1)-q%g(1,3)*q%g(3,1)
    q%g(2,3) = q%g(3,1)*q%g(1,2)-q%g(1,1)*q%g(3,2)


!   initialize the polygon in the facet plane
    call my_cpu_time(time3)
    call ConvertPanelToPolygon(q%g,x,q%xmid,q%poly, &
                               q%num_vertices,q%node,q%edge)
    call my_cpu_time(time4)
    time43 = time43 + (time4-time3)

  end do

  if (showcpu) then
    call my_cpu_time(time2)
    write(0,cpufmt)'cpu time in DefinePanels:',time2-time1
    write(0,cpufmt)'  > time in ConvertPanelToPolygon:',time43
  end if

  return
  end subroutine DefinePanels

  end module DefinePanelsM
