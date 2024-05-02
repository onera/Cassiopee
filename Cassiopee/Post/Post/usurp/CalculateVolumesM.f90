! ================================================================
  module CalculateVolumesM
! ================================================================

  implicit none

  private
  public :: CalculateVolumes

  contains


! ====================================================================
  subroutine CalculateVolumes(ni,nj,nk)
! ====================================================================
! Procedure to calculate the cell volumes.  The logic in this 
! procedure was taken from UNCLE.  The logic appears to be slightly
! different in OVERINT (FVJCC) but the results, so far, look 
! similar enough.

  use IntrType    ,only: rd
  use VolumeArrays,only: volv,xb,yb,zb

  implicit none

  integer,intent(in) :: ni,nj,nk

  real(kind=rd),dimension(ni) :: aix,aiy,aiz
  real(kind=rd),dimension(nj) :: ajx,ajy,ajz
  real(kind=rd),dimension(nk) :: akx,aky,akz

  real(kind=rd) :: dx1,dy1,dz1,dx2,dy2,dz2
  real(kind=rd) :: dxc,dyc,dzc,ax,ay,az
  real(kind=rd) :: xci,yci,zci,xcim1,ycim1,zcim1
  real(kind=rd) :: xcj,ycj,zcj,xcjm1,ycjm1,zcjm1
  real(kind=rd) :: xck,yck,zck,xckm1,yckm1,zckm1

  integer :: i,j,k


  continue


! compute the cell volumes
  do k=2,nk
  do j=2,nj
!   compute i face metrics
    do i = 1,ni
      dx1=xb(i,j  ,k)-xb(i,j-1,k-1)
      dy1=yb(i,j  ,k)-yb(i,j-1,k-1)
      dz1=zb(i,j  ,k)-zb(i,j-1,k-1)
      dx2=xb(i,j-1,k)-xb(i,j  ,k-1)
      dy2=yb(i,j-1,k)-yb(i,j  ,k-1)
      dz2=zb(i,j-1,k)-zb(i,j  ,k-1)
      aix(i)=(dy1*dz2-dz1*dy2)*0.5
      aiy(i)=(dz1*dx2-dx1*dz2)*0.5
      aiz(i)=(dx1*dy2-dy1*dx2)*0.5
    end do
    do i=2,ni
      xci=xb(i,j,k)+xb(i,j,k-1)+xb(i,j-1,k-1)+xb(i,j-1,k)
      yci=yb(i,j,k)+yb(i,j,k-1)+yb(i,j-1,k-1)+yb(i,j-1,k)
      zci=zb(i,j,k)+zb(i,j,k-1)+zb(i,j-1,k-1)+zb(i,j-1,k)
      xcim1=xb(i-1,j,k)+xb(i-1,j,k-1)+xb(i-1,j-1,k-1)+xb(i-1,j-1,k)
      ycim1=yb(i-1,j,k)+yb(i-1,j,k-1)+yb(i-1,j-1,k-1)+yb(i-1,j-1,k)
      zcim1=zb(i-1,j,k)+zb(i-1,j,k-1)+zb(i-1,j-1,k-1)+zb(i-1,j-1,k)
      dxc=xci-xcim1
      dyc=yci-ycim1
      dzc=zci-zcim1
      ax=aix(i)+aix(i-1)
      ay=aiy(i)+aiy(i-1)
      az=aiz(i)+aiz(i-1)
      volv(i,j,k)=ax*dxc+ay*dyc+az*dzc
    end do
  end do
  end do

  do k = 2,nk
  do i = 2,ni
!   calculate j face metrics
    do j = 1,nj
      dx1=xb(i,j,k)-xb(i-1,j,k-1)
      dy1=yb(i,j,k)-yb(i-1,j,k-1)
      dz1=zb(i,j,k)-zb(i-1,j,k-1)
      dx2=xb(i,j,k-1)-xb(i-1,j,k)
      dy2=yb(i,j,k-1)-yb(i-1,j,k)
      dz2=zb(i,j,k-1)-zb(i-1,j,k)
      ajx(j)=(dy1*dz2-dz1*dy2)*0.5
      ajy(j)=(dz1*dx2-dx1*dz2)*0.5
      ajz(j)=(dx1*dy2-dy1*dx2)*0.5
    end do
    do j = 2,nj
      xcj=xb(i,j,k)+xb(i,j,k-1)+xb(i-1,j,k-1)+xb(i-1,j,k)
      ycj=yb(i,j,k)+yb(i,j,k-1)+yb(i-1,j,k-1)+yb(i-1,j,k)
      zcj=zb(i,j,k)+zb(i,j,k-1)+zb(i-1,j,k-1)+zb(i-1,j,k)
      xcjm1=xb(i,j-1,k)+xb(i,j-1,k-1)+xb(i-1,j-1,k-1)+xb(i-1,j-1,k)
      ycjm1=yb(i,j-1,k)+yb(i,j-1,k-1)+yb(i-1,j-1,k-1)+yb(i-1,j-1,k)
      zcjm1=zb(i,j-1,k)+zb(i,j-1,k-1)+zb(i-1,j-1,k-1)+zb(i-1,j-1,k)
      dxc=xcj-xcjm1
      dyc=ycj-ycjm1
      dzc=zcj-zcjm1
      ax=ajx(j)+ajx(j-1)
      ay=ajy(j)+ajy(j-1)
      az=ajz(j)+ajz(j-1)
      volv(i,j,k)=volv(i,j,k)+ax*dxc+ay*dyc+az*dzc
    end do
  end do
  end do

  do j = 2,nj
  do i = 2,ni
!   compute k face metrics
    do k = 1,nk
      dx1=xb(i,j,k)-xb(i-1,j-1,k)
      dy1=yb(i,j,k)-yb(i-1,j-1,k)
      dz1=zb(i,j,k)-zb(i-1,j-1,k)
      dx2=xb(i-1,j,k)-xb(i,j-1,k)
      dy2=yb(i-1,j,k)-yb(i,j-1,k)
      dz2=zb(i-1,j,k)-zb(i,j-1,k)
      akx(k)=(dy1*dz2-dz1*dy2)*0.5
      aky(k)=(dz1*dx2-dx1*dz2)*0.5
      akz(k)=(dx1*dy2-dy1*dx2)*0.5
    end do
    do k=2,nk
      xck=xb(i,j,k)+xb(i,j-1,k)+xb(i-1,j,k)+xb(i-1,j-1,k)
      yck=yb(i,j,k)+yb(i,j-1,k)+yb(i-1,j,k)+yb(i-1,j-1,k)
      zck=zb(i,j,k)+zb(i,j-1,k)+zb(i-1,j,k)+zb(i-1,j-1,k)
      xckm1=xb(i,j,k-1)+xb(i,j-1,k-1)+xb(i-1,j,k-1)+xb(i-1,j-1,k-1)
      yckm1=yb(i,j,k-1)+yb(i,j-1,k-1)+yb(i-1,j,k-1)+yb(i-1,j-1,k-1)
      zckm1=zb(i,j,k-1)+zb(i,j-1,k-1)+zb(i-1,j,k-1)+zb(i-1,j-1,k-1)
      dxc=xck-xckm1
      dyc=yck-yckm1
      dzc=zck-zckm1
      ax=akx(k)+akx(k-1)
      ay=aky(k)+aky(k-1)
      az=akz(k)+akz(k-1)
      volv(i,j,k)=(volv(i,j,k)+ax*dxc+ay*dyc+az*dzc)/24.0
      if (volv(i,j,k) < 0.0) then
        write(unit=*,fmt=*)"negative volume located at ",i,j,k
      end if
    end do
  end do
  end do

  return
  end subroutine CalculateVolumes

  end module CalculateVolumesM
