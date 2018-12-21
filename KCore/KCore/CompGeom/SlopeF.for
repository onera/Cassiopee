C  
C    Copyright 2013-2019 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

C ============================================================================
C Compute the slope
C  ==========================================================================
      SUBROUTINE k6slope(small,big,m,x0,y0,z0,dx,dy,dz)
!*****************************************************************************
!
!     compute slopes for each points of a given edge
!         non symetric formulation at each end of the edge
! 
!        O.P JACQUOTTE  1987
!     
!*****************************************************************************
!
      IMPLICIT NONE
C
      INTEGER_E m,k,kp,km
      REAL_E small,big,dxp,dxm,dyp,dym,dzp,dzm,d2m,d2p   
      REAL_E dm,dp,ambda,cxm,cym,czm,denm,cxp,cyp,czp,denp
      REAL_E cx,cy,cz,den
      REAL_E x0(*),y0(*),z0(*)
      REAL_E dx(*),dy(*),dz(*)
!
!  degenerated edge
!
C      sma = 1000.*small
C      x1m = ABS(x0(1)-x0(m))
C      y1m = ABS(y0(1)-y0(m))
C      z1m = ABS(z0(1)-z0(m))
C      if (x1m.le.sma.and.y1m.le.sma.and.z1m.le.sma)then
C         DO i=1,m
C            dx(i) = big
C            dy(i) = big
C            dz(i) = big
C         ENDDO
C         RETURN
C      endif

!   2 pts case
      IF (m.EQ.2) THEN
         dxp = x0(2)-x0(1)
         dyp = y0(2)-y0(1)
         dzp = z0(2)-z0(1)
         den = 1.D0/DSQRT(dxp**2+dyp**2+dzp**2)
         dx(1) = dxp*den
         dy(1) = dyp*den
         dz(1) = dzp*den
         dx(2) = dx(1)
         dy(2) = dy(1)
         dz(2) = dz(1)
         RETURN
      ENDIF
!
!   first vertex of the edge
!
      dxp = x0(3)-x0(2)
      dxm = x0(2)-x0(1)
      dyp = y0(3)-y0(2)
      dym = y0(2)-y0(1)
      dzp = z0(3)-z0(2)
      dzm = z0(2)-z0(1)
      d2m = dxm*dxm+dym*dym+dzm*dzm
      d2p = dxp*dxp+dyp*dyp+dzp*dzp
      dm = DSQRT(d2m)
      dp = DSQRT(d2p)
      ambda = (1.D0+dp/dm)**2
      cxm = ambda*dxm-dxm-dxp
      cym = ambda*dym-dym-dyp
      czm = ambda*dzm-dzm-dzp
      denm = 1.D0/DSQRT(cxm**2+cym**2+czm**2)
      dx(1) = cxm*denm
      dy(1) = cym*denm
      dz(1) = czm*denm

C
C  last vertex of the edge
C
      dxm=x0(m-1)-x0(m-2)
      dxp=x0(m)-x0(m-1)
      dym=y0(m-1)-y0(m-2)
      dyp=y0(m)-y0(m-1)
      dzm=z0(m-1)-z0(m-2)
      dzp=z0(m)-z0(m-1)
      d2m=dxm*dxm+dym*dym+dzm*dzm
      d2p=dxp*dxp+dyp*dyp+dzp*dzp
      dm=DSQRT(d2m)
      dp=DSQRT(d2p)
      ambda=(1.D0+dm/dp)**2
      cxp=ambda*dxp-dxm-dxp
      cyp=ambda*dyp-dym-dyp
      czp=ambda*dzp-dzm-dzp
      denp=1.D0/DSQRT(cxp**2+cyp**2+czp**2)
      dx(m)=cxp*denp
      dy(m)=cyp*denp
      dz(m)=czp*denp
C
C  other points
C
      DO k=2,m-1
         kp = k+1
         km = k-1
         dxm = x0(k)-x0(km)
         dxp = x0(kp)-x0(k)
         dym = y0(k)-y0(km)
         dyp = y0(kp)-y0(k)
         dzm = z0(k)-z0(km)
         dzp = z0(kp)-z0(k)
         d2m = dxm*dxm+dym*dym+dzm*dzm
         d2p = dxp*dxp+dyp*dyp+dzp*dzp
         ambda = d2p/(d2m+d2p)
         cx=(1.D0-ambda)*dxm+ambda*dxp
         cy=(1.D0-ambda)*dym+ambda*dyp
         cz=(1.D0-ambda)*dzm+ambda*dzp
         den=1.D0/DSQRT(cx**2+cy**2+cz**2)
         dx(k)=cx*den
         dy(k)=cy*den
         dz(k)=cz*den

      ENDDO
      END
C=============================================================================
      SUBROUTINE k6slopebar(small, big, npts, x0, y0, z0, net, cn1, cn2, 
     &                      dx, dy, dz)

      IMPLICIT NONE
C_IN
      INTEGER_E small, big
      INTEGER_E npts            ! dimension de la ligne entrante
      REAL_E x0(0:npts-1)       ! ligne entrante
      REAL_E y0(0:npts-1)       ! ligne entrante
      REAL_E z0(0:npts-1)       ! ligne entrante
      INTEGER_E net             ! nb d elts BAR
      INTEGER_E cn1(0:net-1)    ! connectivite BAR : 1ers sommets
      INTEGER_E cn2(0:net-1)    ! connectivite BAR : 2nds sommets
      
C_OUT
      REAL_E dx(0:npts-1)
      REAL_E dy(0:npts-1)
      REAL_E dz(0:npts-1)
C_LOCAL
      INTEGER_E k, kp, km
      INTEGER_E ind1, ind2, ind3 
      REAL_E dxp, dxm, dyp, dym, dzp, dzm, d2m, d2p   
      REAL_E dm, dp, ambda, cxm, cym, czm, denm, cxp, cyp, czp, denp
      REAL_E cx, cy, cz, den      
      INTEGER_E net1, et
C-----------------------------------------------------------------------------
C
      IF (net.EQ.1) THEN
         ind1 = cn1(0)-1; ind2 = cn2(0)-1
         dxp = x0(ind2)-x0(ind1)
         dyp = y0(ind2)-y0(ind1)
         dzp = z0(ind2)-z0(ind1)
         den = 1.D0/DSQRT(dxp**2+dyp**2+dzp**2)
         dx(0) = dxp*den
         dy(0) = dyp*den
         dz(0) = dzp*den
         dx(1) = dx(0)
         dy(1) = dy(0)
         dz(1) = dz(0)
         RETURN
      ENDIF

      net1 = net-1
C   first vertex of the edge
C
      ind1 = cn1(0)-1; ind2 = cn2(0)-1; ind3 = cn2(1)-1
      dxp = x0(ind3)-x0(ind2)
      dxm = x0(ind2)-x0(ind1)
      dyp = y0(ind3)-y0(ind2)
      dym = y0(ind2)-y0(ind1)
      dzp = z0(ind3)-z0(ind2)
      dzm = z0(ind2)-z0(ind1)
      d2m = dxm*dxm+dym*dym+dzm*dzm
      d2p = dxp*dxp+dyp*dyp+dzp*dzp
      dm = DSQRT(d2m)
      dp = DSQRT(d2p)
      ambda = (1.D0+dp/dm)**2
      cxm = ambda*dxm-dxm-dxp
      cym = ambda*dym-dym-dyp
      czm = ambda*dzm-dzm-dzp
      denm = 1.D0/DSQRT(cxm**2+cym**2+czm**2)
      dx(ind1) = cxm*denm
      dy(ind1) = cym*denm
      dz(ind1) = czm*denm
C
C  last vertex of the edge
C
      ind1 = cn1(net1-1)-1; ind2 = cn1(net1)-1; ind3 = cn2(net1)-1
      dxm = x0(ind2)-x0(ind1)
      dxp = x0(ind3)-x0(ind2)
      dym = y0(ind2)-y0(ind1)
      dyp = y0(ind3)-y0(ind2)
      dzm = z0(ind2)-z0(ind1)
      dzp = z0(ind3)-z0(ind2)
      d2m = dxm*dxm+dym*dym+dzm*dzm
      d2p = dxp*dxp+dyp*dyp+dzp*dzp
      dm = DSQRT(d2m)
      dp = DSQRT(d2p)
      ambda = (1.D0+dm/dp)**2
      cxp = ambda*dxp-dxm-dxp
      cyp = ambda*dyp-dym-dyp
      czp = ambda*dzp-dzm-dzp
      denp = 1.D0/DSQRT(cxp**2+cyp**2+czp**2)
      dx(ind3) = cxp*denp
      dy(ind3) = cyp*denp
      dz(ind3) = czp*denp

C
C  other points
C     
      DO et = 0, net1-1
         km = cn1(et)-1
         k  = cn2(et)-1
         kp = cn2(et+1)-1
         dxm=x0(k)-x0(km)
         dxp=x0(kp)-x0(k)
         dym=y0(k)-y0(km)
         dyp=y0(kp)-y0(k)
         dzm=z0(k)-z0(km)
         dzp=z0(kp)-z0(k)        
         d2m=dxm*dxm+dym*dym+dzm*dzm
         d2p=dxp*dxp+dyp*dyp+dzp*dzp
         ambda=d2p/(d2m+d2p)
         cx=(1.D0-ambda)*dxm+ambda*dxp
         cy=(1.D0-ambda)*dym+ambda*dyp
         cz=(1.D0-ambda)*dzm+ambda*dzp
         den=1.D0/DSQRT(cx**2+cy**2+cz**2)
         dx(k)=cx*den
         dy(k)=cy*den
         dz(k)=cz*den
      ENDDO

      END
c===================   SlopeF.for   ==========================================
