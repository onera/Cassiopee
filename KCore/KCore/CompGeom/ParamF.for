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
      SUBROUTINE k6param (stota, small, m, x0, y0, z0,               
     &                    dx, dy, dz, s0)
!****************************************************************************
!
!     compute the parametrization . F.Desbois.O-P.Jacquotte (1991).
!
! 
!        O.P JACQUOTTE  1987
!        F.MONTIGNY-RANNOU  1991
!     
!****************************************************************************
      IMPLICIT NONE
C      
      INTEGER_E m, k, kp
      REAL_E stota, small, ds
      REAL_E delx,dely,delz,qx,qy,qz,q2,g,f,e
      REAL_E x0(*),y0(*),z0(*),s0(*)
      REAL_E dx(*),dy(*),dz(*)

      REAL_E stotai

      s0(1) = 0.D0
      ds = 0.D0
!
!   degenerated edge
!
C      sma=1000.*small
C      x1m=abs(x0(1)-x0(m))
C      y1m=abs(y0(1)-y0(m))
C      z1m=abs(z0(1)-z0(m))
C      if(x1m.le.sma.and.y1m.le.sma.and.z1m.le.sma)then
C         DO k=2,m
C            s0(k)=sma
C         ENDDO
C         stota=s0(m)
C         return
C      endif
C
      DO k = 1, m-1
         kp = k+1
         delx = x0(kp)-x0(k)
         dely = y0(kp)-y0(k)
         delz = z0(kp)-z0(k)
         qx = dx(kp)+dx(k)
         qy = dy(kp)+dy(k)
         qz = dz(kp)+dz(k)
         q2 = qx**2+qy**2+qz**2
         g = delx**2+dely**2+delz**2
         f = delx*qx+dely*qy+delz*qz
         e = 8.D0-0.5D0*q2
         ds = (3.D0/e)*(DSQRT(f*f+2.D0*g*e)-f)
         s0(kp) = s0(k)+ds
      ENDDO
      stota = s0(m)
      stotai = 1.D0 / stota
      DO k = 1, m-1
         s0(k+1) = s0(k+1)*stotai
      ENDDO
      END
C=============================================================================
      SUBROUTINE k6parambar(stota, small, npts, x0, y0, z0, net,cn1,cn2,      
     &                      dx, dy, dz, s0)
      IMPLICIT NONE
C_IN
      INTEGER_E npts            ! dimension de la ligne entrante
      REAL_E stota, small
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
      REAL_E s0(0:npts-1)
C_LOCAL
      INTEGER_E k, kp, et, ind0
      REAL_E ds, stotai
      REAL_E delx,dely,delz,qx,qy,qz,q2,g,f,e

C-----------------------------------------------------------------------------
      ind0 = cn1(0)-1
      s0(ind0) = 0.D0
      ds = 0.D0
C
      DO et = 0, net-1
         k  = cn1(et)-1
         kp = cn2(et)-1
         delx = x0(kp)-x0(k)
         dely = y0(kp)-y0(k)
         delz = z0(kp)-z0(k)
         qx = dx(kp)+dx(k)
         qy = dy(kp)+dy(k)
         qz = dz(kp)+dz(k)
         q2 = qx**2+qy**2+qz**2
         g = delx**2+dely**2+delz**2
         f = delx*qx+dely*qy+delz*qz
         e = 8.D0-0.5D0*q2
         ds = (3.D0/e)*(DSQRT(f*f+2.D0*g*e)-f)
         s0(kp) = s0(k) + ds
      ENDDO

      stota = s0(npts-1)
      stotai = 1.D0 / stota
      DO k = 1, npts-1
         s0(k) = s0(k)*stotai
      ENDDO
      END
C=============================================================================
