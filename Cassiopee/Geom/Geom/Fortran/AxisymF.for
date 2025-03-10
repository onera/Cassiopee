C  
C    Copyright 2013-2025 Onera.
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
C generation de maillage axisymetrique
C ============================================================================

      SUBROUTINE k6axisym(npts, x, y, z,
     &                    xc, yc, zc, nx, ny, nz,
     &                    teta, rmod, xo, yo, zo)
C
      IMPLICIT NONE
# include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! mesh size
      REAL_E x(0:npts-1)        ! mesh coordinates
      REAL_E y(0:npts-1)
      REAL_E z(0:npts-1)
      REAL_E xc,yc,zc           ! center of rotation
      REAL_E nx,ny,nz           ! rotation vector
      REAL_E teta               ! angle
      REAL_E rmod               ! modulation de r
C_OUT
      REAL_E xo(0:npts-1)       ! rotated mesh
      REAL_E yo(0:npts-1)
      REAL_E zo(0:npts-1)
C_LOCAL
      INTEGER_E ind
      REAL_E unx, uny, unz
      REAL_E norm
      REAL_E px, py, pz
      REAL_E rx, ry, rz
      REAL_E e0, e1, e2, e3
      REAL_E a1, a2, sinteta, sinteta5
      REAL_E a, c1x, c1y, c1z, c2x, c2y, c2z

C==============================================================================
C     nx,ny,nz must be unit vector
      norm = nx*nx+ny*ny+nz*nz
      IF (norm.LE.1.e-12) THEN
         WRITE(*,*) 'nx,ny,nz has null norm'
         RETURN
      ENDIF
         
      norm = 1.D0/SQRT(norm)
      unx = nx*norm
      uny = ny*norm
      unz = nz*norm
         
      sinteta = sin(teta)
      sinteta5 = sin(teta*0.5D0)
C quaternion
      e0 = cos(teta*0.5D0)
      e1 = -unx*sinteta5
      e2 = -uny*sinteta5
      e3 = -unz*sinteta5
      a1 = e0*e0-e1*e1-e2*e2-e3*e3

      DO ind = 0, npts-1
         rx = x(ind)-xc
         ry = y(ind)-yc
         rz = z(ind)-zc

         a = rx*nx+ry*ny+rz*nz
         c1x = a*nx
         c1y = a*ny
         c1z = a*nz
         c2x = rx-c1x
         c2y = ry-c1y
         c2z = rz-c1z
         rx = c1x + rmod*c2x
         ry = c1y + rmod*c2y
         rz = c1z + rmod*c2z

         a2 = e1*rx+e2*ry+e3*rz
         px = a1*rx+2*e1*a2-(ry*unz-rz*uny)*sinteta
         py = a1*ry+2*e2*a2-(rz*unx-rx*unz)*sinteta
         pz = a1*rz+2*e3*a2-(rx*uny-ry*unx)*sinteta
         xo(ind) = xc+px
         yo(ind) = yc+py
         zo(ind) = zc+pz
      ENDDO

      END
