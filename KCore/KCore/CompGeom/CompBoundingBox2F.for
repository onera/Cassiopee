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

C  ============================================================================
C k6boundbox2
C compute the bounding box of a mesh in the absolute frame
C The given coordinates (x,y,z) of mesh are in relative frame.
C (n,r0,xc0) i sthe rotation matrix from relative to absolute frame.
C cf Blk/Compose/BlkBoundingBox2F.for
C =============================================================================
      SUBROUTINE k6boundbox2( im, jm, km, x, y, z,
     &                        m, r0, xc0,
     &                        xmax, ymax, zmax, 
     &                        xmin, ymin, zmin )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E im, jm, km
      REAL_E x(0:im*jm*km-1)
      REAL_E y(0:im*jm*km-1)
      REAL_E z(0:im*jm*km-1)
      REAL_E m(3,3)
      REAL_E r0(3)
      REAL_E xc0(3)
      
C_OUT
      REAL_E xmax, ymax, zmax, xmin, ymin, zmin
C_LOCAL
      REAL_E xs, ys, zs, dx, dy, dz
      INTEGER_E ind, ind2
      INTEGER_E i, j, k
      INTEGER_E imjm, imkm, jmkm, imjmkm
C==============================================================================
      imjm = im*jm
      imkm = im*km
      jmkm = jm*km
      imjmkm = imjm*km

      xmax = -MAXFLOAT
      ymax = -MAXFLOAT
      zmax = -MAXFLOAT
      xmin = +MAXFLOAT
      ymin = +MAXFLOAT
      zmin = +MAXFLOAT

C     i = 1
      ind = 0

      DO ind2 = 1, jmkm

         dx = x(ind)-xc0(1)
         dy = y(ind)-xc0(2)
         dz = z(ind)-xc0(3)

         xs = r0(1) +
     &        m(1,1)*dx +
     &        m(1,2)*dy +
     &        m(1,3)*dz
         ys = r0(2) +
     &        m(2,1)*dx +
     &        m(2,2)*dy +
     &        m(2,3)*dz
         zs = r0(3) +
     &        m(3,1)*dx +
     &        m(3,2)*dy +
     &        m(3,3)*dz

         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)

         ind = ind+im
      ENDDO

C     i = im
      ind = im-1
      DO ind2 = 1, jmkm

         dx = x(ind)-xc0(1)
         dy = y(ind)-xc0(2)
         dz = z(ind)-xc0(3)

         xs = r0(1) +
     &        m(1,1)*dx +
     &        m(1,2)*dy +
     &        m(1,3)*dz
         ys = r0(2) +
     &        m(2,1)*dx +
     &        m(2,2)*dy +
     &        m(2,3)*dz
         zs = r0(3) +
     &        m(3,1)*dx +
     &        m(3,2)*dy +
     &        m(3,3)*dz
         
         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)
         
         ind = ind+im
      ENDDO

      j = 1
      DO k = 1, km
         DO i = 1, im
            
            ind = (i-1)+(j-1)*im+(k-1)*imjm

            dx = x(ind)-xc0(1)
            dy = y(ind)-xc0(2)
            dz = z(ind)-xc0(3)
           
            xs = r0(1) +
     &           m(1,1)*dx +
     &           m(1,2)*dy +
     &           m(1,3)*dz
            ys = r0(2) +
     &           m(2,1)*dx +
     &           m(2,2)*dy +
     &           m(2,3)*dz
            zs = r0(3) +
     &           m(3,1)*dx +
     &           m(3,2)*dy +
     &           m(3,3)*dz         

            xmax = MAX(xmax, xs)
            ymax = MAX(ymax, ys)
            zmax = MAX(zmax, zs)
            xmin = MIN(xmin, xs)
            ymin = MIN(ymin, ys)
            zmin = MIN(zmin, zs)
         ENDDO
      ENDDO
      
      j = jm
      DO k = 1, km
         DO i = 1, im
            
            ind = (i-1)+(j-1)*im+(k-1)*imjm

            dx = x(ind)-xc0(1)
            dy = y(ind)-xc0(2)
            dz = z(ind)-xc0(3)

            xs = r0(1) +
     &           m(1,1)*dx +
     &           m(1,2)*dy +
     &           m(1,3)*dz
            ys = r0(2) +
     &           m(2,1)*dx +
     &           m(2,2)*dy +
     &           m(2,3)*dz
            zs = r0(3) +
     &           m(3,1)*dx +
     &           m(3,2)*dy +
     &           m(3,3)*dz     
            
            xmax = MAX(xmax, xs)
            ymax = MAX(ymax, ys)
            zmax = MAX(zmax, zs)
            xmin = MIN(xmin, xs)
            ymin = MIN(ymin, ys)
            zmin = MIN(zmin, zs)
         ENDDO
      ENDDO

C     k = 1
      DO ind = 0, imjm-1

         dx = x(ind)-xc0(1)
         dy = y(ind)-xc0(2)
         dz = z(ind)-xc0(3)

         xs = r0(1) +
     &        m(1,1)*dx +
     &        m(1,2)*dy +
     &        m(1,3)*dz
         ys = r0(2) +
     &        m(2,1)*dx +
     &        m(2,2)*dy +
     &        m(2,3)*dz
         zs = r0(3) +
     &        m(3,1)*dx +
     &        m(3,2)*dy +
     &        m(3,3)*dz
     
         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)
      ENDDO

C     k = km
      DO ind = (km-1)*imjm, imjmkm-1
         
         xs = r0(1) +
     &        m(1,1)*(x(ind)-xc0(1)) +
     &        m(1,2)*(y(ind)-xc0(2)) +
     &        m(1,3)*(z(ind)-xc0(3))
         ys = r0(2) +
     &        m(2,1)*(x(ind)-xc0(1)) +
     &        m(2,2)*(y(ind)-xc0(2)) +
     &        m(2,3)*(z(ind)-xc0(3))
         zs = r0(3) +
     &        m(3,1)*(x(ind)-xc0(1)) +
     &        m(3,2)*(y(ind)-xc0(2)) +
     &        m(3,3)*(z(ind)-xc0(3))
     
         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)
      ENDDO

      RETURN

      END
