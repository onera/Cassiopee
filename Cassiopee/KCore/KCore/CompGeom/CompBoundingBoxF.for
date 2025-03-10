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

C  ============================================================================
C k6boundbox
C compute the bounding box of a mesh in the frame mesh         
C =============================================================================
      SUBROUTINE k6boundbox( im, jm, km, x, y, z,
     &                       xmax, ymax, zmax, 
     &                       xmin, ymin, zmin )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E im, jm, km
      REAL_E x(0:im*jm*km-1)
      REAL_E y(0:im*jm*km-1)
      REAL_E z(0:im*jm*km-1)
      
C_OUT
      REAL_E xmax, ymax, zmax, xmin, ymin, zmin
C_LOCAL
      REAL_E xs, ys, zs
      INTEGER_E ind, ind2
      INTEGER_E i, j, k
      INTEGER_E imjm, imkm, jmkm, imjmkm
      REAL_E xminl, yminl, zminl, xmaxl, ymaxl, zmaxl
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

!$OMP PARALLEL PRIVATE(ind,ind2,xs,ys,zs,xmaxl,ymaxl,zmaxl,
!$OMP&         xminl,yminl,zminl)
      xmaxl = -MAXFLOAT
      ymaxl = -MAXFLOAT
      zmaxl = -MAXFLOAT
      xminl = +MAXFLOAT
      yminl = +MAXFLOAT
      zminl = +MAXFLOAT

C     i = 1
!$OMP DO
      DO ind2 = 1, jmkm

         ind = (ind2-1)*im
         xs = x(ind)
         ys = y(ind)
         zs = z(ind)
         
         xmaxl = MAX(xmaxl, xs)
         ymaxl = MAX(ymaxl, ys)
         zmaxl = MAX(zmaxl, zs)
         xminl = MIN(xminl, xs)
         yminl = MIN(yminl, ys)
         zminl = MIN(zminl, zs)
         
      ENDDO
!$OMP ENDDO NOWAIT

C     i = im
!$OMP DO
      DO ind2 = 1, jmkm

         ind = ind2*im -1
         xs = x(ind)
         ys = y(ind)
         zs = z(ind)
         
         xmaxl = MAX(xmaxl, xs)
         ymaxl = MAX(ymaxl, ys)
         zmaxl = MAX(zmaxl, zs)
         xminl = MIN(xminl, xs)
         yminl = MIN(yminl, ys)
         zminl = MIN(zminl, zs)
         
      ENDDO
!$OMP ENDDO NOWAIT

      j = 1
!$OMP DO
      DO k = 1, km
         DO i = 1, im
            
            ind = (i-1)+(j-1)*im+(k-1)*imjm

            xs = x(ind)
            ys = y(ind)
            zs = z(ind)

            xmaxl = MAX(xmaxl, xs)
            ymaxl = MAX(ymaxl, ys)
            zmaxl = MAX(zmaxl, zs)
            xminl = MIN(xminl, xs)
            yminl = MIN(yminl, ys)
            zminl = MIN(zminl, zs)
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT

      j = jm
!$OMP DO
      DO k = 1, km
         DO i = 1, im
            
            ind = (i-1)+(j-1)*im+(k-1)*imjm
            
            xs = x(ind)
            ys = y(ind)
            zs = z(ind)

            xmaxl = MAX(xmaxl, xs)
            ymaxl = MAX(ymaxl, ys)
            zmaxl = MAX(zmaxl, zs)
            xminl = MIN(xminl, xs)
            yminl = MIN(yminl, ys)
            zminl = MIN(zminl, zs)
         ENDDO
      ENDDO
!$OMP ENDDO NOWAIT

C     k = 1
!$OMP DO
      DO ind = 0, imjm-1
         xs = x(ind)
         ys = y(ind)
         zs = z(ind)
         
         xmaxl = MAX(xmaxl, xs)
         ymaxl = MAX(ymaxl, ys)
         zmaxl = MAX(zmaxl, zs)
         xminl = MIN(xminl, xs)
         yminl = MIN(yminl, ys)
         zminl = MIN(zminl, zs)
      ENDDO
!$OMP ENDDO NOWAIT

C     k = km
!$OMP DO
      DO ind = (km-1)*imjm, imjmkm-1
         xs = x(ind)
         ys = y(ind)
         zs = z(ind)
         
         xmaxl = MAX(xmaxl, xs)
         ymaxl = MAX(ymaxl, ys)
         zmaxl = MAX(zmaxl, zs)
         xminl = MIN(xminl, xs)
         yminl = MIN(yminl, ys)
         zminl = MIN(zminl, zs)
      ENDDO
!$OMP ENDDO NOWAIT

!$OMP CRITICAL
      xmax = MAX(xmaxl, xmax)
      ymax = MAX(ymaxl, ymax)
      zmax = MAX(zmaxl, zmax)
      xmin = MIN(xminl, xmin)
      ymin = MIN(yminl, ymin)
      zmin = MIN(zminl, zmin)
!$OMP END CRITICAL

!$OMP END PARALLEL
      RETURN

      END

c=============================================================================

