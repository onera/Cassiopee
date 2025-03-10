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
C k6boundboxunstr2
C compute the bounding box of an unstructured mesh in the absolute frame
C The given coordinates (x,y,z) of mesh are in relative frame.
C (n,r0,xc0) is the rotation matrix from relative to absolute frame.
C 
C tous les elts sont parcourus ... cela revient au meme que de determiner
C les facettes externes
C =============================================================================
      SUBROUTINE k6boundboxunstr2( npts, x, y, z, 
     &                             m, r0, xc0,
     &                             xmax, ymax, zmax, 
     &                             xmin, ymin, zmin )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb de noeuds 
      REAL_E x(0:npts-1)
      REAL_E y(0:npts-1)
      REAL_E z(0:npts-1)
      REAL_E m(3,3)
      REAL_E r0(3)
      REAL_E xc0(3)
      
C_OUT
      REAL_E xmax, ymax, zmax, xmin, ymin, zmin

C_LOCAL
      REAL_E xs, ys, zs, dx, dy, dz
      REAL_E xc01, xc02, xc03, r01, r02, r03
      REAL_E m11, m12, m13, m21, m22, m23, m31, m32, m33
      INTEGER_E ind
      
C==============================================================================

      xmax = -MAXFLOAT
      ymax = -MAXFLOAT
      zmax = -MAXFLOAT
      xmin = +MAXFLOAT
      ymin = +MAXFLOAT
      zmin = +MAXFLOAT

      xc01 = xc0(1)
      xc02 = xc0(2)
      xc03 = xc0(3)
      r01 = r0(1)
      r02 = r0(2)
      r03 = r0(3)
      m11 = m(1,1)
      m12 = m(1,2)
      m13 = m(1,3)
      m21 = m(2,1)
      m22 = m(2,2)
      m23 = m(2,3)
      m31 = m(3,1)
      m32 = m(3,2)
      m33 = m(3,3)
      
      DO ind = 0, npts-1
         
         dx = x(ind)-xc01
         dy = y(ind)-xc02
         dz = z(ind)-xc03

         xs = r01 +
     &        m11*dx +
     &        m12*dy +
     &        m13*dz
         ys = r02 +
     &        m21*dx +
     &        m22*dy +
     &        m23*dz
         zs = r03 +
     &        m31*dx +
     &        m32*dy +
     &        m33*dz

         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)
      ENDDO
      END
