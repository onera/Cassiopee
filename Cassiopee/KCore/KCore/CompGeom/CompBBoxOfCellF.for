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
C k6boundboxofstructcell
C compute the bounding box of a structured cell in the mesh frame 
C =============================================================================
      SUBROUTINE k6boundboxofstructcell( ind, npts, x, y, z,
     &                                   xmin,xmax,ymin,ymax,zmin,zmax)

      IMPLICIT NONE

# include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ind(0:7) ! tab of indices of the cell vertices 
      INTEGER_E npts ! size of the tab of coord
      REAL_E x(0:npts-1), y(0:npts-1), z(0:npts-1) ! coordinates of the vertices in the mesh frame
      
C_OUT
      REAL_E xmin, xmax, ymin, ymax, zmin, zmax
C_LOCAL
      INTEGER_E l
C==============================================================================
      xmax = -MAXFLOAT
      ymax = -MAXFLOAT
      zmax = -MAXFLOAT
      xmin = +MAXFLOAT
      ymin = +MAXFLOAT
      zmin = +MAXFLOAT

C     Premier sommet de la cellule (i,j,k)
      l = ind(0)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))

C     Sommet (i+1,j,k)
      l = ind(1)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))

C     Sommet (i+1,j+1,k)
      l = ind(2)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))

C     Sommet (i,j+1,k)
      l = ind(3)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))
      
C     Sommet (i,j,k+1)
      l = ind(4)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))
         
C     Sommet (i+1,j,k+1)
      l = ind(5)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))
    
C    Sommet (i+1,j+1,k+1)
      l = ind(6)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))
    
C    Sommet (i,j+1,k+1)
      l = ind(7)
      xmin = MIN(xmin, x(l))
      xmax = MAX(xmax, x(l))
      ymin = MIN(ymin, y(l))
      ymax = MAX(ymax, y(l))
      zmin = MIN(zmin, z(l))
      zmax = MAX(zmax, z(l))

      END
c====================== CompGeom/CompBBOfCellF.for ==========================
