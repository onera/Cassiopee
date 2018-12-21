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
C k6boundboxunstr
C compute the bounding box of an unstructured mesh in the local frame
C tous les elts sont parcourus ... cela revient au meme que de determiner
C les facettes externes
C =============================================================================
      SUBROUTINE k6boundboxunstr(npts, x, y, z, 
     &                           xmax, ymax, zmax, 
     &                           xmin, ymin, zmin )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb de noeuds 
      REAL_E x(0:npts-1)
      REAL_E y(0:npts-1)
      REAL_E z(0:npts-1)
      
C_OUT
      REAL_E xmax, ymax, zmax, xmin, ymin, zmin

C_LOCAL
      REAL_E xs, ys, zs
      INTEGER_E ind

C==============================================================================

      xmax = -MAXFLOAT
      ymax = -MAXFLOAT
      zmax = -MAXFLOAT
      xmin = +MAXFLOAT
      ymin = +MAXFLOAT
      zmin = +MAXFLOAT
      
      DO ind = 0, npts-1
         
         xs = x(ind)
         ys = y(ind)
         zs = z(ind)

         xmax = MAX(xmax, xs)
         ymax = MAX(ymax, ys)
         zmax = MAX(zmax, zs)
         xmin = MIN(xmin, xs)
         ymin = MIN(ymin, ys)
         zmin = MIN(zmin, zs)

      ENDDO

      END
C ===========================================================================
