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

C==============================================================================
C compute surface integral of the moment M (OM^F.vect(n)), coordinates 
C       and F have the same size
C ============================================================================
      SUBROUTINE k6integmomentnormunstruct(nbt, size, cn, cx, cy, cz, 
     &                                     ratio, xt, yt, zt, sx, sy,sz, 
     &                                     field, result)
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E cx, cy, cz
      REAL_E ratio(0:size-1)
      REAL_E xt(0:size-1)
      REAL_E yt(0:size-1)
      REAL_E zt(0:size-1)
      REAL_E sx(0:nbt-1)
      REAL_E sy(0:nbt-1)
      REAL_E sz(0:nbt-1)
      REAL_E field(0:size-1)
      
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i
      INTEGER_E ind, ind1, ind2, ind3
      REAL_E a, resx, resy, resz
      REAL_E f, f1, f2, f3
      REAL_E mx, my, mz, sx0, sy0, sz0
      REAL_E centerx, centery, centerz
      REAL_E ONE_THIRD
C
      resx = 0.D0
      resy = 0.D0
      resz = 0.D0

      ONE_THIRD = ONE/THREE

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         
         f1 =  ratio(ind1)*field(ind1)
         f2 =  ratio(ind2)*field(ind2)
         f3 =  ratio(ind3)*field(ind3)

         sx0 = sx(i)
         sy0 = sy(i)
         sz0 = sz(i)

         f = ONE_THIRD * (f1+f2+f3)
         centerx = ONE_THIRD * (xt(ind1) + xt(ind2) + xt(ind3) ) - cx
         centery = ONE_THIRD * (yt(ind1) + yt(ind2) + yt(ind3) ) - cy
         centerz = ONE_THIRD * (zt(ind1) + zt(ind2) + zt(ind3) ) - cz
       
         mx = centery*sz0 - centerz*sy0
         my = centerz*sx0 - centerx*sz0
         mz = centerx*sy0 - centery*sx0

         resx = resx + f*mx
         resy = resy + f*my
         resz = resz + f*mz
      ENDDO

      result(0) = resx
      result(1) = resy
      result(2) = resz

      END
C =============================================================================
C compute linear integral of the moment.norm  (OM^ F.n), coordinates 
C       are defined in nodes and F is defined in center, unstructured case
C  ============================================================================
      SUBROUTINE k6integmomentnormunsnodecenter(nbt, size, cn, 
     &                                          cx, cy, cz,
     &                                          ratio, xt, yt, zt, 
     &                                          sx, sy, sz, field,
     &                                          result) 

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E nbt             ! nb de triangles
      INTEGER_E  size           ! nb de noeuds
      INTEGER_E cn(0:nbt-1,3)   ! connect elt->noeud
      REAL_E cx, cy, cz         ! pt d integration des moments
      REAL_E ratio(0:nbt-1)     ! ratio aux centres
      REAL_E xt(0:size-1)       
      REAL_E yt(0:size-1)
      REAL_E zt(0:size-1)
      REAL_E sx(0:nbt-1)
      REAL_E sy(0:nbt-1)
      REAL_E sz(0:nbt-1)      
      REAL_E field(0:nbt-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E f, resx, resy, resz
      REAL_E mx, my, mz, sx0, sy0, sz0
      REAL_E centerx, centery, centerz
      REAL_E ONE_THIRD
C
      resx = 0.D0
      resy = 0.D0
      resz = 0.D0

      ONE_THIRD = ONE/THREE
C
      DO i = 0, nbt-1
         f = ratio(i)*field(i)
         
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1

         sx0 = sx(i)
         sy0 = sy(i)
         sz0 = sz(i)

         centerx = xt(ind1) + xt(ind2) + xt(ind3)
         centerx = ONE_THIRD*centerx - cx
         
         centery = yt(ind1) + yt(ind2) + yt(ind3) 
         centery = ONE_THIRD*centery - cy
         
         centerz = zt(ind1) + zt(ind2) + zt(ind3)
         centerz = ONE_THIRD*centerz -cz

         mx = centery*sz0 - centerz*sy0
         my = centerz*sx0 - centerx*sz0
         mz = centerx*sy0 - centery*sx0
         
         resx = resx + f*mx
         resy = resy + f*my
         resz = resz + f*mz

      ENDDO

      result(0) = resx
      result(1) = resy
      result(2) = resz

      END
C=========================Post/Fortran/IntegMomentNormUnstrF.for==========
