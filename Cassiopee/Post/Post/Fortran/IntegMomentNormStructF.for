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
C compute surface integral of the moment M (OM^F.vect(n)), coordinates 
C       and F have the same size
C =============================================================================
      SUBROUTINE k6integmomentnormstruct(ni, nj, cx, cy, cz, ratio, 
     &                                   xt, yt, zt, sx, sy, sz, field, 
     &                                   result)
      IMPLICIT NONE
      

C==============================================================================
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj
      REAL_E cx, cy, cz
      REAL_E ratio(0:ni*nj-1)
      REAL_E xt(0:ni*nj-1)
      REAL_E yt(0:ni*nj-1)
      REAL_E zt(0:ni*nj-1)      
      REAL_E sx(0:(ni-1)*(nj-1)-1)
      REAL_E sy(0:(ni-1)*(nj-1)-1)
      REAL_E sz(0:(ni-1)*(nj-1)-1)      
      REAL_E field(0:ni*nj-1)
      
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind, ind1, ind2, ind3, ind4
      REAL_E f, f1, f2, f3, f4
      REAL_E mx, my, mz
      REAL_E centerx, centery, centerz
      REAL_E res1, res2, res3
C
C
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ni1 = ni-1
      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i + j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind3 + ni
            ind = i + j*ni1

            f1 = ratio(ind1)*field(ind1)
            f2 = ratio(ind2)*field(ind2)
            f3 = ratio(ind3)*field(ind3)
            f4 = ratio(ind4)*field(ind4)

            f = ONE_FOURTH*(f1+f2+f3+f4)
            
            centerx = ONE_FOURTH * (xt(ind1)+xt(ind2)+xt(ind3)+xt(ind4))
            centerx = centerx - cx

            centery = ONE_FOURTH * (yt(ind1)+yt(ind2)+yt(ind3)+yt(ind4))
            centery = centery - cy

            centerz = ONE_FOURTH * (zt(ind1)+zt(ind2)+zt(ind3)+zt(ind4))
            centerz = centerz - cz

            mx = centery*sz(ind) - centerz*sy(ind)
            my = centerz*sx(ind) - centerx*sz(ind)
            mz = centerx*sy(ind) - centery*sx(ind)

            res1 = res1 + f*mx
            res2 = res2 + f*my
            res3 = res3 + f*mz
         ENDDO
      ENDDO

      result(0) = res1
      result(1) = res2
      result(2) = res3

      END
C==============================================================================
C Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
C       are defined in nodes and F is defined in center
C ============================================================================
      SUBROUTINE k6integmomentnormstructnodecenter(
     &     ni, nj, cx, cy, cz, ratio, 
     &     xt, yt, zt, sx, sy, sz, field, 
     &     result)

      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj
      REAL_E cx, cy, cz
      REAL_E ratio(0:(ni-1)*(nj-1)-1)
      REAL_E xt(0:ni*nj-1)
      REAL_E yt(0:ni*nj-1)
      REAL_E zt(0:ni*nj-1)
      REAL_E sx(0:(ni-1)*(nj-1)-1)
      REAL_E sy(0:(ni-1)*(nj-1)-1)
      REAL_E sz(0:(ni-1)*(nj-1)-1)      
      REAL_E field(0:(ni-1)*(nj-1)-1)
      
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind, ind1, ind2, ind3, ind4
      REAL_E f, sx0, sy0, sz0
      REAL_E mx, my, mz, resx, resy, resz
      REAL_E centerx, centery, centerz
C
      resx = 0.D0
      resy = 0.D0
      resz = 0.D0
      
      ni1 = ni-1
      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i+j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind3 + ni
            ind = i + j*ni1
            
            f =  ratio(ind)*field(ind)
            
            sx0 = sx(ind)
            sy0 = sy(ind)
            sz0 = sz(ind)

            centerx = xt(ind1) + xt(ind2) + xt(ind3) + xt(ind4) 
            centerx = ONE_FOURTH*centerx - cx 

            centery = yt(ind1) + yt(ind2) + yt(ind3) + yt(ind4) 
            centery = ONE_FOURTH*centery - cy 

            centerz = zt(ind1) + zt(ind2) + zt(ind3) + zt(ind4) 
            centerz = ONE_FOURTH*centerz - cz 

            mx = centery * sz0 - centerz * sy0
            my = centerz * sx0 - centerx * sz0
            mz = centerx * sy0 - centery * sx0
            
            resx = resx + f*mx
            resy = resy + f*my
            resz = resz + f*mz
         ENDDO
      ENDDO

      result(0) = resx
      result(1) = resy
      result(2) = resz

      END
C==============Post/Fortran/IntegMomentNormStructF.for=====================
