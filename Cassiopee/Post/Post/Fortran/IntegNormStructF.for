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

C =============================================================================
C  compute surface integral of the field F.vect(n), coordinates 
C        and field have the same size
C        I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
C        Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
C  ============================================================================
      SUBROUTINE k6integnormstruct(ni, nj, ratio, sx, sy, sz, field, 
     &                             result)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj
      REAL_E ratio(0:ni*nj-1)
      REAL_E sx(0:(ni-1)*(nj-1)-1)
      REAL_E sy(0:(ni-1)*(nj-1)-1)
      REAL_E sz(0:(ni-1)*(nj-1)-1)      
      REAL_E field(0:ni*nj-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind1, ind2, ind3, ind4, ind
      REAL_E f1, f2, f3, f4
      REAL_E res1, res2, res3
      
C==============================================================================
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ni1 = ni-1

      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i+j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind3 + ni
            ind = i+j*ni1 
            
            f1 = ratio(ind1)*field(ind1)
            f2 = ratio(ind2)*field(ind2)
            f3 = ratio(ind3)*field(ind3)
            f4 = ratio(ind4)*field(ind4)

            res1 = res1 + sx(ind)*(f1+f2+f3+f4)
            res2 = res2 + sy(ind)*(f1+f2+f3+f4)
            res3 = res3 + sz(ind)*(f1+f2+f3+f4)
         ENDDO
      ENDDO

      result(0) = ONE_FOURTH * res1
      result(1) = ONE_FOURTH * res2
      result(2) = ONE_FOURTH * res3

      END

C==============================================================================
C Compute surface integral of the field F.vect(n), coordinates 
C       are defined in nodes and F is defined in center
C =============================================================================
      SUBROUTINE k6integnormstructnodecenter(ni, nj, ratio, sx, sy, sz,
     &                                       field, result)
      IMPLICIT NONE     

C==============================================================================
C_IN
      INTEGER_E ni, nj
      REAL_E ratio(0:ni*nj-1)
      REAL_E sx(0:ni*nj-1)
      REAL_E sy(0:ni*nj-1)
      REAL_E sz(0:ni*nj-1)
      REAL_E field(0:ni*nj-1)
C_OUT
      REAL_E result(0:2)
      
C_LOCAL
      INTEGER_E i, j, ind
      REAL_E res1, res2, res3, ri

C==============================================================================
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      DO j = 0, nj-1
         DO i = 0, ni-1
            ind = i+j*ni
            ri = ratio(ind) * field(ind)
            res1 = res1 + ri * sx(ind)
            res2 = res2 + ri * sy(ind)
            res3 = res3 + ri * sz(ind)
         ENDDO
      ENDDO
      result(0) = res1
      result(1) = res2
      result(2) = res3
      END
