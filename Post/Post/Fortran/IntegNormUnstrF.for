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
C  compute surface integral of the field F.vect(n), coordinates 
C       and field have the same size
C       I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3
C       Aire(ABCD) = ||AB^AC||/2
C  ============================================================================
      SUBROUTINE k6integnormunstruct(nbt, size, cn, ratio, sx, sy, sz, 
     &                               field, result)
                     
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E ratio(0:size-1)
      REAL_E sx(0:nbt-1)
      REAL_E sy(0:nbt-1)
      REAL_E sz(0:nbt-1)
      REAL_E field(0:size-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E ONE_THIRD
      REAL_E f1, f2, f3, res1, res2, res3, sum
C==============================================================================
      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0

      ONE_THIRD = ONE/THREE

      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         f1 = ratio(ind1)*field(ind1)
         f2 = ratio(ind2)*field(ind2)
         f3 = ratio(ind3)*field(ind3)
         sum = f1 + f2 + f3
         res1 = res1 + sx(i) * sum
         res2 = res2 + sy(i) * sum
         res3 = res3 + sz(i) * sum
      ENDDO
      result(0) = ONE_THIRD * res1
      result(1) = ONE_THIRD * res2
      result(2) = ONE_THIRD * res3
      END
C  ============================================================================
C Compute surface integral of the field F, coordinates are defined
C       in nodes and F is defined in center, unstructured case
C =============================================================================
      SUBROUTINE k6integnormunsnodecenter(nbt, ratio, nsurfx, nsurfy, 
     &                                    nsurfz, field, result) 

      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E nbt
      REAL_E ratio(0:nbt-1)
      REAL_E nsurfx(0:nbt-1)
      REAL_E nsurfy(0:nbt-1)
      REAL_E nsurfz(0:nbt-1)      
      REAL_E field(0:nbt-1)

C_OUT
      REAL_E result(0:2)             

C_LOCAL
      INTEGER_E i
      REAL_E f, res1, res2, res3
C==============================================================================

      res1 = 0.D0
      res2 = 0.D0
      res3 = 0.D0
      DO i = 0, nbt-1
         f = ratio(i)*field(i)
         res1 = res1 + nsurfx(i)*f
         res2 = res2 + nsurfy(i)*f
         res3 = res3 + nsurfz(i)*f
      ENDDO
      result(0) = res1
      result(1) = res2
      result(2) = res3
      END
