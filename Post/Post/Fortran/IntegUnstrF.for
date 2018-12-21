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

C==============================================================================
C Compute surface integral of the field F, coordinates 
C       and field have the same size
C       I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3
C       Aire(ABC) = ||AB^AC||/2
C  ============================================================================
      SUBROUTINE k6integunstruct(nbt, size, cn, ratio, surf, field, 
     &                           result)
      IMPLICIT NONE  

C==============================================================================
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E ratio(0:size-1)
      REAL_E surf(0:nbt-1)
      REAL_E field(0:size-1)

C_OUT
      REAL_E result             

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E f1, f2, f3
C==============================================================================
      
      result = 0.D0
      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         f1 = ratio(ind1)*field(ind1)
         f2 = ratio(ind2)*field(ind2)
         f3 = ratio(ind3)*field(ind3)
         result = result+surf(i)*(f1+f2+f3)
      ENDDO
      result = result/3.D0
      END
C =============================================================================
C Compute linear integral of the field F, coordinates 
C       and field have the same size
C       I(AB) = Length(AB)*(F(A)+F(B))/3
C  ============================================================================
      SUBROUTINE k6integunstruct1d(nbt, size, cn, ratio, length, field, 
     &                             result)
      IMPLICIT NONE  

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,2)
      REAL_E ratio(0:size-1)
      REAL_E length(0:nbt-1)
      REAL_E field(0:size-1)

C_OUT
      REAL_E result             

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2
      REAL_E f1, f2
C==============================================================================
      
      result = 0.D0
      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         f1 = ratio(ind1)*field(ind1)
         f2 = ratio(ind2)*field(ind2)
         result = result+length(i)*(f1+f2)
      ENDDO
      result = ONE_HALF * result
      END

C  ============================================================================
C Compute surface integral of the field F, coordinates are defined
C       in nodes and F is defined in center, unstructured case
C =============================================================================
      SUBROUTINE k6integunsnodecenter(nbt, ratio, surf, field, 
     &                                result)
      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E nbt
      REAL_E ratio(0:nbt-1)
      REAL_E surf(0:nbt-1)
      REAL_E field(0:nbt-1)

C_OUT
      REAL_E result             

C_LOCAL
      INTEGER_E i
      REAL_E f
C------------------------------------------------------------------
      result = 0.D0
      DO i = 0, nbt-1
         f = ratio(i)*field(i)
         result = result + surf(i)*f
      ENDDO
      END

C==========================Post/Fortran/IntegUnstrF.for===================
