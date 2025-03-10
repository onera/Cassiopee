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
C  compute surface integral of the field F, coordinates 
C         and field have the same size
C         I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
C         Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
C  ============================================================================
      SUBROUTINE k6integstruct(ni, nj, ratio, surf, field, result)
      IMPLICIT NONE      

C==============================================================================
C_IN
      INTEGER_E ni, nj
      REAL_E ratio(0:ni*nj-1)
      REAL_E surf(0:(ni-1)*(nj-1)-1)
      REAL_E field(0:ni*nj-1)

C_OUT
      REAL_E result             

C_LOCAL
      INTEGER_E i, j,ni1
      INTEGER_E ind1, ind2, ind3, ind4, ind
      REAL_E f1, f2, f3, f4
C==============================================================================
      ni1 = ni-1
      result = 0.D0

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

            result = result + surf(ind)*(f1+f2+f3+f4)
         ENDDO
      ENDDO
      result = result/4.D0
      END

C  ============================================================================
C  compute linear integral of the field F, coordinates 
C         and field have the same size
C         I(AB) = Length(AB)*(F(A)+F(B))/2
C  ============================================================================
      SUBROUTINE k6integstruct1d(ni, ratio, length, field, result)
      IMPLICIT NONE      
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni
      REAL_E ratio(0:ni-1)
      REAL_E length(0:ni-2)
      REAL_E field(0:ni-1)

C_OUT
      REAL_E result             

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind
      REAL_E f1, f2
C==============================================================================
      
      result = 0.D0
      
      DO i = 0, ni-2
         ind1 = i
         ind2 = i+1
         ind = i
      
         f1 = ratio(ind1)*field(ind1)
         f2 = ratio(ind2)*field(ind2)
         
         result = result + length(ind)*(f1+f2)
      ENDDO

      result = ONE_HALF * result
      END
C  ============================================================================
C Compute surface integral of the field F, coordinates are defined
C         in nodes and F is defined in center, structured case
C         ni1, nj1 : dim en centres
C =============================================================================
      SUBROUTINE k6integstructnodecenter(ni1, nj1, ratio, surf, field, 
     &                                   result)
      IMPLICIT NONE
      
C_IN
      INTEGER_E ni1, nj1
      REAL_E ratio(0:ni1*nj1-1)
      REAL_E surf(0:ni1*nj1-1)
      REAL_E field(0:ni1*nj1-1)
      
C_OUT
      REAL_E result              
      
C_LOCAL
      INTEGER_E i, j, ind
      
      result = 0.D0
      DO j = 0, nj1-1
         DO i = 0, ni1-1
            ind = i + j*ni1
            result = result + ratio(ind)*surf(ind)*field(ind)
         ENDDO
      ENDDO
      
      END    
C =============================================================================
C Compute linear integral of the field F, coordinates are defined
C       in nodes and F is defined in center
C =============================================================================
      SUBROUTINE k6integstructnodecenter1d(ni, ratio, length, field, 
     &                                     result)

      IMPLICIT NONE
      
C_IN
      INTEGER_E ni
      REAL_E ratio(0:ni-1)
      REAL_E length(0:ni-1)
      REAL_E field(0:ni-1)
      
C_OUT
      REAL_E result              
      
C_LOCAL
      INTEGER_E i
      
      result = 0.D0
      DO i = 0, ni-1
         result = result + ratio(i)*length(i)*field(i)
      ENDDO
      
      END  
C ==================== Post/Fortran/IntegStructF.for ======================
