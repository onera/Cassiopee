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
C compute surface integral of the field F.vect(n), coordinates 
C       and field have the same size
C       I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3
C       Aire(ABC) = ||AB^AC||/2
C ============================================================================
      SUBROUTINE k6integnormprodunstruct(nbt, size, cn, ratio, 
     &                                   sx, sy, sz, vx, vy, vz,
     &                                   result)
      IMPLICIT NONE
      
# include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nbt, size
      INTEGER_E cn(0:nbt-1,3)
      REAL_E ratio(0:size-1)
      REAL_E sx(0:nbt-1) !surface ( aux centres ) %x
      REAL_E sy(0:nbt-1) !surface ( aux centres ) %y
      REAL_E sz(0:nbt-1) !surface ( aux centres ) %z
      REAL_E vx(0:size-1) ! vecteur a integrer %x 
      REAL_E vy(0:size-1) ! vecteur a integrer %y
      REAL_E vz(0:size-1) ! vecteur a integrer %z 
C_OUT
      REAL_E result             
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E fx, fy, fz
      REAL_E r1, r2, r3
      REAL_E ONE_THIRD 
C      
      ONE_THIRD = ONE/THREE

      result = 0.D0
      DO i = 0, nbt-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1

         r1 = ratio(ind1)
         r2 = ratio(ind2)
         r3 = ratio(ind3)
         fx = r1*vx(ind1) + r2*vx(ind2) +  r3*vx(ind3) 
         fy = r1*vy(ind1) + r2*vy(ind2) +  r3*vy(ind3)
         fz = r1*vz(ind1) + r2*vz(ind2) +  r3*vz(ind3)
         result = result + sx(i)*fx + sy(i)*fy + sz(i)*fz
         
      ENDDO
      result = ONE_THIRD * result
      END

C  ============================================================================
C Compute surface integral of the field F, coordinates are defined
C       in nodes and F is defined in center, unstructured case
C =============================================================================
      SUBROUTINE k6integnormprodunsnodecenter(nbt, ratio, sx, sy, sz,
     &                                        vx, vy, vz, result)
      IMPLICIT NONE
C_IN
      INTEGER_E nbt
      REAL_E ratio(0:nbt-1)
      REAL_E sx(0:nbt-1) !surface ( aux centres ) %x
      REAL_E sy(0:nbt-1) !surface ( aux centres ) %y
      REAL_E sz(0:nbt-1) !surface ( aux centres ) %z
      REAL_E vx(0:nbt-1) ! vecteur a integrer %x 
      REAL_E vy(0:nbt-1) ! vecteur a integrer %y
      REAL_E vz(0:nbt-1) ! vecteur a integrer %z 

C_OUT
      REAL_E result            
C_LOCAL
      INTEGER_E i
      REAL_E ri, sum
C    
      result = 0.D0
      DO i = 0, nbt-1
         ri = ratio(i)
         sum = sx(i)*vx(i) + sy(i)*vy(i) + sz(i)*vz(i)
         result = result + ri * sum
      ENDDO
      END
C====================Post/Fortran/IntegNormProdUnstrF.for=================
 
