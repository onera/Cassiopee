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

C =============================================================================
C  compute surface integral of the field F.vect(n), coordinates 
C       and field have the same size
C       I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
C       Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
C  ============================================================================
      SUBROUTINE k6integnormprodstruct(ni, nj, ratio, sx, sy, sz,
     &                                 vx, vy, vz, result)
      IMPLICIT NONE
      
# include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj          ! dimensions de la grille surfacique
      REAL_E ratio(0:ni*nj-1)   ! ratio aux noeuds
      REAL_E sx(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %x
      REAL_E sy(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %y
      REAL_E sz(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %z
      REAL_E vx(0:ni*nj-1) ! vecteur a integrer %x 
      REAL_E vy(0:ni*nj-1) ! vecteur a integrer %y
      REAL_E vz(0:ni*nj-1) ! vecteur a integrer %z 

C_OUT
      REAL_E result             
C_LOCAL
      INTEGER_E i, j, ni1
      INTEGER_E ind1, ind2, ind3, ind4, ind
      REAL_E fx, fy, fz
      REAL_E r1, r2, r3, r4
C
      ni1 = ni-1
      result = 0.D0
      
      DO j = 0, nj-2
         DO i = 0, ni-2
            ind1 = i + j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind3 + ni
            
            ind = i + j* ni1

            r1 = ratio(ind1)
            r2 = ratio(ind2)
            r3 = ratio(ind3)
            r4 = ratio(ind4)
            fx = r1*vx(ind1) + r2*vx(ind2) +  r3*vx(ind3) + r4*vx(ind4) 
            fy = r1*vy(ind1) + r2*vy(ind2) +  r3*vy(ind3) + r4*vy(ind4) 
            fz = r1*vz(ind1) + r2*vz(ind2) +  r3*vz(ind3) + r4*vz(ind4) 
            result = result + sx(ind)*fx + sy(ind)*fy + sz(ind)*fz
         ENDDO
      ENDDO

      result = ONE_FOURTH * result

      END
C==============================================================================
C Compute surface integral of the product vect(F).vect(n), coordinates 
C       are defined in nodes and F is defined in center
C  ============================================================================
      SUBROUTINE k6integnormprodstructnodecenter(ni, nj, ratio, 
     &                                           sx, sy, sz, vx, vy, vz, 
     &                                           result)
      IMPLICIT NONE
      
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj
      REAL_E ratio(0:ni*nj-1)
      REAL_E sx(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %x
      REAL_E sy(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %y
      REAL_E sz(0:(ni-1)*(nj-1)-1) !surface ( aux centres ) %z
      REAL_E vx(0:ni*nj-1) ! vecteur a integrer %x 
      REAL_E vy(0:ni*nj-1) ! vecteur a integrer %y
      REAL_E vz(0:ni*nj-1) ! vecteur a integrer %z 
    
C_OUT
      REAL_E result              
      
C_LOCAL
      INTEGER_E i, j, ind
      REAL_E sum
C
      result = 0.D0
      DO j = 0, nj-1
         DO i = 0, ni-1
            ind = i+j*ni
            sum = sx(ind)*vx(ind) + sy(ind)*vy(ind) + sz(ind)*vz(ind)
            result = result + ratio(ind) * sum
         ENDDO
      ENDDO
      
      END
