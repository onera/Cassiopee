C  
C    Copyright 2013-2026 ONERA.
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
C  convert field of floats defined in nodes to field defined in centers,
C  1D 
      SUBROUTINE k6conv2center11d(ni, nfld, fieldnode, fieldcenter)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"      
C==============================================================================
C_IN
      INTEGER_E ni, nfld
      REAL_E fieldnode(0:ni-1, nfld)
      
C_OUT
      REAL_E fieldcenter(0:ni-2, nfld)
      
C_LOCAL
      INTEGER_E i, n
      INTEGER_E ind0, ind1
      REAL_E a
C==============================================================================
      a = ONE_HALF
      DO i = 0, ni-2
         ind0 = i
         ind1 = i+1
         DO n = 1, nfld
            fieldcenter(i,n) = a*(fieldnode(ind0,n)+fieldnode(ind1,n))
         ENDDO
      ENDDO
      
      END

C  ============================================================================
C convert cellnature field defined in nodes to field defined in centers
C cellN=0 (blanked or interpolated) or 1 (normal)         
C Si au moins un 0 ==> 0, sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center211d(ni, fieldnode, fieldcenter)
      IMPLICIT NONE            
C_IN
      INTEGER_E ni
      REAL_E fieldnode(0:ni-1)
      
C_OUT
      REAL_E fieldcenter(0:ni-2)
      
C_LOCAL
      INTEGER_E i, indcenter
      INTEGER_E ind0, ind1
      INTEGER_E nic
      
      nic = ni-1
      DO i = 0, nic-1
         indcenter = i
         ind0 = i
         ind1 = i+1
         fieldcenter(indcenter) = MIN(fieldnode(ind0), fieldnode(ind1))
      ENDDO
      END

C  ============================================================================
C convert cellnature field defined in nodes to field defined in centers
C cellN=0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur: Tous 0 ==> 0
C                                 Si au moins un 2 ==> 2
C                                 Sinon 1
C      Sur la frontiere: Tous 0 ==> 0
C                        Tous 2 ==> 2
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center221d(ni, fieldnode, fieldcenter)
      IMPLICIT NONE
C_IN
      INTEGER_E ni
      REAL_E fieldnode(0:ni-1)      
C_OUT
      REAL_E fieldcenter(0:ni-2)
      
C_LOCAL
      INTEGER_E i, indcenter
      INTEGER_E ind0, ind1, nic
      REAL_E somme

      nic = ni-1
      DO i = 1, nic-1
         indcenter = i
         ind0 = i
         ind1 = ind0+1
         somme = fieldnode(ind0)+fieldnode(ind1)
         IF (somme .EQ. 0.) THEN
            fieldcenter(indcenter) = 0.D0
         ELSE
            IF ((fieldnode(ind0) .EQ. 2.) 
     &              .OR. (fieldnode(ind1) .EQ. 2.))  THEN
               fieldcenter(indcenter) = 2.D0
            ELSE
               fieldcenter(indcenter) = 1.D0
            ENDIF
         ENDIF
      ENDDO
      ! i = 1
      somme = fieldnode(0)+fieldnode(1)
      IF (somme .EQ. 0) THEN
         fieldcenter(0) = 0.D0
      ELSE IF (somme .EQ. 4.) THEN 
         fieldcenter(0) = 2.D0
      ELSE
         fieldcenter(0) = 1.D0
      ENDIF
      ! i = nic-1
      somme = fieldnode(nic)+fieldnode(nic-1)
      IF (somme .EQ. 0) THEN
         fieldcenter(nic-1) = 0.D0
      ELSE IF (somme .EQ. 4.) THEN 
         fieldcenter(nic-1) = 2.D0
      ELSE
         fieldcenter(nic-1) = 1.D0
      ENDIF
      END

C ========================= KCore/Loc/Conv2Center1DF.for ====================
