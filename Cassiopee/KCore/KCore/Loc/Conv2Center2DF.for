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
C  2D 
      SUBROUTINE k6conv2center12d(ni, nj, nfld, fieldnode, fieldcenter)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"      
C==============================================================================
C_IN
      INTEGER_E ni, nj, nfld
      REAL_E fieldnode(0:ni*nj-1, nfld)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)-1, nfld)
      
C_LOCAL
      INTEGER_E i, j, n, indc
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E nic, njc
      REAL_E a
C==============================================================================
      a = ONE_FOURTH
      nic = ni-1
      njc = nj-1
      DO j = 0, nj-2
         DO i = 0, ni-2
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            DO n = 1, nfld
               fieldcenter(indc,n) = a*(fieldnode(ind0,n)+
     &              fieldnode(ind1,n)+
     &              fieldnode(ind2,n)+
     &              fieldnode(ind3,n))
            ENDDO
         ENDDO
      ENDDO
      
      END

C  ============================================================================
C convert cellnature field defined in nodes to field defined in centers
C cellN=0 (blanked or interpolated) or 1 (normal)         
C Si au moins un 0 ==> 0, sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center212d(ni, nj, fieldnode, fieldcenter)
      IMPLICIT NONE
      
C_IN
      INTEGER_E ni, nj
      REAL_E fieldnode(0:ni*nj-1)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)-1)
      
C_LOCAL
      INTEGER_E i, j, indc
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E nic, njc
      REAL_E f

      nic = ni-1
      njc = nj-1
      DO j = 0, nj-2
         DO i = 0, ni-2
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            f = MIN(fieldnode(ind0), fieldnode(ind1))
            f = MIN(f, fieldnode(ind2))
            f = MIN(f, fieldnode(ind3))
            fieldcenter(indc) = f
         ENDDO
      ENDDO

      END

C  ============================================================================
C convert cellnature field defined in nodes to field defined in centers
C cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur: Tous 0 ==> 0
C                                 Si au moins un 2 ==> 2
C                                 Sinon 1
C      Sur la frontiere: Tous 0 ==> 0
C                        Tous 2 ==> 2
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center222d(ni, nj, fieldnode, fieldcenter)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj
      REAL_E fieldnode(0:ni*nj-1)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)-1)
      
C_LOCAL
      INTEGER_E i, j, indc
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E nic, njc
      REAL_E somme

      nic = ni-1
      njc = nj-1
      DO j = 1, nj-3
         DO i = 1, ni-3
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
c$$$  prod = fieldnode(ind0)*fieldnode(ind1)
c$$$  &              *fieldnode(ind2)*fieldnode(ind3)
c$$$     &              *fieldnode(ind4)*fieldnode(ind5)
c$$$     &              *fieldnode(ind6)*fieldnode(ind7)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
c$$$                  IF ( prod == 0.) THEN
c$$$                     fieldcenter(indc) = 2.D0
c$$$  ELSE
c$$$                     fieldcenter(indc) = 1.D0
c$$$                  ENDIF
               IF ((fieldnode(ind0) .EQ. 2.) 
     &              .OR. (fieldnode(ind1) .EQ. 2.) 
     &              .OR. (fieldnode(ind2) .EQ. 2.)
     &              .OR. (fieldnode(ind3) .EQ. 2.)) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     Traitement des fontieres
C     i = 0
      DO j = 0, nj-2
         indc = j*nic
         ind0 = j*ni
         ind1 = ind0+1
         ind2 = ind0+ni
         ind3 = ind2+1
         somme = fieldnode(ind0)+fieldnode(ind1)
     &        +fieldnode(ind2)+fieldnode(ind3)
         IF (somme .EQ. 0.) THEN
            fieldcenter(indc) = 0.D0
         ELSE
            IF (somme .EQ. 8.) THEN
               fieldcenter(indc) = 2.D0
            ELSE
               fieldcenter(indc) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     i = ni-2
      IF( ni > 2) THEN
         i = ni-2
         DO j = 0, nj-2
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
               IF (somme .EQ. 8.) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      
C     j = 0
      DO i = 0, ni-2
         indc = i
         ind0 = i
         ind1 = ind0+1
         ind2 = ind0+ni
         ind3 = ind2+1
         somme = fieldnode(ind0)+fieldnode(ind1)
     &        +fieldnode(ind2)+fieldnode(ind3)
         IF (somme .EQ. 0.) THEN
            fieldcenter(indc) = 0.D0
         ELSE
            IF (somme .EQ. 8.) THEN
               fieldcenter(indc) = 2.D0
            ELSE
               fieldcenter(indc) = 1.D0
            ENDIF
         ENDIF
      ENDDO
      
      IF( nj > 2) THEN
C     j = nj-2
         j = nj-2
         DO i = 0, ni-2
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
               IF (somme .EQ. 8.) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      END
C ===== Post/Fortran/Conv2Center2DF.for =====
