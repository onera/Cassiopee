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
C Converti les champs en centres en noeuds (pas de prise en compte du cellN)
C =============================================================================
      SUBROUTINE k6conv2node11d(ni, nfld, fieldcenter, fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nfld
      REAL_E fieldcenter(0:ni-1, nfld)
      
C_OUT
      REAL_E fieldnode(0:ni, nfld)
      
C_LOCAL
      INTEGER_E i, n
      INTEGER_E i0
      INTEGER_E ind0, ind1
      INTEGER_E beta
      REAL_E a
C==============================================================================
      a = ONE_HALF
      DO i = 0, ni         
         beta = 1
         IF (i.EQ.0 .OR. i.EQ.ni) THEN
            beta = 0
         ENDIF
         i0 = MAX(i, 1)-1         
            
         ind0 = i0
         ind1 = ind0 + beta
         DO n = 1, nfld
            fieldnode(i,n) = a*(fieldcenter(ind0,n)+
     &           fieldcenter(ind1,n))
         ENDDO
      ENDDO

      END

C==============================================================================
C Converti les champ en centres en noeuds (prise en compte du cellN)
C le cellN peut etre mode 0 (0 et 1) ou mode 2 (0, 1, 2)
C =============================================================================
      SUBROUTINE k6conv2node11dp(ni, nfld, fieldcenter, cellN, 
     &     fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nfld
      REAL_E fieldcenter(0:ni-1,nfld)
      REAL_E cellN(0:ni-1)
      
C_OUT
      REAL_E fieldnode(0:ni,nfld)
      
C_LOCAL
      INTEGER_E i, n
      INTEGER_E i0
      INTEGER_E ind0, ind1
      INTEGER_E beta
      REAL_E w, cellN0, cellN1
C==============================================================================
      DO i = 0, ni
         beta = 1
         IF (i.EQ.0 .OR. i.EQ.ni) THEN
            beta = 0
         ENDIF
         i0 = MAX(i, 1)-1  
            
         ind0 = i0
         ind1 = ind0 + beta
         cellN0 = min(cellN(ind0), 1.)
         cellN1 = min(cellN(ind1), 1.)
         w = cellN0 + cellN1
         IF (w.EQ.0) THEN
            w = ONE_HALF
            cellN0 = 1.
            cellN1 = 1.
         ELSE
            w = 1./w
         ENDIF

         DO n = 1, nfld
            fieldnode(i,n) = w*(cellN0*fieldcenter(ind0,n)+
     &           cellN1*fieldcenter(ind1,n))
         ENDDO
      ENDDO

      END

C=============================================================================
C Converti des coords en centres en noeuds
C Extrapolation aux frontieres
C=============================================================================
      SUBROUTINE k6convcoord2node1d(ni, Fcx, Fcy, Fcz, Fnx, Fny, Fnz)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni
      REAL_E Fcx(0:ni-2)
      REAL_E Fcy(0:ni-2)
      REAL_E Fcz(0:ni-2)
C_OUT
      REAL_E Fnx(0:ni-1)
      REAL_E Fny(0:ni-1)
      REAL_E Fnz(0:ni-1)
C_LOCAL
      INTEGER_E i
C==============================================================================
      DO i = 1, ni-2
         Fnx(i) = ONE_HALF*(Fcx(i-1)+Fcx(i))
         Fny(i) = ONE_HALF*(Fcy(i-1)+Fcy(i))
         Fnz(i) = ONE_HALF*(Fcz(i-1)+Fcz(i))
      ENDDO

C  i = 0
      Fnx(0) = 2*Fnx(1)-Fnx(2)
      Fny(0) = 2*Fny(1)-Fny(2)
      Fnz(0) = 2*Fnz(1)-Fnz(2)
         
C  i = ni-1
      Fnx(ni-1) = 2*Fnx(ni-2)-Fnx(ni-3)
      Fny(ni-1) = 2*Fny(ni-2)-Fny(ni-3)
      Fnz(ni-1) = 2*Fnz(ni-2)-Fnz(ni-3)
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node211d(ni, fieldcenter, fieldnode)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni
      REAL_E fieldcenter(0:ni-1)
      
C_OUT
      REAL_E fieldnode(0:ni)
      
C_LOCAL
      INTEGER_E i, i0
      INTEGER_E beta
      INTEGER_E ind0, ind1
      REAL_E prod
C==============================================================================
      DO i = 0, ni
         beta = 1
         IF (i.EQ.0 .OR. i.EQ.ni) THEN
            beta = 0
         ENDIF
         i0 = MAX(i,1)-1        
         
         ind0 = i0
         ind1 = ind0 + beta
         prod = fieldcenter(ind0) * fieldcenter(ind1)
         IF (prod .EQ. 0) THEN
            fieldnode(i) = 0.
         ELSE
            fieldnode(i) = 1.
         ENDIF
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C      Dans le domaine interieur : Si tous 0 ==> 0
C                                  Sinon 1
C      Sur la frontiere : Si tous 0 ==> 0
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node211dp(ni, fieldcenter, fieldnode)
      IMPLICIT NONE
C_IN
      INTEGER_E ni
      REAL_E fieldcenter(0:ni-1)
      
C_OUT
      REAL_E fieldnode(0:ni)
      
C_LOCAL
      INTEGER_E i, i0
      INTEGER_E beta
      INTEGER_E ind0, ind1
      REAL_E sum
      
      DO i = 0, ni
         beta = 1
         IF (i.EQ.0 .OR. i.EQ.ni) THEN
            beta = 0
         ENDIF
         i0 = MAX(i,1)-1      
         
         ind0 = i0
         ind1 = ind0 + beta
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0) THEN
            fieldnode(i) = 0.
         ELSE
            fieldnode(i) = 1.
         ENDIF
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon si au moins un 2 ==> 2
C                                        Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon si tous 2 ==> 2
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node221d(ni, fieldcenter, fieldnode)
      IMPLICIT NONE
C_IN
      INTEGER_E ni
      REAL_E fieldcenter(0:ni-1)
      
C_OUT
      REAL_E fieldnode(0:ni)
      
C_LOCAL
      INTEGER_E i, indcenter
      INTEGER_E ind0, ind1
      REAL_E f, prod

      DO i = 1, ni-1
         ind0 = i-1
         ind1 = ind0+1
         prod = fieldcenter(ind0)*fieldcenter(ind1)
         IF (prod .EQ. 0) THEN
            fieldnode(i) = 0.D0
         ELSE
            IF((fieldcenter(ind0) .EQ. 2.) 
     &           .OR. (fieldcenter(ind1) .EQ. 2.)) THEN
               fieldnode(i) = 2.D0
            ELSE
               fieldnode(i) = 1.D0
            ENDIF
         ENDIF
      ENDDO
      
C     Traitement des fontieres
C     i = 0
      fieldnode(0) = fieldcenter(0)
C     i = ni
      fieldnode(ni-1) = fieldcenter(ni-2)
      END

C==============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon si au moins un 2 ==> 2
C                                        Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon si tous 2 ==> 2
C                         Sinon 1
C==============================================================================
      SUBROUTINE k6conv2node221dp(ni, fieldcenter, fieldnode)
      IMPLICIT NONE
C_IN
      INTEGER_E ni
      REAL_E fieldcenter(0:ni-1)
      
C_OUT
      REAL_E fieldnode(0:ni)
      
C_LOCAL
      INTEGER_E i, indcenter
      INTEGER_E ind0, ind1
      REAL_E sum

      DO i = 1, ni-1
         ind0 = i-1
         ind1 = ind0+1
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0) THEN
            fieldnode(i) = 0.D0
         ELSE
            IF((fieldcenter(ind0) .EQ. 2.) 
     &           .OR. (fieldcenter(ind1) .EQ. 2.)) THEN
               fieldnode(i) = 2.D0
            ELSE
               fieldnode(i) = 1.D0
            ENDIF
         ENDIF
      ENDDO
      
C     Traitement des fontieres
C     i = 0
      fieldnode(0) = fieldcenter(0)
C     i = ni
      fieldnode(ni-1) = fieldcenter(ni-2)

      END

C ===== Post/Fortran/Conv2Node1DF.for =====
