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
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>

C  ============================================================================
C  Convert field defined on centers to field defined on nodes
C  ============================================================================
      SUBROUTINE k6conv2node12d(nj, nk, nfld, fieldcenter, fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nj, nk, nfld
      REAL_E fieldcenter(0:nj*nk-1,nfld)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1,nfld)
      
C_LOCAL
      INTEGER_E j, k, n, indnode
      INTEGER_E j0, k0
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E njn
      INTEGER_E beta, gamma
      REAL_E b
C==============================================================================
      b = ONE_FOURTH
      njn = nj+1
      DO k = 0, nk
         DO j = 0, nj
            indnode = j+k*njn
            
            beta = 1
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = nj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
               
            ind0 = j0 + k0 * nj
            ind1 = ind0 + beta
            ind2 = ind0 + gamma
            ind3 = ind2 + beta
            DO n = 1, nfld
               fieldnode(indnode,n) = b*(fieldcenter(ind0,n)+
     &              fieldcenter(ind1,n)+
     &              fieldcenter(ind2,n)+
     &              fieldcenter(ind3,n))
            ENDDO
         ENDDO
      ENDDO

      END

C  ============================================================================
C  Convert field defined on centers to field defined on nodes
C  ============================================================================
      SUBROUTINE k6conv2node12dp(nj, nk, nfld, fieldcenter, cellN,
     &     fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nj, nk, nfld
      REAL_E fieldcenter(0:nj*nk-1,nfld)
      REAL_E cellN(0:nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1, nfld)
      
C_LOCAL
      INTEGER_E j, k, n, indnode
      INTEGER_E j0, k0
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E njn
      INTEGER_E beta, gamma
      REAL_E b, w, cellN0, cellN1, cellN2, cellN3
C==============================================================================
      b = ONE_FOURTH
      njn = nj+1
      DO k = 0, nk
         DO j = 0, nj
            indnode = j+k*njn
            
            beta = 1
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = nj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
               
            ind0 = j0 + k0 * nj
            ind1 = ind0 + beta
            ind2 = ind0 + gamma
            ind3 = ind2 + beta
            cellN0 = MIN(cellN(ind0),1.)
            cellN1 = MIN(cellN(ind1),1.)
            cellN2 = MIN(cellN(ind2),1.)
            cellN3 = MIN(cellN(ind3),1.)
            w = cellN0 + cellN1 + cellN2 + cellN3
            IF (w.EQ.0) THEN
               w = b
               cellN0 = 1.
               cellN1 = 1.
               cellN2 = 1.
               cellN3 = 1.
            ELSE
               w = 1./w
            ENDIF

            DO n = 1, nfld
               fieldnode(indnode,n) = w*(cellN0*fieldcenter(ind0,n)+
     &              cellN1*fieldcenter(ind1,n)+
     &              cellN2*fieldcenter(ind2,n)+
     &              cellN3*fieldcenter(ind3,n))
            ENDDO
         ENDDO
      ENDDO

      END

C  ============================================================================
C  convert coord defined on centers to coord defined on nodes 
C  ============================================================================
      SUBROUTINE k6convcoord2node2d(nj, nk, Fcx, Fcy, Fcz,
     &     Fnx, Fny, Fnz)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"      
C_IN
      INTEGER_E nj, nk
      REAL_E Fcx(0:(nj-1)*(nk-1)-1)
      REAL_E Fcy(0:(nj-1)*(nk-1)-1)
      REAL_E Fcz(0:(nj-1)*(nk-1)-1)
C_OUT
      REAL_E Fnx(0:nj*nk-1)
      REAL_E Fny(0:nj*nk-1)
      REAL_E Fnz(0:nj*nk-1)
            
C_LOCAL
      INTEGER_E i, j, k, j0, k0, n
      INTEGER_E ind, ind0, ind1, ind2, ind3, indnode
      INTEGER_E nin, njn, ninnjn, ninj
      REAL_E a, b, c
      INTEGER_E gamma, beta

      njn = nj
      DO k = 0, nk-1
         DO j = 0, nj-1
            beta = 1
            IF (j.EQ.0 .OR. j.EQ.nj-1) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = nj-1
            IF (k.EQ.0 .OR. k.EQ.nk-1) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1

            ind = j + k*nj
            ind0 = j0 + k0 * (nj-1)
            ind1 = ind0 + beta
            ind2 = ind0 + gamma
            ind3 = ind2 + beta
            Fnx(ind) = ONE_FOURTH*(Fcx(ind0)+Fcx(ind1)+
     &           Fcx(ind2)+Fcx(ind3))
            Fny(ind) = ONE_FOURTH*(Fcy(ind0)+Fcy(ind1)+
     &           Fcy(ind2)+Fcy(ind3))
            Fnz(ind) = ONE_FOURTH*(Fcz(ind0)+Fcz(ind1)+
     &           Fcz(ind2)+Fcz(ind3))
         ENDDO
      ENDDO

C  j = 0
      DO k = 1, nk-1
         ind0 = k*nj
         ind2 = ind0+1
         ind3 = ind2+1
         Fnx(ind0) = 2*Fnx(ind2)-Fnx(ind3)
         Fny(ind0) = 2*Fny(ind2)-Fny(ind3)
         Fnz(ind0) = 2*Fnz(ind2)-Fnz(ind3)
      ENDDO
         
C  j = nj
      DO k = 1, nk-1
         ind0 = nj-1+k*nj
         ind2 = ind0-1
         ind3 = ind2-1
         Fnx(ind0) = 2*Fnx(ind2)-Fnx(ind3)
         Fny(ind0) = 2*Fny(ind2)-Fny(ind3)
         Fnz(ind0) = 2*Fnz(ind2)-Fnz(ind3)
      ENDDO

C  k = 0
      DO j = 1, nj-1
         ind0 = j
         ind2 = ind0+nj
         ind3 = ind2+nj
         Fnx(ind0) = 2*Fnx(ind2)-Fnx(ind3)
         Fny(ind0) = 2*Fny(ind2)-Fny(ind3)
         Fnz(ind0) = 2*Fnz(ind2)-Fnz(ind3)
      ENDDO
         
C  k = nk
      DO j = 1, nj-1
         ind0 = j+(nk-1)*nj
         ind2 = ind0-nj
         ind3 = ind2-nj
         Fnx(ind0) = 2*Fnx(ind2)-Fnx(ind3)
         Fny(ind0) = 2*Fny(ind2)-Fny(ind3)
         Fnz(ind0) = 2*Fnz(ind2)-Fnz(ind3)
      ENDDO


C ========Les coins ===========================================
C i = 0, j = 0, k = 0
      Fnx(0) = 2*Fnx(1)-Fnx(2)
      Fny(0) = 2*Fny(1)-Fny(2)
      Fnz(0) = 2*Fnz(1)-Fnz(2)
C i = 0, j = nj, k = 0
      Fnx(nj-1) = 2*Fnx(nj-2)-Fnx(nj-3)
      Fny(nj-1) = 2*Fny(nj-2)-Fny(nj-3)
      Fnz(nj-1) = 2*Fnz(nj-2)-Fnz(nj-3)
C i = 0, j = 0, k = nk
      Fnx((nk-1)*nj) = 2*Fnx(1+(nk-1)*nj)
     &                           -Fnx(2+(nk-1)*nj)
      Fny((nk-1)*nj) = 2*Fny(1+(nk-1)*nj)
     &                           -Fny(2+(nk-1)*nj)
      Fnz((nk-1)*nj) = 2*Fnz(1+(nk-1)*nj)
     &                           -Fnz(2+(nk-1)*nj)
C i = 0, j = nj, k = nk
      Fnx((nj-1)+(nk-1)*nj) = 2*Fnx((nj-2)+(nk-1)*nj)
     &                             -Fnx((nj-3)+(nk-1)*nj) 
      Fny((nj-1)+(nk-1)*nj) = 2*Fny((nj-2)+(nk-1)*nj)
     &                             -Fny((nj-3)+(nk-1)*nj) 
      Fnz((nj-1)+(nk-1)*nj) = 2*Fnz((nj-2)+(nk-1)*nj)
     &                             -Fnz((nj-3)+(nk-1)*nj) 

      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node212d(nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE
C_IN
      INTEGER_E nj, nk
      REAL_E fieldcenter(0:nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E j, k, indnode, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E njn, ninnjn
      INTEGER_E beta, gamma
      REAL_E prod
      REAL_E f

      njn = nj+1
      
      DO k = 0, nk
         DO j = 0, nj
            indnode = j+k*njn
            beta = 1
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j,1)-1        
            
            gamma = nj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k,1)-1
            
            ind0 = j0 + k0 * nj
            ind1 = ind0 + beta
            ind2 = ind0 + gamma
            ind3 = ind2 + beta
            prod = fieldcenter(ind0) * fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0) THEN
               fieldnode(indnode) = 0.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDDO
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C      Dans le domaine interieur : Si tous 0 ==> 0
C                                  Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node212dp(nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE
C_IN
      INTEGER_E nj, nk
      REAL_E fieldcenter(0:nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E j, k, indnode, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E njn, ninnjn
      INTEGER_E beta, gamma
      REAL_E sum
      REAL_E f

      njn = nj+1
      
      DO k = 0, nk
         DO j = 0, nj
            indnode = j+k*njn
            beta = 1
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j,1)-1
            
            gamma = nj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k,1)-1
            
            ind0 = j0 + k0 * nj
            ind1 = ind0 + beta
            ind2 = ind0 + gamma
            ind3 = ind2 + beta
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           + fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0) THEN
               fieldnode(indnode) = 0.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDDO
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon si au moins un 2 ==> 2
C                                        Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon si tous 2 ==> 2
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node222d(nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E j, k, indnode, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E njn, beta, gamma
      REAL_E f, prod

      njn = nj+1
      DO k = 1, nk-1
         DO j = 1, nj-1
            indnode = j+k*njn
            ind0 = (j-1)+(k-1)*nj
            ind1 = ind0+1
            ind2 = ind0+nj
            ind3 = ind2+1
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0) THEN
               fieldnode(indnode) = 0.D0
            ELSE
               IF((fieldcenter(ind0) .EQ. 2.) 
     &              .OR. (fieldcenter(ind1) .EQ. 2.) 
     &              .OR. (fieldcenter(ind2) .EQ. 2.)
     &              .OR. (fieldcenter(ind3) .EQ. 2.)) THEN
                  fieldnode(indnode) = 2.D0
               ELSE
                  fieldnode(indnode) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     Traitement des fontieres
C     j = 0
      DO k = 0, nk
         indnode = k*njn
         gamma = nj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         
         ind0 = k0*nj 
         ind1 = ind0+gamma                      
         prod = fieldcenter(ind0)*fieldcenter(ind1)
         IF (prod .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (prod .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     j = nj
      DO k = 0, nk
         indnode = nj+k*njn
         gamma = nj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         
         ind0 = nj-1+k0*nj 
         ind1 = ind0+gamma                      
         prod = fieldcenter(ind0)*fieldcenter(ind1)
         IF (prod .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (prod .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     k = 0
      DO j = 0, nj
         indnode = j
         beta = 1
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            beta = 0
         ENDIF
         j0 = MAX(j, 1)-1
         
         ind0 = j0 
         ind1 = ind0+beta                      
         prod = fieldcenter(ind0)*fieldcenter(ind1)
         IF (prod .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (prod .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     k = nk
      DO j = 0, nj
         indnode = j+nk*njn
         beta = 1
         IF (j.EQ.0 .OR. j.EQ.nj) THEN
            beta = 0
         ENDIF
         j0 = MAX(j, 1)-1
         
         ind0 = j0+(nk-1)*nj 
         ind1 = ind0+beta       
         prod = fieldcenter(ind0)*fieldcenter(ind1)
         IF (prod .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (prod .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur : Si au moins un 0 ==> 0
C                                  Sinon si au moins un 2 ==> 2
C                                        Sinon 1
C      Sur la frontiere : Si au moins un 0 ==> 0
C                         Sinon si tous 2 ==> 2
C                         Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node222dp(nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E j, k, indnode, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3
      INTEGER_E njn, beta, gamma
      REAL_E f, sum

      njn = nj+1
      DO k = 1, nk-1
         DO j = 1, nj-1
            indnode = j+k*njn
            ind0 = (j-1)+(k-1)*nj
            ind1 = ind0+1
            ind2 = ind0+nj
            ind3 = ind2+1
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           +fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0) THEN
               fieldnode(indnode) = 0.D0
            ELSE
               IF((fieldcenter(ind0) .EQ. 2.) 
     &              .OR. (fieldcenter(ind1) .EQ. 2.) 
     &              .OR. (fieldcenter(ind2) .EQ. 2.)
     &              .OR. (fieldcenter(ind3) .EQ. 2.)) THEN
                  fieldnode(indnode) = 2.D0
               ELSE
                  fieldnode(indnode) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     Traitement des fontieres
C     j = 0
      DO k = 0, nk
         indnode = k*njn
         gamma = nj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         
         ind0 = k0*nj 
         ind1 = ind0+gamma                      
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (sum .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     j = nj
      DO k = 0, nk
         indnode = nj+k*njn
         gamma = nj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         
         ind0 = nj-1+k0*nj 
         ind1 = ind0+gamma                      
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (sum .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     k = 0
      DO j = 0, nj
         indnode = j
         beta = 1
         IF (j.EQ.0 .OR. j.EQ.nj) THEN
            beta = 0
         ENDIF
         j0 = MAX(j, 1)-1
         
         ind0 = j0 
         ind1 = ind0+beta                      
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (sum .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

C     k = nk
      DO j = 0, nj
         indnode = j+nk*njn
         beta = 1
         IF (j.EQ.0 .OR. j.EQ.nj) THEN
            beta = 0
         ENDIF
         j0 = MAX(j, 1)-1
         
         ind0 = j0+(nk-1)*nj 
         ind1 = ind0+beta       
         sum = fieldcenter(ind0) + fieldcenter(ind1)
         IF (sum .EQ. 0.) THEN
            fieldnode(indnode) = 0.D0
         ELSE
            IF (sum .EQ. 4.) THEN
               fieldnode(indnode) = 2.D0
            ELSE
               fieldnode(indnode) = 1.D0
            ENDIF
         ENDIF
      ENDDO

      END
C ===== Post/Fortran/Conv2Node2DF.for =====
