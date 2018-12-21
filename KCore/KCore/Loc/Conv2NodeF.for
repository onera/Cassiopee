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
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>

C  ============================================================================
C  convert field defined on centers to field defined on nodes
C  ============================================================================
      SUBROUTINE k6conv2node1(ni, nj, nk, nfld, fieldcenter, fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk, nfld
      REAL_E fieldcenter(0:ni*nj*nk-1,nfld)
      
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1,nfld)
      
C_LOCAL
      INTEGER_E i, j, k, n, indcenter, indn
      INTEGER_E i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nin, njn, ninnjn, ninj
      INTEGER_E alpha, beta, gamma
      REAL_E a
C==============================================================================
      a = ONE_EIGHT
      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i, indn, alpha, i0, ind0, ind1, ind2, ind3, 
!$OMP& ind4, ind5, ind6, ind7, n, gamma, k, k0, beta, j0)
      DO n = 1, nfld
      DO k = 0, nk
         gamma = ninj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         DO j = 0, nj
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1  

!$OMP DO
            DO i = 0, ni
               indn = i+j*nin+k*ninnjn
               alpha = 1
               IF (i .EQ. 0 .OR. i .EQ. ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1)-1               
               
               ind0 = i0 + j0 * ni + k0 * ninj
               ind1 = ind0 + alpha
               ind2 = ind0 + beta
               ind3 = ind2 + alpha
               ind4 = ind0 + gamma
               ind5 = ind4 + alpha
               ind6 = ind4 + beta
               ind7 = ind6 + alpha
               fieldnode(indn,n) = a*(fieldcenter(ind0,n)+
     &              fieldcenter(ind1,n)+
     &              fieldcenter(ind2,n)+
     &              fieldcenter(ind3,n)+
     &              fieldcenter(ind4,n)+
     &              fieldcenter(ind5,n)+
     &              fieldcenter(ind6,n)+
     &              fieldcenter(ind7,n))
            ENDDO !i
!$OMP ENDDO NOWAIT
         ENDDO !j
      ENDDO  !k
      ENDDO  !nfld
!$OMP END PARALLEL
      END

C  ============================================================================
C  convert field defined on centers to field defined on nodes
C  ============================================================================
      SUBROUTINE k6conv2node1p(ni, nj, nk, nfld, fieldcenter, cellN,
     &     fieldnode)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk, nfld
      REAL_E fieldcenter(0:ni*nj*nk-1,nfld)
      REAL_E cellN(0:ni*nj*nk-1)
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1,nfld)
      
C_LOCAL
      INTEGER_E i, j, k, n, indcenter, indn
      INTEGER_E i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nin, njn, ninnjn, ninj
      INTEGER_E alpha, beta, gamma
      REAL_E a, cellN0, cellN1, cellN2, cellN3, cellN4, cellN5
      REAL_E cellN6, cellN7, w
C==============================================================================
      a = ONE_EIGHT
      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i,indn,alpha,i0,ind0,ind1,ind2,ind3,ind4,
!$OMP& ind5,ind6,ind7,cellN0,cellN1,cellN2,cellN3,cellN4,cellN5,
!$OMP& cellN6,cellN7,w)
      DO k = 0, nk
         gamma = ninj
         IF (k.EQ.0 .OR. k.EQ.nk) THEN
            gamma = 0
         ENDIF
         k0 = MAX(k, 1)-1
         DO j = 0, nj
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1
!$OMP DO     
            DO i = 0, ni
               indn = i+j*nin+k*ninnjn
               alpha = 1
               IF (i .EQ. 0 .OR. i .EQ. ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1)-1         

               ind0 = i0 + j0 * ni + k0 * ninj
               ind1 = ind0 + alpha
               ind2 = ind0 + beta
               ind3 = ind2 + alpha
               ind4 = ind0 + gamma
               ind5 = ind4 + alpha
               ind6 = ind4 + beta
               ind7 = ind6 + alpha
               cellN0 = MIN(cellN(ind0), 1.)
               cellN1 = MIN(cellN(ind1), 1.)
               cellN2 = MIN(cellN(ind2), 1.)
               cellN3 = MIN(cellN(ind3), 1.)
               cellN4 = MIN(cellN(ind4), 1.)
               cellN5 = MIN(cellN(ind5), 1.)
               cellN6 = MIN(cellN(ind6), 1.)
               cellN7 = MIN(cellN(ind7), 1.)
               w = cellN0 + cellN1 + cellN2 + cellN3 + cellN4 + cellN5 +
     &              cellN6 + cellN7
               IF (w.EQ.0) THEN
                  w = a
                  cellN0 = 1.
                  cellN1 = 1.
                  cellN2 = 1.
                  cellN3 = 1.
                  cellN4 = 1.
                  cellN5 = 1.
                  cellN6 = 1.
                  cellN7 = 1.
               ELSE
                  w = 1./w
               ENDIF
               DO n = 1, nfld
                  fieldnode(indn,n) = w*(cellN0*fieldcenter(ind0,n)+
     &                 cellN1*fieldcenter(ind1,n)+
     &                 cellN2*fieldcenter(ind2,n)+
     &                 cellN3*fieldcenter(ind3,n)+
     &                 cellN4*fieldcenter(ind4,n)+
     &                 cellN5*fieldcenter(ind5,n)+
     &                 cellN6*fieldcenter(ind6,n)+
     &                 cellN7*fieldcenter(ind7,n))
               ENDDO
            ENDDO
!$OMP ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL
      END
C  ============================================================================
C  convert coord defined on centers to coord defined on nodes 
C  (only on boundaries)
C  ============================================================================
      SUBROUTINE k6convcoord2node(ni, nj, nk, 
     &                            Fcx, Fcy, Fcz, X, Y, Z, Xp, Yp, Zp,
     &                            Fnx, Fny, Fnz)
      IMPLICIT NONE
      
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj, nk
      REAL_E Fcx(0:ni*nj*nk-1)
      REAL_E Fcy(0:ni*nj*nk-1)
      REAL_E Fcz(0:ni*nj*nk-1)
      
C_OUT
      REAL_E Fnx(0:(ni+1)*(nj+1)*(nk+1)-1)
      REAL_E Fny(0:(ni+1)*(nj+1)*(nk+1)-1)
      REAL_E Fnz(0:(ni+1)*(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, i0, j0, k0, n
      INTEGER_E ind, ind0, ind1, ind2, ind3, indn
      INTEGER_E nin, njn, ninnjn, ninj
      REAL_E a, b, c
      REAL_E X(0:ni*nj*nk-1)
      REAL_E Y(0:ni*nj*nk-1)
      REAL_E Z(0:ni*nj*nk-1)
      REAL_E Xp(0:ni+nj+nk-1)
      REAL_E Yp(0:ni+nj+nk-1)
      REAL_E Zp(0:ni+nj+nk-1)
      a = ONE_EIGHT
      b = ONE_HALF
      c = ONE_FOURTH
      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
C     Traitement des fontieres
C=============== DEBUT i = 0 =================================

C     i = 0, interieur
      DO k = 0, nk-1
         DO j = 0, nj-1
            ind = j+k*nj
            ind0 = j*ni+k*ninj
            ind1 = ind0+1
            X(ind) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO
C OK
      DO k = 0, nk-2
         DO j = 0, nj-2
            indn = (j+1)*nin+(k+1)*ninnjn
            ind0 = j+k*nj
            ind1 = ind0+1
            ind2 = ind0+nj
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO
C OK

C i = 0, j = 0

      DO k = 0, nk-1
         ind0 = k*nj
         ind1 = 1+k*nj
         Xp(k) = b*(3*X(ind0)-X(ind1))
         Yp(k) = b*(3*Y(ind0)-Y(ind1))
         Zp(k) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO k = 0, nk-2
         indn = (k+1)*ninnjn
         Fnx(indn) = b*(Xp(k)+Xp(k+1))
         Fny(indn) = b*(Yp(k)+Yp(k+1))
         Fnz(indn) = b*(Zp(k)+Zp(k+1))
       ENDDO
C OK

C i = 0, j = nj

      DO k = 0, nk-1
         ind0 = (nj-1)+k*nj
         ind1 = (nj-2)+k*nj
         Xp(k) = b*(3*X(ind0)-X(ind1))
         Yp(k) = b*(3*Y(ind0)-Y(ind1))
         Zp(k) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO k = 0, nk-2
         indn = nj*nin+(k+1)*ninnjn
         Fnx(indn) = b*(Xp(k)+Xp(k+1))
         Fny(indn) = b*(Yp(k)+Yp(k+1))
         Fnz(indn) = b*(Zp(k)+Zp(k+1))
       ENDDO
C OK

C i = 0, k = 0

      DO j = 0, nj-1
         ind0 = j
         ind1 = j+nj
         Xp(j) = b*(3*X(ind0)-X(ind1))
         Yp(j) = b*(3*Y(ind0)-Y(ind1))
         Zp(j) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO j = 0, nj-2
         indn = (j+1)*nin
         Fnx(indn) = b*(Xp(j)+Xp(j+1))
         Fny(indn) = b*(Yp(j)+Yp(j+1))
         Fnz(indn) = b*(Zp(j)+Zp(j+1))
       ENDDO
C OK

C i = 0, k = nk

      DO j = 0, nj-1
         ind0 = j+(nk-1)*nj
         ind1 = j+(nk-2)*nj
         Xp(j) = b*(3*X(ind0)-X(ind1))
         Yp(j) = b*(3*Y(ind0)-Y(ind1))
         Zp(j) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO j = 0, nj-2
         indn = (j+1)*nin+nk*ninnjn
         Fnx(indn) = b*(Xp(j)+Xp(j+1))
         Fny(indn) = b*(Yp(j)+Yp(j+1))
         Fnz(indn) = b*(Zp(j)+Zp(j+1))
      ENDDO
C OK

C=============== FIN i = 0 =================================
C OK

C=============== DEBUT j = 0 =================================

C     j = 0, interieur

      DO k = 0, nk-1
         DO i = 0, ni-1
            ind = i+k*ni
            ind0 = i+k*ninj
            ind1 = ind0+ni
            X(ind) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO
      DO k = 0, nk-2
         DO i = 0, ni-2
            indn = (i+1)+(k+1)*ninnjn
            ind0 = i+k*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO

C j = 0, i = ni

      DO k = 0, nk-1
         ind0 = (ni-1)+k*ni
         ind1 = (ni-2)+k*ni
         Xp(k) = b*(3*X(ind0)-X(ind1))
         Yp(k) = b*(3*Y(ind0)-Y(ind1))
         Zp(k) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO k = 0, nk-2
         indn = ni+(k+1)*ninnjn
         Fnx(indn) = b*(Xp(k)+Xp(k+1))
         Fny(indn) = b*(Yp(k)+Yp(k+1))
         Fnz(indn) = b*(Zp(k)+Zp(k+1))
      ENDDO

C j = 0, k = 0

      DO i = 0, ni-1
         ind0 = i
         ind1 = i+ni
         Xp(i) = b*(3*X(ind0)-X(ind1))
         Yp(i) = b*(3*Y(ind0)-Y(ind1))
         Zp(i) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO i = 0, ni-2
         indn = (i+1)
         Fnx(indn) = b*(Xp(i)+Xp(i+1))
         Fny(indn) = b*(Yp(i)+Yp(i+1))
         Fnz(indn) = b*(Zp(i)+Zp(i+1))
      ENDDO

C OK 
C j = 0, k = nk

      DO i = 0, ni-1
         ind0 = i+(nk-1)*ni
         ind1 = i+(nk-2)*ni
         Xp(i) = b*(3*X(ind0)-X(ind1))
         Yp(i) = b*(3*Y(ind0)-Y(ind1))
         Zp(i) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO i = 0, ni-2
         indn = (i+1)+nk*ninnjn
         Fnx(indn) = b*(Xp(i)+Xp(i+1))
         Fny(indn) = b*(Yp(i)+Yp(i+1))
         Fnz(indn) = b*(Zp(i)+Zp(i+1))
      ENDDO
C OK

C=============== FIN j = 0 =================================

C=============== DEBUT k = 0 =================================

C     k = 0, interieur

      DO j = 0, nj-1
         DO i = 0, ni-1
            ind0 = i+j*ni
            ind1 = ind0+ninj
            X(ind0) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind0) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind0) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO
      DO j = 0, nj-2
         DO i = 0, ni-2
            indn = i+1+(j+1)*nin
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO
C OK

C k = 0, i = ni

      DO j = 0, nj-1
         ind0 = (ni-1)+j*ni
         ind1 = (ni-2)+j*ni
         Xp(j) = b*(3*X(ind0)-X(ind1))
         Yp(j) = b*(3*Y(ind0)-Y(ind1))
         Zp(j) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO j = 0, nj-2
         indn = ni+(j+1)*nin
         Fnx(indn) = b*(Xp(j)+Xp(j+1))
         Fny(indn) = b*(Yp(j)+Yp(j+1))
         Fnz(indn) = b*(Zp(j)+Zp(j+1))
       ENDDO
C OK

C k = 0, j = nj

      DO i = 0, ni-1
         ind0 = i+(nj-1)*ni
         ind1 = i+(nj-2)*ni
         Xp(i) = b*(3*X(ind0)-X(ind1))
         Yp(i) = b*(3*Y(ind0)-Y(ind1))
         Zp(i) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO i = 0, ni-2
         indn = i+1+nj*nin
         Fnx(indn) = b*(Xp(i)+Xp(i+1))
         Fny(indn) = b*(Yp(i)+Yp(i+1))
         Fnz(indn) = b*(Zp(i)+Zp(i+1))
      ENDDO


C=============== FIN k = 0 =================================

C OK

C=============== DEBUT i = ni =================================

C     i = 0, interieur

      DO k = 0, nk-1
         DO j = 0, nj-1
            ind = j+k*nj
            ind0 = ni-1+j*ni+k*ninj
            ind1 = ind0-1
            X(ind) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO
C OK
      DO k = 0, nk-2
         DO j = 0, nj-2
            indn = ni+(j+1)*nin+(k+1)*ninnjn
            ind0 = j+k*nj
            ind1 = ind0+1
            ind2 = ind0+nj
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO

C i = ni, j = nj

      DO k = 0, nk-1
         ind0 = (nj-1)+k*nj
         ind1 = (nj-2)+k*nj
         Xp(k) = b*(3*X(ind0)-X(ind1))
         Yp(k) = b*(3*Y(ind0)-Y(ind1))
         Zp(k) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO k = 0, nk-2
         indn = ni+nj*(ni+1)+(k+1)*(ni+1)*(nj+1)
         Fnx(indn) = b*(Xp(k)+Xp(k+1))
         Fny(indn) = b*(Yp(k)+Yp(k+1))
         Fnz(indn) = b*(Zp(k)+Zp(k+1))
       ENDDO

C i = ni, k = nk

      DO j = 0, nj-1
         ind0 = j+(nk-1)*nj
         ind1 = j+(nk-2)*nj
         Xp(j) = b*(3*X(ind0)-X(ind1))
         Yp(j) = b*(3*Y(ind0)-Y(ind1))
         Zp(j) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO j = 0, nj-2
         indn = ni+(j+1)*nin+nk*ninnjn
         Fnx(indn) = b*(Xp(j)+Xp(j+1))
         Fny(indn) = b*(Yp(j)+Yp(j+1))
         Fnz(indn) = b*(Zp(j)+Zp(j+1))
      ENDDO


C=============== FIN i = ni =================================

C=============== DEBUT j = nj =================================

C     j = nj, interieur

      DO k = 0, nk-1
         DO i = 0, ni-1
            ind = i+k*ni
            ind0 = i+(nj-1)*ni+k*ninj
            ind1 = ind0-ni
            X(ind) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO

      DO k = 0, nk-2
         DO i = 0, ni-2
            indn = (i+1)+nj*nin+(k+1)*ninnjn
            ind0 = i+k*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO


C j = nj, k = nk

      DO i = 0, ni-1
         ind0 = i+(nk-1)*ni
         ind1 = i+(nk-2)*ni
         Xp(i) = b*(3*X(ind0)-X(ind1))
         Yp(i) = b*(3*Y(ind0)-Y(ind1))
         Zp(i) = b*(3*Z(ind0)-Z(ind1))        
      ENDDO

      DO i = 0, ni-2
         indn = (i+1)+nj*nin+nk*ninnjn
         Fnx(indn) = b*(Xp(i)+Xp(i+1))
         Fny(indn) = b*(Yp(i)+Yp(i+1))
         Fnz(indn) = b*(Zp(i)+Zp(i+1))
      ENDDO

C=============== FIN j = nj =================================


C=============== DEBUT k = nk =================================

C     k = nk, interieur

      DO j = 0, nj-1
         DO i = 0, ni-1
            ind = i+j*ni
            ind0 = i+j*ni+(nk-1)*ninj
            ind1 = ind0-ninj
            X(ind) = b*(3*Fcx(ind0)-Fcx(ind1))
            Y(ind) = b*(3*Fcy(ind0)-Fcy(ind1))
            Z(ind) = b*(3*Fcz(ind0)-Fcz(ind1))
         ENDDO
      ENDDO

      DO j = 0, nj-2
         DO i = 0, ni-2
            indn = i+1+(j+1)*nin+nk*ninnjn
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            Fnx(indn) = c*(X(ind0)+X(ind1)+X(ind2)+X(ind3))
            Fny(indn) = c*(Y(ind0)+Y(ind1)+Y(ind2)+Y(ind3))
            Fnz(indn) = c*(Z(ind0)+Z(ind1)+Z(ind2)+Z(ind3))
         ENDDO
      ENDDO
C=============== FIN k = nk =================================

C ========Les coins ===========================================
C i = 0, j = 0, k = 0
      Fnx(0) = 2*Fnx(ni+1)-Fnx(2*(ni+1))
      Fny(0) = 2*Fny(ni+1)-Fny(2*(ni+1))
      Fnz(0) = 2*Fnz(ni+1)-Fnz(2*(ni+1))
C i = 0, j = nj, k = 0
      Fnx(nj*nin) = 2*Fnx((nj-1)*nin)-Fnx((nj-2)*nin)
      Fny(nj*nin) = 2*Fny((nj-1)*nin)-Fny((nj-2)*nin)
      Fnz(nj*nin) = 2*Fnz((nj-1)*nin)-Fnz((nj-2)*nin)
C i = 0, j = 0, k = nk
      Fnx(nk*ninnjn) = 2*Fnx(nin+nk*ninnjn)
     &                           -Fnx(2*nin+nk*ninnjn)
      Fny(nk*ninnjn) = 2*Fny(nin+nk*ninnjn)
     &                           -Fny(2*nin+nk*ninnjn)
      Fnz(nk*ninnjn) = 2*Fnz(nin+nk*ninnjn)
     &                           -Fnz(2*nin+nk*ninnjn)
C i = 0, j = nj, k = nk
      Fnx(nj*nin+nk*ninnjn) = 2*Fnx((nj-1)*nin+nk*ninnjn)
     &                             -Fnx((nj-2)*nin+nk*ninnjn) 
      Fny(nj*nin+nk*ninnjn) = 2*Fny((nj-1)*nin+nk*ninnjn)
     &                             -Fny((nj-2)*nin+nk*ninnjn) 
      Fnz(nj*nin+nk*ninnjn) = 2*Fnz((nj-1)*nin+nk*ninnjn)
     &                             -Fnz((nj-2)*nin+nk*ninnjn) 

C i = ni, j = 0, k = 0
      Fnx(ni) = 2*Fnx(2*ni+1)-Fnx(ni+2*nin)
      Fny(ni) = 2*Fny(2*ni+1)-Fny(ni+2*nin)
      Fnz(ni) = 2*Fnz(2*ni+1)-Fnz(ni+2*nin)

C i = ni, j = nj, k = 0
      Fnx(ni+(nj)*nin) = 
     &                   2*Fnx(ni+(nj-1)*nin)-Fnx(ni+(nj-2)*nin)
      Fny(ni+(nj)*nin) = 
     &                   2*Fny(ni+(nj-1)*nin)-Fny(ni+(nj-2)*nin)
      Fnz(ni+(nj)*nin) = 
     &                   2*Fnz(ni+(nj-1)*nin)-Fnz(ni+(nj-2)*nin)


C i = ni, j = 0, k = nk
      Fnx(ni+nk*ninnjn) = 2*Fnx(ni+nin+nk*ninnjn)
     &                           -Fnx(ni+2*nin+nk*ninnjn)
      Fny(ni+nk*ninnjn) = 2*Fny(ni+nin+nk*ninnjn)
     &                           -Fny(ni+2*nin+nk*ninnjn)
      Fnz(ni+nk*ninnjn) = 2*Fnz(ni+nin+nk*ninnjn)
     &                           -Fnz(ni+2*nin+nk*ninnjn)
C i = ni, j = nj, k = nk
      Fnx(ni+nj*nin+nk*ninnjn) = 
     &                       2*Fnx(ni+(nj-1)*nin+nk*ninnjn)
     &                       -Fnx(ni+(nj-2)*nin+nk*ninnjn) 
      Fny(ni+nj*nin+nk*ninnjn) = 
     &                       2*Fny(ni+(nj-1)*nin+nk*ninnjn)
     &                       -Fny(ni+(nj-2)*nin+nk*ninnjn) 
      Fnz(ni+nj*nin+nk*ninnjn) = 
     &                       2*Fnz(ni+(nj-1)*nin+nk*ninnjn)
     &                       -Fnz(ni+(nj-2)*nin+nk*ninnjn) 


      END

C  ============================================================================
C  Convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C  Dans le domaine interieur: Si au moins un 0 ==> 0
C                             Sinon 1
C  Sur la frontiere: Si au moins un 0 ==> 0
C                    Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node21( ni, nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indcenter, indn, i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E ninj, nin, njn, ninnjn
      INTEGER_E alpha, beta, gamma
      REAL_E prod
      REAL_E f

      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
      DO k = 0, nk
         DO j = 0, nj
            DO i = 0, ni
               indn = i+j*nin+k*ninnjn
               alpha = 1
               IF (i.EQ.0 .OR. i.EQ.ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1)-1         
               
               beta = ni
               IF (j.EQ.0 .OR. j.EQ.nj) THEN
                  beta = 0
               ENDIF
               j0 = MAX(j, 1)-1        
               
               gamma = ninj
               IF (k.EQ.0 .OR. k.EQ.nk) THEN
                  gamma = 0
               ENDIF
               k0 = MAX(k, 1)-1
               
               ind0 = i0 + j0 * ni + k0 * ninj
               ind1 = ind0+alpha
               ind2 = ind0+beta
               ind3 = ind2+alpha
               ind4 = ind0+gamma
               ind5 = ind4+alpha
               ind6 = ind4+beta
               ind7 = ind6+alpha
               prod = fieldcenter(ind0)*fieldcenter(ind1)
     &              *fieldcenter(ind2)*fieldcenter(ind3)
     &              *fieldcenter(ind4)*fieldcenter(ind5)
     &              *fieldcenter(ind6)*fieldcenter(ind7)
               IF (prod .EQ. 0) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C      Dans le domaine interieur: Si au moins un 0 ==> 0
C                                 Sinon 1
C      Sur la frontiere: Si au moins un 0 ==> 0
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node21p( ni, nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indcenter, indn, i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E ninj, nin, njn, ninnjn
      INTEGER_E alpha, beta, gamma
      REAL_E sum, f

      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
      DO k = 0, nk
         DO j = 0, nj
            DO i = 0, ni
               indn = i+j*nin+k*ninnjn
               alpha = 1
               IF (i.EQ.0 .OR. i.EQ.ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1) -1         
               
               beta = ni
               IF (j.EQ.0 .OR. j.EQ.nj) THEN
                  beta = 0
               ENDIF
               j0 = MAX(j, 1) -1        
               
               gamma = ninj
               IF (k.EQ.0 .OR. k.EQ.nk) THEN
                  gamma = 0
               ENDIF
               k0 = MAX(k, 1) -1
               
               ind0 = i0 + j0 * ni + k0 * ninj
               ind1 = ind0+alpha
               ind2 = ind0+beta
               ind3 = ind2+alpha
               ind4 = ind0+gamma
               ind5 = ind4+alpha
               ind6 = ind4+beta
               ind7 = ind6+alpha
               sum = fieldcenter(ind0) + fieldcenter(ind1)
     &              + fieldcenter(ind2) + fieldcenter(ind3)
     &              + fieldcenter(ind4) + fieldcenter(ind5)
     &              + fieldcenter(ind6) + fieldcenter(ind7)
               IF (sum .EQ. 0) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      END

C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur: Si au moins un 0 ==> 0
C                                 Sinon si au moins un 2 ==> 2
C                                 Sinon 1
C      Sur la frontiere: Si au moins un 0 ==> 0
C                        Sinon si tous 2 ==> 2
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node22(ni, nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indcenter, indn, i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nin, njn, ninnjn, ninj
      REAL_E f, prod, alpha, beta, gamma, somme2

      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i, j, indn, ind0, ind1, ind2, ind3, ind4, ind5,
!$OMP& ind6, ind7, prod, somme2)
      DO k = 1, nk-1
!$OMP DO
         DO j = 1, nj-1
            DO i = 1, ni-1
               indn = i+j*nin+k*ninnjn
               ind0 = i-1+(j-1)*ni+(k-1)*ninj
               ind1 = ind0+1
               ind2 = ind0+ni
               ind3 = ind2+1
               ind4 = ind0+ninj
               ind5 = ind4+1
               ind6 = ind4+ni
               ind7 = ind6+1
               prod = fieldcenter(ind0)*fieldcenter(ind1)
     &              *fieldcenter(ind2)*fieldcenter(ind3)
     &              *fieldcenter(ind4)*fieldcenter(ind5)
     &              *fieldcenter(ind6)*fieldcenter(ind7)
               IF (prod .EQ. 0) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  somme2 = MAX(fieldcenter(ind0),1.)
     &                 + MAX(fieldcenter(ind1),1.)
     &                 + MAX(fieldcenter(ind2),1.)
     &                 + MAX(fieldcenter(ind3),1.)
     &                 + MAX(fieldcenter(ind4),1.)
     &                 + MAX(fieldcenter(ind5),1.)
     &                 + MAX(fieldcenter(ind6),1.)
     &                 + MAX(fieldcenter(ind7),1.)
                  IF (somme2 .GT. 8.5) THEN
                     fieldnode(indn) = 2.D0
                  ELSE
                     fieldnode(indn) = 1.D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!$OMP ENDDO
      ENDDO
!$OMP END PARALLEL

C     Traitement des fontieres
C     i = 0
      DO k = 0, nk
         DO j = 0, nj
            indn = j*nin+k*ninnjn
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = j0*ni + k0*ninj 
            ind1 = ind0+beta              
            ind2 = ind0+gamma                      
            ind3 = ind2+beta
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (prod .EQ. 16.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C
      i = ni
      i0 = i-1
      DO k = 0, nk
         DO j = 0, nj
            indn = i+j*nin+k*ninnjn
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = i0 + j0*ni + k0*ninj 
            ind1 = ind0+beta              
            ind2 = ind0+gamma                      
            ind3 = ind2+beta
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (prod .EQ. 16.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     j = 0
      DO k = 0, nk
         DO i = 0, ni
            indn = i+k*ninnjn
               alpha = 1
               IF (i.EQ.0 .OR. i.EQ.ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1)-1         
               
               gamma = ninj
               IF (k.EQ.0 .OR. k.EQ.nk) THEN
                  gamma = 0
               ENDIF
               k0 = MAX(k, 1)-1
               
               ind0 = i0 + k0 * ninj 
               ind1 = ind0+alpha              
               ind2 = ind0+gamma                      
               ind3 = ind2+alpha
               prod = fieldcenter(ind0)*fieldcenter(ind1)
     &              *fieldcenter(ind2)*fieldcenter(ind3)
               IF (prod .EQ. 0.) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  IF (prod .EQ. 16.) THEN
                     fieldnode(indn) = 2.D0
                  ELSE
                     fieldnode(indn) = 1.D0
                  ENDIF
               ENDIF
         ENDDO
      ENDDO
      
      j = nj
      j0 = j-1
      DO k = 0, nk
         DO i = 0, ni
            indn = i+j*nin+k*ninnjn
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = i0 +j0 * ni +k0 * ninj
            ind1 = ind0+alpha  
            ind4 = ind0+gamma                      
            ind5 = ind4+alpha
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (prod .EQ. 16.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     k = 0
      DO j = 0, nj
         DO i = 0, ni
            indn = i+j*nin
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
         
            ind0 = i0 + j0 * ni
            ind1 = ind0+alpha              
            ind2 = ind0+beta              
            ind3 = ind2+alpha              
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (prod .EQ. 16.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C     k = nk
      k = nk
      k0 = k-1
      DO j = 0, nj
         DO i = 0, ni
            indn = i+j*nin+k*ninnjn
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1         
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1
         
            ind0 = i0 + j0 * ni+k0*ninj 
            ind1 = ind0+alpha              
            ind2 = ind0+beta              
            ind3 = ind2+alpha              
            prod = fieldcenter(ind0)*fieldcenter(ind1)
     &           *fieldcenter(ind2)*fieldcenter(ind3)
            IF (prod .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (prod .EQ. 16.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      END
C  ============================================================================
C  convert cellnature field defined on centers to field defined on nodes
C  cellN = 0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur: Si au moins un 0 ==> 0
C                                 Sinon si au moins un 2 ==> 2
C                                 Sinon 1
C      Sur la frontiere: Si au moins un 0 ==> 0
C                        Sinon si tous 2 ==> 2
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2node22p(ni, nj, nk, fieldcenter, fieldnode)
      IMPLICIT NONE

C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldcenter(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldnode(0:(ni+1)*(nj+1)*(nk+1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indcenter, indn, i0, j0, k0
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nin, njn, ninnjn, ninj
      REAL_E f, sum, alpha, beta, gamma, somme2

      nin = ni+1
      njn = nj+1
      ninnjn = nin*njn
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i, j, indn, ind0, ind1, ind2, ind3, ind4, ind5,
!$OMP& ind6, ind7, sum, somme2)
      DO k = 1, nk-1
!$OMP DO
         DO j = 1, nj-1
            DO i = 1, ni-1
               indn = i+j*nin+k*ninnjn
               ind0 = i-1+(j-1)*ni+(k-1)*ninj
               ind1 = ind0+1
               ind2 = ind0+ni
               ind3 = ind2+1
               ind4 = ind0+ninj
               ind5 = ind4+1
               ind6 = ind4+ni
               ind7 = ind6+1
               sum = fieldcenter(ind0) + fieldcenter(ind1)
     &              + fieldcenter(ind2) + fieldcenter(ind3)
     &              + fieldcenter(ind4) + fieldcenter(ind5)
     &              + fieldcenter(ind6) + fieldcenter(ind7)
               IF (sum .EQ. 0) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  somme2 = MAX(fieldcenter(ind0),1.)
     &                 + MAX(fieldcenter(ind1),1.)
     &                 + MAX(fieldcenter(ind2),1.)
     &                 + MAX(fieldcenter(ind3),1.)
     &                 + MAX(fieldcenter(ind4),1.)
     &                 + MAX(fieldcenter(ind5),1.)
     &                 + MAX(fieldcenter(ind6),1.)
     &                 + MAX(fieldcenter(ind7),1.)
                  IF (somme2. GT. 8.5) THEN
                     fieldnode(indn) = 2.D0
                  ELSE
                     fieldnode(indn) = 1.D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!$OMP ENDDO
      ENDDO
!$OMP END PARALLEL

C     Traitement des fontieres
C     i = 0
      DO k = 0, nk
         DO j = 0, nj
            indn = j*nin+k*ninnjn
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = j0*ni + k0*ninj 
            ind1 = ind0+beta              
            ind2 = ind0+gamma                      
            ind3 = ind2+beta
            sum = fieldcenter(ind0)+ fieldcenter(ind1)
     &           + fieldcenter(ind2)+ fieldcenter(ind3)
            IF (sum .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (sum .EQ. 8.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C
      i = ni
      i0 = i-1
      DO k = 0, nk
         DO j = 0, nj
            indn = i+j*nin+k*ninnjn
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = i0 + j0*ni + k0*ninj 
            ind1 = ind0+beta              
            ind2 = ind0+gamma                      
            ind3 = ind2+beta
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           + fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (sum .EQ. 8.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     j = 0
      DO k = 0, nk
         DO i = 0, ni
            indn = i+k*ninnjn
               alpha = 1
               IF (i.EQ.0 .OR. i.EQ.ni) THEN
                  alpha = 0
               ENDIF
               i0 = MAX(i, 1)-1         
               
               gamma = ninj
               IF (k.EQ.0 .OR. k.EQ.nk) THEN
                  gamma = 0
               ENDIF
               k0 = MAX(k, 1)-1
               
               ind0 = i0 + k0 * ninj 
               ind1 = ind0+alpha              
               ind2 = ind0+gamma                      
               ind3 = ind2+alpha
               sum = fieldcenter(ind0) + fieldcenter(ind1)
     &              + fieldcenter(ind2) + fieldcenter(ind3)
               IF (sum .EQ. 0.) THEN
                  fieldnode(indn) = 0.D0
               ELSE
                  IF (sum .EQ. 8.) THEN
                     fieldnode(indn) = 2.D0
                  ELSE
                     fieldnode(indn) = 1.D0
                  ENDIF
               ENDIF
         ENDDO
      ENDDO
      
      j = nj
      j0 = j-1
      DO k = 0, nk
         DO i = 0, ni
            indn = i+j*nin+k*ninnjn
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1
            
            gamma = ninj
            IF (k.EQ.0 .OR. k.EQ.nk) THEN
               gamma = 0
            ENDIF
            k0 = MAX(k, 1)-1
            
            ind0 = i0 +j0 * ni +k0 * ninj 
            ind1 = ind0+alpha              
            ind4 = ind0+gamma                      
            ind5 = ind4+alpha
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           + fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (sum .EQ. 8.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
C     k = 0
      DO j = 0, nj
         DO i = 0, ni
            indn = i+j*nin
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1         
         
            ind0 = i0 + j0 * ni
            ind1 = ind0+alpha              
            ind2 = ind0+beta              
            ind3 = ind2+alpha              
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           + fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (sum .EQ. 8.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C     k = nk
      k = nk
      k0 = k-1
      DO j = 0, nj
         DO i = 0, ni
            indn = i+j*nin+k*ninnjn
            alpha = 1
            IF (i.EQ.0 .OR. i.EQ.ni) THEN
               alpha = 0
            ENDIF
            i0 = MAX(i, 1)-1         
            
            beta = ni
            IF (j.EQ.0 .OR. j.EQ.nj) THEN
               beta = 0
            ENDIF
            j0 = MAX(j, 1)-1
         
            ind0 = i0 + j0 * ni+k0*ninj 
            ind1 = ind0+alpha              
            ind2 = ind0+beta              
            ind3 = ind2+alpha              
            sum = fieldcenter(ind0) + fieldcenter(ind1)
     &           + fieldcenter(ind2) + fieldcenter(ind3)
            IF (sum .EQ. 0.) THEN
               fieldnode(indn) = 0.D0
            ELSE
               IF (sum .EQ. 8.) THEN
                  fieldnode(indn) = 2.D0
               ELSE
                  fieldnode(indn) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      END
C ===== Post/Fortran/Conv2NodeF.for =====
