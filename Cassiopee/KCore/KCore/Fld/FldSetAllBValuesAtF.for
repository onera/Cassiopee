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
C Set all border values of a an interface defined field.
C==============================================================================
      SUBROUTINE k6setallbvaluesatf(length, nfld, x, val, ni, nj, nk, 
     &                              border)

      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E  length
      INTEGER_E  nfld
      REAL_E     val
      INTEGER_E  ni, nj, nk
      INTEGER_E  border

C_OUT
      REAL_E     x(0:length-1, 1:nfld)

C_LOCAL
      INTEGER_E i, j, k, n, ni1, nj1, nk1, int
      INTEGER_E nbinti, nbintj, nbintij, ninj1, ni1nj1, ni1nj
C==============================================================================

      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
      nbinti = ni*nj1*nk1
      nbintj = ni1*nj*nk1
      nbintij = nbinti + nbintj
      ninj1 = ni*nj1
      ni1nj = ni1*nj
      ni1nj1 = ni1*nj1

      DO n = 1, nfld

C     I interfaces      
         DO k = 1, nk1
            DO j = 1, nj1
               DO i = 1, border
                  int = i-1 + (j-1)*ni + (k-1)*ninj1
                  x(int, n) = val
               END DO
            ENDDO
         ENDDO
      
         DO k = 1, nk1
            DO j = 1, nj1
               DO i = ni-border+1, ni
                  int = i-1 + (j-1)*ni + (k-1)*ninj1
                  x(int, n) = val
               END DO
            ENDDO
         ENDDO

C     J interfaces
         DO k = 1, nk1
            DO i = nbinti, nbinti+border*ni1-1
               int = i + (k-1)*ni1nj
               x(int, n) = val
            END DO
         ENDDO

         DO k = 1, nk1
            DO i = (nj-border)*ni1+nbinti, ni1nj-1+nbinti
               int = i + (k-1)*ni1nj
               x(int, n) = val
            END DO
         ENDDO
      ENDDO

C     K interfaces
      IF (nk.GT.2) THEN
         DO n = 1, nfld
            DO int = nbintij, nbintij+border*ni1nj1-1
               x(int, n) = val
            END DO
         
            DO int = (nk-border)*ni1nj1+nbintij, nbintij+nk*ni1nj1-1
               x(int, n) = val
            END DO
         ENDDO
      ENDIF

      END
