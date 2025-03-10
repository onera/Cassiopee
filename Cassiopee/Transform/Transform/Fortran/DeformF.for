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
      SUBROUTINE k6deformpoint( size, dx, dy, dz,
     &                          amort,
     &                          ni, nj, nk,
     &                          x, y, z)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E size                   ! Number of points of the deformed window
      REAL_E dx(0:size-1)              ! Delta vector of deformed window
      REAL_E dy(0:size-1) 
      REAL_E dz(0:size-1) 
      INTEGER_E ni, nj, nk             ! Size of mesh to be deformed
      REAL_E amort                     ! Damping in j
C_OUT
      REAL_E x(0:ni*nj*nk-1)           ! Mesh coordinates
      REAL_E y(0:ni*nj*nk-1)
      REAL_E z(0:ni*nj*nk-1)
C_LOCAL
      INTEGER_E i, j, k, ind, indd
      REAL_E alpha, deltax, deltay, deltaz
      REAL_E a, b
C==============================================================================
      IF (nj.EQ.1) THEN
         a = 0.D0
         b = 1.D0
      ELSE
         a = -amort/(nj-1)
         b = 1.+amort/(nj-1)
      ENDIF
     
      DO k = 1, nk
         DO i = 1, ni
        
            indd = i-1 + (k-1)*ni
            deltax = dx(indd)
            deltay = dy(indd)
            deltaz = dz(indd)

            DO j = 1, nj
               alpha = a*j+b
               ind = i-1 + (j-1)*ni + (k-1)*ni*nj
               x(ind) = x(ind) + alpha*deltax
               y(ind) = y(ind) + alpha*deltay
               z(ind) = z(ind) + alpha*deltaz
            ENDDO
         ENDDO
      ENDDO

      END
C ===== Transform/Fortran/DeformF.for ======
