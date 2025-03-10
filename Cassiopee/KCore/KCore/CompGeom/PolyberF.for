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
C Generate a demi-profile using bernstein polynomials
C ===========================================================================
      SUBROUTINE k6polyber(np, nb, xp, yp, zp,
     &                     xb, yb, zb)

      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E np  ! nombre de points de controle
      INTEGER_E nb  ! nombre de points dans la ligne finale
      REAL_E xp(np), yp(np), zp(np) ! coord. points de controle
      
C_OUT
      REAL_E xb(nb), yb(nb), zb(nb) ! coord. ligne de sortie

C_LOCAL 
      INTEGER_E i
      INTEGER m 
      REAL_E brst, u, du
C==============================================================================
c Approximation de bernstein
c  u varie entre 0 et 1
c 
      m = np-1
      du = 1.D0/(nb-1)
      DO i = 1, nb
         u = (i-1)*du
         xb(i) = brst(m, u, xp)
         yb(i) = brst(m, u, yp)
         zb(i) = brst(m, u, zp)
      ENDDO   
      
      RETURN
      END

C==============================================================================
      FUNCTION brst(m, u, s) 
c Approximation de bernstein
c  u varie entre 0 et 1

      IMPLICIT NONE

      INTEGER m, j
      REAL_E brst, s(m+1), u
      REAL_E cb
      brst = 0.D0

      DO j = 0, m
         brst = brst
     &        +cb(m,j)*s(j+1)*u**j*(1.-u)**(m-j)
      ENDDO
      RETURN
      END

C==============================================================================
      FUNCTION cb(m,i)

      IMPLICIT NONE

      INTEGER m,i
      INTEGER j
      REAL_E cb

      cb = 1.D0
      DO j = 1, i
       cb = cb*(m-i+j)/j
      ENDDO

      RETURN
      END
