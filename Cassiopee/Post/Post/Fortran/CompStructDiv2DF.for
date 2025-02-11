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

C Calcul de la divergence d'un champ defini aux noeuds d une grille surfacique
C structuree
C retourne la divergence defini aux centres des cellules
C ============================================================================
      SUBROUTINE k6compstructdiv2d(ni, nj, nbcell, xt, yt, zt,
     &                             fieldX, fieldY, fieldZ,
     &                             surf, nxt, nyt, nzt,
     &                             div)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ni, nj          ! dimensions de la grille aux noeuds
      INTEGER_E nbcell          ! nb de cellules
      REAL_E xt(0:ni*nj-1)      ! coordonnees des noeuds de la grille
      REAL_E yt(0:ni*nj-1)      ! coordonnees des noeuds de la grille
      REAL_E zt(0:ni*nj-1)      ! coordonnees des noeuds de la grille
      REAL_E fieldX(0:ni*nj-1)  ! champ (x-comp) defini aux noeuds auquel on applique div
      REAL_E fieldY(0:ni*nj-1)  ! champ (y-comp) defini aux noeuds auquel on applique div
      REAL_E fieldZ(0:ni*nj-1)  ! champ (z-comp) defini aux noeuds auquel on applique div
      REAL_E surf(0:nbcell-1)   ! surface des cellules 2d
      REAL_E nxt(0:nbcell-1) ! normale aux cellules
      REAL_E nyt(0:nbcell-1) ! normale aux cellules
      REAL_E nzt(0:nbcell-1) ! normale aux cellules
C_OUT
      REAL_E div(0:nbcell-1) ! divergence du champ vectoriel

C_LOCAL
      INTEGER_E indcell, i, j, ni1, nj1, nk
      INTEGER_E indA, indB, indC, indD
      REAL_E xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA
      REAL_E nx, n1x, n2x, n3x, n4x
      REAL_E ny, n1y, n2y, n3y, n4y
      REAL_E nz, n1z, n2z, n3z, n4z
      REAL_E fxAB, fxBC, fxCD, fxDA
      REAL_E fyAB, fyBC, fyCD, fyDA
      REAL_E fzAB, fzBC, fzCD, fzDA
      REAL_E vinv, nn
      REAL_E gradxx, gradyy, gradzz
C ------------------------------------------------------------------
C_CONSTANTES
      ni1 = ni-1
      nj1 = nj-1
C
C ------------------------------------------------------------------
C     calcul de la surface totale des cellules
      nk = 1
      call k6structsurft(ni, nj, nk, nbcell, xt, yt, zt, surf)
C
      call k6normstructsurft(ni, nj, ni*nj, xt, yt, zt, nxt, nyt, nzt)

      DO j = 0, nj1-1
      DO i = 0, ni1-1
         indA = i + j*ni
         indB = indA + 1
         indC = indB + ni
         indD = indA + ni

         indcell = i + j * ni1

         nx = nxt(indcell)
         ny = nyt(indcell)
         nz = nzt(indcell)

         nn = sqrt(nx*nx+ny*ny+nz*nz)
         vinv = 2D0 * surf(indcell) * nn
         vinv = ONE/MAX(vinv, E_MIN_VOL)

         xAB = xt(indB)-xt(indA)
         yAB = yt(indB)-yt(indA)
         zAB = zt(indB)-zt(indA)

         xBC = xt(indC)-xt(indB)
         yBC = yt(indC)-yt(indB)
         zBC = zt(indC)-zt(indB)

         xCD = xt(indD)-xt(indC)
         yCD = yt(indD)-yt(indC)
         zCD = zt(indD)-zt(indC)

         xDA = xt(indA)-xt(indD)
         yDA = yt(indA)-yt(indD)
         zDA = zt(indA)-zt(indD)

         n1x = yAB*nz - zAB*ny
         n1y = zAB*nx - xAB*nz
         n1z = xAB*ny - yAB*nx

         n2x = yBC*nz - zBC*ny
         n2y = zBC*nx - xBC*nz
         n2z = xBC*ny - yBC*nx

         n3x = yCD*nz - zCD*ny
         n3y = zCD*nx - xCD*nz
         n3z = xCD*ny - yCD*nx

         n4x = yDA*nz - zDA*ny
         n4y = zDA*nx - xDA*nz
         n4z = xDA*ny - yDA*nx

         fxAB = fieldX(indA) + fieldX(indB)
         fxBC = fieldX(indB) + fieldX(indC)
         fxCD = fieldX(indC) + fieldX(indD)
         fxDA = fieldX(indD) + fieldX(indA)

         fyAB = fieldY(indA) + fieldY(indB)
         fyBC = fieldY(indB) + fieldY(indC)
         fyCD = fieldY(indC) + fieldY(indD)
         fyDA = fieldY(indD) + fieldY(indA)

         fzAB = fieldZ(indA) + fieldZ(indB)
         fzBC = fieldZ(indB) + fieldZ(indC)
         fzCD = fieldZ(indC) + fieldZ(indD)
         fzDA = fieldZ(indD) + fieldZ(indA)

         gradxx = fxAB*n1x + fxBC*n2x + fxCD*n3x + fxDA*n4x
         gradyy = fyAB*n1y + fyBC*n2y + fyCD*n3y + fyDA*n4y
         gradzz = fzAB*n1z + fzBC*n2z + fzCD*n3z + fzDA*n4z

         div(indcell) = vinv *(gradxx + gradyy + gradzz)

      ENDDO
      ENDDO
      END

C=============================================================================

