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

C Calcul du rotationnel d'un vecteur defini aux noeuds d une grille surfacique
C retourne le rotationnel defini aux centres des cellules
C ============================================================================
      SUBROUTINE k6compstructcurl2dt(ni, nj, nbcell, 
     &     xt, yt, zt, ux, uy, uz, 
     &     rotux, rotuy, rotuz) 
      
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C=============================================================================
C_IN
      INTEGER_E ni, nj          ! dimensions de la grille
      INTEGER_E nbcell          ! nb de cellules
      REAL_E xt(0:ni*nj-1)
      REAL_E yt(0:ni*nj-1)
      REAL_E zt(0:ni*nj-1)
      REAL_E ux(0:ni*nj-1)     ! vecteur aux noeuds
      REAL_E uy(0:ni*nj-1)     ! vecteur aux noeuds
      REAL_E uz(0:ni*nj-1)     ! vecteur aux noeuds

C_OUT 
      REAL_E rotux(0:nbcell-1) ! rotationnel de u au centre des cellules
      REAL_E rotuy(0:nbcell-1) ! rotationnel de u au centre des cellules
      REAL_E rotuz(0:nbcell-1) ! rotationnel de u au centre des cellules

C_LOCAL
      INTEGER_E i, j, ni1, nj1
      INTEGER_E indA, indB, indC, indD, indcell
      REAL_E vinv, nn 
      REAL_E vx, vy, vz
      REAL_E xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA
      REAL_E nx, ny, nz, n1x,n1y, n1z
      REAL_E curlx, curly, curlz
      REAL_E surf(0:nbcell-1)   ! surface des cellules 2d
      REAL_E nxt(0:nbcell-1)
      REAL_E nyt(0:nbcell-1)
      REAL_E nzt(0:nbcell-1)

C=============================================================================
C_CONSTANTES
      ni1 = ni-1
      nj1 = nj-1
      
      call k6structsurft(ni, nj, 1, nbcell, xt, yt, zt, surf)
C
      call k6normstructsurft(ni, nj, ni*nj, xt, yt, zt, nxt,nyt,nzt)
     
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
         vinv = -ONE/MAX(vinv, E_MIN_VOL)

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
C
         curlx = 0.D0
         curly = 0.D0
         curlz = 0.D0
C
         n1x = yAB*nz - zAB*ny
         n1y = zAB*nx - xAB*nz
         n1z = xAB*ny - yAB*nx

         vx = ux(indA) + ux(indB)
         vy = uy(indA) + uy(indB)
         vz = uz(indA) + uz(indB)

         curlx = curlx + vy*n1z - vz*n1y
         curly = curly + vz*n1x - vx*n1z
         curlz = curlz + vx*n1y - vy*n1x
C
         n1x = yBC*nz - zBC*ny
         n1y = zBC*nx - xBC*nz
         n1z = xBC*ny - yBC*nx

         vx = ux(indC) + ux(indB)
         vy = uy(indC) + uy(indB)
         vz = uz(indC) + uz(indB)
         curlx = curlx + vy*n1z - vz*n1y
         curly = curly + vz*n1x - vx*n1z
         curlz = curlz + vx*n1y - vy*n1x
C
         n1x = yCD*nz - zCD*ny
         n1y = zCD*nx - xCD*nz
         n1z = xCD*ny - yCD*nx

         vx = ux(indC) + ux(indD)
         vy = uy(indC) + uy(indD)
         vz = uz(indC) + uz(indD)
         curlx = curlx + vy*n1z - vz*n1y
         curly = curly + vz*n1x - vx*n1z
         curlz = curlz + vx*n1y - vy*n1x
C
         n1x = yDA*nz - zDA*ny
         n1y = zDA*nx - xDA*nz
         n1z = xDA*ny - yDA*nx

         vx = ux(indD) + ux(indA)
         vy = uy(indD) + uy(indA)
         vz = uz(indD) + uz(indA)

         curlx = curlx + vy*n1z - vz*n1y
         curly = curly + vz*n1x - vx*n1z
         curlz = curlz + vx*n1y - vy*n1x
C
         rotux(indcell) = vinv * curlx
         rotuy(indcell) = vinv * curly
         rotuz(indcell) = vinv * curlz
      ENDDO
      ENDDO
      END
C====================Post/Fortran/CompStructCurl2DF.for=======================
