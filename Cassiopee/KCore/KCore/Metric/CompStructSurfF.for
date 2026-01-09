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
C  Calcul de la surface pour une grille surfacique structuree
C  ============================================================================
      SUBROUTINE  k6structsurft(ni, nj, nk, ncells, xt, yt, zt, surface)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      INTEGER_E ncells
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E surface(0:ncells-1)              

C_LOCAL
      INTEGER_E i, j, k, ninj, ni1, nj1, nk1
      INTEGER_E ind1, ind2, ind3, ind4, indcell
      REAL_E l1x, l1y, l1z
      REAL_E l2x, l2y, l2z
      REAL_E l3x, l3y, l3z
      REAL_E l4x, l4y, l4z
      REAL_E surf1x, surf1y, surf1z
      REAL_E surf2x, surf2y, surf2z
      REAL_E surface1, surface2, ps
C
C==============================================================================
C 
      ninj = ni*nj
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1

      IF (ni .eq. 1) THEN
         DO k = 0, nk-2
         DO j = 0, nj-2
            indcell = j + k * nj1
            ind1 = j + k*nj
            ind2 = ind1 + nj
            ind3 = ind1 + 1
            ind4 = ind2 + 1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            l2x = xt(ind1)-xt(ind3)
            l2y = yt(ind1)-yt(ind3)
            l2z = zt(ind1)-zt(ind3)
            surf1x = (l1y*l2z-l1z*l2y)
            surf1y = (l1z*l2x-l1x*l2z)
            surf1z = (l1x*l2y-l1y*l2x)
            surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
            l3x = xt(ind4)-xt(ind3)
            l3y = yt(ind4)-yt(ind3)
            l3z = zt(ind4)-zt(ind3)
            l4x = xt(ind4)-xt(ind2)
            l4y = yt(ind4)-yt(ind2)
            l4z = zt(ind4)-zt(ind2)
            surf2x = (l3y*l4z-l3z*l4y)
            surf2y = (l3z*l4x-l3x*l4z)
            surf2z = (l3x*l4y-l3y*l4x)
            surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
            ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
            surface(indcell) = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
         ENDDO
         ENDDO
      ELSE IF (nj .eq. 1) THEN 
         DO i = 0, ni-2
         DO k = 0, nk-2
            indcell = i + k *ni1
            ind1 = i + k*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind2 + 1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            l2x = xt(ind1)-xt(ind3)
            l2y = yt(ind1)-yt(ind3)
            l2z = zt(ind1)-zt(ind3)  
            surf1x = (l1y*l2z-l1z*l2y)
            surf1y = (l1z*l2x-l1x*l2z)
            surf1z = (l1x*l2y-l1y*l2x)
            surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
            l3x = xt(ind4)-xt(ind3)
            l3y = yt(ind4)-yt(ind3)
            l3z = zt(ind4)-zt(ind3)
            l4x = xt(ind4)-xt(ind2)
            l4y = yt(ind4)-yt(ind2)
            l4z = zt(ind4)-zt(ind2)
            surf2x = (l3y*l4z-l3z*l4y)
            surf2y = (l3z*l4x-l3x*l4z)
            surf2z = (l3x*l4y-l3y*l4x)
            surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
            ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
            surface(indcell) = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
         ENDDO
         ENDDO
      ELSE IF (nk .eq. 1) THEN 
         DO j = 0, nj-2
         DO i = 0, ni-2
            indcell = i + j * ni1
            ind1 = i + j*ni
            ind2 = ind1 + ni
            ind3 = ind1 + 1
            ind4 = ind2 + 1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            l2x = xt(ind1)-xt(ind3)
            l2y = yt(ind1)-yt(ind3)
            l2z = zt(ind1)-zt(ind3)            
         
            surf1x = (l1y*l2z-l1z*l2y)
            surf1y = (l1z*l2x-l1x*l2z)
            surf1z = (l1x*l2y-l1y*l2x)
            surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
            l3x = xt(ind4)-xt(ind3)
            l3y = yt(ind4)-yt(ind3)
            l3z = zt(ind4)-zt(ind3)
            l4x = xt(ind4)-xt(ind2)
            l4y = yt(ind4)-yt(ind2)
            l4z = zt(ind4)-zt(ind2)     
            surf2x = (l3y*l4z-l3z*l4y)
            surf2y = (l3z*l4x-l3x*l4z)
            surf2z = (l3x*l4y-l3y*l4x)
            surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
            ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
            surface(indcell) = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
         ENDDO
         ENDDO
      ELSE 
         write(*,*) 'In k6structsurf: 2D field is required.'
         STOP
      ENDIF

      END

C  ============================================================================
C  Calcul de la longueur entre chaque sommet pour une ligne structuree
C  ============================================================================
      SUBROUTINE  k6structsurf1dt(ni, nj, nk, xt, yt, zt, length)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E length(0:ni*nj*nk-2)              

C_LOCAL
      INTEGER_E i, j, k
      INTEGER_E ind1, ind2, indcell
      REAL_E l1x, l1y, l1z
C
C==============================================================================
C 
      IF ((nj .eq. 1) .AND. (nk .eq. 1)) THEN
         DO i = 0, ni-2
            indcell = i 
            ind1 = i 
            ind2 = i+1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            length(indcell) = sqrt(l1x*l1x+l1y*l1y+l1z*l1z)
         ENDDO
      ELSE IF ((ni .eq. 1) .AND. (nk .eq. 1)) THEN
         DO j = 0, nj-2
            indcell = j 
            ind1 = j 
            ind2 = j+1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            length(indcell) = sqrt(l1x*l1x+l1y*l1y+l1z*l1z)
         ENDDO
      ELSE IF ((ni .eq. 1) .AND. (nj .eq. 1)) THEN
         DO k = 0, nk-2
            indcell = k 
            ind1 = k 
            ind2 = k+1
            l1x = xt(ind1)-xt(ind2)
            l1y = yt(ind1)-yt(ind2)
            l1z = zt(ind1)-zt(ind2)
            length(indcell) = sqrt(l1x*l1x+l1y*l1y+l1z*l1z)
         ENDDO
      ELSE 
         write(*,*) 'In k6structsurf1d: 1D field is required.'
         STOP
      ENDIF

      END


C  ============================================================================
C  Calcul de la surface pour une grille surfacique structuree
C  ============================================================================
      SUBROUTINE  k6compsurfofstructcell(ni, nj, nk, indcell, 
     &     xt, yt, zt, surface)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      INTEGER_E indcell
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E surface              

C_LOCAL
      INTEGER_E i, j, k, ninj, ni1, nj1, nk1
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E l1x, l1y, l1z
      REAL_E l2x, l2y, l2z
      REAL_E l3x, l3y, l3z
      REAL_E l4x, l4y, l4z
      REAL_E surf1x, surf1y, surf1z
      REAL_E surf2x, surf2y, surf2z
      REAL_E surface1, surface2, ps
C
C==============================================================================
C 
      ninj = ni*nj
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
      
C     indcell = i + j*ni + k*ni*nj
      i = MOD(MOD(indcell, (ni*nj)), ni);
      j = MOD(indcell, (ni*nj))/ni;
      k = indcell/(ni*nj);

      IF (i .eq. ni1) THEN 
         i = i-1
      ENDIF
      IF (j .eq. nj1) THEN 
         j = j-1
      ENDIF  
      IF (k .eq. nk1) THEN 
         k = k-1
      ENDIF

      IF (ni .eq. 1) THEN
         ind1 = j + k*nj
         ind2 = ind1 + nj
         ind3 = ind1 + 1
         ind4 = ind2 + 1
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
         surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
         l3x = xt(ind4)-xt(ind3)
         l3y = yt(ind4)-yt(ind3)
         l3z = zt(ind4)-zt(ind3)
         l4x = xt(ind4)-xt(ind2)
         l4y = yt(ind4)-yt(ind2)
         l4z = zt(ind4)-zt(ind2)
         surf2x = (l3y*l4z-l3z*l4y)
         surf2y = (l3z*l4x-l3x*l4z)
         surf2z = (l3x*l4y-l3y*l4x)
         surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
         ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
         surface = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
      ELSE IF (nj .eq. 1) THEN 
         indcell = i + k *ni1
         ind1 = i + k*ni
         ind2 = ind1 + ni
         ind3 = ind1 + 1
         ind4 = ind2 + 1
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)  
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
         surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
         l3x = xt(ind4)-xt(ind3)
         l3y = yt(ind4)-yt(ind3)
         l3z = zt(ind4)-zt(ind3)
         l4x = xt(ind4)-xt(ind2)
         l4y = yt(ind4)-yt(ind2)
         l4z = zt(ind4)-zt(ind2)
         surf2x = (l3y*l4z-l3z*l4y)
         surf2y = (l3z*l4x-l3x*l4z)
         surf2z = (l3x*l4y-l3y*l4x)
         surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
         ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
         surface = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
      ELSE IF (nk .eq. 1) THEN 
         indcell = i + j * ni1
         ind1 = i + j*ni
         ind2 = ind1 + ni
         ind3 = ind1 + 1
         ind4 = ind2 + 1
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)            
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
         surface1 = sqrt(surf1x*surf1x+surf1y*surf1y+surf1z*surf1z)
         l3x = xt(ind4)-xt(ind3)
         l3y = yt(ind4)-yt(ind3)
         l3z = zt(ind4)-zt(ind3)
         l4x = xt(ind4)-xt(ind2)
         l4y = yt(ind4)-yt(ind2)
         l4z = zt(ind4)-zt(ind2)     
         surf2x = (l3y*l4z-l3z*l4y)
         surf2y = (l3z*l4x-l3x*l4z)
         surf2z = (l3x*l4y-l3y*l4x)
         surface2 = sqrt(surf2x*surf2x+surf2y*surf2y+surf2z*surf2z)
         ps = surf1x*surf2x+surf1y*surf2y+surf1z*surf2z
         surface = SIGN(ONE,ps)*ONE_HALF*(surface1+surface2)
      ELSE 
         write(*,*) 'In k6compsurfofstructcell: 2D field is required.'
         STOP
      ENDIF

      END
