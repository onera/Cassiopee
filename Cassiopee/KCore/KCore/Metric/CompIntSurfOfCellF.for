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

C  ============================================================================
C  Calcul des vecteurs surfaces des 6 interfaces d'une cellule 3D 
C  structuree 
C  IN: ind: indice du premier sommet de la cellule
C  IN: ni,nj,nk: dimensions du maillage en noeuds
C  OUT: surf: vecteur surface de la cellule
C  ============================================================================
      SUBROUTINE k6compintsurfofcell(ind, ni, nj, nk, xt, yt, zt, surf) 
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ind             ! indice du premier noeud de la cellule
      INTEGER_E ni, nj, nk      ! dimensions de la grille contenant la cellule
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
      
C_OUT
      REAL_E   surf(0:5,3)      ! tableau des vecteurs surfaces sur chaque interface 

C_LOCAL
      INTEGER_E i,j,k
      INTEGER_E ninj
      INTEGER_E l1, l2, l3, l4 
      REAL_E   d13x, d13y, d13z, d24x, d24y, d24z
C
C-----------------------------------------------------------------------------
C
      ninj = ni*nj
C
      k = ind / ninj           
      j = ind/ni - k * nj 
      i = ind - j*ni - k*ninj
      
      i = i+1
      j = j+1
      k = k+1

C     interface en i
      l1 = ind                  ! i,j,k
      l2 = l1 + ni              ! i,j+1,k
      l3 = l2 + ninj            ! i,j+1,k+1
      l4 = l1 + ninj            ! i,j,k+1
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(0,1) = -ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(0,2) = -ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(0,3) = -ONE_HALF * (d13x*d24y - d13y*d24x)
      
C     interface en i+1
      l1 = ind + 1              ! i+1,j,k
      l2 = l1 + ni              ! i+1,j+1,k
      l3 = l2 + ninj            ! i+1,j+1,k+1
      l4 = l1 + ninj            ! i+1,j,k+1
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(1,1) = ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(1,2) = ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(1,3) = ONE_HALF * (d13x*d24y - d13y*d24x)    

C     interface en j
      l1 = ind + 1              ! i+1,j,k
      l2 = l1-1                 ! i,j,k
      l3 = l2 + ninj            ! i,j,k+1
      l4 = l1 + ninj            ! i+1,j,k+1
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(2,1) = -ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(2,2) = -ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(2,3) = -ONE_HALF * (d13x*d24y - d13y*d24x)

C     interface en j+1
      l1 = ind + ni + 1         ! i+1,j+1,k
      l2 = l1-1                 ! i,j+1,k
      l3 = l2 + ninj            ! i,j+1,k+1
      l4 = l1 + ninj            ! i,j+1,k+1
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(3,1) = ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(3,2) = ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(3,3) = ONE_HALF * (d13x*d24y - d13y*d24x)

C     interface en k
      l1 = ind                  ! i,j,k
      l2 = l1 + 1               ! i+1,j,k
      l3 = l2 + ni              ! i+1,j+1,k
      l4 = l1 + ni              ! i,j+1,k
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(4,1) = -ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(4,2) = -ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(4,3) = -ONE_HALF * (d13x*d24y - d13y*d24x)

C     interface en k+1      
      l1 = ind + ninj           ! i,j,k+1
      l2 = l1 + 1               ! i+1,j,k+1
      l3 = l2 + ni              ! i+1,j+1,k+1
      l4 = l1 + ni              ! i,j+1,k+1
      d13x = xt(l3) - xt(l1)
      d13y = yt(l3) - yt(l1)
      d13z = zt(l3) - zt(l1)

      d24x = xt(l4) - xt(l2)
      d24y = yt(l4) - yt(l2)
      d24z = zt(l4) - zt(l2)

      surf(5,1) = ONE_HALF * (d13y*d24z - d13z*d24y)
      surf(5,2) = ONE_HALF * (d13z*d24x - d13x*d24z)
      surf(5,3) = ONE_HALF * (d13x*d24y - d13y*d24x)

      END
C======== KCore/Metric/CompIntSurfOfCellF.for ====
