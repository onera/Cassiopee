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
C  Calcul du rotationnel moyen d un champ (u,v,w) sur une cellule
C  attention cette routine est uniquement 3d
C  IN: ind: indice du premier sommet de la cellule
C  IN: ni,nj,nk: dimensions du maillage en noeuds
C  IN: velo: vecteur dont le rotationnel est a calculer. Defini sur le maillage
C  OUT: rotu,rotv,rotw: rotationnel moyen de (u,v,w) sur la cellule
C  ============================================================================
      SUBROUTINE k6compmeancurlofstructcell(ind, ni, nj, nk, 
     &                                      velox, veloy, veloz,  
     &                                      xt, yt, zt,rotu,rotv,rotw)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E ind             ! indice du premier noeud de la cellule
      INTEGER_E ni, nj, nk      ! dimensions de la grille contenant la cellule
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
      REAL_E velox(0:ni*nj*nk-1) !vecteur dont le rotationnel est a calculer ds la cellule
      REAL_E veloy(0:ni*nj*nk-1)
      REAL_E veloz(0:ni*nj*nk-1)
C_OUT
      REAL_E  rotu, rotv, rotw ! curl of u,v,w 

C_LOCAL
      INTEGER_E i,j,k,l
      INTEGER_E ninj
      INTEGER_E l1, l2, l3, l4 
      REAL_E   d13x, d13y, d13z, d24x, d24y, d24z
      REAL_E   surf(0:5,3)      ! tableau des vecteurs surfaces sur chaque interface 
      REAL_E   uint(0:5,3)      ! vitesses aux interfaces
      REAL_E   curlx, curly, curlz
      REAL_E   cellvol
C
C-----------------------------------------------------------------------------
C
C     Calcul des surfaces aux interfaces
      call k6compintsurfofcell(ind, ni, nj, nk, xt, yt, zt, surf) 

C     Calcul du volume de la cellule 
      call k6compvolofstructcell(ni,nj,nk,-1,ind,xt,yt,zt,cellvol)

C-----------------------------------------------------------------------------
C     Calcul des vecteurs vitesse aux interfaces
      ninj = ni*nj
C     interface en i
      l1 = ind                  ! i,j,k
      l2 = l1 + ni              ! i,j+1,k
      l3 = l2 + ninj            ! i,j+1,k+1
      l4 = l1 + ninj            ! i,j,k+1
C
      uint(0,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(0,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(0,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4)
C
C     interface en i+1
      l1 = ind + 1              ! i+1,j,k
      l2 = l1 + ni              ! i+1,j+1,k
      l3 = l2 + ninj            ! i+1,j+1,k+1
      l4 = l1 + ninj            ! i+1,j,k+1
c
      uint(1,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(1,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(1,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4)
C
C     interface en j
      l1 = ind + 1              ! i+1,j,k
      l2 = l1-1                 ! i,j,k
      l3 = l2 + ninj            ! i,j,k+1
      l4 = l1 + ninj            ! i+1,j,k+1
c
      uint(2,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(2,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(2,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4)  
C
C     interface en j+1
      l1 = ind + ni + 1         ! i+1,j+1,k
      l2 = l1-1                 ! i,j+1,k
      l3 = l2 + ninj            ! i,j+1,k+1
      l4 = l1 + ninj            ! i,j+1,k+1
c
      uint(3,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(3,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(3,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4) 
C
C     interface en k
      l1 = ind                  ! i,j,k
      l2 = l1 + 1               ! i+1,j,k
      l3 = l2 + ni              ! i+1,j+1,k
      l4 = l1 + ni              ! i,j+1,k
c
      uint(4,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(4,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(4,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4) 
C
C     interface en k+1      
      l1 = ind + ninj           ! i,j,k+1
      l2 = l1 + 1               ! i+1,j,k+1
      l3 = l2 + ni              ! i+1,j+1,k+1
      l4 = l1 + ni              ! i,j+1,k+1
c
      uint(5,1) = velox(l1) + velox(l2) + velox(l3) + velox(l4)
      uint(5,2) = veloy(l1) + veloy(l2) + veloy(l3) + veloy(l4)
      uint(5,3) = veloz(l1) + veloz(l2) + veloz(l3) + veloz(l4) 
C
C-----------------------------------------------------------------------------
      curlx = 0.
      curly = 0.
      curlz = 0.
      DO l = 0, 5
         curlx = curlx + uint(l,2)*surf(l,3) - uint(l,3)*surf(l,2)
         curly = curly + uint(l,3)*surf(l,1) - uint(l,1)*surf(l,3)
         curlz = curlz + uint(l,1)*surf(l,2) - uint(l,2)*surf(l,1)
      ENDDO
     
      cellvol = -ONE_FOURTH / MAX(cellvol, E_MIN_VOL)
      
      rotu = cellvol * curlx
      rotv = cellvol * curly
      rotw = cellvol * curlz

      END

C ========= Post/Fortran/CompMeanCurlOfStructCellF.for ===================
