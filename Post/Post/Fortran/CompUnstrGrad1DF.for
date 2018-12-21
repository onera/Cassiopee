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
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

C Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
C retourne le gradient defini aux centres des elts
C CAS 1D
C ============================================================================
      SUBROUTINE k6compunstrgrad1d( npts, nelts, nnodes,   
     &                              cn, xt, yt, zt, field, 
     &                              gradx, grady, gradz )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de points ds le maillage
      INTEGER_E nelts           ! nb d elts dans le maillage
      INTEGER_E nnodes          ! nb de noeuds par elt 
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E field(0:npts-1)    ! champ defini aux noeuds 
      REAL_E xt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des noeuds de la grille

C_OUT
      REAL_E gradx(0:nelts-1) ! gradient de field %x
      REAL_E grady(0:nelts-1) ! gradient de field %y
      REAL_E gradz(0:nelts-1) ! gradient de field %z

C_LOCAL
      INTEGER_E nfaces, et
      INTEGER_E indA, indB
      REAL_E xAB, yAB, zAB
      REAL_E fAB
      REAL_E lengthAB2
      REAL_E EPS

C------------------------------------------------------------------------- 
      nfaces = 1
      EPS = 1.e-12

      DO et = 0, nelts-1
         indA = cn(et,1)-1
         indB = cn(et,2)-1
                  
         xAB = xt(indB)- xt(indA)
         yAB = yt(indB)- yt(indA)
         zAB = zt(indB)- zt(indA)
         lengthAB2 = xAB*xAB + yAB*yAB + zAB*zAB

         fAB = field(indB) - field(indA)

         gradx(et) = fAB*xAB/MAX(lengthAB2,EPS)
         grady(et) = fAB*yAB/MAX(lengthAB2,EPS)
         gradz(et) = fAB*zAB/MAX(lengthAB2,EPS)
         
      ENDDO

      END

C=====================Post/Fortran/CompUnstrGrad1DF.for=====================

