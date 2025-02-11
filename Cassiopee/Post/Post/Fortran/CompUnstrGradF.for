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

C Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
C retourne le gradient defini aux centres des elts
C CAS 3D
C ============================================================================
      SUBROUTINE k6compunstrgrad(dim, npts, nelts, nedges, nnodes,   
     &                           cn, xt, yt, zt, field, fieldf, 
     &                           snx, sny, snz, surf, vol,
     &                           xint, yint, zint,
     &                           gradx, grady, gradz )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E dim             ! dimension du pb
      INTEGER_E npts            ! nb de points ds le maillage
      INTEGER_E nelts           ! nb d elts dans le maillage
      INTEGER_E nedges          ! nb de facettes par elt
      INTEGER_E nnodes          ! nb de noeuds par elt 
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E field(0:npts-1)    ! champ defini aux noeuds auquel on applique grad
      REAL_E xt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des noeuds de la grille

C_OUT
      REAL_E gradx(0:nelts-1) ! gradient de field %x
      REAL_E grady(0:nelts-1) ! gradient de field %y
      REAL_E gradz(0:nelts-1) ! gradient de field %z

C_LOCAL
      REAL_E vol(0:nelts-1)
      REAL_E surf(0:nelts-1,nedges) ! surfaces des facettes
      REAL_E snx(0:nelts-1,nedges) ! normales aux facettes %x
      REAL_E sny(0:nelts-1,nedges) ! normales aux facettes %y
      REAL_E snz(0:nelts-1,nedges) ! normales aux facettes %z
      REAL_E fieldf(0:nelts-1,nedges) ! ieme champ defini aux facettes des elts
      REAL_E xint(0:nelts-1, nedges) ! coordonnees du centre des facettes 
      REAL_E yint(0:nelts-1, nedges)
      REAL_E zint(0:nelts-1, nedges)
      INTEGER_E i, eti, fi
      REAL_E sumx, sumy, sumz
      REAL_E invvol
      
C ------------------------------------------------------------------
C     calcul de la surface + volume des elts
      call k6compunstrmetric(npts, nelts, nedges, nnodes, cn, 
     &     xt, yt, zt, xint, yint, zint, snx, sny, snz, surf, vol)
c     calcul du champ aux centres des facettes
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,  
     &     cn, field, fieldf)
      
c     calcul du gradient au centre des elts

!$OMP PARALLEL PRIVATE(eti, sumx, sumy, sumz, invvol, fi)
!$OMP DO
      DO eti = 0, nelts-1
         sumx = ZERO
         sumy = ZERO
         sumz = ZERO
         invvol = ONE/MAX(vol(eti), E_MIN_VOL)
         
         DO fi = 1, nedges
            sumx = sumx + snx(eti, fi) * fieldf(eti,fi)
            sumy = sumy + sny(eti, fi) * fieldf(eti,fi)
            sumz = sumz + snz(eti, fi) * fieldf(eti,fi)
         ENDDO
         
         gradx(eti) = sumx * invvol
         grady(eti) = sumy * invvol
         gradz(eti) = sumz * invvol
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C=====================Post/Fortran/CompUnstrGradF.for=====================

