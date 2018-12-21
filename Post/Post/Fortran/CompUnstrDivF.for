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

C Calcul de la divergence d'un champ defini aux noeuds d une grille non structuree
C retourne la divergence defini aux centres des elts
C CAS 3D
C ============================================================================
      SUBROUTINE k6compunstrdiv(dim, npts, nelts, nedges, nnodes,
     &                          cn, xt, yt, zt, fieldX, fieldY, fieldZ,
     &                          fieldfx, fieldfy, fieldfz,
     &                          snx, sny, snz, surf, vol,
     &                          xint, yint, zint,
     &                          div)

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
      REAL_E fieldX(0:npts-1)    ! champ (x-comp) defini aux noeuds auquel on applique div
      REAL_E fieldY(0:npts-1)    ! champ (y-comp) defini aux noeuds auquel on applique div
      REAL_E fieldZ(0:npts-1)    ! champ (z-comp) defini aux noeuds auquel on applique div
      REAL_E xt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des noeuds de la grille

C_OUT
      REAL_E div(0:nelts-1) ! divergence du champ

C_LOCAL
      REAL_E vol(0:nelts-1)
      REAL_E surf(0:nelts-1,nedges) ! surfaces des facettes
      REAL_E snx(0:nelts-1,nedges) ! normales aux facettes %x
      REAL_E sny(0:nelts-1,nedges) ! normales aux facettes %y
      REAL_E snz(0:nelts-1,nedges) ! normales aux facettes %z
      REAL_E fieldfx(0:nelts-1,nedges) ! x-champ defini aux facettes des elts
      REAL_E fieldfy(0:nelts-1,nedges) ! y-champ defini aux facettes des elts
      REAL_E fieldfz(0:nelts-1,nedges) ! z-champ defini aux facettes des elts
      REAL_E xint(0:nelts-1, nedges) ! coordonnees du centre des facettes
      REAL_E yint(0:nelts-1, nedges)
      REAL_E zint(0:nelts-1, nedges)
      INTEGER_E i, eti, fi
      REAL_E gradxx, gradyy, gradzz
      REAL_E invvol

C ------------------------------------------------------------------
C     calcul de la surface + volume des elts
      call k6compunstrmetric(npts, nelts, nedges, nnodes, cn,
     &     xt, yt, zt, xint, yint, zint, snx, sny, snz, surf, vol)

c     calcul des champs aux centres des facettes
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,
     &     cn, fieldX, fieldfx)
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,
     &     cn, fieldY, fieldfy)
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,
     &     cn, fieldZ, fieldfz)

c     calcul de la divergence au centre des elts

!$OMP PARALLEL PRIVATE(eti, gradxx, gradyy, gradzz, invvol, fi)
!$OMP DO
      DO eti = 0, nelts-1
         gradxx = ZERO
         gradyy = ZERO
         gradzz = ZERO
         DO fi = 1, nedges
            gradxx = gradxx + snx(eti, fi) * fieldfx(eti,fi)
            gradyy = gradyy + sny(eti, fi) * fieldfy(eti,fi)
            gradzz = gradzz + snz(eti, fi) * fieldfz(eti,fi)
         ENDDO

         invvol = ONE/vol(eti)
         div(eti) = (gradxx + gradyy + gradzz)* invvol
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C=====================Post/Fortran/CompUnstrDivF.for======================

