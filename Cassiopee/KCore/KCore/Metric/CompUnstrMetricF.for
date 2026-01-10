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
C  Calcul du volume de toutes les cellules  et des surface des interfaces
C  CAS NON STRUCTURE
C  ============================================================================
      SUBROUTINE  k6compunstrmetric(npts, nelts, nedges, nnodes, cn, 
     &     coordx, coordy, coordz, xint, yint, zint,
     &     snx, sny, snz, surf, vol)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb de pts dans le maillage
      INTEGER_E nelts           ! nb d elements
      INTEGER_E nedges          ! nb de facettes par elemt
      INTEGER_E nnodes          ! nb de noeuds par elemt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E coordx(0:npts-1)   !coordonnees des pts de la grille    
      REAL_E coordy(0:npts-1)   !coordonnees des pts de la grille    
      REAL_E coordz(0:npts-1)   !coordonnees des pts de la grille    

C_OUT
      REAL_E vol(0:nelts-1)
      REAL_E surf(0:nelts-1,nedges) ! surfaces des facettes
      REAL_E snx(0:nelts-1,nedges) ! normales aux facettes %x
      REAL_E sny(0:nelts-1,nedges) ! normales aux facettes %y
      REAL_E snz(0:nelts-1,nedges) ! normales aux facettes %z
      REAL_E xint(0:nelts-1, nedges) ! coordonnees du centre des facettes 
      REAL_E yint(0:nelts-1, nedges) ! Tableau local
      REAL_E zint(0:nelts-1, nedges)
C_LOCAL
      INTEGER_E i, edge
      REAL_E voli
      REAL_E ONE_THIRD

C==============================================================================
      ONE_THIRD = ONE/ THREE

C     calcul des centres des facettes
      call k6unstructcenterint(npts, nelts, nedges, nnodes, cn, 
     &     coordx, coordy, coordz, xint, yint, zint)

C     calcul des normales / des surfaces
      call k6unstructsurf(npts, nelts, nedges, nnodes, cn,
     &     coordx, coordy, coordz, snx, sny, snz, surf) 
!$OMP PARALLEL PRIVATE(i,edge,voli)
!$OMP DO
      DO i = 0, nelts-1
         voli = ZERO
         DO edge = 1, nedges
            voli = voli + xint(i,edge)*snx(i,edge) + 
     &           yint(i,edge)*sny(i,edge) + zint(i,edge)*snz(i,edge)       
         ENDDO
         
         vol(i) = ONE_THIRD * voli
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END
C========================= Metric/CompUnstrMetricF.for ======================
