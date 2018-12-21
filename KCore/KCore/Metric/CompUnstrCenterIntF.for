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

C  ============================================================================
C  Calcul des centres des interfaces pour des mailles non structurees
C  ============================================================================
      SUBROUTINE  k6unstructcenterint(npts, nelts, nedges, nnodes, cn, 
     &                                xt, yt, zt, xint, yint, zint)
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de pts dans le maillage
      INTEGER_E nelts           ! nb d elements
      INTEGER_E nedges          ! nb de facettes par elemt
      INTEGER_E nnodes          ! nb de noeuds par elemt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E xt(0:npts-1)       !coordonnees des pts de la grille 
      REAL_E yt(0:npts-1)       !coordonnees des pts de la grille 
      REAL_E zt(0:npts-1)       !coordonnees des pts de la grille 

C_OUT 
      REAL_E xint(0:nelts-1, nedges) ! coordonnees du centre des facettes 
      REAL_E yint(0:nelts-1, nedges)
      REAL_E zint(0:nelts-1, nedges)
C_LOCAL 
      INTEGER_E i, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
      REAL_E ONE_THIRD
      REAL_E x1, x2, x3, x4, x5, x6, x7, x8
      REAL_E y1, y2, y3, y4, y5, y6, y7, y8
      REAL_E z1, z2, z3, z4, z5, z6, z7, z8

C ============================================================================
      ONE_THIRD = ONE/ THREE

      IF (nedges .eq. 1 .and. nnodes .eq. 3) THEN ! TRI
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,x1,x2,x3,y1,y2,y3,z1,z2,z3)
!$OMP DO
         DO i = 0, nelts-1
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            
            xint(i,1) = ONE_THIRD * ( x1 + x2 + x3 )
            yint(i,1) = ONE_THIRD * ( y1 + y2 + y3 )
            zint(i,1) = ONE_THIRD * ( z1 + z2 + z3 )
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE IF (nedges .eq. 1 .and. nnodes .eq. 4) THEN ! QUAD
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,ind4,x1,x2,x3,x4)
!$OMP& PRIVATE(y1,y2,y3,y4,z1,z2,z3,z4)
!$OMP DO
         DO i = 0, nelts-1
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            ind4 = cn(i,4)-1 
            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            x4 = xt(ind4)

            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            y4 = yt(ind4)

            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            z4 = zt(ind4)

            xint(i,1) = ONE_FOURTH * ( x1 + x2 + x3 + x4 )
            yint(i,1) = ONE_FOURTH * ( y1 + y2 + y3 + y4 )
            zint(i,1) = ONE_FOURTH * ( z1 + z2 + z3 + z4 )
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE IF (nedges .eq. 4 .and. nnodes .eq. 4) THEN ! TETRA
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,ind4,x1,x2,x3,x4)
!$OMP& PRIVATE(y1,y2,y3,y4,z1,z2,z3,z4)
!$OMP DO
         DO i = 0, nelts-1
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            ind4 = cn(i,4)-1
            
            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            x4 = xt(ind4)

            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            y4 = yt(ind4)

            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            z4 = zt(ind4)

C     facette 123
            xint(i,1) = ONE_THIRD * ( x1 + x2 + x3 )
            yint(i,1) = ONE_THIRD * ( y1 + y2 + y3 )
            zint(i,1) = ONE_THIRD * ( z1 + z2 + z3 )
C     facette 124
            xint(i,2) = ONE_THIRD * ( x1 + x2 + x4 )
            yint(i,2) = ONE_THIRD * ( y1 + y2 + y4 )
            zint(i,2) = ONE_THIRD * ( z1 + z2 + z4 )
C     facette 234
            xint(i,3) = ONE_THIRD * ( x2 + x3 + x4 )
            yint(i,3) = ONE_THIRD * ( y2 + y3 + y4 )
            zint(i,3) = ONE_THIRD * ( z2 + z3 + z4 )  
C     facette 134
            xint(i,4) = ONE_THIRD * ( x1 + x3 + x4 )
            yint(i,4) = ONE_THIRD * ( y1 + y3 + y4 )
            zint(i,4) = ONE_THIRD * ( z1 + z3 + z4 )
         ENDDO            
!$OMP END DO
!$OMP END PARALLEL
      ELSE IF (nedges .eq. 5 .and. nnodes .eq. 5) THEN ! PYRA
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,ind4,ind5,x1,x2,x3,x4,x5)
!$OMP& PRIVATE(y1,y2,y3,y4,y5,z1,z2,z3,z4,z5)
!$OMP DO
         DO i = 0, nelts-1
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            ind4 = cn(i,4)-1
            ind5 = cn(i,5)-1

            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            x4 = xt(ind4)
            x5 = xt(ind5)

            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            y4 = yt(ind4)
            y5 = yt(ind5)

            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            z4 = zt(ind4)
            z5 = zt(ind5)

C     facette 1234 : quad
            xint(i,1) = ONE_FOURTH * ( x1 + x2 + x3 + x4 )
            yint(i,1) = ONE_FOURTH * ( y1 + y2 + y3 + y4 )
            zint(i,1) = ONE_FOURTH * ( z1 + z2 + z3 + z4 )
C     facette 125
            xint(i,2) = ONE_THIRD * ( x1 + x2 + x5 )
            yint(i,2) = ONE_THIRD * ( y1 + y2 + y5 )
            zint(i,2) = ONE_THIRD * ( z1 + z2 + z5 )
C     facette 235
            xint(i,3) = ONE_THIRD * ( x2 + x3 + x5 )
            yint(i,3) = ONE_THIRD * ( y2 + y3 + y5 )
            zint(i,3) = ONE_THIRD * ( z2 + z3 + z5 )  
C     facette 345
            xint(i,4) = ONE_THIRD * ( x5 + x3 + x4 )
            yint(i,4) = ONE_THIRD * ( y5 + y3 + y4 )
            zint(i,4) = ONE_THIRD * ( z5 + z3 + z4 )
C     facette 415
            xint(i,5) = ONE_THIRD * ( x1 + x5 + x4 )
            yint(i,5) = ONE_THIRD * ( y1 + y5 + y4 )
            zint(i,5) = ONE_THIRD * ( z1 + z5 + z4 )
         ENDDO            
!$OMP END DO
!$OMP END PARALLEL
      ELSE IF (nedges .eq. 6 .and. nnodes .eq. 8) THEN ! HEXA
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8)
!$OMP& PRIVATE(x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8)
!$OMP& PRIVATE(z1,z2,z3,z4,z5,z6,z7,z8)
!$OMP DO
         DO i = 0, nelts-1
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            ind4 = cn(i,4)-1
            ind5 = cn(i,5)-1
            ind6 = cn(i,6)-1
            ind7 = cn(i,7)-1
            ind8 = cn(i,8)-1

            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            x4 = xt(ind4)
            x5 = xt(ind5)
            x6 = xt(ind6)
            x7 = xt(ind7)
            x8 = xt(ind8)

            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            y4 = yt(ind4)
            y5 = yt(ind5)
            y6 = yt(ind6)
            y7 = yt(ind7)
            y8 = yt(ind8)

            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            z4 = zt(ind4)
            z5 = zt(ind5)
            z6 = zt(ind6)
            z7 = zt(ind7)
            z8 = zt(ind8)
C     facette 1234
            xint(i,1) = ONE_FOURTH * ( x1 + x2 + x3 + x4 )
            yint(i,1) = ONE_FOURTH * ( y1 + y2 + y3 + y4 )
            zint(i,1) = ONE_FOURTH * ( z1 + z2 + z3 + z4 )
            
C     facette 5678
            xint(i,2) = ONE_FOURTH * ( x5 + x6 + x7 + x8 )
            yint(i,2) = ONE_FOURTH * ( y5 + y6 + y7 + y8 )
            zint(i,2) = ONE_FOURTH * ( z5 + z6 + z7 + z8 )
C     facette 1485
            xint(i,3) = ONE_FOURTH * ( x1 + x4 + x8 + x5 )
            yint(i,3) = ONE_FOURTH * ( y1 + y4 + y8 + y5 )
            zint(i,3) = ONE_FOURTH * ( z1 + z4 + z8 + z5 )            
C     facette 2376
            xint(i,4) = ONE_FOURTH * ( x2 + x3 + x7 + x6 )
            yint(i,4) = ONE_FOURTH * ( y2 + y3 + y7 + y6 )
            zint(i,4) = ONE_FOURTH * ( z2 + z3 + z7 + z6 ) 
C     facette 1265 
            xint(i,5) = ONE_FOURTH * ( x1 + x2 + x6 + x5 )
            yint(i,5) = ONE_FOURTH * ( y1 + y2 + y6 + y5 )
            zint(i,5) = ONE_FOURTH * ( z1 + z2 + z6 + z5 )             
C     facette 4378
            xint(i,6) = ONE_FOURTH * ( x4 + x3 + x7 + x8 )
            yint(i,6) = ONE_FOURTH * ( y4 + y3 + y7 + y8 )
            zint(i,6) = ONE_FOURTH * ( z4 + z3 + z7 + z8 ) 
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE IF (nedges .eq. 5 .and. nnodes .eq. 6) THEN ! PENTA
!$OMP PARALLEL PRIVATE(i,ind1,ind2,ind3,ind4,ind5,ind6)
!$OMP& PRIVATE(x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6)
!$OMP& PRIVATE(z1,z2,z3,z4,z5,z6)
!$OMP DO
         DO i = 0, nelts-1
            
            ind1 = cn(i,1)-1
            ind2 = cn(i,2)-1
            ind3 = cn(i,3)-1
            ind4 = cn(i,4)-1
            ind5 = cn(i,5)-1
            ind6 = cn(i,6)-1

            x1 = xt(ind1)
            x2 = xt(ind2)
            x3 = xt(ind3)
            x4 = xt(ind4)
            x5 = xt(ind5)
            x6 = xt(ind6)

            y1 = yt(ind1)
            y2 = yt(ind2)
            y3 = yt(ind3)
            y4 = yt(ind4)
            y5 = yt(ind5)
            y6 = yt(ind6)

            z1 = zt(ind1)
            z2 = zt(ind2)
            z3 = zt(ind3)
            z4 = zt(ind4)
            z5 = zt(ind5)
            z6 = zt(ind6)

C     facette triangle 1 : 123
            xint(i,1) = ONE_THIRD * ( x1 + x2 + x3 )
            yint(i,1) = ONE_THIRD * ( y1 + y2 + y3 )
            zint(i,1) = ONE_THIRD * ( z1 + z2 + z3 )
C     facette triangle 2 : 456
            xint(i,2) = ONE_THIRD * ( x4 + x5 + x6 )
            yint(i,2) = ONE_THIRD * ( y4 + y5 + y6 )
            zint(i,2) = ONE_THIRD * ( z4 + z5 + z6 )
C     troisieme facette quad 1254 
            xint(i,3) = ONE_FOURTH * ( x1 + x2 + x5 + x4 )
            yint(i,3) = ONE_FOURTH * ( y1 + y2 + y5 + y4 )
            zint(i,3) = ONE_FOURTH * ( z1 + z2 + z5 + z4 )
C     quatrieme facette quad 2365
            xint(i,4) = ONE_FOURTH * ( x2 + x3 + x6 + x5 )
            yint(i,4) = ONE_FOURTH * ( y2 + y3 + y6 + y5 )
            zint(i,4) = ONE_FOURTH * ( z2 + z3 + z6 + z5 )
C     cinquieme facette quad 3146 
            xint(i,5) = ONE_FOURTH * ( x3 + x1 + x4 + x6 )
            yint(i,5) = ONE_FOURTH * ( y3 + y1 + y4 + y6 )
            zint(i,5) = ONE_FOURTH * ( z3 + z1 + z4 + z6 )
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ELSE 
         write(*,*) 'Error in KCore/Metric/CompUnstrSurfF.for '
         write(*,*) 'Unknown type of elements.'
         STOP
      ENDIF
      END
C =========================== CompUnstructCenterIntF.for =====================
