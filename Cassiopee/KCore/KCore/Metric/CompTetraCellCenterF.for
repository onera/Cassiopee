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
C  Compute barycenter of tetra cells.
C =============================================================================
      SUBROUTINE  k6comptetracellcenter(npts, nelts, cn, xt, yt, zt, 
     &                                  bary)
      
      IMPLICIT NONE
      
#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de noeuds
      INTEGER_E nelts           ! nb d elements
      INTEGER_E cn(0:nelts-1,4) ! connectivite Elt->noeuds
      REAL_E    xt(0:npts-1)    ! coord %x des noeuds
      REAL_E    yt(0:npts-1)    ! coord %y des noeuds
      REAL_E    zt(0:npts-1)    ! coord %z des noeuds

C_OUT
      REAL_E    bary(0:nelts-1,3) ! cell centers

C_LOCAL
      INTEGER_E et, ind1, ind2, ind3, ind4

C =============================================================================
C
!$OMP PARALLEL PRIVATE(et,ind1,ind2,ind3,ind4)
!$OMP DO
      DO et = 0, nelts-1
         ind1 = cn(et,1)-1
         ind2 = cn(et,2)-1
         ind3 = cn(et,3)-1
         ind4 = cn(et,4)-1
                  
         bary(et,1) = ONE_FOURTH * ( 
     &        xt(ind1) + xt(ind2) + xt(ind3) + xt(ind4) )
         bary(et,2) = ONE_FOURTH * ( 
     &        yt(ind1) + yt(ind2) + yt(ind3) + yt(ind4) )
         bary(et,3) = ONE_FOURTH * ( 
     &        zt(ind1) + zt(ind2) + zt(ind3) + zt(ind4) )         
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
C ======================= Metric/CompTetraCellCenterF.for ====================
