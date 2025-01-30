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

      SUBROUTINE k6adjustcellnaturefield(
     &     ncells, blankedcells, cellnaturefield )
C
      IMPLICIT NONE
C
#include "Def/DefFortranConst.h"

C==============================================================================
C_IN 
      INTEGER_E ncells          ! size of grid
      INTEGER_E blankedcells(0:ncells-1)
C_OUT
      INTEGER_E cellNatureField(0:ncells-1)
C_LOCAL
      INTEGER_E i
C==============================================================================
!$OMP PARALLEL PRIVATE(i)
!$OMP DO     
      DO i = 0, ncells-1
         cellNatureField(i) = MIN(cellNatureField(i),
     &                            blankedcells(i))
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
C==============================================================================
