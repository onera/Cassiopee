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

C ============================================================================
C 2D hyperbolic mesh generator from Chan and Steger
C      What is new regarding hyper2D :
C      - delta implementation
C
C      Input :
C      - 1D line specifying geometry.
C      - 2D distribution :
C      In x a [0,1] distribution specifying the stretching on geometry.
C      In y, the real height of mesh.
C      - Boundary condition on the left
C      - Boundary condition on the right
C
C      BC can be : 
C      - 
C
C      Output :
C      - The generated mesh.
C  ===========================================================================
      SUBROUTINE k6hgcs2D(ni, nj, d)

      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"

C_IN
      INTEGER_E ni              ! nbre de points sur une ligne eta=cte
      INTEGER_E nj              ! nbre de points sur une ligne xi=cte
      REAL_E d(ni*nj)           ! distribution 2D
      REAL_E xi(ni)             ! Ligne 1D, geometrie
      REAL_E yi(ni)
      REAL_E zi(ni)
      INTEGER_E typeLeft        ! Type de CL a gauche
      INTEGER_E typeRight       ! Type de CL a droite
      
C_OUT
      REAL_E xd(ni*nj)
      REAL_E yd(ni*nj)
      REAL_E zd(ni*nj)          ! Maillage final


      END
