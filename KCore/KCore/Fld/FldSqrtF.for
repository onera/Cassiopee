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
C  Sqrt de chaque element
C  ============================================================================
      SUBROUTINE k6fldsqrt(f, size)
C
      IMPLICIT NONE
C
C==============================================================================
C_IN
      INTEGER_E size
      REAL_E    f(0:size-1)
C_LOCAL
      INTEGER_E i
C==============================================================================
!$OMP PARALLEL PRIVATE(i)
!$OMP DO     
      DO i = 0, size-1
         f(i) = SQRT(f(i))
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
C =============== Fld/FldSqrtF.for ===============
