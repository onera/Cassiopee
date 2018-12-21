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
C  Copy a field as a subpart of another field.
C  E_Float arrays
C  ============================================================================
      SUBROUTINE k6fldcopyfrom(deb, sizerhs, sizelhs,
     &                         nfld, rhs, lhs)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E deb, sizerhs, sizelhs, nfld
      REAL_E    rhs(0:sizerhs-1, 1:nfld)
C_OUT
      REAL_E    lhs(0:sizelhs-1, 1:nfld)
C_LOCAL
      INTEGER ifld, no
C==============================================================================
      DO ifld = 1, nfld
         DO no = 0, sizerhs-1
            lhs(no+deb, ifld) = rhs(no, ifld)
         END DO
      END DO
      END
C =============== Fld/FldCopyFromF.for ===============
