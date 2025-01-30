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
C  ============================================================================
C  Add sections of fields
C  Vectorized Fortran version of FldFieldF::addNewArrayInPart
      SUBROUTINE k6fldaddpart(fldIn, size, nb, deb,
     &                        fldOut, nfld)
      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E size, nfld, nb, deb
      REAL_E fldIn(0:size*nb-1)
C_INOUT
      REAL_E fldOut(0:size*nfld-1)
C_LOCAL
      INTEGER_E ifld, nt, inb
C==============================================================================
      ifld = size*(deb-1)
      inb = size*nb-1 + ifld
      DO nt = ifld, inb
         fldOut(nt) = fldOut(nt) + fldIn(nt - ifld)
      ENDDO
      END
C =============== Fld/FldAddArrayPartF.for ===============
