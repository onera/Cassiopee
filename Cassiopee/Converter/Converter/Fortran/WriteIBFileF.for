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
C Write the generic.ib file according to USURP         
C iblank data is added here for each block
C         assume the right order is respected in writing 
C  ============================================================================
      SUBROUTINE k6writeibfile(nzone, ncellmax, ncelltot, iblank)


      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E nzone           !number of blocks
      INTEGER_E ncellmax        !number of cells max of blocks (for allocation)
      INTEGER_E ncelltot(nzone) !nombre de cellules totales pour chaque bloc
      INTEGER_E iblank(ncellmax,nzone) !tableau des valeurs de iblank aux centres des cellules

C_LOCAL
      INTEGER_E i, ncell, noblk
C==============================================================================
      open(unit=3,file="generic.ib",status="unknown")
      do noblk = 1,nzone
         ncell = ncelltot(noblk)        
         write(3,*) (iblank(i,noblk),i=1,ncell)
      enddo
      close(3)
      return
      END
