C  ============================================================================
C  THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
C  ============================================================================
C  File   : WriteIBFileF.for
C  SVN    : $Rev$ $Date$
C  Cksum  : 
C  ============================================================================
#ifdef ELSA_DOCUMENTATION
// ============================================================================
// @Name writeIBFile
// @Memo Write the generic.ib file according to USURP         
/* @Text iblank data is added here for each block
         assume the right order is respected in writing. 
*/
// ============================================================================
Fortran writeibfile();
#endif
      SUBROUTINE writeibfile(nzone, ncellmax, ncelltot, iblank)


      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nzone           !number of blocks
      INTEGER_E ncellmax        !number of cells max of blocks (for allocation)
      INTEGER_E ncelltot(nzone) !nombre de cellules totales pour chaque bloc
      INTEGER_E iblank(ncellmax,nzone) !tableau des valeurs de iblank aux centres des cellules

C_LOCAL
      INTEGER_E i, ncell, noblk
C==============================================================================

      open(unit=3, file="generic.ib", status="unknown")
      do noblk = 1, nzone
         ncell = ncelltot(noblk)        
         write(3,*) (iblank(i,noblk),i=1,ncell)
      enddo
      close(3)
      return
      END
