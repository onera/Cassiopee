C =============================================================================
C THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
C =============================================================================
C File   : Convert/WriteBCFileF.for
C SVN    : $Rev$ $Date$
C Cksum  : 
C =============================================================================
#ifdef ELSA_DOCUMENTATION
// ============================================================================
// @Name writeBCFile
// @Memo Write the generic.bc file according to USURP         
/* @Text
*/
// ============================================================================
Fortran writebcfile();
#endif
      SUBROUTINE writebcfile( nzone, nit, njt, nkt)


      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E nzone
      INTEGER_E nit(nzone)
      INTEGER_E njt(nzone)
      INTEGER_E nkt(nzone)
C_LOCAL
      INTEGER_E n
      INTEGER_E ierr
C==============================================================================
      open(unit=3, file="generic.bc", status="unknown")     
      if (ierr .ne. 0) stop
      do n = 1, nzone
         write(unit=3,fmt=*,iostat=ierr) n,1,nit(n),1,njt(n),1,nkt(n)
      enddo
      close(3)
      return
      END
