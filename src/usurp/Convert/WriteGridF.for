C  ============================================================================
C  THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
C  ============================================================================
C  File   : WriteGridF.for
C  SVN    : $Rev$ $Date$
C  Cksum  : 
C  ============================================================================
#ifdef ELSA_DOCUMENTATION
// ============================================================================
// @Name writegrid
// @Memo Write the generic.grd file according to USURP         
/* @Text this contains information about mesh coordinates for a given block
*/
// ============================================================================
Fortran writegrid();
#endif
      SUBROUTINE writegrid(nzone, nptsmax, nit, njt, nkt, xb, yb, zb)


      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nzone
      INTEGER_E nptsmax
      INTEGER_E nit(nzone)
      INTEGER_E njt(nzone)
      INTEGER_E nkt(nzone)
      REAL_E xb(nptsmax,nzone)
      REAL_E yb(nptsmax,nzone)
      REAL_E zb(nptsmax,nzone)
C_LOCAL
      INTEGER_E i,n, npts 
C==============================================================================
      open(unit=3,file="generic.grd",status="unknown",
     &     action="write")

      write(unit=3,fmt=*) nzone
      write(unit=3,fmt=*) (nit(n),njt(n),nkt(n),n=1,nzone)

      do n=1,nzone
         npts = nit(n)*njt(n)*nkt(n)
         write(unit=3,fmt=*) (xb(i,n),i=1,npts), 
     &        (yb(i,n),i=1,npts), 
     &        (zb(i,n),i=1,npts)
      enddo
      close(3)
      return
      END
