C 
      SUBROUTINE myfunction(n, a)
      
      IMPLICIT NONE
#include "Def/DefFortranConst.h"
      INTEGER_E n
      REAL_E a(n)

      WRITE(*,*) n, a(1)
      END
