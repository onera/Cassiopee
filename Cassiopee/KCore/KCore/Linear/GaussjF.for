C  ============================================================================
C Matrix inversion by Gauss pivoting
C Resoud A x = b pour plusieurs seconds membres
C IN: a(n,n), b(n,m)
C OUT: a, b, ok 
C ok = 1 if matrix is not singular, 0 otherwise
C a(n,n): inverse
C b(n,m) -> solution (x) 
C==============================================================================
      SUBROUTINE kgaussj(a, n, b, m, ok, indxc, indxr, ipiv)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"  
C==============================================================================
C_IN_OUT
      INTEGER_E m, n
      REAL_E a(n, n), b(n, m)
      INTEGER_E ok
C_WORK
      INTEGER_E i, icol, irow, j, k, l, ll
      INTEGER_E indxc(n), indxr(n), ipiv(n)
      REAL_E big, dum, pivinv
C==============================================================================
      ok = 1
      icol = 0
      irow = 0

      DO j = 1, n
         ipiv(j) = 0
      ENDDO
      DO 22 i = 1, n
         big = 0.
         DO 13 j = 1, n
            IF(ipiv(j).NE.1)THEN
               DO 12 k = 1, n
                  if (ipiv(k).EQ.0) THEN
                     if (ABS(a(j,k)).GE.big)THEN
                        big = ABS(a(j,k))
                        irow = j
                        icol = k
                     ENDIF
                  ELSE IF (ipiv(k).GT.1) THEN
                     PRINT *, 'singular matrix in gaussj, 1'
                     ok = 0
                     RETURN
                  ENDIF
 12            CONTINUE
            ENDIF
 13      CONTINUE
         ipiv(icol)=ipiv(icol)+1
         if (irow.NE.icol) THEN
            DO 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         CONTINUE
            DO 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         CONTINUE
         ENDIF
         indxr(i)=irow
         indxc(i)=icol
         IF (a(icol,icol).EQ.0.) THEN
            PRINT *, 'singular matrix in gaussj, 2'
            ok = 0
            RETURN
         END IF
         pivinv = 1./a(icol,icol)
         a(icol,icol) = 1.
         DO l=1,n
            a(icol,l)=a(icol,l)*pivinv
         ENDDO
         DO l=1,m
            b(icol,l)=b(icol,l)*pivinv
         ENDDO
         DO 21 ll = 1, n
            IF (ll.NE.icol) THEN
               dum = a(ll,icol)
               a(ll,icol) = 0.
               DO l = 1, n
                  a(ll,l) = a(ll,l)-a(icol,l)*dum
               ENDDO
               DO l = 1,m
                  b(ll,l) = b(ll,l)-b(icol,l)*dum
               ENDDO
            ENDIF
 21      CONTINUE
 22   CONTINUE
      DO l = n,1,-1
         IF (indxr(l).NE.indxc(l)) THEN
            DO k = 1, n
               dum = a(k,indxr(l))
               a(k,indxr(l)) = a(k,indxc(l))
               a(k,indxc(l)) = dum
            ENDDO
         ENDIF
      ENDDO
      
      END
