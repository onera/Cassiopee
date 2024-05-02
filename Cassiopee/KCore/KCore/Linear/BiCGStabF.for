     SUBROUTINE bicgstab(np, ia, ianp, ja, aa, b, x)
      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E np              ! Implicit number nodes
      INTEGER_E ianp            ! ia(np) 
      INTEGER_E ia(0:np)
      INTEGER_E ja(0:ianp-1)    ! CSR array
      REAL_E aa(0:ianp-1)       ! matrix (to be inversed)
      REAL_E b(0:np-1)          ! Source term of interpolation system

C_OUT
      REAL_E x(0:np-1)          ! Implicit point field

C_LOCAL
      INTEGER_E i, it
      INTEGER_E itmax           ! Max number of iterations accepted for convergence
      REAL_E invnormer0, normer 
      REAL_E prodscal1, prodscal2
      REAL_E rAp, sAs, AsAs
      REAL_E alpha, beta , gamma
      REAL_E erreur
      REAL_E eps
      REAL_E r(0:np-1), rt(0:np-1), p(0:np-1)
      REAL_E s(0:np-1)
      REAL_E As(0:np-1), Ap(0:np-1)  
      REAL_E normit,unsurnormit
C==============================================================================
      itmax = 200
      eps = 1.e-6

      DO i = 0, np-1
         x(i) = 0.D0
         r(i) = b(i) 
         rt(i) = r(i)
         p(i) = r(i)
      END DO

      CALL prodscal(r, r, np, prodscal1)

C NORMALISATION
      normit = 0.
      IF ( (prodscal1.LT.1.e-4) .AND. (prodscal1.GT.1.e-20) ) THEN
         normit=prodscal1
         unsurnormit=1./normit
         DO i=0, np-1
            b(i)=b(i)*unsurnormit
         ENDDO
      ENDIF

      IF (ABS(prodscal1).LT.1.e-20) THEN
         RETURN
      ENDIF

      invnormer0 = 1.D0/SQRT(prodscal1)

      iterations: DO it = 1, itmax
      IF (it == itmax) THEN
         PRINT *,'BiCGStab: nb iteration >= itmax=',itmax
         EXIT iterations
      END IF
      
      CALL multmatvec(aa, p, Ap, ia, ia(np), ja, np)
      CALL prodscal(Ap, rt, np, rAp)
      alpha = prodscal1/rAp
      
      DO i = 0, np-1       
         s(i) = r(i)-alpha*Ap(i)
      END DO
      CALL multmatvec(aa, s, As, ia, ia(np), ja, np)

      CALL prodscal(As, s, np, sAs)
      CALL prodscal(As, As, np, AsAs)
      beta = sAs/AsAs
      DO i = 0, np-1
         x(i) = x(i)+alpha*p(i)+beta*s(i)
         r(i) = s(i)-beta*As(i)
      END DO
      CALL prodscal(r, r, np, normer)
      normer = SQRT(normer)
      
      erreur = normer*invnormer0
      IF (erreur <= eps) THEN
C         PRINT *, 'it  = ', it

         IF (normit.NE.0) THEN
            DO i=0,np-1
               x(i)=x(i)*normit
            ENDDO
         ENDIF

         RETURN
      END IF

      CALL prodscal(r,rt, np, prodscal2)
      gamma = (prodscal2/prodscal1)*(alpha/beta)

      prodscal1 = prodscal2

      DO i = 0, np-1
         p(i) = r(i)+gamma*(p(i)-beta*Ap(i))
      END DO

      END DO iterations

      IF (normit.NE.0.) THEN
         DO i=0,np-1
            x(i)=x(i)*normit
         ENDDO
      ENDIF

      END 

C==============================================================================
      SUBROUTINE multmatvec(aa, x, r, ia, ianp, ja, np)
      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
C     Matrix*vector product 
C==============================================================================
C_IN
      INTEGER_E np                  ! Implicit number nodes
      INTEGER_E ia(0:np)           
      INTEGER_E ianp                ! ia(np)
      INTEGER_E ja(0:ianp-1)        ! CSR array
      REAL_E aa(0:ianp-1)           ! Interpolation coefficient matrix 
      REAL_E x(0:np-1)
C_OUT
      REAL_E r(0:np-1)
C_LOCAL
      INTEGER_E i, j, k
      REAL_E somme
C==============================================================================
      DO i = 0, np-1
         somme = 0.       
         DO k = ia(i), ia(i+1)-1
            j = ja(k)-1
            somme = somme + aa(k)*x(j)
         END DO
         r(i) = somme
      END DO
      
      END 

C==============================================================================
      SUBROUTINE multmatvecpara(aa, x, r, ia, ianp, ja, np, npglobal)
      IMPLICIT NONE

C     SUBROUTINE multmatvecf(lengthi, lenghtj, x, array) 

#include "Def/Global/DefFortran.h"
C     Matrix*vector product 
C==============================================================================
C_IN
      INTEGER_E np                  ! Implicit number nodes
      INTEGER_E npGlobal                  ! Implicit number nodes
      INTEGER_E ia(0:np)           
      INTEGER_E ianp                ! ia(np)
      INTEGER_E ja(0:ianp-1)        ! CSR array
      REAL_E aa(0:ianp-1)           ! Interpolation coefficient matrix 
      REAL_E x(0:npglobal-1)
C_OUT
      REAL_E r(0:np-1)
C_LOCAL
      INTEGER_E i, j, k
      REAL_E somme
C==============================================================================
      
      DO i = 0, np-1
         somme = 0.D0   
         DO k = ia(i), ia(i+1)-1
            j = ja(k)-1
            somme = somme + aa(k)*x(j)
         END DO
         r(i) = somme
      END DO
      
      END 

C==============================================================================
C Calcul a.b
C IN : a : vecteur a(n)
C IN : b : vecteur b(n)
C IN : n : taille des vecteurs
C OUT : prod : a.b
C==============================================================================
      SUBROUTINE prodscal(a, b, n, prod) 
      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
C==============================================================================
C_IN
      INTEGER_E n
      REAL_E a(n)
      REAL_E b(n)

C_OUT
      REAL_E prod

C_LOCAL
      INTEGER_E i
C==============================================================================
      prod = 0.D0
      DO i = 1, n
         prod = prod + a(i)*b(i)
      END DO

      END 
C ===== BiCGStabF.for =====
