C  ============================================================================
C Compute a least-square matrix, make a QR decomposition and 
C         solve an upper triangular system.
C==============================================================================
      SUBROUTINE householder(N, M, x, y, z, F, weight, dim, type, coeff,
     &                       invdiagR, invBeta, A, diag, perm, b2, FS)

      IMPLICIT NONE

#include "Def/Global/DefFortran.h"
#include "Def/Global/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E N, M, dim, type
      REAL_E x(N), y(N), z(N), F(N), weight(N)
C_WORK
      REAL_E invDiagR(M), invBeta(M), A(N,M),diag(M),b2(M),FS(N)
      INTEGER_E perm(M)

C_OUT
      REAL_E coeff(M)

C_LOCAL
      INTEGER_E i, j, k, l
      REAL_E s, ak, alpha, betak, gamma, invBetak, t
      INTEGER_E M2
      REAL_E x1,y1,z1,F1
      REAL_E    smax, normF
C==============================================================================
C Express x in delta
      x1 = x(1)
      y1 = y(1)
      z1 = z(1)
      DO i = 1, N
         x(i) = x(i)-x1
         y(i) = y(i)-y1
         z(i) = z(i)-z1
      ENDDO

C Express F in delta
      F1 = F(1)
      DO i=1, N
         FS(i) = F(i)
         F(i) = (F(i)-F1)*weight(i)
      ENDDO

C     A matrix computation
      IF (dim.EQ.2) THEN
         IF (type.EQ.1) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.2) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = x(i)*x(i)*weight(i)
               A(i,5) = y(i)*y(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.3) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = x(i)*x(i)*weight(i)
               A(i,5) = y(i)*y(i)*weight(i)
               A(i,6) = x(i)*x(i)*x(i)*weight(i)
               A(i,7) = y(i)*y(i)*y(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.4) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = x(i)*x(i)*weight(i)
               A(i,5) = y(i)*y(i)*weight(i)
               A(i,6) = x(i)*x(i)*x(i)*weight(i)
               A(i,7) = y(i)*y(i)*y(i)*weight(i)  
               A(i,8) = x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,9) = y(i)*y(i)*y(i)*y(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.5) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = x(i)*x(i)*weight(i)
               A(i,5) = y(i)*y(i)*weight(i)
               A(i,6) = x(i)*x(i)*x(i)*weight(i)
               A(i,7) = y(i)*y(i)*y(i)*weight(i)  
               A(i,8) = x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,9) = y(i)*y(i)*y(i)*y(i)*weight(i)
               A(i,10) = x(i)*x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,11) = y(i)*y(i)*y(i)*y(i)*y(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.20) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = x(i)*x(i)*weight(i)
               A(i,5) = y(i)*y(i)*weight(i)
               A(i,6) = x(i)*y(i)*weight(i)
            ENDDO
         ENDIF
      ELSE
         IF (type.EQ.1) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.2) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
               A(i,5) = x(i)*x(i)*weight(i)
               A(i,6) = y(i)*y(i)*weight(i)
               A(i,7) = z(i)*z(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.3) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
               A(i,5) = x(i)*x(i)*weight(i)
               A(i,6) = y(i)*y(i)*weight(i)
               A(i,7) = z(i)*z(i)*weight(i) 
               A(i,8) = x(i)*x(i)*x(i)*weight(i)
               A(i,9) = y(i)*y(i)*y(i)*weight(i)
               A(i,10) = z(i)*z(i)*z(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.4) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
               A(i,5) = x(i)*x(i)*weight(i)
               A(i,6) = y(i)*y(i)*weight(i)
               A(i,7) = z(i)*z(i)*weight(i) 
               A(i,8) = x(i)*x(i)*x(i)*weight(i)
               A(i,9) = y(i)*y(i)*y(i)*weight(i)
               A(i,10) = z(i)*z(i)*z(i)*weight(i)  
               A(i,11) = x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,12) = y(i)*y(i)*y(i)*y(i)*weight(i)
               A(i,13) = z(i)*z(i)*z(i)*z(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.5) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
               A(i,5) = x(i)*x(i)*weight(i)
               A(i,6) = y(i)*y(i)*weight(i)
               A(i,7) = z(i)*z(i)*weight(i) 
               A(i,8) = x(i)*x(i)*x(i)*weight(i)
               A(i,9) = y(i)*y(i)*y(i)*weight(i)
               A(i,10) = z(i)*z(i)*z(i)*weight(i)  
               A(i,11) = x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,12) = y(i)*y(i)*y(i)*y(i)*weight(i)
               A(i,13) = z(i)*z(i)*z(i)*z(i)*weight(i)
               A(i,14) = x(i)*x(i)*x(i)*x(i)*x(i)*weight(i)
               A(i,15) = y(i)*y(i)*y(i)*y(i)*y(i)*weight(i)
               A(i,16) = z(i)*z(i)*z(i)*z(i)*z(i)*weight(i)
            ENDDO
         ELSE IF (type.EQ.20) THEN
            DO i = 1, N
               A(i,1) = weight(i)
               A(i,2) = x(i)*weight(i)
               A(i,3) = y(i)*weight(i)
               A(i,4) = z(i)*weight(i)
               A(i,5) = x(i)*x(i)*weight(i)
               A(i,6) = y(i)*y(i)*weight(i)
               A(i,7) = z(i)*z(i)*weight(i)
               A(i,8) = x(i)*y(i)*weight(i)
               A(i,9) = x(i)*z(i)*weight(i)
               A(i,10) = y(i)*z(i)*weight(i)
            ENDDO
         ENDIF
      ENDIF

C Normalisation droite
      normF = 0.D0
      DO i=1,N
         normF = MAX(normF, ABS(F(i)))
      ENDDO
      normF = MAX(normF, 1.e-30)

C Preconditionnement matrice
      DO j=1,M
         diag(j) = 0.D0
         DO i=1,N
            diag(j) = diag(j) + A(i,j)*A(i,j)
         ENDDO
         diag(j) = 1.D0/SQRT(diag(j))
C         WRITE(*,*) 'diag',diag(j)
      ENDDO
C
      DO i=1,N
         DO j=1,M
            A(i,j) = A(i,j)*diag(j)
         ENDDO
      ENDDO
C fin normalisation

C      WRITE(*,*) 'mat A'
C      DO i=1,N
C         WRITE(*,*) (A(i,j),j=1,M)
C      ENDDO
C      WRITE(*,*) 'Flux norm'
C      DO i=1,N
C         WRITE(*,*) F(i)
C      ENDDO

C     QR Decompostion
      M2 = M
      DO k = 1, M
         perm(k) = k
      END DO
      DO k = 1, M
         smax = 0.D0
         DO j = k, M
            s = 0.
            DO i = k, N
               s = s + A(i,perm(k))*A(i,perm(k))
            ENDDO
            IF ( s.GT.smax ) THEN
               i = perm(k)
               perm(j) = i
               smax = s
            END IF
         END DO
         s = smax
         ak = A(k,perm(k))
         alpha = -SIGN(ONE,ak)*SQRT(s)
         invDiagR(k) = 1.D0/alpha
C         WRITE(*,*) 'inv ',invDiagR(k)
         betak = s - alpha*ak
         IF (betak<E_CUTOFF) THEN
C            WRITE(*,*) 'betak = ', betak
            M2 = k-1
            GOTO 100
         ENDIF
         invBetak = 1.D0/betak
         A(k,perm(k))= ak - alpha

         DO j = k+1, M
            s = 0.
            DO i = k, N
               s = s + A(i,perm(k))*A(i,perm(j))
            ENDDO
            gamma = s*invBetak
            DO i = k, N
               A(i,perm(j)) = A(i,perm(j))-gamma*A(i,perm(k))
            ENDDO
         ENDDO
         invBeta(k) = invBetak
      ENDDO
      


C     Source term computation QT F
 100  CONTINUE
            
      DO k = 1, M2
         invBetak = invBeta(k)
         s = 0.D0
         DO j = k, N
            s = s + A(j,perm(k))*F(j)
         ENDDO
         DO i = k, N
            F(i) = F(i)-invBetak*A(i,perm(k))*s
         ENDDO
      ENDDO

C     Solver
      coeff(M2) = F(M2)*invDiagR(M2)

      DO k = M2-1, 1, -1
         s = 0.D0
         DO l = k+1, M2
            s = s + A(k,perm(l))*coeff(l)
         ENDDO
         t = F(k)-s
         coeff(k) = invDiagR(k)*t
      ENDDO
C
      DO i = M2+1, M
         coeff(i) = 0.D0
      END DO
      DO i = 1, M
         b2(perm(i)) = coeff(i)
      END DO
      DO i = 1, M
         coeff(i) = b2(i)
      END DO

C      WRITE(*,*) 'mesure de l ecart '
C      norme=0.
C      somme=0.
C      normeFS = 0.
C      DO i=1,N
C         IF (M.EQ.3) THEN
C            WRITE(*,*) i,coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)-FS(i)  
C            norme=norme+weight(i)*
C     &        (coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)-FS(i))*
C     &        (coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)-FS(i))
C            somme=somme+weight(i)
C         ELSE
C            WRITE(*,*) i,coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)+
C     &   coeff(4)*x(i)*x(i)+coeff(5)*y(i)*y(i)+coeff(6)*x(i)*y(i)-FS(i)
C            norme=norme+weight(i)*
C     &           (coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)+
C     &           coeff(4)*x(i)*x(i)+coeff(5)*y(i)*y(i)+
C     &           coeff(6)*x(i)*y(i)-FS(i))*
C     &           (coeff(1)+coeff(2)*x(i)+coeff(3)*y(i)+
C     &           coeff(4)*x(i)*x(i)+coeff(5)*y(i)*y(i)+
C     &           coeff(6)*x(i)*y(i)-FS(i))
C            somme=somme+weight(i)
C         ENDIF
C         normeFS = normeFS + FS(i)*FS(i)
C      ENDDO
      
C      WRITE(*,*) 'somme ',somme
C      WRITE(*,*) 'norme ',norme/normeFS/somme
C      IF (M.EQ.3) THEN
C         WRITE(*,*) 'ecart moyen ',coeff(1)+coeff(2)*x(N)+
C     &        coeff(3)*y(N)-FS(N)
C      ELSE
C         WRITE(*,*) 'ecart moyen ',coeff(1)+coeff(2)*x(N)+
C     &        coeff(3)*y(N)-FS(N)+
C     &        coeff(4)*x(N)*x(N)+coeff(5)*y(N)*y(N)+
C     &        coeff(6)*x(N)*y(N)
C      ENDIF

C      IF (M.EQ.6) THEN
C         WRITE(*,*) 'coeff norm ',coeff(1),coeff(2),coeff(3),coeff(4),
C     &   coeff(5), coeff(6)
C      ELSE
C         WRITE(*,*) 'coeff norm ',coeff(1),coeff(2),coeff(3)
C      ENDIF

C      IF (M.EQ.6) THEN
C         WRITE(*,*) 'coeff ',coeff(1),coeff(2),coeff(3),coeff(4),
C     &   coeff(5), coeff(6)
C      ELSE
C         WRITE(*,*) 'coeff ',coeff(1),coeff(2),coeff(3)
C      ENDIF

      DO i = 1,M
         coeff(i) = coeff(i) * diag(i)
      ENDDO

      DO i=1,N
         F(i) = FS(i)
      ENDDO
      DO i = 1, N
         x(i) = x(i) + x1
         y(i) = y(i) + y1
         z(i) = z(i) + z1
      ENDDO

      RETURN
 
      END
