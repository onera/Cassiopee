C ============================================================================
C THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
C ============================================================================
C File  : Generator/Fortran/SolF.for
C SVN   : $Rev$ $Date$
C Cksum : 
C ============================================================================
C Solution of linear system
C  ==========================================================================
      SUBROUTINE k6SOL(N,NDIM,A,B,IP)
C
      IMPLICIT NONE
C
C***********************************************************************
C
C  SOLUTION OF LINEAR SYSTEME A=X =B USING OUTPUT OF DEC
C  INPUT..
C     NDIM=DECLARED FIRST DIMENSION OF ARRAY A.
C     N=ORDER OF MATRIX.
C     A=TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C     B=RIGHT HAND SIDE VECTOR.
C     IP=PIVOT INFORMATION VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER.NE.O.
C  OUTPUT..
C     B=SOLUTION VECTOR,X.
C
C***********************************************************************
C
C_IN
      INTEGER_E N, NDIM
      REAL_E A(NDIM,N), B(N)
C_OUT
      INTEGER_E IP(N)
C_LOCAL
      INTEGER_E NM1,K,KP1,M,I,KB,KM1
      REAL_E T
C
      IF (N.EQ.1) GO TO 50
      NM1 = N-1
C  APPLY ROW PERMUTATIONS AND MULTIPLIERS TO B.-------------------------
      DO K = 1,NM1
         KP1 = K+1
         M = IP(K)
         T = B(M)
         B(M) = B(K)
         B(K) = T
         DO I = KP1,N
            B(I) = B(I)+A(I,K)*T
         ENDDO
      ENDDO
C  BACK SOLVE.----------------------------------------------------------
      DO KB=1,NM1
         KM1 = N-KB
         K= KM1+1
         B(K)=B(K)/A(K,K)
         T=-B(K)
         DO I=1,KM1
            B(I)=B(I)+A(I,K)*T
         ENDDO
      ENDDO
 50   B(1) = B(1)/A(1,1)
      RETURN
      END
C
