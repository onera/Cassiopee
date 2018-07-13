C ============================================================================
C Triangulate a matrix by gauss pivoting
C  ===========================================================================
      SUBROUTINE k6DEC(N,NDIM,A,IP,IER)
C
      IMPLICIT NONE
C
C
C******************************************************************************
C
C  MATRIX TRIANGULARISATION BY GAUSS ELIMINATION WITH PARTIAL PIVOTING
C  INPUT..
C     N=ORDER OF MATRIX
C     NDIM=DECLARED FIRST DIMENSION OF ARRAY A.
C     A=MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J),I.LE.J=UPPER TRIANGULAR FACTOR,U.
C     A(I,J),I.GT.J=MULTIPLIERS=LOWER TRIANGULAR FACTOR, I-L.
C     IP(K),K.LT.N=INDEX OF  K-TH PIVOT ROW.
C     IER=0 IF MATRIX A IS NONSINGULAR,OR  K IF FOUND TO BE
C             SINGULAR AT STAGE K.
C  ROW INTERCHANGES ARE FINISHED IN U,ONLY PARTLY IN L.
C  USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM
C  IF IER .NE. 0,A IS SINGULAR, SOL WILL DIVIDE BY ZERO
C
C***********************************************************************
C
C_IN
      INTEGER_E N, NDIM
      REAL_E A(NDIM,N)
C_OUT
      INTEGER_E IP(N)
      INTEGER_E IER
C_LOCAL
      INTEGER_E I,J,K,NM1,M,KP1
      REAL_E T

      IER = 0
      IF (N.EQ.1) GO TO 70
      NM1 = N-1
      DO 60 K = 1, NM1
        KP1 = K+1
C  FIND THE PIVOT IN COLUMN K.SEARCH ROWS K TO N.-----------------------
        M = K
      DO 10 I = KP1, N
 10      IF (ABS(A(I,K)).GT.ABS(A(M,K))) M = I
        IP(K) = M
C  INTERCHANGE ELEMENTS IN ROWS K AND M.--------------------------------
        T = A(M,K)
        IF (M.EQ.K) GO TO 20
        A(M,K) = A(K,K)
        A(K,K) = T
 20      IF (T.EQ.0.) GO TO 80
C  STORE MULTIPLIERS IN A(I,K),I=K+1,...,N.-----------------------------
        T = 1.D0/T
        DO 30 I = KP1, N
 30       A(I,K) = -A(I,K)*T
C APPLY MULTIPLIERS TO OTHER COLUMNS OF A.------------------------------
        DO 50 J = KP1, N
          T = A(M,J)
          A(M ,J) = A(K,J)
          A(K,J) = T
          IF (T.EQ.0) GO TO 50
          DO 40 I = KP1,N
 40       A(I,J) = A(I,J)+A(I,K)*T
 50        CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N).EQ.0.) GO TO 80
      RETURN
 80   IER = K
      RETURN
      END
C
