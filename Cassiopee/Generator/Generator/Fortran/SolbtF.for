C ============================================================================
C Solution of block diagonal linear system
C  ==========================================================================
      SUBROUTINE k6SOLBT (M,N,A,B,C,Y,IP)
C
      IMPLICIT NONE
C
C***********************************************************************
C
C  SOLUTION OF BLOCK-TRIDIAGONAL LINEAR SYSTEM
C COEFFICIENT MATRIX MUST HAVE BEEN PREVIOUSLY PROCESSED BY DECBT.
C INPUT
C     M=ORDER OF EACH BLOCK.
C     N=NUMBER OF BLOCKS IN EACH DIRECTION OF MATRIX
C A,B,C=M BY M BY N ARRAYS CONTAINING BLOCK LU DECOMPOSITION
C       OF COEFFICIENT MATRIX FROM DECBT.
C     IP=M BY N INTEGER ARRAY OF PIVOT INFORMATION FROM DECBT .
C     Y=ARRAY OF LENGHT M*N CONTAINING THE RIGHT-HAND SIDE VECTOR
C        (TREATED AS AN M BY N ARRAY HERE).
C OUTPUT..
C     Y=SOLUTION VECTOR OF LENGTH M*N.
C SOLBT MAKES CALLS TO SUBROUTINE SOL(M,MO,A,Y,IP)
C FOR SOLUTION OF M BY M LINEAR SYSTEMS.
C
C***********************************************************************
C
C_IN
      INTEGER_E M,N
      REAL_E A(M,M,N),B(M,M,N),C(M,M,N),Y(M,N)
      INTEGER_E IP(M,N)
C_LOCAL
      INTEGER_E NM1,NM2,KM1,KP1,KB,I,J,K
      REAL_E DP
C
      NM1 = N-1
      NM2 = N-2
C FORWARD SOLUTION SWEEP
      CALL k6SOL(M,M,A,Y,IP)
      DO 30 K=2,NM1
        KM1=K-1
        DO 20 I=1,M
          DP=0.D0
          DO 10 J=1,M
            DP=DP+C(I,J,K)*Y(J,KM1)
 10       CONTINUE
          Y(I,K)=Y(I,K)-DP
 20     CONTINUE
        CALL k6SOL (M,M,A(1,1,K),Y(1,K),IP(1,K))
 30   CONTINUE
      DO 50 I=1,M
        DP=0.D0
        DO 40 J=1,M
          DP=DP+C(I,J,N)*Y(J,NM1)+B(I,J,N)*Y(J,NM2)
 40     CONTINUE
        Y(I,N)=Y(I,N)-DP
 50   CONTINUE
      CALL k6SOL (M,M,A(1,1,N),Y(1,N),IP(1,N))
C BACKWARD SOLUTION SWEEP.
      DO 80 KB=1,NM1
        K=N-KB
        KP1=K+1
        DO 70 I=1,M
          DP=0.D0
          DO 60 J=1,M
            DP=DP+B(I,J,K)*Y(J,KP1)
 60       CONTINUE
          Y(I,K)=Y(I,K)-DP
 70     CONTINUE
 80   CONTINUE
      DO 100 I=1,M
        DP=0.D0
        DO 90 J=1,M
          DP=DP+C(I,J,1)*Y(J,3)
 90     CONTINUE
        Y(I,1)=Y(I,1)-DP
 100  CONTINUE
      RETURN
      END
C
