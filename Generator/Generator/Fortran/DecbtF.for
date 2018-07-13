C ============================================================================
C block tridiagonal block decomposition
C ============================================================================
      SUBROUTINE k6DECBT(M,N,A,B,C,IP,IER)
C
      IMPLICIT NONE
C
C*****************************************************************************
C
C BLOCK-TRIDIAGONAL MATRIX DECOMPOSITION ROUTINE.
C WRITTEN BY A.C. HINDMARSH,
C THE INPUT MATRIX CONTAINS THREE BLOCKS OF ELEMENTS IN EACH BLOCK ROW,
C INCLUDING BLOCKS IN THE (1,3) AND(N,N-2) BLOCK POSITIONS.
C DECBT USES BLOCK GAUSS ELIMINATION AND SUBROUTINES DEC AND SOL
C FOR SOLUTION OF BLOCKS,PARTIAL PIVOTING IS DONE WITHIN
C BLOCK ROWS ONLY,
C INPUT..
C     M=ORDER OF EACH BLOCK.
C     N=NUMBER OF BLOCKS IN EACH DIRECTION OF THE MATRIX.
C        N MUST BE 4 OR MORE. THE COMPLETE MATRIX HAS ORDER M*N.
C     A=M BY M BY N ARRAY CONTAINING DIAGONAL BLOCKS.
C         A(I,J,K) CONTAINS THE (I,J) ELEMENT OF THE K-TH BLOCK.
C     B=M BY M BY N ARRAY CONTAINING THE SUPER-DIAGONAL BLOCKS
C         (IN B(,,K) FOR K=1,...,N-1) AND THE BLOCK IN THE (N,N-2)
C         BLOCK POSITION (IN B(,,N)).
C     C=M BY M BY N ARRAY  CONTAINING THE SUBDIAGONAL BLOCKS
C         (IN C(,,K) FOR K=2,3,....N) AND THE BLOCK IN THE
C         (1,3) BLOCK POSITION (IN C(,,1)).
C     IP=INTEGER ARRAY OF LENGTH M*N FOR WORKING  STORAGE.
C OUTPUT..
C A,B,C=M BY M BYN ARRAYS CONTAINING THE BLOCK LU DECOMPOSITION
C         OF THE INPUT MATRIX.
C     IP=M BY N ARRAY OF PIVOT INFORMATION.IP(,,K) CONTAINS
C         INFORMATION FOR THE K-TH DIAGONAL BLOCK.
C     IER=0 IF NO TROUBLE OCCURRED,OR
C        =-1 IF THE INPUT VALUE OF M OR N WAS ILLEGAL,OR
C        =K IF A SINGULAR MATRIX WAS FOUND IN THE K-TH DIAGONAL BLOCK.
C DECBT CALLS SUBROUTINES  DEC(M,MO,A,IP,IER) AND SOL(M,MO,A,Y,IP)
C FOR SOLUTION OF M BY M LINEAR SYSTEMS.
C
C****************************************************************************
C
C_IN     
      INTEGER_E M,N
      REAL_E A(M,M,N), B(M,M,N), C(M,M,N)
C_OUT
      INTEGER_E IER
      INTEGER_E IP(M,N)
C_LOCAL
      INTEGER_E NM1, NM2, KM1
      INTEGER_E I,J,K,L
      REAL_E DP

      IF (M.LT.1.OR.N.LT.4) GO TO 210
      NM1 = N-1
      NM2 = N-2
C PROCESS THE FIRST BLOCK-ROW.------------------------------------------
      
      CALL k6DEC (M,M,A,IP,IER)
      K = 1
      IF (IER .NE.0) GO TO 200
      DO 10 J = 1,M
        CALL k6SOL (M,M,A,B(1,J,1),IP)
        CALL k6SOL (M,M,A,C(1,J,1),IP)
 10     CONTINUE
C ADJUST B(,,2).--------------------------------------------------------
      DO 40 J=1,M
        DO 30 I=1,M
          DP = 0.D0
          DO 20 L=1,M
 20         DP=DP+C(I,L,2)*C(L,J,1)
 30         B(I,J,2)=B(I,J,2)-DP
40      CONTINUE
C MAIN LOOP. PROCESS BLOCK-ROWS 2 TO N-1.-------------------------------
      DO 100 K=2,NM1
        KM1=K-1
        DO 70 J=1,M
          DO 60 I=1,M
            DP=0.D0
            DO 50 L=1,M
 50           DP=DP+C(I,L,K)*B(L,J,KM1)
 60         A(I,J,K)=A(I,J,K)-DP
 70       CONTINUE
 
      CALL k6DEC (M,M,A(1,1,K),IP(1,K),IER)
      
      IF (IER.NE.0) GO TO 200
      DO 80 J=1,M
 80     CALL k6SOL (M,M,A(1,1,K),B(1,J,K),IP(1,K))
 100    CONTINUE
C PROCESS LAST BLOCK-ROW AND RETURN.------------------------------------
      DO 130 J=1,M
        DO 120 I=1,M
          DP=0.D0
          DO 110 L=1,M
 110        DP=DP+B(I,L,N)*B(L,J,NM2)
 120      C(I,J,N)=C(I,J,N)-DP
 130    CONTINUE
      DO 160 J=1,M
        DO 150 I=1,M
          DP=0.D0
          DO 140 L=1,M
 140        DP=DP+C(I,L,N)*B(L,J,NM1)
 150      A(I,J,N)=A(I,J,N)-DP
 160    CONTINUE
 	K = N
      CALL k6DEC (M,M,A(1,1,N),IP(1,N),IER)
      K = N
      IF (IER.NE.0) GOTO 200
      RETURN
C ERROR RETURNS-.-------------------------------------------------------
 200  IER = K
      RETURN
 210  IER = -1
      RETURN
      END
C
