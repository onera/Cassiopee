C ============================================================================
C Solution of periodic linear system
C  ==========================================================================
	SUBROUTINE k6PTRID (M,N,A,B,C,D,Z,ZA,IP)
C
	IMPLICIT NONE
C
C***********************************************************************
C
C	SOLUTION OF PERIODIC BLOCK TRIDIAGONAL SYSTEM
C
C	THIS SUBROUTINE SOLVES A SYSTEM OF THE TYPE
C
C	--                      --   --    --   --    --
C	| B1  C1             A1  |   |  X1  |   |  D1  |
C	| A2  B2  C2             |   |  X2  |   |  D2  | 
C	|                        |   |      |   |      |
C	|     ... ... ...        | * |  ... | = |  ... |
C     |                        |   |      |   |      |
C	|         AN-1 BN-1 CN-1 |   |  XN-1|   |  DN-1|
C	| CN            AN   BN  |   |  XN  |   |  DN  |
C	--                      --   --    --   --    --
C
C       INPUT:
C	     M ORDER OF BLOCKS IN EACH DIRECTION OF THE MATRIX
C	     N ORDER OF LINES OF THE MATRIX
C	     A,B,C  MATRICES MxMxN (DESTROYED DURING COMPUTATION)
C	     X,D    MATRICES MxN (STORED IN THE SAME MATRIX,D AT
C		    INPUT,X AT THE OUTPUT)
C            Z      WORKING MATRIX  MxMxN
C	       ZA     WORKING MATRIX Mx(N-1)
C            IP     WORKING MATRIX FOR PIVOT INFORMATION Mx(N-1)
C	THE SUBROUTINES DECBT,SOLBT,DEC,SOL ARE USED FOR SOLVING THE
C	ASSOCIATED BLOCK TRIDIAGONAL SYSTEMS.
C	SUBROUTINE GELG  IS USED FOR THE CALCULATION OF THE MATRIX XN
C	THE EXECUTION IS STOPPED IF THE BLOCK DECOMPOSITION OR THE 
C	SOLUTION OF THE SYSTEM FOR XN ARE FAILED
C
C***********************************************************************
C
C_IN
	INTEGER_E M,N
	REAL_E A(M,M,N),B(M,M,N),C(M,M,N),D(M,N)
	REAL_E Z(M,M,N),ZA(M,N-1)
	INTEGER_E IP(M,N-1)
C_LOCAL
	REAL_E EPS,S1,S2,T
	INTEGER_E I,J,L,IER,JK,IR
C
	DO 10 I=1,M
	DO 10 J=1,M
	Z(I,J,1)=A(I,J,1)
	A(I,J,1)=0.D0
	Z(I,J,N-1)=C(I,J,N-1)
   10   C(I,J,N-1)=0.D0
C
C
	EPS=1E-10
	
	CALL k6DECBT (M,N-1,B,C,A,IP,IER)
	IF (IER.NE.0.) THEN
	WRITE (6,99) 
   99   FORMAT (1X,'BLOCK DECOMPOSITION FAILED')
	STOP
	ENDIF
C
	DO 11 L=1,N-1
	DO 11 I=1,M
   11   ZA(I,L)=D(I,L)
	CALL k6SOLBT (M,N-1,B,C,A,ZA,IP)
	DO 12 L=1,N-1
	DO 12 I=1,M
   12   D(I,L)=ZA(I,L)
	DO 20 JK=1,M
	DO 21 L=1,N-1
	DO 21 I=1,M
   21   ZA(I,L)=Z(I,JK,L)
	CALL k6SOLBT (M,N-1,B,C,A,ZA,IP)
	DO 22 L=1,N-1
	DO 22 I=1,M
   22   Z(I,JK,L)=ZA(I,L)
   20   CONTINUE
	DO 23 I=1,M
   23   Z(I,1,N)=D(I,N)
C
C
	DO 30 I=1,M
	DO 30 J=1,M
	S1=0.D0
	S2=0.D0
	DO 31 IR=1,M
	S1=S1+C(I,IR,N)*Z(IR,J,1)
   31   S2=S2+A(I,IR,N)*Z(IR,J,N-1)
   30   ZA(I,J)=B(I,J,N)-S1-S2
	DO 40 I=1,M
	S1=0.D0
	S2=0.D0
	DO 41 IR=1,M
	S1=S1+C(I,IR,N)*D(IR,1)
   41   S2=S2+A(I,IR,N)*D(IR,N-1)
   40   D(I,N)=D(I,N)-S1-S2
	CALL k6GELG (D(1,N),ZA,M,1,EPS,IER)
	IF (IER.NE.0.) THEN
	WRITE (6,100)
  100   FORMAT (1X,'SOLUTION OF SYSTEM A*(XN)=B FAILED')
	ENDIF
C
	DO 60 L=1,N-1
	DO 60 I=1,M
	T=0.D0
	DO 61 IR=1,M
   61   T=T+Z(I,IR,L)*D(IR,N)
   60   D(I,L)=D(I,L)-T
	RETURN
	END
C
