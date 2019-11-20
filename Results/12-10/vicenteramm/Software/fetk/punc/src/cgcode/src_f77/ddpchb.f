      SUBROUTINE DDPCHB(IOUNIT,AA,BB,N,CL,CU,COND,EREPS,IADAPT)
C***BEGIN PROLOGUE  DDPCHB
C***REFER TO  DPPCG,DCGCHB
C***ROUTINES CALLED  D1MACH
C***REVISION DATE  860715   (YYMMDD)
C***END PROLOGUE  DDPCHB
C
C     THIS SUBROUTINE DYNAMICALLY DETERMINES THE NEW PRECONDITIONING 
C     POLYNOMIAL, C(A).  THIS IS DONE BY DETERMINING THE NEW INTERVAL
C     [AA,BB] OVER WHICH TO DEFINE THE SCALED AND TRANSLATED CHEBYSHEV
C     POLYNOMIAL, P(A).  ONCE WE HAVE P(A), WE DEFINE C(A)A = I-P(A).
C     AA AND BB ARE UPDATED VIA INFORMATION OBTAINED BY DONEST, WHICH
C     RETURNS EGVAL ESTIMATES FOR THE CURRENT PRECONDITIONED SYSTEM,
C     C(A)A.  THESE ESTIMATES ARE MAPPED TO EGVAL ESTIMATES FOR A, AND
C     THEN [AA,BB] IS EXPANDED IF NECESSARY.  IF THE NEW CG CONVERGENCE
C     FACTOR (CF) IS LARGER THAN THE CURRENT CF, THEN THE OLD POLY IS
C     RETAINED AND THE ITERATION RESUMED.  OTHERWISE, THE CG ITERATION
C     IS RESTARTED WITH THE NEW PRECONDITIONER.  (NOTE: N, THE DEGREE
C     OF THE CHEBYSHEV POLY, MUST BE ODD FOR THE ADAPTIVE PROCEDURE TO
C     WORK CORRECTLY.  IF N IS EVEN, THE ROUTINE RETURNS IMMEDIATELY.)
C
      IMPLICIT  DOUBLE PRECISION(A-H,O-Z)
C
C***FIRST EXECUTABLE STATEMENT  DDPCHB
 1    CONTINUE
C
C     *** CHECK THAT N IS POSITIVE AND ODD ***
      IF ((N .LE. 0) .OR. (MOD(N,2) .EQ. 0)) THEN
         IADAPT = -1
         RETURN
      END IF
C
C     *** COMPUTE D, C, RECIPN, SQTMEP ***
      D = 0.5D0*(BB+AA)
      C = 0.5D0*(BB-AA)
      RECIPN = 1.0D0/DBLE(N)
      SQTMEP = DSQRT(D1MACH(4))
C
C     *** CHECK FOR SMALL C (RELATIVE TO D) ***
      IF (C .LE. D*SQTMEP) THEN
C        *** C IS SMALL ***
         TAU = D**N
         IF (CL .GT. 1.0D0) CL = 1.0D0
         IF (CU .LT. 1.0D0) CU = 1.0D0
         AN = D - ((1.0D0-CL)*TAU)**RECIPN
         BN = D + ((CU-1.0D0)*TAU)**RECIPN
         GOTO 20
      END IF
C
C     *** COMPUTE POLY DEVIATION FROM 1 OVER (AA,BB) ***
      G = D/C
      TAU = COSH(N*DLOG(G+DSQRT(G*G-1.0D0)))
C
C     *** DETERMINE NEW LEFT ENDPOINT AN ***
      ETA = (1.0D0-CL)*TAU
      IF (ETA .GT. 1.0D0) THEN
C        *** COMPUTE NEW AA ***
         CSHINV = DLOG(ETA + DSQRT(ETA*ETA - 1.0D0))
         AN = D - C*COSH(RECIPN*CSHINV)
      ELSE
C        *** USE OLD AA ***
         AN = AA
         CL = (TAU-1.0D0)/TAU
      END IF
C
C     *** DETERMINE NEW RIGHT ENDPOINT BN ***
      ETA = (CU-1.0D0)*TAU
      IF (ETA .GT. 1.0D0) THEN
C        *** COMPUTE NEW BB ***
         CSHINV = DLOG(ETA + DSQRT(ETA*ETA - 1.0D0))
         BN = D + C*COSH(RECIPN*CSHINV)
      ELSE
C        *** USE OLD BB ***
         BN = BB
         CU = (TAU+1.0D0)/TAU
      END IF
C
C     *** CHECK FOR NO CHANGE IN ENDPOINTS ***  
 20   IF ((AN .EQ. AA) .AND. (BN .EQ. BB)) THEN
C        *** NO NEW ENDPOINTS; SET IADAPT ***
         IADAPT = 0 
         COND = CU/CL
         CF = (DSQRT(COND)-1.0D0) / (DSQRT(COND)+1.0D0)
         IF (IOUNIT .GT. 0) WRITE(IOUNIT,30) CF, COND
 30      FORMAT(' THESE CA EIGENVALUES YIELD NO NEW AA OR BB', /,
     2          ' REFINED CF = ', D12.5, ' AND CONDCA = ', D12.5, /)
         RETURN
      END IF
C
C     *** NEW ENDPOINTS: COMPARE CONVERGENCE FACTORS ***
      CONDO = CU/CL
      GN = (BN+AN) / (BN-AN)
      TAUN = COSH(N*DLOG(GN + DSQRT(GN*GN -1.0D0)))
      CONDN = (TAUN+1.0D0) / (TAUN-1.0D0)
      CFO = (DSQRT(CONDO)-1.0D0) / (DSQRT(CONDO)+1.0D0)
      CFN = (DSQRT(CONDN)-1.0D0) / (DSQRT(CONDN)+1.0D0)
      EPZ = DMAX1(EREPS, SQTMEP)
      TEST = DLOG(EPZ) * (1.0D0/DLOG(CFO) - 1.0D0/DLOG(CFN))
      IF (TEST .LT. 1.0D0) THEN
C        *** RESUME ITERATION WITH OLD AA,BB ***
         IADAPT = 1
         COND = CONDO
         CF = CFO
         IF (IOUNIT .GT. 0) WRITE(IOUNIT,40) AN, BN, CL, CU, CF, COND
 40      FORMAT(' NEW AA = ', D12.5, ' AND NEW BB = ', D12.5, /,
     2          ' CURRENT POLYNOMIAL IS SUPERIOR; RESUME CG' ,/,
     3          ' CURRENT CL = ', D12.5, ' AND CURRENT CU = ', D12.5,/,
     4          ' CURRENT CF = ', D12.5, ' AND CONDCA     = ', D12.5,/)
      ELSE
C        *** RESTART ITERATION WITH NEW AA,BB ***
         IADAPT = 2
         AA = AN
         BB = BN
         CL = (TAUN-1.0D0)/TAUN
         CU = (TAUN+1.0D0)/TAUN
         COND = CONDN
         CF = CFN
         IF (IOUNIT .GT. 0) WRITE(IOUNIT,50) AN, BN, CL, CU, CF, COND
 50      FORMAT(' NEW AA = ', D12.5, ' AND NEW BB = ', D12.5, /,
     2          ' NEW POLYNOMIAL IS SUPERIOR; RESTART CG ', /,
     3          ' NEW CL = ', D12.5, ' AND NEW CU = ', D12.5, /,
     4          ' NEW CF = ', D12.5, ' AND CONDCA = ', D12.5, /)
      END IF    
C
      RETURN
      END
