       REAL*8 FUNCTION PCHQA(N,X,F,D,A,B,IERR)
C***BEGIN PROLOGUE  PCHQA
C***DATE WRITTEN   870829   (YYMMDD)
C***REVISION DATE  870829   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  EASY TO USE CUBIC HERMITE OR SPLINE INTEGRATION
C             NUMERICAL INTEGRATION, QUADRATURE
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the definite integral of a piecewise cubic Hermite
C            or spline function over an arbitrary interval, easy to use.
C***DESCRIPTION
C
C          PCHQA:  Piecewise Cubic Hermite or Spline Integrator, 
C                  Arbitrary limits, Easy to Use.
C
C          From the book "Numerical Methods and Software"
C                  by  D. Kahaner, C. Moler, S. Nash
C                          Prentice Hall 1988
C
C     Evaluates the definite integral of the cubic Hermite or spline
C     function defined by  N, X, F, D  over the interval [A, B].  This 
C     is an easy to use driver for the routine PCHIA by F.N. Fritsch 
C     described in reference (2) below. That routine also has other
C     capabilities.
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C           VALUE = PCHQA (N, X, F, D, A, B, IERR)
C
C     INTEGER  N, IERR 
C     REAL*8  X(N), F(N), D(N), A, B
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The 
C           elements of X must be strictly increasing: 
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.) 
C
C     F -- (input) real array of function values.  F(I) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(I) is
C           the value corresponding to X(I).
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag. 
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 . 
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION 
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHIA
C***END PROLOGUE  PCHQA
      INTEGER  N, IERR 
      REAL*8  X(N), F(N), D(N), A, B
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  INCFD
      REAL*8  PCHIA
      LOGICAL SKIP
C
C  INITIALIZE.
C
      DATA  INCFD /1/
      DATA  SKIP /.TRUE./
C
C
C***FIRST EXECUTABLE STATEMENT  PCHQA

      PCHQA  =  PCHIA( N, X, F, D, INCFD, SKIP, A, B, IERR )
C
C ERROR MESSAGES ARE FROM LOWER LEVEL ROUTINES
      RETURN
C
C------------- LAST LINE OF PCHQA FOLLOWS ------------------------------
      END 
      REAL*8 FUNCTION PCHIA(N,X,F,D,INCFD,SKIP,A,B,IERR)
C***BEGIN PROLOGUE  PCHIA
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  PCHIA
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      REAL*8  X(N), F(INCFD,N), D(INCFD,N), A, B
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IA, IB, IERD, IERV, IL, IR
      REAL*8  VALUE, XA, XB, ZERO
      REAL*8  CHFIV, PCHID
C
C  INITIALIZE.
C
      DATA  ZERO /0./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIA
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (A .EQ. B)  THEN
         VALUE = ZERO
      ELSE
         XA = AMIN1 (A, B)
         XB = AMAX1 (A, B)
         IF (XB .LE. X(2))  THEN
C           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
C                   --------------------------------------------
            VALUE = CHFIV (X(1),X(2), F(1,1),F(1,2),
     *                                D(1,1),D(1,2), A, B, IERV)
C                   --------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE IF (XA .GE. X(N-1))  THEN
C           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
C                   -----------------------------------------------
            VALUE = CHFIV(X(N-1),X(N), F(1,N-1),F(1,N),
     *                                 D(1,N-1),D(1,N), A, B, IERV)
C                   -----------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE
C           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
C      ......LOCATE IA AND IB SUCH THAT
C               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
C             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
C             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
C
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
C             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
C             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
C
C     ......COMPUTE THE INTEGRAL.
            IERV = 0
            IF (IB .LT. IA)  THEN
C              THIS MEANS IB = IA-1 AND
C                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
C                      ------------------------------------------------
               VALUE = CHFIV (X(IB),X(IA), F(1,IB),F(1,IA),
     *                                     D(1,IB),D(1,IA), A, B, IERV)
C                      ------------------------------------------------
               IF (IERV .LT. 0)  GO TO 5004
            ELSE
C
C              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
               IF (IB .EQ. IA)  THEN
                  VALUE = ZERO
               ELSE
C                         ---------------------------------------------
                  VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
C                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5005
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX0 (1, IA-1)
                  IR = IL + 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), XA, X(IA), IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN0 (IB+1, N)
                  IL = IR - 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), X(IB), XB, IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
      PCHIA = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHIA -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHIA -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHIA -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     TROUBLE IN CHFIV.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERROR ('PCHIA -- TROUBLE IN CHFIV'
     *           , 25, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
      IERR = -5
      CALL XERROR ('PCHIA -- TROUBLE IN PCHID'
     *           , 25, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHIA FOLLOWS ------------------------------
      END
      REAL*8 FUNCTION PCHID(N,X,F,D,INCFD,SKIP,IA,IB,IERR)
C***BEGIN PROLOGUE  PCHID
C     THIS PROLOGUE HAS BEEN REMOVED FOR REASONS OF SPACE
C     FOR A COMPLETE COPY OF THIS ROUTINE CONTACT THE AUTHORS
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C***END PROLOGUE  PCHID
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IA, IB, IERR
      REAL*8  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IUP, LOW
      REAL*8  H, HALF, SIX, SUM, VALUE, ZERO
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  HALF /0.5/,  SIX /6./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHID
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (IA .EQ. IB)  THEN
         VALUE = ZERO
      ELSE
         LOW = MIN0(IA, IB)
         IUP = MAX0(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +
     *                      (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
C
C  NORMAL RETURN.
C
      PCHID = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHID -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHID -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHID -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERROR ('PCHID -- IA OR IB OUT OF RANGE'
     *           , 30, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHID FOLLOWS ------------------------------
      END
      REAL*8 FUNCTION CHFIV(X1,X2,F1,F2,D1,D2,A,B,IERR)
C***BEGIN PROLOGUE  CHFIV
C***REFER TO  PCHIA
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***END PROLOGUE  CHFIV
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IERR
      REAL*8  X1, X2, F1, F2, D1, D2, A, B
C
C  DECLARE LOCAL VARIABLES.
C
      REAL*8  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2, PHIB1, PHIB2,
     *      PSIA1, PSIA2, PSIB1, PSIB2, TA1, TA2, TB1, TB2, THREE, TWO,
     *      UA1, UA2, UB1, UB2
C
C  INITIALIZE.
C
      DATA  HALF /0.5/,  TWO /2./,  THREE /3./,  FOUR /4./,  SIX /6./
C
C  VALIDITY CHECK INPUT.
C
C***FIRST EXECUTABLE STATEMENT  CHFIV
      IF (X1 .EQ. X2)  GO TO 5001
      IERR = 0
C
C  COMPUTE INTEGRAL.
C
      H = X2 - X1
      TA1 = (A - X1) / H
      TA2 = (X2 - A) / H
      TB1 = (B - X1) / H
      TB2 = (X2 - B) / H
C
      UA1 = TA1**3
      PHIA1 = UA1 * (TWO - TA1)
      PSIA1 = UA1 * (THREE*TA1 - FOUR)
      UA2 = TA2**3
      PHIA2 =  UA2 * (TWO - TA2)
      PSIA2 = -UA2 * (THREE*TA2 - FOUR)
C
      UB1 = TB1**3
      PHIB1 = UB1 * (TWO - TB1)
      PSIB1 = UB1 * (THREE*TB1 - FOUR)
      UB2 = TB2**3
      PHIB2 =  UB2 * (TWO - TB2)
      PSIB2 = -UB2 * (THREE*TB2 - FOUR)
C
      FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
      DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
C
C  RETURN VALUE.
C
      CHFIV = (HALF*H) * (FTERM + DTERM)
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('CHFIV -- X1 EQUAL TO X2'
     *           , 23, IERR, 1)
      RETURN
C------------- LAST LINE OF CHFIV FOLLOWS ------------------------------
      END
