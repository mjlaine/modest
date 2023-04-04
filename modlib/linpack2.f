      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END

      INTEGER FUNCTION IDMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDMAX = 0
      IF( N .LT. 1 ) RETURN
      IDMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DX(1)
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DX(IX).LE.DMAX) GO TO 5
         IDMAX = I
         DMAX = DX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DX(1)
      DO 30 I = 2,N
         IF(DX(I).LE.DMAX) GO TO 30
         IDMAX = I
         DMAX = DX(I)
   30 CONTINUE
      RETURN
      END

      INTEGER FUNCTION IDMIN(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMIN
      INTEGER I,INCX,IX,N
C
      IDMIN = 0
      IF( N .LT. 1 ) RETURN
      IDMIN = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMIN = DX(1)
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DX(IX).GE.DMIN) GO TO 5
         IDMIN = I
         DMIN = DX(IX)
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = DX(1)
      DO 30 I = 2,N
         IF(DX(I).GE.DMIN) GO TO 30
         IDMIN = I
         DMIN = DX(I)
   30 CONTINUE
      RETURN
      END


      SUBROUTINE  DSET(N,DA,DX,INCX)
C
C     SETS A VECTOR TO A CONSTANT VALUE.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,IX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX.LT.0) IX = (-N + 1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA
        DX(I + 1) = DA
        DX(I + 2) = DA
        DX(I + 3) = DA
        DX(I + 4) = DA
   50 CONTINUE
      RETURN
      END



      SUBROUTINE  ICOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END


      SUBROUTINE  DHPROD(N,DX,INCX,DY,INCY,DZ,INCZ)
C
C     SETS Z TO HADAMARD PRODUCT OF X AND Y.                            
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DZ(1)
      INTEGER I,INCX,INCY,INCZ,IX,IY,IZ,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1.AND.INCZ.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IZ = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      IF(INCZ.LT.0)IZ = (-N+1)*INCZ + 1
      DO 10 I = 1,N
        DZ(IZ) = DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
        IZ = IZ + INCZ
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DZ(I) = DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DZ(I) = DX(I)*DY(I)
        DZ(I + 1) = DX(I + 1)*DY(I + 1)
        DZ(I + 2) = DX(I + 2)*DY(I + 2)
        DZ(I + 3) = DX(I + 3)*DY(I + 3)
        DZ(I + 4) = DX(I + 4)*DY(I + 4)
        DZ(I + 5) = DX(I + 5)*DY(I + 5)
        DZ(I + 6) = DX(I + 6)*DY(I + 6)
   50 CONTINUE
      RETURN
      END


      DOUBLE PRECISION FUNCTION DSUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,IX
C
      DSUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX.LT.0) IX = (-N + 1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)
        IX = IX + INCX
   10 CONTINUE
      DSUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DX(I) + DX(I + 1) + DX(I + 2)
     *  + DX(I + 3) + DX(I + 4) + DX(I + 5)
   50 CONTINUE
   60 DSUM = DTEMP
      RETURN
      END



