
      SUBROUTINE JAC (NEQ, T, Y, P, ML, MU, PD, NROWPD)
      INTEGER NEQ, ML, MU, NROWPD, MITER
      DOUBLE PRECISION T, Y, P, PD, PARAM1, PARAM2
      DIMENSION NEQ(*), Y(*), P(*), PD(NROWPD,*)

c
      stop 'you should provide your own JAC'

      MITER = NEQ(3)
      PARAM1 = P(4)*Y(1)
      PARAM2 = P(5)*Y(2)
      IF (MITER .EQ. 4) GO TO 100
c      PD(1,1) = -PARAM1
c      PD(2,1) =  PARAM1
c      PD(2,2) = -PARAM2
c      PD(3,2) =  PARAM2
      RETURN
 100  continue
c     PD(1,1) = -PARAM1
c      PD(2,1) =  PARAM1
c      PD(1,2) = -PARAM2
c      PD(2,2) =  PARAM2
      RETURN
      END


