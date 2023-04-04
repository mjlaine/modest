
      SUBROUTINE DF (NEQ, T, Y, P, DFDP, JPAR)
      INTEGER NEQ, JPAR
      DOUBLE PRECISION T, Y, P, DFDP
      DIMENSION NEQ(*), Y(*), P(*), DFDP(*)
      GO TO (1,2,3), JPAR
  1   continue
c      DFDP(1) = -Y(1)*Y(1)
c      DFDP(2) = -DFDP(1)
      RETURN
  2   continue
c     DFDP(2) = -Y(2)*Y(2)
c      DFDP(3) = -DFDP(2)
      RETURN
  3   RETURN
      END

