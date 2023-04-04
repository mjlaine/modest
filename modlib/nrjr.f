************************************************************************
*                                                                      *
*     NLE-SOLVER BY NEWTON-RAPHSON METHOD.                             *
*                                                                      *
*     Last modifications in effect        2.2.1995                     *
*                                                                      *
************************************************************************

          Subroutine NRJR(FUN,N,X,Xmin,Xmax,Tol,Res,MaxLL,WA,LWA,IPR,LL,
     &                ipar,IWA,LIWA,DERFAC,ldim)
C
C     solve set of non linear equations by Newton-Raphson method
C
C     FUN         subroutine for equations to be solved
C     IPR         debugging parameter
C     N           number of variables
C     X           variable vector
C     XMIN        lower bounds
C     XMAX        upper bounds
C     TOL         error tolerance input
C     RES         final error after calculations
C     DERFAC      derivation factor in jacobian computation
C     LL          iteration counter
C     MAXLL       max allowed iterations
C     ERROR       error indicator true if solution not found
C     WA          work area for computations
C     ipar        Integer parameter passed to the routine
C     ipar(1) =   info  =     return status indicator
C                 (-1=User stop from FUN,0=OK,1=Improper input,
C                   2=MAXLL reached,3=No solution,4=S<1.d-6)
C     ipar(2) =   MAXGRP=     see subroutine FDJS
C     ipar(3) =   MINGRP=     see subroutine FDJS
C     ipar(4) =   NNZ   =     see subroutine FDJS
C     ipar(5) =   (0=absolute differences,1=relative differences)
C     ipar(6..10) vapaita
C     IWA         Integer work area for the computations
C     LIWA        IWA dimension
C     ldim        Dimension reserved for the sparse Jacobian INDROW and INDCOL
C                 (max = N*N, min = sparsity*N*N) Must be given in the calling
C                 routine/program!
C
C     WA    N     F
C     WA    N*N   DF
C     WA    N     XOLD
C     WA    N     DX
C     WA    N     D
C     WA    N     XD
C     WA    N     FJACD
C     WA    N*N   FJAC
C     WA    N     RIP (in DEC and SOL)
C
C     space required totally N*N+7*N elements from WA vector
C
C     IWA   ldim  INDCOL
C     IWA   ldim  INDROW
C     IWA   N+1   IPNTR
C     IWA   N+1   JPNTR
C     IWA   N     NGRP
C     IWA   6*N   free workspace (here)
C
C     space required totally 2*ldim+9*N+2 elements from IWA vector
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      EXTERNAL    FUN
	INTEGER*4   LWA, N, LIWA
C
C
      DIMENSION   X(N), XMIN(N), XMAX(N), WA(LWA), IPAR(*), IWA(LIWA)
C
C
      LOGICAL*4   LOG, ERROR		   
      
      ipar(1) = 3       ! info=3
      LOG  = .FALSE.
C
      IF(LWA.LT.(N*N+7*N)) THEN
          WRITE( * , 106 ) N*N+7*N, LWA
          LOG = .TRUE.
      END IF
C
c      IF(LIWA.LT.(2*ldim+9*N+2)) THEN
c          WRITE( * , 106 )  2*ldim+9*N+2, LIWA
c          LOG = .TRUE.
c      END IF
C
      DO 1030 I = 1, N
          If( X(I) .LT. XMIN(I) ) THEN
c              WRITE( * , 101) I, X(I), XMIN(I)
              X(I) = XMIN(I)
          ElseIf( X(I) .GT. XMAX(I) ) THEN
c              WRITE( * , 102) I, X(I), XMAX(I)
              X(I) = XMAX(I)
          EndIf
          If(XMIN(I) .GT. XMAX(I)) THEN
            LOG = .TRUE.
            WRITE ( * , 105) XMIN(I), XMAX(I)
          EndIf

1030  CONTINUE
C
C
C
      IF (LOG) Then
          ipar(1) = 1   ! info=1
          Return
      EndIf
C
C
      IFF   = 1
      IDF   = 1 + N
      IXO   = IDF + N*N
      IDX   = IXO + N
      IDW   = IDX + N
      IXD   = IDW + N
      IFD   = IXD + N
      IFJ   = IFD + N*N
      IWX   = IFJ + N
C
      IWCOL = 1
      IWROW = ldim+1
      IWIP  = 2*ldim+1
      IWJP  = 2*ldim+N+2
      IWNG  = 2*ldim+2*N+3
      IWIWA = 2*ldim+3*N+3
C
      CALL NEWRAQ (FUN, IPR, N, X, XMIN, XMAX, TOL, RES, DERFAC,
     &  LL, MAXLL, WA(IFF), WA(IDF), WA(IXO), WA(IDX),
     &  WA(IDW), WA(IXD), WA(IFD), WA(IFJ), WA(IWX),
     &  ERROR, ipar, IWA(IWCOL), IWA(IWROW), IWA(IWIP), IWA(IWJP),
     &  IWA(IWNG), IWA(IWIWA), ldim)
C
      RETURN
C
c101   FORMAT (' New-Rap variable ', I4, ' initial value ', G12.4,
c     &  ' set to minimum bound ', G12.4)
c102   FORMAT (' New-Rap variable ', I4, ' initial value ', G12.4,
c     &  ' set to maximum bound ', G12.4)
105   FORMAT (' New-Rap lower bound ', G12.3,
     &        ' exceeds upper bound'
     &        ,G12.3)
106   FORMAT (' New-Rap requires workspace ', G12.3,
     &' ,more than defined',G12.3)
C
      END


      SUBROUTINE NEWRAQ (FUN, IPR, N, X, XMIN, XMAX,
     &  TOL, RES, DERFAC, LL, MAXLL, F, DF, XOLD, DX,
     &  D, XD, FJACD, FJAC, WA, ERROR, ipar,
     &  INDCOL, INDROW, IPNTR, JPNTR, NGRP, IWA, ldim)
C
      IMPLICIT REAL*8 ( A-H, O-Z )
C
      EXTERNAL    FUN
C
      LOGICAL*4   ERROR
C
C
      DIMENSION   X(N),XMIN(N),XMAX(N), F(N), DF(N,N)
      DIMENSION   XOLD(N), DX(N), IPAR(*)
      DIMENSION   D(N), XD(N), FJACD(N), FJAC(N*N), WA(N)
      DIMENSION   INDCOL(ldim),INDROW(ldim),IPNTR(N+1),JPNTR(N+1)
      DIMENSION   NGRP(N), IWA(6*N)
C
      S  = 0.5D0
      LL = 0
C
      If (IPR .GT. 1) THEN
          Write(*,*) 'N, N*N, ldim, N+1, 6*N'
          Write(*,*)  N, N*N, ldim, N+1, 6*N
      EndIf

      If (IPR .GT. 3) THEN
        do i = 1,n
          WRITE ( * , 111) XMIN(I), XMAX(I)
        enddo
      EndIf
C
2500  CONTINUE
C
C
      If (IPR .GT. 1) THEN
          WRITE( * , 106) LL, N, S, (X(I), I=1, N)
      EndIf
C
      iflg = 1
      CALL FUN (N, X, F, iflg)
C
      IF (IPR .GT. 1) WRITE ( * , 107) (F(I),I = 1,N)
C
      RES = VECNOR (F, N, IPR)
C
      IF (IPR .GT. 0) WRITE ( * , 108) RES
C
      If(iflg.lt.0) Then
          ipar(1) = iflg
          Return
      EndIf
C
      IF (RES .LT. TOL .AND. LL.GT.0) THEN
          ERROR = .FALSE.
          ipar(1) = 0
          RETURN
      ElseIf (LL .GT. MAXLL) THEN
          ERROR = .TRUE.
          ipar(1) = 2
          RETURN
      EndIf
C
      LL = LL + 1
C
C     save X to XOLD and put F to DX
C
      DO 2610 I = 1, N
          DX(I)   = F(I)
2610      XOLD(I) = X(I)
C
C-------------------------- NEW ----------------------------------------
C
C     APPROXIMATE THE JACOBIAN MATRIX.
C
      Do i=1,N
          Do j=1,N
              DF(i,j) = 0.d0
          Enddo
      Enddo

******* Dense jacobian*******
      iflg = 2
      DO 2750 J = 1, N
          XAVER    =   0.5D0*(XMIN(J)+XMAX(J))
          SX   = X(J)
          DXX  = DERFAC * DABS (SX)
          IF (DXX .LT. 1.D-10) DXX = 1.0D-10
          IF (XOLD(J) .GT. XAVER) THEN
              DXX = - DXX
          ELSE
              DXX =  DXX
          EndIF
          X(J) = SX + DXX
C
          CALL FUN (N, X, DF(1,J), iflg)
          X(J) = SX
C
          DO 2730 I = 1, N
2730          DF (I, J) = (DF(I,J) - F(I)) / DXX
C
2750  CONTINUE
C

C
c      MAXGRP = ipar(2)
c      MINGRP = ipar(3)
c      iflg=2
c      DO 2730 NUMGRP = 1, MAXGRP
c          DO 2731 J = 1, N
c              D(J) = 0.d0
c              IF (NGRP(J) .EQ. NUMGRP) THEN
c                  IF(IPAR(5).eq.0) THEN
c                      D(J) = DERFAC
c                  ELSE
c                      D(J) = DERFAC*X(J)
c                  ENDIF
c              ENDIF
c              XD(J) = X(J) + D(J)
c2731      CONTINUE
c          CALL FUN(N,XD,FJACD,iflg)
c          DO 2732 I = 1, N
c              FJACD(I) = FJACD(I) - DX(I)
c2732      CONTINUE
c          CALL FDJS(N,N,.TRUE.,INDROW,JPNTR,NGRP,NUMGRP,D,FJACD,FJAC)
c2730  CONTINUE
c
c      DO 2733 J = 1, N
c          DO 2733 JP = JPNTR(J), JPNTR(J+1)-1
c              DF(INDROW(JP),J) = FJAC(JP)
c2733  CONTINUE

      IF (IPR .GT. 3) THEN
          WRITE ( * , 101)
          DO 2770 I = 1, N
2770          WRITE ( * , 102) I, (DF(I,J), J=1, N)
      EndIF
C
C     LU-decomposition
C
      CALL DEC(N, N, DF, WA, IER)

      IF (IER .GT. 0) THEN
          WRITE ( * , * ) 'Error in New-Rap LU-decomposition. Ier =',IER
          STOP
      EndIF

C     Backsubstitution
C
      CALL SOL(N, N, DF, DX, WA)


      IF (IPR .GT. 2) WRITE ( * , 104) (DX(I), I=1, N)
C
C
      IIS = 0
      IF (S .GT. 1.5D0) S =   1.5D0
C
C     compute new X
C
3500  CONTINUE
C
C
C
C     create new x values
C
      DO 3520 I = 1, N
          XDEV =   - S * DX(I)
          CALL STEPSI (XMIN(I), XMAX(I), XOLD(I), XDEV, IPR, I)
          X(I) =   XOLD(I) + XDEV
3520  CONTINUE
C
C
      IIS = IIS + 1
C
      IF (S .LT. 1.0D-6) THEN
          ERROR   =   .TRUE.
          ipar(1) = 4
          RETURN
      END IF
C
      IF (IPR .GT. 0) WRITE ( * , 109) IIS, S
      IF (IPR .GT. 2) WRITE ( * , 105) (X(I), I=1, N)
      IF (IPR .GT. 3) THEN
        do i = 1,n
          WRITE ( * , 111) XMIN(I), XMAX(I)
        enddo
      EndIF
C
C
C     compute new function value
C
      iflg = 1
      CALL FUN (N, X, F, iflg)
C
      IF ( IPR .GT. 2) WRITE ( * , 107) (F(I),I = 1,N)
C
      RESN = VECNOR  (F, N, IPR)
C
      IF (IPR .GT. 0) WRITE ( * , 110) RESN, RES, TOL
C
      If(iflg.lt.0) Then
          ipar(1) = iflg
          Return
      EndIf
C
      IF (RESN .LT. TOL) THEN
          RES   = RESN
          ERROR = .FALSE.
          ipar(1) = 0
          RETURN
      EndIF
C
C
      IF (IIS .EQ. 1) THEN
          RES1 = RESN
          S1   = S
          IF (RES1 .LT. RES) THEN
              S = 1.2D0 * S
          ELSE
              S = 0.62D0 * S
          EndIF
          GOTO 3500
      ELSE IF (IIS .EQ. 2) THEN
          RES2  = RESN
          S2    = S
          Y1M   = RES1 - RES
          Y2M   = RES2 - RES
          YDEV  = Y1M*S2 - Y2M*S1
          IF (DABS(YDEV) .GT. 1.0D-30) THEN
              S = 0.5D0*(Y1M*S2*S2 - Y2M*S1*S1) / YDEV
          ELSE
              S = 0.5D0 * S
          EndIF
          SMIN  = 0.2*DMIN1(S1, S2, 2.0D0 )
          SMAX  = 2.0*DMAX1(S1, S2, 1.0D-5)
          S     = DMIN1(S, SMAX)
          S     = DMAX1(S, SMIN)
          GOTO 3500
      EndIF
C
      RES3 = RESN
      S3   = S
C
      REST = DMIN1(RES1, RES2, RES3)
C
      IF (RES .LT. REST) THEN
          S = 0.6D0 * DMIN1(S1, S2, S3)
          IF (S .LT. 1.D-7) GO TO 2500
          GOTO 3500
      ELSE
          IF (RES1.LT.RES2 .AND. RES1.LT.RES3) THEN
              S   =   S1
          ElseIF (RES2.LT.RES1 .AND. RES2.LT.RES3) THEN
              S   =   S2
          ELSE
              S   =   S3
          EndIF
C
          DO 3840 I = 1, N
              XDEV =   - S * DX(I)
              CALL STEPSI (XMIN(I), XMAX(I), XOLD(I), XDEV, IPR, I)
              X(I) =   XOLD(I) + XDEV
3840      CONTINUE
C
          GOTO 2500
      EndIF
C
      RETURN
C
101   FORMAT (' Jacobian')
102   FORMAT (' I', I4, (10G12.4) )
104   FORMAT (' dx ', (1X, 8G13.6) )
105   FORMAT (' x  ', (1X, 8G13.6) )
106   FORMAT ('0New-Rap  ll iter ', I4, ' N ', I4, ' S  ', G15.5,
     1   /, (1X, 8G13.6))
107   FORMAT (' New-Rap  f     ',  (8G13.6))
108   FORMAT (' New-Rap  res   ',  (8G13.6))
109   FORMAT (' New-Rap  iis ', I4, '  S ', G15.5)
110   FORMAT (' New-Rap  resn  res tol ',  (3G18.10))
111   FORMAT (' New-Rap xmin  ', (G14.7) ,' New-Rap xmax  ', (G14.7))
C
      END

      REAL*8 FUNCTION VECNOR (F, N, IPR)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
      DIMENSION F(N)
C
C
      VECNOR = 0.0D0
C
      DO 1010 I = 1, N
          VECNOR = VECNOR + F(I)*F(I)
          IF( VECNOR .GT. 0.D0 .AND. F(I) .GT. 0.D0 ) THEN
              IF((2.D0*DLOG10(DABS(F(I))).GT.DLOG10(DABS(D1MACH(2)))-2)
     &        .OR.(DLOG10(DABS(VECNOR)).GT.DLOG10(DABS(D1MACH(2)))-2 ))
     &        WRITE(*,*) 'Vecnor Overflow !',VECNOR,'(',D1MACH(2),')'
          ENDIF
1010  CONTINUE
      vecnor = dsqrt(vecnor)
C
      RETURN
C
      END


      SUBROUTINE STEPSI (XMIN, XMAX, X, DX, IPR, IVAR)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
*******      Call LIMITS(IVAR,XMIN,XMAX,XSTEP)
C      XSTEP = 0.60D0
C
      IF (DX .GT. 0.0D0) THEN
          XINC = DX+1.d0
          If(X+DX.GT.XMAX) THEN
              XSTEP = (1.d0-(XMAX/(X+DX))**4.d0)
              XSTEP = DMAX1(XSTEP,0.6d0)
              XSTEP = DMIN1(XSTEP,0.99d0)
              XINC  = XSTEP * (XMAX - X)
          EndIf
          DX  =   DMIN1 (DX, XINC)
          If(IPR.GT.0.AND.DX.EQ.XINC) Write(*,*) 'XMAX hit with ',IVAR,
     &' while X,DX,XSTEP were ',X,DX,XSTEP
      ELSE
          XDEC = DX-1.d0
          If(X+DX.LT.XMIN) THEN
              XSTEP = (1.d0-((X-XMIN)/(-DX))**4.d0)
              XSTEP = DMAX1(XSTEP,0.6d0)
              XSTEP = DMIN1(XSTEP,0.99d0)
              XDEC  = XSTEP * (XMIN - X)
          EndIf
          DX  =   DMAX1 (DX, XDEC)
          If(IPR.GT.0.AND.DX.EQ.XDEC) Write(*,*) 'XMIN hit with ',IVAR,
     &' while X,DX,XSTEP were ',X,DX,XSTEP
      END IF
C
      RETURN
      END
C
C     LU-DECOMPOSITION OF A QUADRATIC MATRIX
C     **************************************
 
      SUBROUTINE DEC(N,NDIM,A,RIP,IER)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDIM,N),RIP(N)
C

C     ERRORS: IER=1, ONLY ONE EQUATION; IER=2,  ZERO PIVOT
C     ELEMENT

C     GAUSSIAN ELIMINATION
C     --------------------
      IER=0
      RIP(N)=DFLOAT(1)
      IF(N.EQ.1) GOTO 170
      NM1=N-1
      DO 60 K=1,NM1
      KP1=K+1
      M=K
      DO 10 I=KP1,N
      AIK=A(I,K)
      AMK=A(M,K)
   10 IF(DABS(AIK).GT.DABS(AMK))M=I
      RIP(K)=DFLOAT(M)
      T=A(M,K)
      IF(M.EQ.K)GOTO 20
      RIP(N)=-RIP(N)
      A(M,K)=A(K,K)
      A(K,K)=T
   20 IF(T.EQ.0.0) GOTO 160
      T=1.d0/T
      DO 30 I=KP1,N
   30 A(I,K)=-A(I,K)*T
      DO 50 J=KP1,N
      T=A(M,J)
      A(M,J)=A(K,J)
      A(K,J)=T
      IF(T.EQ.0.0) GOTO 50
      DO 40 I=KP1,N
   40 A(I,J)=A(I,J)+A(I,K)*T
   50 CONTINUE
   60 CONTINUE
   70 K=N
      RETURN
  160 IER=2
      RETURN
  170 IER=3
      END



C     BACK SUBSTITUTION (SOLUTION OF A*X=B)
C     *************************************
 
      SUBROUTINE SOL(N,NDIM,A,B,RIP)

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDIM,N),B(N),RIP(N)
C

      IF(N.EQ.1) GOTO 50
      NM1=N-1
      DO 20 K=1,NM1
      KP1=K+1
      M=IDINT(RIP(K))
      T=B(M)
      B(M)=B(K)
      B(K)=T
      DO 10 I=KP1,N
   10 B(I)=B(I)+A(I,K)*T
   20 CONTINUE
      DO 40 KB=1,NM1
      KM1=N-KB
      K=KM1+1
      IF(B(K).NE.0.0)B(K)=B(K)/A(K,K)
      T=-B(K)
      DO 30 I=1,KM1
   30 B(I)=B(I)+A(I,K)*T
   40 CONTINUE
   50 IF(B(1).NE.0.0)B(1)=B(1)/A(1,1)
      RETURN
      END
