      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  870930   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERROR processes a diagnostic message. It is a stub routine
C        written for the book above. Actually, XERROR is a sophisticated
C        error handling package with many options, and is described
C        in the reference below. Our version has the same calling sequence
C        but only prints an error message and either returns (if the
C        input value of ABS(LEVEL) is less than 2) or stops (if the
C        input value of ABS(LEVEL) equals 2).
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C                (this is ignored in this stub routine)
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C                (this is ignored in this stub routine)
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C                 (in this stub routine
C                       LEVEL=2 causes a message to be printed and then a 
C                                         stop.
C                       LEVEL<2 causes a message to be printed and then a 
C                                         return.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL nXERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END 
c
c  Muutin nimen 28.10.1999 ML
c
c

      SUBROUTINE nXERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C   FOR REASONS OF SPACE THIS PROLOGUE HAS BEEN OMITTED
C   FOR A COMPLETE COPY CONTACT THE AUTHORS
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERRWV
      WRITE(*,*) MESSG
      IF(NI.EQ.2)THEN
        WRITE(*,*) I1,I2
      ELSEIF(NI.EQ.1) THEN
        WRITE(*,*) I1
      ENDIF
      IF(NR.EQ.2) THEN
        WRITE(*,*) R1,R2
      ELSEIF(NR.EQ.1) THEN
        WRITE(*,*) R1
      ENDIF
      IF(ABS(LEVEL).LT.2)RETURN
      STOP 
      END
